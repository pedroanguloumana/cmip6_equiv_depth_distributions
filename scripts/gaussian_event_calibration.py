#!/usr/bin/env python3
"""Prepare, simulate and summarize a reproducible Gaussian packet calibration.

Use submit_gaussian_event_calibration.sh to run the dependent Slurm stages.
All inputs/outputs enter through data/. Checkpoints preserve complete simulation
chunks, full null candidate scores and the exact source/recipe/seed provenance.
"""

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import inspect
import json
from pathlib import Path
import platform
import socket
import time

import numpy as np
import pandas as pd
import scipy
from scipy.signal import hilbert
from scipy.stats import beta, kurtosis, skew
import xarray as xr
import xr_ccew as tw
from xr_ccew.matsuno import profile_curve_frequency

from src import wave_events as we, wave_profiles as wp, wave_significance as sg
from src.cmip import inspect_cmip_rlut, load_cmip_rlut
from src.paths import DATA_DIR, NOAA_OLR_FILE, PROJECT_ROOT, WAVE_SIGNIFICANCE_DIR


EVENT_COLUMNS = ["event_id", "peak_time", "peak_lon_deg", "duration_days", "span_deg", "peak_amplitude_W_m2"]
SIGNED_EDGES = np.r_[-np.inf, np.linspace(-8,8,161), np.inf]
ENVELOPE_EDGES = np.r_[np.linspace(0,12,241), np.inf]
FREQUENCY_EDGES = np.array([0,1/90,1/30,.1,.2,.4,.500001])
K_EDGES = np.array([0,2,5,10,20,60,181])
# Calibration follows the explicitly selected boxes, independently of the
# legacy sinusoid-fitting worker's profile choices.
PROFILES = wp.BOX_PROFILES


def code_paths():
    return dict(workflow=Path(__file__), significance=Path(sg.__file__), wave_events=Path(we.__file__),
                wave_profiles=Path(wp.__file__),
                xr_filtering=Path(inspect.getfile(tw.filter_field)),
                xr_spectral=Path(inspect.getfile(tw.space_time_fft)),
                xr_matsuno=Path(inspect.getfile(profile_curve_frequency)))


def simulation_hashes(background_method="gaussian", detection_sigma=we.DETECTION_SIGMA):
    # Analysis/plotting changes can reuse saved scores. Changes to the generator,
    # detector, worker or preparation cannot silently mix a simulation ensemble.
    functions = dict(prepare=prepare,simulate=simulate,project=sg.project_olr,
        background=sg.estimate_background_power,draw=sg.draw_background,
        filter=sg.filter_projected_olr,detector=we.detect_events)
    # Archived curved-profile ensembles retain their own snapshots. These box
    # definitions and their source are part of every new ensemble's provenance.
    if background_method=="wk99_121":
        functions.update(prepare_wk99=prepare_wk99,background_wk99=sg.estimate_wk99_background_power)
    elif background_method!="gaussian":
        raise ValueError(f"Unknown background method: {background_method}")
    hashes = {key:hashlib.sha256(inspect.getsource(value).encode()).hexdigest() for key,value in functions.items()}
    hashes.update({key:hashlib.sha256(value.read_bytes()).hexdigest() for key,value in code_paths().items() if key.startswith("xr_")})
    recipe = dict(profiles={wave:profile.as_dict() for wave,profile in PROFILES.items()},
        lat_range=we.LAT_RANGE,average_lat_range=we.AVERAGE_LAT_RANGE,harmonics=we.SEASONAL_HARMONICS,
        taper=we.TAPER_FRACTION,reference_days=we.REFERENCE_EXCLUDE_DAYS,detection_sigma=detection_sigma,
        peak_sigma=detection_sigma,min_days=we.MIN_DURATION_DAYS,min_span=we.MIN_LONGITUDE_SPAN_DEG)
    hashes["recipe"] = hashlib.sha256(json.dumps(recipe,sort_keys=True).encode()).hexdigest()
    hashes["diagnostic_bins"] = hashlib.sha256(SIGNED_EDGES.tobytes()+ENVELOPE_EDGES.tobytes()+FREQUENCY_EDGES.tobytes()+K_EDGES.tobytes()).hexdigest()
    return hashes


def write_csv(frame, path):
    temporary = path.with_name(path.name+".tmp")
    frame.to_csv(temporary,index=False,compression="gzip" if path.name.endswith(".gz") else None)
    temporary.replace(path)


def read_manifest(path):
    path = Path(path).absolute()
    path.relative_to(DATA_DIR)
    manifest = json.loads(path.read_text())
    cfg = manifest["configuration"]
    if simulation_hashes(cfg.get("background_method","gaussian"),
                         cfg["detection"]["detection_sigma"])!=cfg["simulation_sha256"]:
        raise RuntimeError("Simulation recipe/code changed. Restore its snapshot or plan a new ensemble.")
    versions = dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        pandas=pd.__version__,xarray=xr.__version__,xr_ccew=version("xr-ccew"),cftime=version("cftime"))
    if versions!=manifest["configuration"]["versions"]:
        raise RuntimeError("Numerical environment changed since planning.")
    return path, manifest, path.parent


def open_template(root, manifest, wave):
    symmetry = manifest["configuration"]["profiles"][wave]["symmetry"]
    coder = xr.coders.CFDatetimeCoder(use_cftime=manifest["use_cftime"])
    with xr.open_dataarray(root/f"template.{symmetry}.nc",decode_times=coder) as source:
        return source.load()


def plan(args):
    if args.n_calibration<99 or args.n_validation<100 or args.chunk_size<1:
        raise ValueError("Use >=99 calibration records (1% resolution), >=100 validation records and a positive chunk size.")
    if args.start_year>args.end_year or args.seed<0:
        raise ValueError("Require ordered years and a nonnegative seed.")
    if args.detection_sigma not in (1.0, 1.5, 2.0):
        raise ValueError("Use detection_sigma of 1.0, 1.5, or 2.0.")
    if args.model=="NOAA_OLR":
        if args.experiment!="observations" or args.member is not None:
            raise ValueError("NOAA_OLR requires observations and no member.")
        if not NOAA_OLR_FILE.parent.is_symlink() or not NOAA_OLR_FILE.is_file():
            raise FileNotFoundError("Configure data/noaa_olr.")
        path = NOAA_OLR_FILE
        with xr.open_dataset(path) as source:
            calendar = str(source.time.dt.calendar)
        identity = dict(dataset="NOAA_OLR",source_id="NOAA_CDR",experiment_id="observations",
                        member_id="not_applicable",calendar=calendar,source_grid_label="NOAA_1deg",analysis_grid="NOAA_1deg")
    else:
        if args.experiment not in ("historical","amip"):
            raise ValueError("CMIP6 requires historical or amip.")
        info = inspect_cmip_rlut(args.model,experiment=args.experiment)
        if args.member is not None and args.member!=info["member_id"]:
            raise ValueError("Requested member does not match the selected input.")
        path = info["path"]
        identity = dict(dataset="CMIP6",**{key:info[key] for key in ("source_id","experiment_id","member_id","calendar","source_grid_label","analysis_grid")})
    stat = path.stat()
    config = dict(schema_version=2,model=args.model,identity=identity,start_year=args.start_year,end_year=args.end_year,
        profile_collection="BOX_PROFILES",
        source=dict(path=str(path.relative_to(PROJECT_ROOT)),size_bytes=stat.st_size,mtime_ns=stat.st_mtime_ns),
        n_calibration=args.n_calibration,n_validation=args.n_validation,chunk_size=args.chunk_size,seed=args.seed,
        alpha=.01,diagnostic_alphas=[.10,.05,.01],bootstrap_replicates=200,injection_realizations=8,
        frequency_smoothing_cpd=1/96,wavenumber_smoothing=2.0,
        null_model="stationary Gaussian projected anomalies; taper-corrected smoothed spectrum; variance matched in expectation",
        score="candidate maximum longitudinal Hilbert envelope, W m-2",
        testing_scope="one source/experiment/member/wave/native-calendar record; not joint across waves/models; not catalog FDR",
        profiles={wave:profile.as_dict() for wave,profile in PROFILES.items()},
        preprocessing=dict(lat_range=we.LAT_RANGE,average_lat_range=we.AVERAGE_LAT_RANGE,
            latitude_weights="unweighted",seasonal_harmonics=we.SEASONAL_HARMONICS,taper_fraction=we.TAPER_FRACTION,
            remove_mean_and_trend=True,reference_exclude_days=we.REFERENCE_EXCLUDE_DAYS),
        detection=dict(detection_sigma=args.detection_sigma,peak_sigma=args.detection_sigma,
            min_duration_days=we.MIN_DURATION_DAYS,min_longitude_span_deg=we.MIN_LONGITUDE_SPAN_DEG,
            connectivity=8,periodic_longitude=True),
        injection=dict(depths_m=[8.,25.,90.],amplitudes_sigma=[2.,4.,8.],
            wavenumbers={"Kelvin":[4,10],"n1ERW":[-3,-8],"MRG":[-3,-8],"n0EIG":[3,10],"n1WIG":[-3,-8]},
            widths={"compact":[.5,.2],"extended":[1.,.5]}),
        code_sha256={key:hashlib.sha256(value.read_bytes()).hexdigest() for key,value in code_paths().items()},
        simulation_sha256=simulation_hashes(args.background_method,args.detection_sigma),
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                      pandas=pd.__version__,xarray=xr.__version__,xr_ccew=version("xr-ccew"),cftime=version("cftime")))
    if args.background_method=="wk99_121":
        config.pop("frequency_smoothing_cpd")
        config.pop("wavenumber_smoothing")
        config.update(background_method="wk99_121",
            null_model="stationary Gaussian; shared WK99-style 1-2-1 background of projected OLR; pooled segment variance matched in expectation",
            background_estimation=dict(segment_days=96,overlap_days=60,passes_frequency=10,
                wavenumber_passes=[dict(f_min=0.,f_max=.2,passes=10),dict(f_min=.2,f_max=.3,passes=20),
                                   dict(f_min=.3,f_max=.5,passes=40)],
                pooling="mean of symmetric and antisymmetric projected-field segment powers",
                adaptations="catalog meridional means; zero frequency included; daily Nyquist; interpolation to full-record FFT conserving pooled segment variance",
                reference="https://atmos.washington.edu/~breth/classes/AS552/papers/WheelerKiladis1999.pdf"))
    config = json.loads(json.dumps(config))
    digest = hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    name = f"{args.model}.{args.experiment}.{identity['member_id']}.{args.start_year}-{args.end_year}.{digest[:12]}"
    folder = "gaussian_calibration_wk99_121" if args.background_method=="wk99_121" else "gaussian_calibration"
    if args.background_method=="wk99_121" and args.model!="NOAA_OLR":
        # CMIP ensembles are larger than the home quota; configure this data/
        # entry as a machine-local symlink to the intended output filesystem.
        folder = "gaussian_calibration_wk99_121_cmip6"
    root = WAVE_SIGNIFICANCE_DIR/folder/name
    root.mkdir(parents=True,exist_ok=True)
    manifest_path = root/"manifest.json"
    if manifest_path.exists():
        if json.loads(manifest_path.read_text())["configuration_sha256"]!=digest:
            raise RuntimeError("Existing ensemble configuration differs.")
    else:
        tasks = []
        n_total = args.n_calibration+args.n_validation
        for wave in PROFILES:
            for start in range(0,n_total,args.chunk_size):
                tasks.append(dict(task_id=len(tasks),wave=wave,start=start,stop=min(start+args.chunk_size,n_total)))
        we.write_json(manifest_path,dict(status="planned",ensemble_id=digest[:12],configuration_sha256=digest,
            configuration=config,tasks=tasks,created_utc=datetime.now(timezone.utc).isoformat()))
        snapshot = root/"code_snapshot"
        snapshot.mkdir()
        for key,value in code_paths().items():
            (snapshot/f"{key}.py").write_bytes(value.read_bytes())
    print(manifest_path)


def prepare(args):
    path, manifest, root = read_manifest(args.manifest)
    cfg = manifest["configuration"]
    if manifest["status"] in ("prepared","calibrated","complete"):
        print("Already prepared",path)
        return
    source_path = PROJECT_ROOT/cfg["source"]["path"]
    stat = source_path.stat()
    if (stat.st_size,stat.st_mtime_ns)!=(cfg["source"]["size_bytes"],cfg["source"]["mtime_ns"]):
        raise RuntimeError("Input changed since planning.")
    if cfg["model"]=="NOAA_OLR":
        with xr.open_dataset(source_path,chunks={"time":365}) as source:
            olr = source.olr.sel(time=slice(str(cfg["start_year"]),str(cfg["end_year"])),lat=slice(*we.LAT_RANGE))
            if (str(olr.time.values[0])[:10]!=f"{cfg['start_year']}-01-01"
                    or str(olr.time.values[-1])[:10]!=f"{cfg['end_year']}-12-31"):
                raise ValueError("NOAA input does not cover the requested complete years.")
            components = sg.project_olr(olr)
    else:
        olr = load_cmip_rlut(cfg["model"],experiment=cfg["identity"]["experiment_id"],
            member_id=cfg["identity"]["member_id"],year_range=(cfg["start_year"],cfg["end_year"]))
        components = sg.project_olr(olr)
        del olr
    diagnostics = []
    for symmetry in ("symmetric","antisymmetric"):
        field = components[symmetry]
        field.to_netcdf(root/f"template.{symmetry}.nc")
        power = sg.estimate_background_power(field,frequency_smoothing_cpd=cfg["frequency_smoothing_cpd"],
                                              wavenumber_smoothing=cfg["wavenumber_smoothing"])
        anomaly = tw.remove_mean_and_linear_trend(tw.remove_harmonics_of_seasonal_cycle(field,num_harmonics=we.SEASONAL_HARMONICS))
        tapered = tw.apply_window(anomaly,window="tukey",pct=we.TAPER_FRACTION)
        taper = tw.apply_window(xr.ones_like(field.isel(lon=0)),window="tukey",pct=we.TAPER_FRACTION)
        raw_power = power.copy(data=abs(np.fft.fft2(tapered.values))**2/(field.size*float((taper**2).mean())))
        raw_power.attrs = dict(units="(W m-2)^2",normalization="abs(FFT(tapered anomaly))**2 / N / mean(taper**2)")
        xr.Dataset(dict(background_power=power,observed_taper_corrected_power=raw_power)).to_netcdf(root/f"spectrum.{symmetry}.nc")
        for row in sg.distribution_diagnostics(anomaly.isel(time=slice(we.REFERENCE_EXCLUDE_DAYS,-we.REFERENCE_EXCLUDE_DAYS))):
            diagnostics.append(dict(stage=symmetry+"_before_filter",**row))
    for wave in PROFILES:
        field = components[PROFILES[wave].symmetry]
        filtered = sg.filter_projected_olr(field,wave,profile=PROFILES[wave])
        events,sigma = we.detect_events(filtered,detection_sigma=cfg["detection"]["detection_sigma"])
        frame = pd.DataFrame(events,columns=EVENT_COLUMNS)
        frame["reference_sigma_W_m2"] = sigma
        write_csv(frame,root/f"observed_candidates.{wave}.csv")
        filtered.rename("filtered_olr").to_netcdf(root/f"observed_filtered.{wave}.nc")
        for row in sg.distribution_diagnostics(filtered.isel(time=slice(we.REFERENCE_EXCLUDE_DAYS,-we.REFERENCE_EXCLUDE_DAYS))):
            diagnostics.append(dict(stage=wave,**row))
    write_csv(pd.DataFrame(diagnostics),root/"observed_distributions.csv")
    manifest.update(status="prepared",use_cftime=not np.issubdtype(components.time.dtype,np.datetime64),
        record_samples=components.sizes["time"],longitude_samples=components.sizes["lon"],
        start=str(components.time.values[0]),end=str(components.time.values[-1]),
        prepared_utc=datetime.now(timezone.utc).isoformat())
    we.write_json(path,manifest)
    print("Prepared",path,flush=True)


def prepare_wk99(args):
    """Prepare the shared 96-day, 1-2-1 background without changing old runs."""
    path,manifest,root = read_manifest(args.manifest)
    cfg = manifest["configuration"]
    if manifest["status"] in ("prepared","calibrated","complete"):
        print("Already prepared",path)
        return
    source_path = PROJECT_ROOT/cfg["source"]["path"]
    stat = source_path.stat()
    if (stat.st_size,stat.st_mtime_ns)!=(cfg["source"]["size_bytes"],cfg["source"]["mtime_ns"]):
        raise RuntimeError("Input changed since planning.")
    if cfg["model"]=="NOAA_OLR":
        with xr.open_dataset(source_path,chunks={"time":365}) as source:
            olr = source.olr.sel(time=slice(str(cfg["start_year"]),str(cfg["end_year"])),lat=slice(*we.LAT_RANGE))
            if (str(olr.time.values[0])[:10]!=f"{cfg['start_year']}-01-01"
                    or str(olr.time.values[-1])[:10]!=f"{cfg['end_year']}-12-31"):
                raise ValueError("NOAA input does not cover the requested complete years.")
            components = sg.project_olr(olr)
    else:
        olr = load_cmip_rlut(cfg["model"],experiment=cfg["identity"]["experiment_id"],
            member_id=cfg["identity"]["member_id"],year_range=(cfg["start_year"],cfg["end_year"]))
        components = sg.project_olr(olr)
        del olr
    background = sg.estimate_wk99_background_power(components)
    background.to_netcdf(root/"wk99_background.nc")
    we.write_json(root/"background_diagnostics.json",dict(background.attrs))
    diagnostics = []
    for symmetry in ("symmetric","antisymmetric"):
        field = components[symmetry]
        field.to_netcdf(root/f"template.{symmetry}.nc")
        anomaly = tw.remove_mean_and_linear_trend(tw.remove_harmonics_of_seasonal_cycle(field,num_harmonics=we.SEASONAL_HARMONICS))
        tapered = tw.apply_window(anomaly,window="tukey",pct=we.TAPER_FRACTION)
        taper = tw.apply_window(xr.ones_like(field.isel(lon=0)),window="tukey",pct=we.TAPER_FRACTION)
        raw = background.background_power.copy(data=abs(np.fft.fft2(tapered.values))**2/(field.size*float((taper**2).mean())))
        raw.attrs = dict(units="(W m-2)^2",normalization="abs(FFT(tapered anomaly))**2 / N / mean(taper**2)")
        xr.Dataset(dict(background_power=background.background_power,observed_taper_corrected_power=raw)).to_netcdf(root/f"spectrum.{symmetry}.nc")
        for row in sg.distribution_diagnostics(anomaly.isel(time=slice(we.REFERENCE_EXCLUDE_DAYS,-we.REFERENCE_EXCLUDE_DAYS))):
            diagnostics.append(dict(stage=symmetry+"_before_filter",**row))
    for wave in PROFILES:
        field = components[PROFILES[wave].symmetry]
        filtered = sg.filter_projected_olr(field,wave,profile=PROFILES[wave])
        events,sigma = we.detect_events(filtered,detection_sigma=cfg["detection"]["detection_sigma"])
        frame = pd.DataFrame(events,columns=EVENT_COLUMNS)
        frame["reference_sigma_W_m2"] = sigma
        write_csv(frame,root/f"observed_candidates.{wave}.csv")
        filtered.rename("filtered_olr").to_netcdf(root/f"observed_filtered.{wave}.nc")
        for row in sg.distribution_diagnostics(filtered.isel(time=slice(we.REFERENCE_EXCLUDE_DAYS,-we.REFERENCE_EXCLUDE_DAYS))):
            diagnostics.append(dict(stage=wave,**row))
    write_csv(pd.DataFrame(diagnostics),root/"observed_distributions.csv")
    manifest.update(status="prepared",use_cftime=not np.issubdtype(components.time.dtype,np.datetime64),
        record_samples=components.sizes["time"],longitude_samples=components.sizes["lon"],
        start=str(components.time.values[0]),end=str(components.time.values[-1]),
        prepared_utc=datetime.now(timezone.utc).isoformat())
    we.write_json(path,manifest)
    print("Prepared WK99",path,flush=True)


def simulate(args):
    _, manifest, root = read_manifest(args.manifest)
    if manifest["status"] not in ("prepared","calibrated","complete"):
        raise RuntimeError("Preparation has not completed.")
    task = manifest["tasks"][args.task]
    if task["task_id"]!=args.task:
        raise ValueError("Unexpected task index.")
    cfg, wave = manifest["configuration"], task["wave"]
    chunks = root/"chunks"
    chunks.mkdir(exist_ok=True)
    prefix = chunks/f"task_{args.task:05d}"
    done = prefix.with_suffix(".json")
    if done.exists():
        old = json.loads(done.read_text())
        if old.get("status")=="complete" and old.get("ensemble_id")==manifest["ensemble_id"]:
            print("Already complete",args.task)
            return
        raise RuntimeError("Unexpected existing chunk status.")
    template = open_template(root,manifest,wave)
    symmetry = cfg["profiles"][wave]["symmetry"]
    with xr.open_dataset(root/f"spectrum.{symmetry}.nc") as source:
        power = source.background_power.load()
    fbin = np.searchsorted(FREQUENCY_EDGES,abs(np.fft.fftfreq(template.sizes["time"])),side="right")-1
    kbin = np.searchsorted(K_EDGES,abs(np.fft.fftfreq(template.sizes["lon"])*template.sizes["lon"]),side="right")-1
    band_index = (fbin[:,None]*(len(K_EDGES)-1)+kbin[None,:]).ravel()
    target_bands = np.bincount(band_index,weights=power.values.ravel(),minlength=36)/template.size
    event_frames, records, bands, signed_histograms, envelope_histograms = [], [], [], [], []
    started = time.monotonic()
    for realization in range(task["start"],task["stop"]):
        tick = time.monotonic()
        role = "calibration" if realization<cfg["n_calibration"] else "validation"
        seed = np.random.SeedSequence([cfg["seed"],0 if symmetry=="symmetric" else 1,realization])
        background = sg.draw_background(power,template,np.random.default_rng(seed))
        actual_bands = np.bincount(band_index,weights=(abs(np.fft.fft2(background.values))**2/template.size).ravel(),minlength=36)/template.size
        for fi in range(len(FREQUENCY_EDGES)-1):
            for ki in range(len(K_EDGES)-1):
                index = fi*(len(K_EDGES)-1)+ki
                bands.append(dict(realization=realization,role=role,wave=wave,f_min_cpd=FREQUENCY_EDGES[fi],
                    f_max_cpd=FREQUENCY_EDGES[fi+1],abs_K_min=K_EDGES[ki],abs_K_max=K_EDGES[ki+1],
                    target_variance_W2_m4=target_bands[index],realized_variance_W2_m4=actual_bands[index]))
        filtered = sg.filter_projected_olr(background,wave,profile=PROFILES[wave])
        events,sigma = we.detect_events(filtered,detection_sigma=cfg["detection"]["detection_sigma"])
        frame = pd.DataFrame(events,columns=EVENT_COLUMNS)
        frame["realization"],frame["role"],frame["wave"] = realization,role,wave
        event_frames.append(frame)
        values = filtered.isel(time=slice(we.REFERENCE_EXCLUDE_DAYS,-we.REFERENCE_EXCLUDE_DAYS)).values
        envelope = abs(hilbert(filtered.values,axis=1))
        records.append(dict(realization=realization,role=role,wave=wave,candidates=len(events),
            maximum_candidate_peak_W_m2=max((e["peak_amplitude_W_m2"] for e in events),default=0.0),
            reference_sigma_W_m2=sigma,filtered_mean_W_m2=float(values.mean()),
            filtered_skewness=float(skew(values,axis=None)),filtered_excess_kurtosis=float(kurtosis(values,axis=None)),
            elapsed_seconds=time.monotonic()-tick))
        signed_histograms.append(np.histogram(values/sigma,bins=SIGNED_EDGES)[0])
        envelope_histograms.append(np.histogram(envelope/sigma,bins=ENVELOPE_EDGES)[0])
        if realization==0:
            xr.Dataset(dict(background=background,filtered_olr=filtered,envelope=filtered.copy(data=envelope))).to_netcdf(root/f"example_null.{wave}.nc")
        if (realization-task["start"]+1)%10==0:
            print(f"{wave} task {args.task}: {realization-task['start']+1}/{task['stop']-task['start']} records",flush=True)
    write_csv(pd.concat(event_frames,ignore_index=True),prefix.with_suffix(".events.csv.gz"))
    write_csv(pd.DataFrame(records),prefix.with_suffix(".records.csv"))
    write_csv(pd.DataFrame(bands),prefix.with_suffix(".bands.csv.gz"))
    hist = xr.Dataset(dict(signed_count=(("realization","signed_bin"),np.asarray(signed_histograms)),
        envelope_count=(("realization","envelope_bin"),np.asarray(envelope_histograms)),
        signed_edges=("signed_edge",SIGNED_EDGES),envelope_edges=("envelope_edge",ENVELOPE_EDGES)),
        coords=dict(realization=np.arange(task["start"],task["stop"])),
        attrs=dict(wave=wave,signed_units="filtered OLR / per-record reference sigma; reference time interval only",
                   envelope_units="Hilbert envelope / per-record reference sigma; complete record"))
    hist.to_netcdf(prefix.with_suffix(".histograms.nc"))
    we.write_json(done,dict(status="complete",ensemble_id=manifest["ensemble_id"],task=task,
        elapsed_seconds=time.monotonic()-started,host=socket.gethostname(),completed_utc=datetime.now(timezone.utc).isoformat()))
    print("Complete",args.task,flush=True)


def summarize(args):
    path, manifest, root = read_manifest(args.manifest)
    cfg = manifest["configuration"]
    n_total = cfg["n_calibration"]+cfg["n_validation"]
    threshold_rows, convergence_rows, all_records = [], [], []
    for wave in cfg["profiles"]:
        vectors = [np.empty(0) for _ in range(n_total)]
        record_frames, band_frames, histograms = [], [], []
        for task in (item for item in manifest["tasks"] if item["wave"]==wave):
            prefix = root/"chunks"/f"task_{task['task_id']:05d}"
            status = json.loads(prefix.with_suffix(".json").read_text())
            if status.get("status")!="complete" or status.get("ensemble_id")!=manifest["ensemble_id"] or status["task"]!=task:
                raise RuntimeError(f"Incomplete or mismatched simulation chunk: {prefix}")
            record_frames.append(pd.read_csv(prefix.with_suffix(".records.csv")))
            scores = pd.read_csv(prefix.with_suffix(".events.csv.gz"),usecols=["realization","peak_amplitude_W_m2"])
            for index, group in scores.groupby("realization"):
                if not task["start"]<=index<task["stop"]:
                    raise ValueError("Candidate realization lies outside its chunk.")
                vectors[int(index)] = np.sort(group.peak_amplitude_W_m2.to_numpy())
            band_frames.append(pd.read_csv(prefix.with_suffix(".bands.csv.gz")))
            with xr.open_dataset(prefix.with_suffix(".histograms.nc")) as source:
                histograms.append(source.load())
        records = pd.concat(record_frames,ignore_index=True).sort_values("realization").reset_index(drop=True)
        np.testing.assert_array_equal(records.realization,np.arange(n_total))
        expected_roles = np.where(np.arange(n_total)<cfg["n_calibration"],"calibration","validation")
        np.testing.assert_array_equal(records.role,expected_roles)
        np.testing.assert_array_equal(records.candidates,[len(values) for values in vectors])
        np.testing.assert_allclose(records.maximum_candidate_peak_W_m2,
                                   [values[-1] if len(values) else 0 for values in vectors])
        all_records.append(records)
        observed = pd.read_csv(root/f"observed_candidates.{wave}.csv")
        observed_peaks = observed.peak_amplitude_W_m2.to_numpy()
        maxima = records.maximum_candidate_peak_W_m2.to_numpy()
        cal_maxima = maxima[:cfg["n_calibration"]]
        calibration_peaks = np.sort(np.concatenate(vectors[:cfg["n_calibration"]]))
        if not calibration_peaks.size:
            raise RuntimeError("No calibration candidates; inspect the source/background/filter.")
        observed["record_max_p_mc"] = sg.record_max_pvalues(cal_maxima,observed_peaks)
        observed["null_candidate_peak_tail_fraction"] = (len(calibration_peaks)-np.searchsorted(calibration_peaks,observed_peaks,side="left"))/len(calibration_peaks)
        for key,value in cfg["identity"].items(): observed[key] = value
        observed["wave_type"],observed["ensemble_id"] = wave,manifest["ensemble_id"]
        cuts = [sg.record_max_threshold(cal_maxima,alpha) for alpha in cfg["diagnostic_alphas"]]
        pooled_cuts = np.quantile(calibration_peaks,[.9,.95,.99,.995,.999],method="higher")
        top = max(1.,float(maxima.max()),float(observed_peaks.max()) if observed_peaks.size else 0.)
        grid = np.unique(np.r_[np.linspace(0,top,513),observed_peaks,pooled_cuts,
                               [cut["threshold_W_m2"] for cut in cuts if cut["available"]]])
        null_counts = np.asarray([len(values)-np.searchsorted(values,grid,side="right") for values in vectors],dtype=np.int32)
        validation_counts = null_counts[cfg["n_calibration"]:]
        n_validation = cfg["n_validation"]
        hits = (validation_counts>0).sum(axis=0)
        probability_low = np.where(hits==0,0.,beta.ppf(.025,hits,n_validation-hits+1))
        probability_high = np.where(hits==n_validation,1.,beta.ppf(.975,hits+1,n_validation-hits))
        rng = np.random.default_rng(cfg["seed"]+100)
        weights = rng.multinomial(n_validation,np.full(n_validation,1/n_validation),size=cfg["bootstrap_replicates"])/n_validation
        boot_means = weights@validation_counts.astype(float)
        mean_count_low,mean_count_high = np.quantile(boot_means,[.025,.975],axis=0)
        mean_count_high[hits==0] = np.nan  # No claim of a zero upper bound for unobserved rare events.
        sorted_observed = np.sort(observed_peaks)
        observed_counts = len(sorted_observed)-np.searchsorted(sorted_observed,grid,side="right")
        curve = pd.DataFrame(dict(wave_type=wave,threshold_W_m2=grid,observed_candidates_gt=observed_counts,
            validation_mean_false_candidates=validation_counts.mean(axis=0),mean_false_count_MC_low=mean_count_low,
            mean_false_count_MC_high=mean_count_high,validation_probability_any=hits/n_validation,
            probability_any_CI_low=probability_low,probability_any_CI_high=probability_high,
            no_validation_exceedances=hits==0))
        curve["model_estimated_false_fraction"] = np.divide(curve.validation_mean_false_candidates,observed_counts,
            out=np.full(len(grid),np.nan),where=observed_counts>0)
        write_csv(curve,root/f"threshold_curve.{wave}.csv")
        count_dataset = xr.Dataset(dict(null_candidates_gt=(("realization","threshold_W_m2"),null_counts),
            role=("realization",expected_roles)),coords=dict(realization=np.arange(n_total),threshold_W_m2=grid),
            attrs=dict(wave_type=wave,comparison="strict peak > threshold",scope=cfg["testing_scope"]))
        count_dataset.to_netcdf(root/f"null_threshold_counts.{wave}.nc")
        for cut in cuts:
            if not cut["available"]: continue
            threshold = cut["threshold_W_m2"]
            index = np.searchsorted(grid,threshold)
            flag = observed.record_max_p_mc<=cut["alpha"]
            np.testing.assert_array_equal(flag,observed_peaks>threshold)
            observed[f"significant_record_{cut['alpha']:g}"] = flag
            row = dict(cfg["identity"],wave_type=wave,ensemble_id=manifest["ensemble_id"],**cut,
                observed_candidates=len(observed),observed_above_threshold=int(flag.sum()),
                validation_records=n_validation,validation_records_with_false_event=int(hits[index]),
                validation_probability_any=float(hits[index]/n_validation),probability_any_CI_low=float(probability_low[index]),
                probability_any_CI_high=float(probability_high[index]),
                validation_mean_false_candidates=float(validation_counts[:,index].mean()),
                mean_false_count_MC_low=float(mean_count_low[index]),
                mean_false_count_MC_high=float(mean_count_high[index]) if np.isfinite(mean_count_high[index]) else None)
            threshold_rows.append(row)
        observed["significant"] = observed.record_max_p_mc<=cfg["alpha"]
        write_csv(observed,root/f"observed_significance.{wave}.csv")
        ecdf = []
        for role,values in (("calibration",cal_maxima),("validation",maxima[cfg["n_calibration"]:])):
            ordered = np.sort(values)
            ecdf.extend(dict(wave_type=wave,role=role,maximum_candidate_peak_W_m2=value,ecdf=(i+1)/len(ordered)) for i,value in enumerate(ordered))
        write_csv(pd.DataFrame(ecdf),root/f"record_max_ecdf.{wave}.csv")
        for n in sorted({99,100,200,400,600,800,cfg["n_calibration"]}):
            if n>cfg["n_calibration"]: continue
            for alpha in cfg["diagnostic_alphas"]:
                cut = sg.record_max_threshold(cal_maxima[:n],alpha)
                if cut["available"]:
                    convergence_rows.append(dict(wave_type=wave,**cut,
                        observed_above_threshold=int(np.count_nonzero(observed_peaks>cut["threshold_W_m2"])),
                        validation_probability_any=float(np.mean(maxima[cfg["n_calibration"]:]>cut["threshold_W_m2"]))))
        xr.concat(histograms,dim="realization",data_vars="minimal",coords="minimal",compat="equals").sortby("realization").to_netcdf(root/f"null_histograms.{wave}.nc")
        band_data = pd.concat(band_frames,ignore_index=True)
        write_csv(band_data,root/f"spectral_bands.{wave}.csv.gz")
        group_columns = ["role","wave","f_min_cpd","f_max_cpd","abs_K_min","abs_K_max"]
        summary = band_data.groupby(group_columns,as_index=False).agg(
            target_variance_W2_m4=("target_variance_W2_m4","first"),
            realized_variance_mean_W2_m4=("realized_variance_W2_m4","mean"),
            realized_variance_std_W2_m4=("realized_variance_W2_m4","std"),records=("realization","size"))
        summary["mean_over_target"] = summary.realized_variance_mean_W2_m4/summary.target_variance_W2_m4
        write_csv(summary,root/f"spectral_summary.{wave}.csv")
        print("Calibrated",wave,flush=True)
    write_csv(pd.concat(all_records,ignore_index=True),root/"null_records.csv")
    write_csv(pd.DataFrame(threshold_rows),root/"thresholds.csv")
    write_csv(pd.DataFrame(convergence_rows),root/"threshold_convergence.csv")
    analysis_snapshot = root/"analysis_code_snapshot"
    analysis_snapshot.mkdir(exist_ok=True)
    for key,value in code_paths().items(): (analysis_snapshot/f"{key}.py").write_bytes(value.read_bytes())
    we.write_json(root/"calibration_report.json",dict(status="calibrated",ensemble_id=manifest["ensemble_id"],
        configuration=cfg,thresholds=threshold_rows,
        analysis_code_sha256={key:hashlib.sha256(value.read_bytes()).hexdigest() for key,value in code_paths().items()},
        interpretation="Significance of a candidate envelope maximum under the fitted Gaussian full-record null; not physical wave identity or catalog FDR.",
        uncertainty="Quantile intervals and validation binomial intervals condition on the fitted null; count intervals bootstrap entire validation records. No upper mean-count interval is claimed when validation has no exceedances."))
    manifest.update(status="calibrated",calibrated_utc=datetime.now(timezone.utc).isoformat())
    we.write_json(path,manifest)


def inject(args):
    _, manifest, root = read_manifest(args.manifest)
    if manifest["status"] not in ("calibrated","complete"):
        raise RuntimeError("Threshold reduction has not completed.")
    cfg = manifest["configuration"]
    thresholds = pd.read_csv(root/"thresholds.csv")
    records = pd.read_csv(root/"null_records.csv")
    identity = cfg["identity"]
    results = []
    waves = [wave for wave,profile in cfg["profiles"].items() if profile["symmetry"]==args.symmetry]
    template = open_template(root,manifest,waves[0])
    with xr.open_dataset(root/f"spectrum.{args.symmetry}.nc") as source:
        power = source.background_power.load()
    backgrounds = {}
    for realization in range(cfg["n_calibration"],cfg["n_calibration"]+cfg["injection_realizations"]):
        rng = np.random.default_rng(np.random.SeedSequence([cfg["seed"],0 if args.symmetry=="symmetric" else 1,realization]))
        field = sg.draw_background(power,template,rng)
        for wave in waves:
            filtered = sg.filter_projected_olr(field,wave,profile=PROFILES[wave])
            events,_ = we.detect_events(filtered,detection_sigma=cfg["detection"]["detection_sigma"])
            backgrounds[wave,realization] = (filtered,events)
    for wave in waves:
        cut = thresholds.loc[(thresholds.wave_type==wave)&np.isclose(thresholds.alpha,cfg["alpha"])].iloc[0]
        calibration = records.loc[(records.wave==wave)&(records.role=="calibration")]
        maxima = calibration.sort_values("realization").maximum_candidate_peak_W_m2.to_numpy()
        sigma = float(calibration.reference_sigma_W_m2.mean())
        case = 0
        for depth in cfg["injection"]["depths_m"]:
            for K in cfg["injection"]["wavenumbers"][wave]:
                for width,(time_width,lon_width) in cfg["injection"]["widths"].items():
                    packet,truth = sg.synthetic_wave_packet(template,wave,K,depth,
                        center_index=template.sizes["time"]//2+(case%5-2)*120,
                        center_lon=[60.5,180.5,300.5,359.5][case%4],
                        time_sigma_periods=time_width,lon_sigma_wavelengths=lon_width,profile=PROFILES[wave])
                    filtered_packet = sg.filter_projected_olr(packet,wave,profile=PROFILES[wave])
                    energy = float((filtered_packet**2).sum()/(packet**2).sum())
                    peak = float(abs(hilbert(filtered_packet.values,axis=1)).max())
                    for amplitude_sigma in cfg["injection"]["amplitudes_sigma"]:
                        amplitude = amplitude_sigma*sigma
                        for (background_wave,realization),(background,events) in backgrounds.items():
                            if background_wave!=wave: continue
                            combined,_ = we.detect_events(background+amplitude*filtered_packet,
                                detection_sigma=cfg["detection"]["detection_sigma"])
                            before = sg.match_packet_peaks(events,template,truth,cut.threshold_W_m2)
                            after = sg.match_packet_peaks(combined,template,truth,cut.threshold_W_m2)
                            results.append(dict(truth,case=case,width=width,realization=realization,
                                amplitude_sigma=amplitude_sigma,input_amplitude_W_m2=amplitude,
                                carrier_below_daily_nyquist=truth["frequency_cpd"]<.5,
                                retained_energy_fraction=energy,noiseless_filtered_peak_W_m2=amplitude*peak,
                                threshold_W_m2=float(cut.threshold_W_m2),
                                baseline_p_record=float(sg.record_max_pvalues(maxima,[before["nearby_maximum_W_m2"]])[0]),
                                injected_p_record=float(sg.record_max_pvalues(maxima,[after["nearby_maximum_W_m2"]])[0]),
                                **{"baseline_"+key:value for key,value in before.items()},
                                **{"injected_"+key:value for key,value in after.items()}))
                    case += 1
                    print(f"Injection {wave} {case}/12",flush=True)
                    write_csv(pd.DataFrame(results),root/f"injections.{args.symmetry}.csv.gz")
    we.write_json(root/f"injections.{args.symmetry}.json",dict(status="complete",ensemble_id=manifest["ensemble_id"],
        trials=len(results),identity=identity,replication=cfg["injection_realizations"],
        interpretation="Recovery of controlled modulated carriers, not fitted k/omega or depth accuracy; baseline matches identify chance associations."))


def finalize(args):
    path, manifest, root = read_manifest(args.manifest)
    cfg = manifest["configuration"]
    frames = []
    for symmetry in ("symmetric","antisymmetric"):
        status = json.loads((root/f"injections.{symmetry}.json").read_text())
        if status["status"]!="complete" or status["ensemble_id"]!=manifest["ensemble_id"]:
            raise RuntimeError("Injection stage incomplete or mismatched.")
        frames.append(pd.read_csv(root/f"injections.{symmetry}.csv.gz"))
    injections = pd.concat(frames,ignore_index=True)
    write_csv(injections,root/"injections.csv.gz")
    injections["in_band_and_resolved"] = injections.carrier_inside_rectangular_bounds & injections.carrier_below_daily_nyquist
    rows = []
    grouping = ["wave_type","h_eff_m","width","amplitude_sigma","in_band_and_resolved"]
    for keys,group in injections.groupby(grouping):
        for alpha in cfg["diagnostic_alphas"]:
            before,after = group.baseline_p_record<=alpha,group.injected_p_record<=alpha
            rows.append(dict(zip(grouping,keys),alpha=alpha,trials=len(group),
                significant_baseline=int(before.sum()),significant_injected=int(after.sum()),
                detection_fraction=float(group.injected_nearby_candidates.gt(0).mean()),
                significance_recovery_fraction=float(after.mean()),baseline_significance_fraction=float(before.mean()),
                paired_change=float(after.astype(int).sub(before.astype(int)).mean()),
                mean_retained_packet_energy=float(group.retained_energy_fraction.mean()),
                mean_noiseless_filtered_peak_W_m2=float(group.noiseless_filtered_peak_W_m2.mean())))
    recovery = pd.DataFrame(rows)
    write_csv(recovery,root/"injection_recovery.csv")
    report = json.loads((root/"calibration_report.json").read_text())
    report.update(status="complete",completed_utc=datetime.now(timezone.utc).isoformat(),
        injection_trials=len(injections),injection_recovery_table="injection_recovery.csv",
        saved_products=dict(thresholds="thresholds.csv",observed_flags="observed_significance.<wave>.csv",
            null_scores="chunks/task_*.events.csv.gz",null_records="null_records.csv",convergence="threshold_convergence.csv",
            threshold_curves="threshold_curve.<wave>.csv",null_counts="null_threshold_counts.<wave>.nc",
            distributions="observed_distributions.csv; null_histograms.<wave>.nc",spectra="spectrum.<symmetry>.nc; spectral_summary.<wave>.csv",
            representative_fields="observed_filtered.<wave>.nc; example_null.<wave>.nc",injections="injections.csv.gz; injection_recovery.csv"),
        limitations="Conditional on a stationary Gaussian smoothed-spectrum null and these filter/detector settings; full-record significance is per source/wave, not catalog FDR or joint across all models/waves. Null parameter uncertainty, nonstationary variance, physical wave identity and fitted propagation/depth errors are separate issues.")
    we.write_json(root/"report.json",report)
    manifest.update(status="complete",completed_utc=report["completed_utc"])
    we.write_json(path,manifest)
    print(pd.read_csv(root/"thresholds.csv").to_string(index=False),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="stage",required=True)
    planning = commands.add_parser("plan")
    planning.add_argument("--model",required=True)
    planning.add_argument("--experiment",required=True,choices=("observations","historical","amip"))
    planning.add_argument("--member")
    planning.add_argument("--start-year",type=int,default=1990)
    planning.add_argument("--end-year",type=int,default=2010)
    planning.add_argument("--n-calibration",type=int,default=1000)
    planning.add_argument("--n-validation",type=int,default=1000)
    planning.add_argument("--chunk-size",type=int,default=100)
    planning.add_argument("--seed",type=int,default=20261031)
    planning.add_argument("--detection-sigma",type=float,choices=(1.0,1.5,2.0),default=1.5,
                          help="Envelope activity cutoff in signed-field standard deviations")
    planning.add_argument("--background-method",choices=("gaussian","wk99_121"),default="gaussian")
    for name in ("prepare","simulate","summarize","inject","finalize"):
        command = commands.add_parser(name)
        command.add_argument("--manifest",required=True)
        if name=="simulate": command.add_argument("--task",type=int,required=True)
        if name=="inject": command.add_argument("--symmetry",choices=("symmetric","antisymmetric"),required=True)
    args = parser.parse_args()
    if PROJECT_ROOT!=Path(__file__).resolve().parents[1]:
        raise RuntimeError("Install this checkout with pip -e .")
    stages = dict(plan=plan,prepare=prepare,simulate=simulate,summarize=summarize,inject=inject,finalize=finalize)
    if args.stage=="prepare" and json.loads(Path(args.manifest).read_text())["configuration"].get("background_method")=="wk99_121":
        stages["prepare"] = prepare_wk99
    stages[args.stage](args)


if __name__=="__main__":
    main()
