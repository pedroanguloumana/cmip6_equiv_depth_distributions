"""Shared NOAA/CMIP OLR preprocessing, packet detection and event-catalog workflow.

Scientific settings are explicit below and identical for both datasets. Workers
supply the dataset identity and loading function; this module owns the common
calculation, quality screens, atomic CSV checkpoints and resource measurements.
"""

import csv
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import time

import numpy as np
import scipy
from scipy import ndimage
from scipy.signal import hilbert
import xarray as xr
import xr_ccew as tw

from .paths import DATA_DIR, PROJECT_ROOT
from . import wave_fitting
from .wave_fitting import fit_event_sinusoid, resting_wave_depth
from .wave_profiles import BOX_PROFILES


WAVE_PROFILES = BOX_PROFILES

FIT_OPTIONS = {
    "Kelvin": dict(
        k_grid=np.arange(1.0, 15.5, 0.5),
        frequency_grid=np.linspace(1/14, 1/2.5, 17),
        half_window_days=10, half_width_deg=90, time_sigma_days=5, lon_sigma_deg=45,
    ),
    "n1ERW": dict(
        k_grid=np.arange(-15.0, -0.5, 0.5),
        frequency_grid=np.linspace(1/90, 1/10, 17),
        half_window_days=90, half_width_deg=90, time_sigma_days=45, lon_sigma_deg=45,
    ),
    "n0EIG": dict(
        k_grid=np.arange(1.0, 14.5, 0.5),
        frequency_grid=np.linspace(1/6, 1/2.25, 17),
        half_window_days=10, half_width_deg=90, time_sigma_days=5, lon_sigma_deg=45,
    ),
    "MRG": dict(
        k_grid=np.arange(-10.0, -0.5, 0.5),
        frequency_grid=np.linspace(1/10, 1/2.25, 17),
        half_window_days=10, half_width_deg=90, time_sigma_days=5, lon_sigma_deg=45,
    ),
    "n1WIG": dict(
        k_grid=np.arange(-15.0, -0.5, 0.5),
        frequency_grid=np.linspace(1/4, 1/2, 17),
        half_window_days=10, half_width_deg=90, time_sigma_days=5, lon_sigma_deg=45,
    ),
}

LAT_RANGE = (-20, 20)
AVERAGE_LAT_RANGE = (0, 20)
SEASONAL_HARMONICS = 3
TAPER_FRACTION = 0.1
REFERENCE_EXCLUDE_DAYS = 150  # Preserve the NOAA notebook's sigma reference selection.
DETECTION_SIGMA = 1.5
PEAK_SIGMA = 1.5
MIN_DURATION_DAYS = 2
MIN_LONGITUDE_SPAN_DEG = 10
MIN_WEIGHTED_R2 = 0.3
MIN_SPATIAL_CYCLES = 0.5  # Global |K|=1 spans half a cycle in the 180-degree fit window.
MIN_TEMPORAL_CYCLES = 1.0

CSV_COLUMNS = [
    "dataset", "source_id", "experiment_id", "member_id", "source_grid_label", "analysis_grid", "calendar",
    "wave_type", "event_id", "peak_time", "peak_lon_deg", "duration_days", "span_deg",
    "peak_amplitude_W_m2", "K_global", "frequency_cpd", "k_rad_m", "omega_rad_s",
    "phase_speed_m_s", "c0_m_s", "h_eff_m", "period_days", "amplitude_W_m2",
    "offset_W_m2", "phase_rad", "weighted_r2", "spatial_cycles", "temporal_cycles",
    "fit_success", "at_search_edge", "dispersion_valid", "fit_quality_ok",
    "keep_for_histogram", "fit_status", "fit_error",
    "fit_elapsed_seconds",
]


def filter_olr_for_wave(olr, wave_type):
    """Return a (time, lon) wave field from daily OLR in W m-2 on (time, lat, lon).

    The input must span matching north/south latitudes within LAT_RANGE and all
    longitudes. Seasonal harmonics use the input calendar; elapsed daily samples
    drive the FFT. No calendar interpolation or rescaling of model days is applied.
    """
    if set(olr.dims) != {"time", "lat", "lon"}:
        raise ValueError("Expected OLR dimensions (time, lat, lon).")
    olr = olr.transpose("time", "lat", "lon")
    if olr.sizes["time"] <= 2 * REFERENCE_EXCLUDE_DAYS or not olr.sizes["lat"]:
        raise ValueError("The selected record must contain more than 300 daily samples and paired latitudes.")
    if olr.attrs.get("units") != "W m-2" or not np.isfinite(olr.values).all():
        raise ValueError("Expected finite daily OLR in W m-2.")
    times = olr.time.values
    if np.issubdtype(times.dtype, np.datetime64):
        intervals = np.diff(times) / np.timedelta64(1, "D")
    else:
        intervals = np.array([(b-a).total_seconds()/86400 for a, b in zip(times[:-1], times[1:])])
    if not np.allclose(intervals, 1, rtol=0, atol=1e-9):
        raise ValueError("Time must be regularly sampled daily in its native calendar.")
    dlon = np.diff(olr.lon.values)
    if not (len(dlon) and np.allclose(dlon, dlon[0]) and dlon[0] > 0
            and np.isclose(olr.sizes["lon"] * dlon[0], 360)):
        raise ValueError("Filtering requires an increasing, regular, global longitude grid.")
    if (not np.all(np.diff(olr.lat.values) > 0)
            or not np.allclose(olr.lat.values, -olr.lat.values[::-1])
            or olr.lat.values[0] < LAT_RANGE[0] or olr.lat.values[-1] > LAT_RANGE[1]):
        raise ValueError("Latitudes must be increasing north/south pairs within LAT_RANGE.")

    anomaly = tw.apply_window(
        tw.remove_mean_and_linear_trend(
            tw.remove_harmonics_of_seasonal_cycle(olr, num_harmonics=SEASONAL_HARMONICS)
        ), window="tukey", pct=TAPER_FRACTION,
    )
    del olr
    # Project symmetry before averaging one hemisphere: a full equatorial mean
    # would cancel antisymmetric MRG/EIG signals. Preserve the notebook's order.
    filtered = tw.filter_field(anomaly, WAVE_PROFILES[wave_type])
    del anomaly
    return filtered.sel(lat=slice(*AVERAGE_LAT_RANGE)).mean("lat").transpose("time", "lon")


def detect_events(field, *, detection_sigma=DETECTION_SIGMA, return_support=False):
    """Return envelope components above detection_sigma times the field sigma.

    The same value seeds active pixels and supplies the candidate peak minimum.
    Optionally map event IDs to their active time/longitude pixel indices.
    """
    if not np.isfinite(detection_sigma) or detection_sigma <= 0:
        raise ValueError("detection_sigma must be finite and positive.")
    sigma = float(field.isel(time=slice(REFERENCE_EXCLUDE_DAYS, -REFERENCE_EXCLUDE_DAYS)).std())
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("The filtered field has no finite positive reference standard deviation.")
    envelope = np.abs(hilbert(field.values, axis=1))
    labels, _ = ndimage.label(envelope > detection_sigma*sigma, structure=np.ones((3, 3), dtype=int))
    # Preserve the notebook's same-day and adjacent-day periodic connections.
    for row in range(labels.shape[0]):
        for neighbor in range(max(0, row-1), min(labels.shape[0], row+2)):
            left, right = labels[row, 0], labels[neighbor, -1]
            if left and right and left != right:
                labels[labels == max(left, right)] = min(left, right)

    events = []
    support = {} if return_support else None
    # Bounding boxes avoid rescanning the entire 21-year field for every label.
    # The candidate definition and label IDs remain those of the notebook.
    for component, box in enumerate(ndimage.find_objects(labels), start=1):
        if box is None:
            continue
        local_t, local_x = np.where(labels[box] == component)
        ti, xi = local_t + box[0].start, local_x + box[1].start
        peak = np.argmax(envelope[ti, xi])
        amplitude = float(envelope[ti[peak], xi[peak]])
        duration = int(ti.max() - ti.min() + 1)
        longitudes = np.unique(field.lon.values[xi])
        span = float(360 - np.diff(np.r_[longitudes, longitudes[0]+360]).max())
        if (ti.min() == 0 or ti.max() == field.sizes["time"]-1
                or duration < MIN_DURATION_DAYS or span < MIN_LONGITUDE_SPAN_DEG
                or amplitude <= detection_sigma*sigma):
            continue
        events.append(dict(
            event_id=int(component), peak_time=field.time.values[ti[peak]],
            peak_lon_deg=float(field.lon.values[xi[peak]]), duration_days=duration,
            span_deg=span, peak_amplitude_W_m2=amplitude,
        ))
        if return_support:
            support[int(component)] = (ti, xi)
    events.sort(key=lambda event: event["peak_time"])
    return (events, sigma, support) if return_support else (events, sigma)


def make_catalog_configuration(dataset_id, wave_type, start, end, *, identity,
                               input_paths, code_paths, max_events=None,
                               detection_sigma=DETECTION_SIGMA, fit_support="window"):
    """Record exact scientific settings, source identity and code/input versions.

    dataset_id is the readable filename prefix; identity supplies the row-level
    dataset/model/member/grid/calendar labels. Paths must enter through data/.
    """
    if set(identity) != {"dataset", "source_id", "experiment_id", "member_id",
                         "source_grid_label", "analysis_grid", "calendar"}:
        raise ValueError("Supply the complete dataset/model/member/grid/calendar identity.")
    if not np.isfinite(detection_sigma) or detection_sigma <= 0:
        raise ValueError("detection_sigma must be finite and positive.")
    if fit_support not in ("window", "component"):
        raise ValueError("fit_support must be 'window' or 'component'.")
    fit_keys = ("k_grid", "frequency_grid") if fit_support == "component" else FIT_OPTIONS[wave_type]
    options = {key: value.tolist() if isinstance(value, np.ndarray) else value
               for key, value in FIT_OPTIONS[wave_type].items() if key in fit_keys}
    inputs = []
    for path in input_paths:
        path = Path(path)
        stat = path.stat()
        inputs.append(dict(file=str(path.relative_to(PROJECT_ROOT)),
                           size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns))
    code_paths = dict(code_paths, wave_events=Path(__file__), wave_fitting=Path(wave_fitting.__file__))
    return {
        "schema_version": 2, "dataset_id": dataset_id, "identity": identity,
        "wave_type": wave_type, "start": start, "end": end,
        "fit_support": fit_support,
        "max_events": max_events, "profile": WAVE_PROFILES[wave_type].as_dict(),
        "fit_options": options,
        "preprocessing": dict(lat_range=LAT_RANGE, average_lat_range=AVERAGE_LAT_RANGE,
                              latitude_weights="unweighted", seasonal_harmonics=SEASONAL_HARMONICS,
                              remove_mean_and_trend=True, taper_fraction=TAPER_FRACTION,
                              reference_exclude_days=REFERENCE_EXCLUDE_DAYS,
                              calendar_handling="native; one model day = 86400 s"),
        "detection": dict(detection_sigma=float(detection_sigma), peak_sigma=float(detection_sigma),
                          min_duration_days=MIN_DURATION_DAYS, min_longitude_span_deg=MIN_LONGITUDE_SPAN_DEG,
                          connectivity=8, periodic_longitude=True, envelope_axis="lon"),
        "fit_screens": dict(min_weighted_r2=MIN_WEIGHTED_R2, min_spatial_cycles=MIN_SPATIAL_CYCLES,
                            min_temporal_cycles=MIN_TEMPORAL_CYCLES),
        "background_u_m_s": 0.0, "input": inputs,
        "code_sha256": {name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                        for name, path in code_paths.items()},
        "versions": dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                         xarray=xr.__version__, xr_ccew=version("xr-ccew"), cftime=version("cftime")),
    }


def catalog_path(configuration, output_dir=DATA_DIR):
    """Use readable fitting choices plus a hash covering every recorded setting."""
    options = configuration["fit_options"]
    number = lambda value: f"{value:g}".replace(".", "p")
    identity = json.dumps(configuration, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(identity.encode()).hexdigest()[:12]
    k, f = options["k_grid"], options["frequency_grid"]
    name = (
        f"{configuration['dataset_id']}.{configuration['wave_type']}."
        f"{configuration['start'].replace('-', '')}-{configuration['end'].replace('-', '')}."
        f"k{number(k[0])}to{number(k[-1])}_dk{number(k[1]-k[0])}."
        f"f{number(f[0])}to{number(f[-1])}_nf{len(f)}."
        f"ds{number(configuration['detection']['detection_sigma'])}."
        f"{'component' if configuration['fit_support'] == 'component' else 'window'}."
    )
    if configuration["fit_support"] == "window":
        name += (f"tw{number(options['half_window_days'])}_lw{number(options['half_width_deg'])}."
                 f"ts{number(options['time_sigma_days'])}_ls{number(options['lon_sigma_deg'])}.")
    name += digest
    if configuration["max_events"] is not None:
        name += f".pilot{configuration['max_events']}"
    return Path(output_dir) / (name + ".csv")


def write_json(path, contents):
    """Replace the small provenance/status document atomically."""
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(contents, indent=2) + "\n")
    temporary.replace(path)


def run_event_catalog(configuration, load_field, *, output_dir=DATA_DIR, overwrite=False):
    """Run one source/wave catalog, checkpointing each completed row atomically.

    load_field is a zero-argument loader returning the filtered (time, lon) field.
    It is called only after checking for completed outputs. Loading/preprocessing
    remains inside the timing and error-reporting boundary, without holding fields
    from multiple models in memory. CSV peak dates retain their native calendar.
    """
    wave_type = configuration["wave_type"]
    output = catalog_path(configuration, output_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    parts = output.with_name(output.name + ".parts")
    metadata_path = output.with_suffix(".json")
    if output.exists() and not overwrite:
        print(f"SKIP (complete): {output}", flush=True)
        return
    if metadata_path.exists():
        previous = json.loads(metadata_path.read_text())
        # JSON normalizes tuples to lists, so compare serialized configurations.
        if json.dumps(previous["configuration"], sort_keys=True) != json.dumps(configuration, sort_keys=True):
            raise RuntimeError("Existing provenance does not match; refusing to reuse its checkpoints.")
    if overwrite and parts.exists():
        shutil.rmtree(parts)  # Only this exact run's generated checkpoint directory.
    parts.mkdir(parents=True, exist_ok=True)
    metadata = dict(configuration=configuration, status="running", csv=str(output),
                    started_utc=datetime.now(timezone.utc).isoformat(),
                    slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                    slurm_array_task_id=os.environ.get("SLURM_ARRAY_TASK_ID"),
                    slurm_partition=os.environ.get("SLURM_JOB_PARTITION"),
                    slurm_cpus_per_task=os.environ.get("SLURM_CPUS_PER_TASK"))
    write_json(metadata_path, metadata)
    started = time.monotonic()
    cpu_started = time.process_time()
    print(f"{configuration['dataset_id']}.{wave_type}: {configuration['start']} through {configuration['end']}\nOutput: {output}", flush=True)

    try:
        print("Loading and filtering the full selected record...", flush=True)
        field = load_field()
        preprocessing_seconds = time.monotonic() - started
        print(f"Load/filter: {preprocessing_seconds:.1f} s", flush=True)
        detection_started = time.monotonic()
        component_fit = configuration.get("fit_support", "window") == "component"
        detected = detect_events(field,
            detection_sigma=configuration["detection"]["detection_sigma"],
            return_support=component_fit)
        if component_fit:
            events, sigma, support = detected
        else:
            events, sigma = detected
        detection_seconds = time.monotonic() - detection_started
        options = configuration["fit_options"]
        half_window = int(options["half_window_days"]) if not component_fit else 0
        time_index = field.get_index("time")
        # Count all eligible windows before applying a pilot limit. This gives
        # the correct multiplier when extrapolating a pilot to the full catalog.
        full_record_fit_candidates = sum(
            half_window <= time_index.get_loc(event["peak_time"]) < field.sizes["time"]-half_window
            for event in events
        )
        metadata.update(reference_sigma_W_m2=sigma, detected_candidates=len(events),
                        full_record_fit_candidates=int(full_record_fit_candidates),
                        timings_seconds=dict(load_filter=preprocessing_seconds, detection=detection_seconds))
        if configuration["max_events"] is not None:
            events = events[:configuration["max_events"]]
        metadata["catalog_candidates"] = len(events)
        write_json(metadata_path, metadata)
        print(f"{len(events)} candidates to fit; reference sigma={sigma:g} W m-2", flush=True)

        resumed = 0
        for index, event in enumerate(events, start=1):
            part = parts / f"event_{event['event_id']:08d}.csv"
            if part.exists():
                resumed += 1
                continue
            row = dict(configuration["identity"], wave_type=wave_type, **event, fit_success=False, at_search_edge=False,
                       dispersion_valid=False, fit_quality_ok=False, keep_for_histogram=False,
                       fit_status="incomplete_window", fit_error="", fit_elapsed_seconds=0.0)
            peak = event["peak_time"]
            full_window = half_window <= time_index.get_loc(peak) < field.sizes["time"]-half_window
            if full_window:
                fit_started = time.monotonic()
                try:
                    fit_options = dict(options)
                    if component_fit:
                        fit_options["support_indices"] = support[event["event_id"]]
                    fit, diagnostic = fit_event_sinusoid(
                        field, peak, event["peak_lon_deg"], **fit_options,
                    )
                    del diagnostic  # Do not accumulate patches/models for thousands of events.
                    fit.update(resting_wave_depth(wave_type, fit["k_rad_m"], fit["omega_rad_s"]))
                    quality = bool(fit["fit_success"] and not fit["at_search_edge"]
                                   and fit["weighted_r2"] >= MIN_WEIGHTED_R2
                                   and fit["spatial_cycles"] >= MIN_SPATIAL_CYCLES
                                   and fit["temporal_cycles"] >= MIN_TEMPORAL_CYCLES)
                    row.update(fit, fit_quality_ok=quality,
                               keep_for_histogram=quality and fit["dispersion_valid"],
                               fit_status="ok" if fit["fit_success"] else "optimizer_failed")
                except (ValueError, np.linalg.LinAlgError) as error:
                    row.update(fit_status="fit_error", fit_error=f"{type(error).__name__}: {error}")
                    print(f"Event {event['event_id']}: {row['fit_error']}", flush=True)
                finally:
                    row["fit_elapsed_seconds"] = time.monotonic() - fit_started

            # A part exists only after its single row has been completely written.
            temporary = part.with_name(part.name + ".tmp")
            with temporary.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, restval="nan")
                writer.writeheader()
                writer.writerow(row)
                handle.flush()
                os.fsync(handle.fileno())
            temporary.replace(part)
            if index % 25 == 0 or index == len(events):
                print(f"[{index}/{len(events)}] {time.monotonic()-started:.1f} s elapsed", flush=True)

        # Stream completed rows into one final CSV; existence of that CSV marks completion.
        temporary_output = output.with_name(output.name + ".tmp")
        counts = dict(fit_success=0, fit_errors=0, incomplete_windows=0, retained=0)
        fit_seconds = []
        with temporary_output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            for event in events:
                part = parts / f"event_{event['event_id']:08d}.csv"
                with part.open(newline="") as source:
                    reader = csv.DictReader(source)
                    if reader.fieldnames != CSV_COLUMNS:
                        raise RuntimeError(f"Unexpected checkpoint schema: {part}")
                    row = next(reader)
                    if int(row["event_id"]) != event["event_id"] or next(reader, None) is not None:
                        raise RuntimeError(f"Unexpected checkpoint contents: {part}")
                writer.writerow(row)
                counts["fit_success"] += row["fit_success"] == "True"
                counts["fit_errors"] += row["fit_status"] == "fit_error"
                counts["incomplete_windows"] += row["fit_status"] == "incomplete_window"
                counts["retained"] += row["keep_for_histogram"] == "True"
                if row["fit_status"] in ("ok", "optimizer_failed"):
                    fit_seconds.append(float(row["fit_elapsed_seconds"]))
            handle.flush()
            os.fsync(handle.fileno())
        temporary_output.replace(output)
        metadata.update(status="complete", completed_utc=datetime.now(timezone.utc).isoformat(),
                        elapsed_seconds=time.monotonic()-started, resumed_events=resumed, counts=counts,
                        cpu_seconds=time.process_time()-cpu_started)
        # getrusage records the process high-water mark, including short FFT
        # peaks that periodic scheduler accounting can miss. Linux reports KiB;
        # macOS reports bytes. This measures the worker, not filesystem cache.
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        metadata["peak_rss_gib"] = float(rss / (1024**3 if platform.system() == "Darwin" else 1024**2))
        if fit_seconds:
            metadata["fit_seconds_summary"] = dict(
                count=len(fit_seconds), mean=float(np.mean(fit_seconds)),
                median=float(np.median(fit_seconds)), p90=float(np.percentile(fit_seconds, 90)),
                maximum=float(np.max(fit_seconds)), total=float(np.sum(fit_seconds)),
            )
        write_json(metadata_path, metadata)
        print(f"OK {output}: {len(events)} rows; {counts}; resumed={resumed}", flush=True)
        print(f"Peak worker RSS: {metadata['peak_rss_gib']:.2f} GiB; "
              f"CPU time: {metadata['cpu_seconds']:.1f} s; "
              f"fit seconds: {metadata.get('fit_seconds_summary')}", flush=True)
    except BaseException as error:
        metadata.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                        error=f"{type(error).__name__}: {error}")
        write_json(metadata_path, metadata)
        raise


def load_event_catalog(model, experiment, wave_type, *, member_id=None,
                       year_range=(1990, 2010), accepted_only=False, run_id=None):
    """Return a completed NOAA/CMIP6 event catalog as a pandas DataFrame.

    Parameters
    ----------
    model, experiment, wave_type : str
        Use model="NOAA_OLR" (alias "NOAA") and experiment="observations"
        for NOAA. For CMIP6 use the source_id and "historical" or "amip".
        Wave labels are Kelvin, n1ERW, MRG, n0EIG and n1WIG.
    member_id : str, optional
        Select a particular CMIP member; otherwise the match must be unique.
    year_range : pair of int or None
        Inclusive complete record years, not a filter on individual peak dates.
        None leaves the record period unconstrained. Native calendar year ends
        are respected, including December 30 for 360_day.
    accepted_only : bool
        Default False returns all candidates, including unsuccessful/poor fits.
        True selects the existing keep_for_histogram quality flag.
    run_id : str, optional
        The final 12-character filename hash, used to resolve multiple runs.

    Returns
    -------
    pandas.DataFrame
        Event rows with common dataset/simulation/calendar labels. h_eff_m is
        in meters, peak_lon_deg is degrees east, and peak_time stays a string:
        Gregorian datetime parsing would lose valid dates such as February 30.
        Boolean quality columns use pandas' nullable boolean dtype. The actual
        CSV path and complete JSON sidecar are in attrs["catalog_path"] and
        attrs["catalog_metadata"]. No raw OLR files are required.

    Only active, completed full catalogs are searched, never archived files,
    pilots or checkpoint directories. Missing and ambiguous matches raise
    errors; an arbitrary newest run or a different ensemble is never selected.
    """
    import re

    import pandas as pd

    from .paths import CMIP6_EVENT_DIR

    if wave_type not in WAVE_PROFILES:
        raise ValueError(f"Unknown wave_type {wave_type!r}; choose {tuple(WAVE_PROFILES)}")
    noaa = model in ("NOAA_OLR", "NOAA")
    if noaa:
        if experiment != "observations":
            raise ValueError("NOAA catalogs use experiment='observations'.")
        directory, pattern = DATA_DIR, f"NOAA_OLR.{wave_type}.*.csv"
    else:
        if experiment not in ("historical", "amip"):
            raise ValueError("CMIP6 experiment must be 'historical' or 'amip'.")
        if not isinstance(model, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", model):
            raise ValueError("model must be an exact CMIP6 source_id or 'NOAA_OLR'.")
        directory = CMIP6_EVENT_DIR
        pattern = f"CMIP6.{model}.{experiment}.*.rlut.{wave_type}.*.csv"
    if year_range is not None:
        if (len(year_range) != 2 or any(not isinstance(year, int) for year in year_range)
                or not 1 <= year_range[0] <= year_range[1] <= 9999):
            raise ValueError("year_range must be an increasing pair of integer years or None.")
    if run_id is not None and not re.fullmatch(r"[0-9a-f]{12}", run_id):
        raise ValueError("run_id must be the final 12-character hexadecimal filename hash.")

    identity_columns = ("dataset", "source_id", "experiment_id", "member_id",
                        "source_grid_label", "analysis_grid", "calendar")
    matches = []
    for path in sorted(directory.glob(pattern)):
        if re.search(r"\.pilot\d+\.csv$", path.name):
            continue
        if run_id is not None and path.stem.rsplit(".", 1)[-1] != run_id:
            continue
        sidecar = path.with_suffix(".json")
        if not sidecar.is_file():
            raise FileNotFoundError(f"Missing catalog provenance: {sidecar}")
        metadata = json.loads(sidecar.read_text())
        configuration = metadata.get("configuration", {})
        if metadata.get("status") != "complete" or configuration.get("max_events") is not None:
            continue
        if configuration.get("wave_type") != wave_type:
            raise ValueError(f"Filename and metadata wave_type disagree: {path}")
        if noaa and configuration.get("schema_version") == 1:
            # The original NOAA writer saved NumPy Gregorian dates but no
            # row-level identity. Supply its known labels without opening OLR.
            identity = dict(dataset="NOAA_OLR", source_id="NOAA_CDR", experiment_id="observations",
                            member_id="not_applicable", source_grid_label="NOAA_1deg",
                            analysis_grid="NOAA_1deg", calendar="proleptic_gregorian")
        else:
            identity = configuration.get("identity", {})
        if not all(identity.get(key) for key in identity_columns):
            raise ValueError(f"Catalog metadata lacks the complete simulation/calendar identity: {sidecar}")
        expected = ("NOAA_OLR", "NOAA_CDR") if noaa else ("CMIP6", model)
        if ((identity["dataset"], identity["source_id"]) != expected
                or identity["experiment_id"] != experiment):
            raise ValueError(f"Filename and metadata simulation identity disagree: {path}")
        if not noaa and identity["member_id"] != path.name.split(".")[3]:
            raise ValueError(f"Filename and metadata ensemble member disagree: {path}")
        if member_id is not None and identity["member_id"] != member_id:
            continue
        if year_range is not None:
            last_day = 30 if identity["calendar"] == "360_day" else 31
            if (configuration.get("start") != f"{year_range[0]:04d}-01-01"
                    or configuration.get("end") != f"{year_range[1]:04d}-12-{last_day}"):
                continue
        matches.append((path, metadata, identity))
    if not matches:
        raise FileNotFoundError(
            f"No completed full event catalog for {model}/{experiment}/{wave_type}, "
            f"member={member_id}, years={year_range}, run_id={run_id} in {directory}. "
            "The requested catalog may still be running."
        )
    if len(matches) > 1:
        names = "\n".join(str(item[0]) for item in matches)
        raise ValueError(f"Multiple event catalogs match; specify member_id, year_range or run_id:\n{names}")
    path, metadata, identity = matches[0]
    flags = ("fit_success", "at_search_edge", "dispersion_valid", "fit_quality_ok", "keep_for_histogram")
    frame = pd.read_csv(path, dtype={"peak_time": "string", "event_id": "Int64", "h_eff_m": "float64",
                                    **{key: "boolean" for key in flags}})
    required = {"event_id", "peak_time", "peak_lon_deg", "h_eff_m", "keep_for_histogram"}
    if missing := required.difference(frame.columns):
        raise ValueError(f"Missing catalog columns {sorted(missing)}: {path}")
    if "catalog_candidates" in metadata and len(frame) != metadata["catalog_candidates"]:
        raise ValueError(f"CSV row count disagrees with completed catalog metadata: {path}")
    for key, value in {**identity, "wave_type": wave_type}.items():
        if key in frame and not frame[key].eq(value).fillna(False).all():
            raise ValueError(f"CSV {key} disagrees with catalog metadata: {path}")
        frame[key] = value
    # Put legacy NOAA and current CMIP columns in the same order, retaining any
    # additional producer columns rather than silently dropping information.
    columns = [key for key in CSV_COLUMNS if key in frame]
    columns.extend(key for key in frame if key not in CSV_COLUMNS)
    frame = frame.loc[:, columns]
    if accepted_only:
        frame = frame.loc[frame["keep_for_histogram"].fillna(False)].copy()
    frame.attrs.update(catalog_path=str(path), catalog_metadata=metadata)
    return frame


def load_noaa_event_catalog(wave_type, *, year_range=(1990, 2010), accepted_only=False, run_id=None):
    """Load NOAA events with the same columns/options as load_event_catalog."""
    return load_event_catalog("NOAA_OLR", "observations", wave_type, year_range=year_range,
                              accepted_only=accepted_only, run_id=run_id)


def load_cmip6_event_catalog(model, experiment, wave_type, *, member_id=None,
                            year_range=(1990, 2010), accepted_only=False, run_id=None):
    """Load one explicitly identified CMIP6 simulation/wave catalog."""
    if model in ("NOAA_OLR", "NOAA"):
        raise ValueError("Use load_noaa_event_catalog for NOAA observations.")
    return load_event_catalog(model, experiment, wave_type, member_id=member_id,
                              year_range=year_range, accepted_only=accepted_only, run_id=run_id)
