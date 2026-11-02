"""Monte Carlo rejection semantics and an archive-free calibration integration."""

import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import xarray as xr

from src import wave_events as we, wave_significance as sg


class RecordMaxTests(unittest.TestCase):
    def test_strict_threshold_equals_monte_carlo_test_including_ties(self):
        for maxima in (np.arange(1000,dtype=float),np.repeat(np.arange(100,dtype=float),10)):
            scores = np.unique(np.r_[maxima,maxima+.01,maxima.max()+1])
            for alpha in (.01,.05,.1):
                cut = sg.record_max_threshold(maxima,alpha)
                self.assertTrue(cut["available"])
                np.testing.assert_array_equal(scores>cut["threshold_W_m2"],sg.record_max_pvalues(maxima,scores)<=alpha)
        self.assertFalse(sg.record_max_threshold(np.arange(98),.01)["available"])
        cut = sg.record_max_threshold(np.arange(99),.01)
        self.assertEqual(cut["threshold_W_m2"],98)
        self.assertTrue(cut["interval_upper_unbounded"])
        self.assertEqual(float(sg.record_max_pvalues(np.arange(99),[99])[0]),.01)
        with self.assertRaises(ValueError): sg.record_max_pvalues([], [1])

    def test_order_statistic_interval_contains_cut(self):
        cut = sg.record_max_threshold(np.arange(1000),.01)
        self.assertLess(cut["quantile_CI_low_W_m2"],cut["threshold_W_m2"])
        self.assertGreater(cut["quantile_CI_high_W_m2"],cut["threshold_W_m2"])


class CalibrationIntegrationTests(unittest.TestCase):
    background_method = "gaussian"

    def test_chunking_reduction_flags_validation_and_injections(self):
        script = Path(__file__).resolve().parents[1]/"scripts/gaussian_event_calibration.py"
        spec = importlib.util.spec_from_file_location("gaussian_calibration_test",script)
        workflow = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(workflow)
        with tempfile.TemporaryDirectory() as temporary:
            data = Path(temporary)/"data"
            root = data/"fixture"
            root.mkdir(parents=True)
            times = xr.date_range("2000-01-01",periods=730,freq="D")
            field = xr.DataArray(np.random.default_rng(12).normal(size=(730,48)),dims=("time","lon"),
                coords=dict(time=times,lon=np.arange(48)*7.5+.5),attrs=dict(units="W m-2"))
            wk99 = (sg.estimate_wk99_background_power(xr.Dataset(dict(symmetric=field,antisymmetric=field)))
                    if self.background_method=="wk99_121" else None)
            for symmetry in ("symmetric","antisymmetric"):
                field.rename(symmetry).to_netcdf(root/f"template.{symmetry}.nc")
                power = wk99.background_power if wk99 is not None else sg.estimate_background_power(field)
                xr.Dataset(dict(background_power=power)).to_netcdf(root/f"spectrum.{symmetry}.nc")
            tasks=[]
            for wave in workflow.PROFILES:
                filtered = sg.filter_projected_olr(field,wave,profile=workflow.PROFILES[wave])
                events,sigma = we.detect_events(filtered)
                observed = pd.DataFrame(events,columns=workflow.EVENT_COLUMNS)
                observed["reference_sigma_W_m2"] = sigma
                observed.to_csv(root/f"observed_candidates.{wave}.csv",index=False)
                for start,stop in ((0,100),(100,199)):
                    tasks.append(dict(task_id=len(tasks),wave=wave,start=start,stop=stop))
            cfg = dict(n_calibration=99,n_validation=100,seed=6,alpha=.01,diagnostic_alphas=[.10,.05,.01],
                bootstrap_replicates=20,injection_realizations=1,testing_scope="fixture",
                detection=dict(detection_sigma=1.5,peak_sigma=1.5),
                profiles={wave:profile.as_dict() for wave,profile in workflow.PROFILES.items()},
                identity=dict(dataset="test",source_id="fixture",experiment_id="observations",member_id="none",calendar="proleptic_gregorian"),
                injection=dict(depths_m=[8.,25.,90.],amplitudes_sigma=[2.,4.,8.],
                    wavenumbers={"Kelvin":[4,10],"n1ERW":[-3,-8],"MRG":[-3,-8],"n0EIG":[3,10],"n1WIG":[-3,-8]},
                    widths=dict(compact=[.5,.2],extended=[1.,.5])),
                simulation_sha256=workflow.simulation_hashes(self.background_method),
                versions=dict(python=platform.python_version(),numpy=workflow.np.__version__,scipy=workflow.scipy.__version__,
                    pandas=pd.__version__,xarray=xr.__version__,xr_ccew=workflow.version("xr-ccew"),cftime=workflow.version("cftime")))
            if self.background_method=="wk99_121": cfg["background_method"] = self.background_method
            path=root/"manifest.json"
            path.write_text(json.dumps(dict(status="prepared",ensemble_id="fixture",configuration=cfg,tasks=tasks,use_cftime=False)))
            with patch.object(workflow,"DATA_DIR",data):
                for task in tasks:
                    workflow.simulate(SimpleNamespace(manifest=str(path),task=task["task_id"]))
                # Completed chunks can be reused without changing a score file.
                sample=root/"chunks/task_00000.events.csv.gz"
                before=hashlib.sha256(sample.read_bytes()).hexdigest()
                workflow.simulate(SimpleNamespace(manifest=str(path),task=0))
                self.assertEqual(before,hashlib.sha256(sample.read_bytes()).hexdigest())
                workflow.summarize(SimpleNamespace(manifest=str(path)))
                for symmetry in ("symmetric","antisymmetric"):
                    workflow.inject(SimpleNamespace(manifest=str(path),symmetry=symmetry))
                workflow.finalize(SimpleNamespace(manifest=str(path)))
            report=json.loads((root/"report.json").read_text())
            self.assertEqual(report["status"],"complete")
            self.assertEqual(report["injection_trials"],36*len(workflow.PROFILES))
            thresholds=pd.read_csv(root/"thresholds.csv")
            self.assertEqual(len(thresholds),3*len(workflow.PROFILES))
            for wave in workflow.PROFILES:
                observed=pd.read_csv(root/f"observed_significance.{wave}.csv")
                threshold=thresholds.loc[(thresholds.wave_type==wave)&np.isclose(thresholds.alpha,.01)].threshold_W_m2.iloc[0]
                np.testing.assert_array_equal(observed.significant,observed.peak_amplitude_W_m2>threshold)
                self.assertTrue((observed.record_max_p_mc>=.01).all())
                with xr.open_dataset(root/f"null_threshold_counts.{wave}.nc") as counts:
                    self.assertEqual(counts.sizes["realization"],199)
                with xr.open_dataset(root/f"null_histograms.{wave}.nc") as hist:
                    self.assertEqual(hist.signed_edges.ndim,1)


class WK99CalibrationIntegrationTests(CalibrationIntegrationTests):
    background_method = "wk99_121"


if __name__=="__main__":
    unittest.main()
