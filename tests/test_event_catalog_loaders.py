"""Catalog selection/normalization checks, using only generated CSV/JSON files."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from src import paths, wave_events


class EventCatalogLoaderTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.data = Path(temporary.name) / "data"
        self.data.mkdir()
        self.cmip = self.data / "cmip6_event_catalogs"
        self.cmip.mkdir()
        self.noaa_patch = patch.object(wave_events, "DATA_DIR", self.data)
        self.cmip_patch = patch.object(paths, "CMIP6_EVENT_DIR", self.cmip)
        self.noaa_patch.start()
        self.cmip_patch.start()
        self.addCleanup(self.noaa_patch.stop)
        self.addCleanup(self.cmip_patch.stop)

    def write_catalog(self, *, noaa=False, legacy=False, experiment="historical", wave="Kelvin",
                      member="r1i1p1f1", run_id="aaaaaaaaaaaa", years=(1990, 2010),
                      calendar="proleptic_gregorian", status="complete", pilot=None, empty=False):
        identity = dict(dataset="NOAA_OLR" if noaa else "CMIP6",
                        source_id="NOAA_CDR" if noaa else "TestModel",
                        experiment_id="observations" if noaa else experiment,
                        member_id="not_applicable" if noaa else member,
                        source_grid_label="NOAA_1deg" if noaa else "gn",
                        analysis_grid="NOAA_1deg" if noaa else "noaa_1deg", calendar=calendar)
        rows = [dict(wave_type=wave, event_id=1, peak_time="2000-02-30 12:00:00" if calendar=="360_day" else "2000-02-29 12:00:00",
                     peak_lon_deg=359.5, h_eff_m=25.0, fit_success=True, at_search_edge=False,
                     dispersion_valid=True, fit_quality_ok=True, keep_for_histogram=True),
                dict(wave_type=wave, event_id=2, peak_time="2000-03-01 12:00:00", peak_lon_deg=0.5,
                     h_eff_m=80.0, fit_success=False, at_search_edge=False,
                     dispersion_valid=True, fit_quality_ok=False, keep_for_histogram=False)]
        frame = pd.DataFrame(rows)
        if not legacy:
            for key, value in identity.items():
                frame[key] = value
        if empty:
            frame = frame.iloc[:0]
        prefix = f"NOAA_OLR.{wave}" if noaa else f"CMIP6.TestModel.{experiment}.{member}.rlut.{wave}"
        last_day = 30 if calendar=="360_day" else 31
        name = f"{prefix}.{years[0]}0101-{years[1]}12{last_day}.options.{run_id}"
        if pilot is not None:
            name += f".pilot{pilot}"
        path = (self.data if noaa else self.cmip) / (name + ".csv")
        frame.to_csv(path, index=False)
        configuration = dict(schema_version=1 if legacy else 2, wave_type=wave,
                             start=f"{years[0]}-01-01", end=f"{years[1]}-12-{last_day}", max_events=pilot)
        if not legacy:
            configuration["identity"] = identity
        metadata = dict(status=status, configuration=configuration, catalog_candidates=len(frame))
        path.with_suffix(".json").write_text(json.dumps(metadata))
        return path

    def test_common_loop_normalizes_legacy_noaa_and_both_experiments(self):
        self.write_catalog(noaa=True, legacy=True)
        self.write_catalog()
        self.write_catalog(experiment="amip", calendar="360_day")
        frames = {}
        for model, experiment in (("NOAA_OLR", "observations"), ("TestModel", "historical"), ("TestModel", "amip")):
            frame = wave_events.load_event_catalog(model, experiment, "Kelvin")
            self.assertEqual(len(frame), 2)  # Poor fits remain present by default.
            self.assertEqual(frame.experiment_id.unique().tolist(), [experiment])
            self.assertEqual(str(frame.keep_for_histogram.dtype), "boolean")
            self.assertIsInstance(frame.peak_time.iloc[0], str)
            self.assertEqual(frame.attrs["catalog_metadata"]["status"], "complete")
            self.assertTrue(Path(frame.attrs["catalog_path"]).is_file())
            frames[experiment] = frame
        self.assertEqual(list(frames["observations"]), list(frames["historical"]))
        self.assertEqual(list(frames["observations"]), list(frames["amip"]))
        self.assertEqual(frames["amip"].peak_time.iloc[0], "2000-02-30 12:00:00")
        self.assertEqual(frames["observations"].calendar.iloc[0], "proleptic_gregorian")
        self.assertEqual(frames["observations"].source_id.iloc[0], "NOAA_CDR")

    def test_named_loaders_and_quality_selection(self):
        self.write_catalog(noaa=True, legacy=True)
        self.write_catalog(experiment="amip")
        noaa = wave_events.load_noaa_event_catalog("Kelvin", accepted_only=True)
        cmip = wave_events.load_cmip6_event_catalog("TestModel", "amip", "Kelvin", accepted_only=True)
        self.assertEqual(noaa.event_id.tolist(), [1])
        self.assertEqual(cmip.event_id.tolist(), [1])
        pd.testing.assert_frame_equal(noaa, wave_events.load_event_catalog("NOAA", "observations", "Kelvin", accepted_only=True))

    def test_pilots_running_records_other_periods_and_archives_are_excluded(self):
        expected = self.write_catalog()
        self.write_catalog(pilot=5, run_id="bbbbbbbbbbbb")
        self.write_catalog(status="running", run_id="cccccccccccc")
        self.write_catalog(years=(2000, 2010), run_id="dddddddddddd")
        archived = self.write_catalog(run_id="eeeeeeeeeeee")
        destination = self.data/"archived"
        destination.mkdir()
        archived.rename(destination/archived.name)
        archived.with_suffix(".json").rename(destination/archived.with_suffix(".json").name)
        (self.cmip/(expected.name+".parts")).mkdir()
        result = wave_events.load_event_catalog("TestModel", "historical", "Kelvin")
        self.assertEqual(Path(result.attrs["catalog_path"]), expected)

    def test_ambiguous_members_and_runs_require_explicit_selection(self):
        self.write_catalog()
        second = self.write_catalog(member="r2i1p1f1", run_id="bbbbbbbbbbbb")
        with self.assertRaisesRegex(ValueError, "Multiple event catalogs"):
            wave_events.load_event_catalog("TestModel", "historical", "Kelvin")
        result = wave_events.load_event_catalog("TestModel", "historical", "Kelvin", member_id="r2i1p1f1")
        self.assertEqual(Path(result.attrs["catalog_path"]), second)
        third = self.write_catalog(member="r2i1p1f1", run_id="cccccccccccc")
        with self.assertRaisesRegex(ValueError, "Multiple event catalogs"):
            wave_events.load_event_catalog("TestModel", "historical", "Kelvin", member_id="r2i1p1f1")
        result = wave_events.load_event_catalog("TestModel", "historical", "Kelvin", run_id="cccccccccccc")
        self.assertEqual(Path(result.attrs["catalog_path"]), third)

    def test_empty_catalog_is_valid_and_missing_catalog_is_explicit(self):
        self.write_catalog(empty=True)
        result = wave_events.load_event_catalog("TestModel", "historical", "Kelvin", accepted_only=True)
        self.assertTrue(result.empty)
        self.assertIn("calendar", result)
        with self.assertRaisesRegex(FileNotFoundError, "may still be running"):
            wave_events.load_event_catalog("TestModel", "historical", "MRG")

    def test_missing_sidecar_and_mismatched_row_identity_fail(self):
        path = self.write_catalog()
        sidecar = path.with_suffix(".json")
        contents = sidecar.read_text()
        sidecar.unlink()
        with self.assertRaisesRegex(FileNotFoundError, "Missing catalog provenance"):
            wave_events.load_event_catalog("TestModel", "historical", "Kelvin")
        sidecar.write_text(contents)
        frame = pd.read_csv(path)
        frame.loc[0, "member_id"] = "r9i1p1f1"
        frame.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, "CSV member_id"):
            wave_events.load_event_catalog("TestModel", "historical", "Kelvin")


if __name__ == "__main__":
    unittest.main()
