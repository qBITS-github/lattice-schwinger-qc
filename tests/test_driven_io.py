"""Check output isolation and provenance without generating physics datasets.

All test writes go to temporary directories. The reproduction and the real
driven-output directories are never used for test results.
"""

import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from driven import io as driven_io


class DrivenIOTests(unittest.TestCase):
    def test_actual_source_provenance(self):
        source = driven_io._source_metadata()
        self.assertEqual(len(source["commit"]), 40)
        self.assertIn("src/driven/io.py", source["file_sha256"])
        self.assertIn("configs/driven/case_a_pilot_n8.json", source["file_sha256"])

    def test_paths_cannot_escape_new_output_roots(self):
        for run_id in ("", "../ED", "/tmp/result", "a/b", "a\\b", ".", ".."):
            with self.subTest(run_id=run_id), self.assertRaises(ValueError):
                driven_io.run_path(run_id)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "data").mkdir()
            (root / "baseline").mkdir()
            (root / "data/driven").symlink_to(root / "baseline", target_is_directory=True)
            with patch.object(driven_io, "ROOT", root), self.assertRaises(ValueError):
                driven_io.run_path("test_002")

    def test_scan_resume_cache_and_corruption_checks(self):
        with tempfile.TemporaryDirectory() as temporary:
            source={"commit":"fixture","dirty":True,"file_sha256":{}}
            with patch.object(driven_io,"ROOT",Path(temporary)),patch.object(
                driven_io,"_source_metadata",return_value=source):
                config={"model":{"N":8}}
                store=driven_io.StudyStore("scan",config,"B")
                metadata=json.loads((store.data/"metadata.json").read_text())
                self.assertEqual(metadata["signature"]["config"],config)
                self.assertIn("numpy",metadata["packages"])
                self.assertTrue((store.data/"source.tar.gz").is_file())
                self.assertFalse((Path(temporary)/"figures").exists())
                with self.assertRaises(FileExistsError):
                    driven_io.StudyStore("scan",config,"B")
                specification={"solver":{"rtol":1e-10},"pulse":{"amplitude":.01}}
                def compute(): return {"values":np.arange(3)}, {"Y":.1}
                arrays,info,identifier=store.point(specification,compute)
                resumed=driven_io.StudyStore("scan",config,"B",resume=True)
                def must_not_run(): raise AssertionError("completed point was recomputed")
                reused,info,_=resumed.point(specification,must_not_run)
                np.testing.assert_array_equal(arrays["values"],reused["values"])
                self.assertEqual(info,{"Y":.1})
                with self.assertRaises(ValueError):
                    driven_io.StudyStore("scan",{"model":{"N":12}},"B",resume=True)
                other=driven_io.StudyStore("other",config,"B")
                other.point(specification,must_not_run)
                numeric=store.data/"points"/(identifier+".npz")
                numeric.write_bytes(b"corrupt")
                with self.assertRaises(ValueError): resumed.point(specification,must_not_run)


if __name__ == "__main__":
    unittest.main()
