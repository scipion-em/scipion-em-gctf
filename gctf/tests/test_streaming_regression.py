# ***************************************************************************
# *
# * Regression tests for GCTF streaming/Continue behaviour.
# *
# ***************************************************************************

import os
import tempfile
import unittest
from unittest.mock import patch

from gctf.protocols.protocol_gctf import ProtGctf


class _Value:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _Mic:
    def __init__(self, obj_id, file_name):
        self._obj_id = obj_id
        self._file_name = file_name

    def getObjId(self):
        return self._obj_id

    def getFileName(self):
        return self._file_name


class _Program:
    def getCommand(self, **kwargs):
        return "gctf", "args"

    def getExt(self):
        return ".epa"


class _BatchHarness:
    def __init__(self, tmp_dir):
        self.tmp_dir = tmp_dir
        self.work_dir = os.path.join(tmp_dir, "work-0002")
        self.extra_dir = os.path.join(tmp_dir, "extra")
        os.makedirs(self.extra_dir, exist_ok=True)

        self.ctfDownFactor = _Value(1)
        self._params = {"scannedPixelSize": 1.0}
        self._gctfProgram = _Program()
        self.errors = []

    def _getMicrographDir(self, mic):
        return os.path.join(self.tmp_dir, "work")

    def _getPsdPath(self, mic_fn):
        base = os.path.splitext(os.path.basename(mic_fn))[0]
        return os.path.join(self.extra_dir, base + "_ctf.mrc")

    def _getCtfOutPath(self, mic_fn):
        base = os.path.splitext(os.path.basename(mic_fn))[0]
        return os.path.join(self.extra_dir, base + "_ctf.log")

    def _getCtfFitOutPath(self, mic_fn):
        base = os.path.splitext(os.path.basename(mic_fn))[0]
        return os.path.join(self.extra_dir, base + "_ctf_EPA.log")

    def runJob(self, *args, **kwargs):
        os.makedirs(self.work_dir, exist_ok=True)

        # mic_001 simulates a partial GCTF failure: its main CTF output is
        # missing. mic_002 has a complete, valid result and must not be lost.
        for base in ("mic_001", "mic_002"):
            with open(os.path.join(self.work_dir, base + "_gctf.log"), "w") as handle:
                handle.write("log\n")
            with open(os.path.join(self.work_dir, base + "_EPA.log"), "w") as handle:
                handle.write("fit\n")

        with open(os.path.join(self.work_dir, "mic_002.epa"), "w") as handle:
            handle.write("valid ctf output\n")

    def error(self, message):
        self.errors.append(message)


class TestGctfStreamingRegression(unittest.TestCase):
    def testBatchFailureDoesNotDiscardLaterValidOutputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            mic1_fn = os.path.join(tmp, "mic_001.mrc")
            mic2_fn = os.path.join(tmp, "mic_002.mrc")

            for mic_fn in (mic1_fn, mic2_fn):
                with open(mic_fn, "w") as handle:
                    handle.write("micrograph\n")

            protocol = _BatchHarness(tmp)
            mic1 = _Mic(1, mic1_fn)
            mic2 = _Mic(2, mic2_fn)

            with patch(
                "gctf.protocols.protocol_gctf.Plugin.getEnviron",
                return_value={},
            ):
                ProtGctf._estimateCtfList(protocol, [mic1, mic2])

            self.assertTrue(
                os.path.exists(protocol._getPsdPath(mic2_fn)),
                "A failed micrograph in a GCTF batch must not prevent a later "
                "micrograph with valid outputs from being preserved.",
            )
            self.assertTrue(os.path.exists(protocol._getCtfOutPath(mic2_fn)))
            self.assertTrue(os.path.exists(protocol._getCtfFitOutPath(mic2_fn)))


if __name__ == "__main__":
    unittest.main()
