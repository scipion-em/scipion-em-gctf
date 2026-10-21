# **************************************************************************
# *
# * A micrograph whose CTF could not be estimated is never published (no
# * -999 placeholder): Continue estimates it again, and the run fails only
# * when no micrograph at all could be estimated.
# *
# **************************************************************************

import unittest

from pwem.objects import CTFModel

from gctf.protocols import ProtGctf

from .test_streaming_architecture import (
    _CtfStep,
    _GctfOutputCompletionHarness,
    _GctfResumeHarness,
    _Mic,
)


def _step(micName, finished=True):
    return _CtfStep("estimateCtfStep", '["%s"]' % micName, finished)


class _Value:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _Program:
    def __init__(self, ctf):
        self.ctf = ctf

    def parseOutputAsCtf(self, ctfFn, psdFile=None):
        return self.ctf


class _ModelHarness:
    def __init__(self, ctf):
        self._gctfProgram = _Program(ctf)
        self.ctfDownFactor = _Value(1)

    def _getPsdPath(self, mic):
        return "psd.mrc"

    def _getCtfOutPath(self, mic):
        return "gctf.log"


class _SampledMic(_Mic):
    def getSamplingRate(self):
        return 1.0

    def setSamplingRate(self, value):
        pass


class TestAFailedEstimationIsNeverPublished(unittest.TestCase):
    def test_AnUnusableResultIsRefusedInsteadOfPublishedAsMinus999(self):
        failed = CTFModel()
        failed.setWrongDefocus()

        with self.assertRaises(ValueError):
            ProtGctf._createCtfModel(_ModelHarness(failed),
                                     _SampledMic(1, "mic_001"))

    def test_AGoodResultIsPublished(self):
        good = CTFModel()
        good.setStandardDefocus(12000.0, 11000.0, 45.0)

        ctf = ProtGctf._createCtfModel(_ModelHarness(good),
                                       _SampledMic(1, "mic_001"))

        self.assertEqual(12000.0, ctf.getDefocusU())


class _UsabilityResumeHarness(_GctfResumeHarness):
    def __init__(self, usable):
        super().__init__()
        self._usable = usable

    def _ctfResultIsUsable(self, mic):
        return self._usable.get(mic.getMicName(), True)


class TestContinueEstimatesFailedMicrographsAgain(unittest.TestCase):
    def test_OnlyMicrographsWithAnUnusableResultAreEstimatedAgain(self):
        # mic_001 is published; mic_002 and mic_003 finished in the earlier
        # run, with an unusable and a usable result.
        protocol = _UsabilityResumeHarness({"mic_002": False})
        protocol._steps = []
        protocol._prevSteps = [_step("mic_002"), _step("mic_003")]

        ProtGctf._restoreProcessedMicsFromPersistentState(protocol)
        ProtGctf._checkNewInput(protocol)

        self.assertEqual(["mic_002"], protocol.scheduledMicNames)
        self.assertIn("mic_003", protocol.micDict)

    def test_ARetriedMicrographWaitsForItsNewStep(self):
        protocol = _GctfOutputCompletionHarness()
        protocol._prevSteps, protocol._steps = protocol._steps, []
        protocol._retriedMicNames = {"mic_001"}

        ProtGctf._checkNewOutput(protocol)
        self.assertEqual([], protocol.publishedBatches)

        protocol._steps = [_step("mic_001")]
        ProtGctf._checkNewOutput(protocol)
        self.assertEqual([["mic_001"]], protocol.publishedBatches)


class _ClosingHarness(_GctfOutputCompletionHarness):
    def __init__(self, knownIds, publishes=True):
        super().__init__()
        self.streamClosed = True
        self._knownMicIds = set(knownIds)
        self._publishes = publishes

    def _updateOutputCTFSet(self, micList, streamMode):
        if self._publishes:
            return super()._updateOutputCTFSet(micList, streamMode)
        return list(micList)

    def _updateStreamState(self, streamMode):
        pass

    def _failIfNoCtfWasEstimated(self):
        return ProtGctf._failIfNoCtfWasEstimated(self)


class TestTheRunFailsOnlyWhenEveryMicrographFailed(unittest.TestCase):
    def test_NothingPublishedAtTheEndFailsTheRun(self):
        with self.assertRaisesRegex(Exception, "No CTF could be estimated"):
            ProtGctf._checkNewOutput(_ClosingHarness({1}, publishes=False))

    def test_SomeFailuresDoNotFailTheRun(self):
        protocol = _ClosingHarness({1, 2})

        ProtGctf._checkNewOutput(protocol)

        self.assertTrue(protocol.finished)


if __name__ == '__main__':
    unittest.main()
