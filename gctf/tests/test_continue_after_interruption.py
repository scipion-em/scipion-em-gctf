# **************************************************************************
# *
# * Continue after an interrupted run: every micrograph ends published or
# * estimated again, and a stream that never completes is reported.
# *
# **************************************************************************

import unittest

from gctf.protocols import ProtGctf
from gctf.protocols.protocol_streaming_base import GctfStreamingBase

from .logical_set_fakes import LogicalSetFake
from .test_streaming_architecture import _CtfStep, _GctfResumeHarness, _Mic


class TestStepsOfAnEarlierRunThatDidNotFinish(unittest.TestCase):
    """pyworkflow only reruns the steps _insertAllSteps creates: the ones the
    generator inserted in an earlier run stay in _prevSteps and never run
    again, so an unfinished one there is not "scheduled" any more."""

    def test_OnlyFinishedStepsOfAnEarlierRunCountAsScheduled(self):
        protocol = GctfStreamingBase()
        protocol._prevSteps = [
            _CtfStep("estimateCtfStep", '["mic_001"]', False),
            _CtfStep("estimateCtfStep", '["mic_002"]', True)]
        protocol._steps = [_CtfStep("estimateCtfStep", '["mic_003"]', False)]

        scheduled = protocol._collectStepArgKeys(("estimateCtfStep",),
                                                 onlyFinished=False)

        # mic_003 belongs to this run and will still run.
        self.assertEqual({"mic_002", "mic_003"}, scheduled)

    def test_AMicrographWhoseStepWasInterruptedIsEstimatedAgain(self):
        protocol = _GctfResumeHarness()
        # mic_002's step ran in the earlier run and was cut off.
        protocol._prevSteps, protocol._steps = protocol._steps, []

        ProtGctf._restoreProcessedMicsFromPersistentState(protocol)
        ProtGctf._checkNewInput(protocol)

        self.assertEqual(["mic_002", "mic_003"], protocol.scheduledMicNames)


class _StalledStream(LogicalSetFake):
    """A closed producer that declares more items than it ever shows."""

    def __init__(self, visibleIds, declared):
        super().__init__([_Mic(i, "mic_%03d" % i) for i in visibleIds],
                         streamClosed=True)
        self._declared = declared

    def getSize(self):
        return self._declared


class _Idle(GctfStreamingBase):
    def __init__(self, micDict, pendingMics):
        self.micDict, self._pendingMics = micDict, pendingMics

    def _hasActiveStreamingWork(self):
        return ProtGctf._hasActiveStreamingWork(self)


class TestAStreamThatNeverCompletesIsReported(unittest.TestCase):
    def _poll(self, protocol, times):
        stream = _StalledStream([1], declared=3)
        for _ in range(times):
            protocol._reconcileClosedStreamIds(stream, [], {1}, True,
                                               '_lastInputId')

    def test_WithNothingInFlightItFailsInsteadOfWaitingForever(self):
        with self.assertRaisesRegex(RuntimeError, "declaring 3 items"):
            self._poll(_Idle({}, {}), GctfStreamingBase.TERMINAL_STALL_POLLS)

    def test_WorkInFlightIsNeverCountedAsAStall(self):
        busy = _Idle({"mic_001": _Mic(1, "mic_001")}, {})
        self._poll(busy, GctfStreamingBase.TERMINAL_STALL_POLLS * 3)


if __name__ == '__main__':
    unittest.main()
