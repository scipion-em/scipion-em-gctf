# **************************************************************************
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""The streaming generator must abandon its loop once the run is over.

A failed step makes pyworkflow mark the protocol as FAILED and the step
executor break out of its own loop - and then join every thread it
started, the steps generator among them. A generator that keeps polling
is never joined, so the run hangs for good with nothing left to do.
"""
import unittest

import pyworkflow.protocol.constants as cons

from gctf.protocols.protocol_gctf import ProtGctf
from gctf.protocols.protocol_streaming_base import GctfStreamingBase


# A real hang cannot be asserted on, so the harness caps its own polling:
# reaching the cap is the failure signal.
MAX_POLLS = 50


class _Status:
    """Mimics the String attribute pyworkflow keeps the status in."""

    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _GeneratorHarness:
    """The collaborators the generator leans on, stubbed out."""

    _streamingMustStop = GctfStreamingBase._streamingMustStop

    def __init__(self, status=cons.STATUS_RUNNING, closeAfter=None):
        self.status = _Status(status)
        self.finished = False
        self.polls = 0
        self.sleeps = 0
        self._closeAfter = closeAfter

    def _defineCtfParamsDict(self):
        pass

    def _insertInitialSteps(self):
        return []

    def _restoreProcessedMicsFromPersistentState(self):
        pass

    def _checkNewInput(self):
        self.polls += 1
        if self.polls > MAX_POLLS:
            raise AssertionError(
                "The generator polled %d times without leaving its loop: "
                "the executor's join() would never return." % self.polls
            )

    def _checkNewOutput(self):
        if self._closeAfter is not None and self.polls >= self._closeAfter:
            self.finished = True

    def _getStreamingSleepOnWait(self):
        return 1

    def _streamingSleepOnWait(self):
        self.sleeps += 1


class TestGeneratorStopsWhenTheRunIsOver(unittest.TestCase):

    def _runGenerator(self, harness):
        ProtGctf.stepsGeneratorStep(harness)

    def testFailedRunLeavesTheLoopImmediately(self):
        harness = _GeneratorHarness(status=cons.STATUS_FAILED)

        self._runGenerator(harness)

        self.assertEqual(
            harness.polls,
            0,
            "A FAILED protocol must not poll even once: the executor is "
            "already joining this very thread.",
        )

    def testAbortedRunLeavesTheLoopImmediately(self):
        harness = _GeneratorHarness(status=cons.STATUS_ABORTED)

        self._runGenerator(harness)

        self.assertEqual(harness.polls, 0)

    def testFailureMidStreamStopsTheLoop(self):
        harness = _GeneratorHarness()
        realCheck = harness._checkNewInput

        def failOnThirdPoll():
            realCheck()
            if harness.polls == 3:
                harness.status = _Status(cons.STATUS_FAILED)

        harness._checkNewInput = failOnThirdPoll

        self._runGenerator(harness)

        self.assertEqual(
            harness.polls,
            3,
            "The loop must stop on the poll right after the failure, not "
            "keep spinning until the producer closes the stream.",
        )

    def testHealthyRunStillPollsUntilTheStreamCloses(self):
        harness = _GeneratorHarness(closeAfter=4)

        self._runGenerator(harness)

        self.assertEqual(
            harness.polls,
            4,
            "A running protocol must keep polling: the stop guard must "
            "not short-circuit a healthy stream.",
        )

    def testAPlainStringStatusIsUnderstood(self):
        """status is not always wrapped: accept the bare value too."""
        harness = _GeneratorHarness()
        harness.status = cons.STATUS_FAILED

        self._runGenerator(harness)

        self.assertEqual(harness.polls, 0)

    def testAMissingStatusDoesNotBreakTheLoop(self):
        harness = _GeneratorHarness(closeAfter=1)
        del harness.status

        self._runGenerator(harness)

        self.assertEqual(harness.polls, 1)


if __name__ == '__main__':
    unittest.main()
