# **************************************************************************
# *
# * Focused regression tests for Gctf streaming architecture.
# *
# **************************************************************************

import tempfile
import os
import unittest
from collections import OrderedDict

import pyworkflow.utils as pwutils

from gctf.protocols import ProtGctf
from gctf.protocols.program_gctf import ProgramGctf
from gctf.protocols.protocol_streaming_base import GctfStreamingBase

from .logical_set_fakes import LogicalSetFake


class _Mic:
    def __init__(self, objId, micName):
        self._objId = objId
        self._micName = micName

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return self._micName

    def clone(self):
        return _Mic(self._objId, self._micName)


class _LogicalMicrographSet(LogicalSetFake):
    pass


class _GctfLogicalInputHarness(GctfStreamingBase):
    def __init__(self, inputSet):
        self._inputSet = inputSet
        self.micDict = OrderedDict()
        self._pendingMics = OrderedDict()
        self._knownMicIds = set()
        self._lastInputId = 0
        self._steps = []
        self.debugMessages = []

    def getInputMicrographs(self):
        return self._inputSet

    def debug(self, message):
        self.debugMessages.append(message)

    def _getScheduledCtfMicNames(self):
        return ProtGctf._getScheduledCtfMicNames(self)

    def _loadSet(self, inputSet, SetClass, getKeyFunc):
        return ProtGctf._loadSet(
            self,
            inputSet,
            SetClass,
            getKeyFunc,
        )


class TestGctfStreamingArchitecture(unittest.TestCase):
    def test_GctfStreamingLoadsLogicalMicrographsWithoutStorageFilename(self):
        inputSet = _LogicalMicrographSet(
            [
                _Mic(1, "mic_001"),
                _Mic(2, "mic_002"),
            ],
            streamClosed=False,
        )
        protocol = _GctfLogicalInputHarness(inputSet)

        newMics, streamClosed = ProtGctf._loadInputList(protocol)

        self.assertEqual(
            ["mic_001", "mic_002"],
            list(newMics),
        )
        self.assertFalse(streamClosed)
        self.assertEqual(1, inputSet.loadCalls)


if __name__ == "__main__":
    unittest.main()


class _NoJoinStepGuard:
    """Fail loudly if the protocol goes back to the legacy join-step pattern.

    Streaming must not create a WAITING output step that the generator
    unlocks with setStatus(STATUS_NEW).
    """

    def _getFirstJoinStep(self):
        raise AssertionError(
            "Streaming must not depend on a waiting join step."
        )


class _GctfInputCheckHarness(_NoJoinStepGuard, _GctfLogicalInputHarness):
    def __init__(self, inputSet):
        super().__init__(inputSet)
        self.streamClosed = False
        self.insertedMicNames = []
        self.updateStepsCalls = 0

    def _loadInputList(self):
        return ProtGctf._loadInputList(self)

    def _insertNewMicsSteps(self, newMics):
        newMics = list(newMics)
        self.insertedMicNames.extend(
            mic.getMicName()
            for mic in newMics
        )
        for mic in newMics:
            self.micDict[mic.getMicName()] = mic

        return [
            101 + index
            for index, _ in enumerate(newMics)
        ]

    def updateSteps(self):
        self.updateStepsCalls += 1


class TestGctfStreamingInputChecks(unittest.TestCase):
    def test_GctfStreamingChecksLogicalInputWithoutFilesystemMtime(self):
        inputSet = _LogicalMicrographSet(
            [
                _Mic(1, "mic_001"),
                _Mic(2, "mic_002"),
            ],
            streamClosed=False,
        )
        protocol = _GctfInputCheckHarness(inputSet)

        ProtGctf._checkNewInput(protocol)

        self.assertEqual(
            ["mic_001", "mic_002"],
            protocol.insertedMicNames,
        )
        self.assertEqual(1, protocol.updateStepsCalls)
        self.assertFalse(protocol.streamClosed)


class _StoredValue:
    def __init__(self, value):
        self._value = value

    def get(self, default=None):
        return self._value if self._value is not None else default


class _CtfStep:
    def __init__(self, funcName, argsStr, finished):
        self.funcName = _StoredValue(funcName)
        self.argsStr = _StoredValue(argsStr)
        self._finished = finished

    def isFinished(self):
        return self._finished


class _GctfFinishedStepsHarness(GctfStreamingBase):
    def __init__(self):
        self._steps = [
            _CtfStep(
                "estimateCtfStep",
                '["mic_001"]',
                True,
            ),
            _CtfStep(
                "estimateCtfListStep",
                '[["mic_002", "mic_003"]]',
                True,
            ),
            _CtfStep(
                "estimateCtfStep",
                '["mic_pending"]',
                False,
            ),
            _CtfStep(
                "createOutputStep",
                '[]',
                True,
            ),
        ]


class TestGctfFinishedStepState(unittest.TestCase):
    def test_GctfGetsFinishedMicrographsFromPersistedSteps(self):
        protocol = _GctfFinishedStepsHarness()

        finishedNames = ProtGctf._getFinishedCtfMicNames(protocol)

        self.assertEqual(
            {"mic_001", "mic_002", "mic_003"},
            finishedNames,
        )


class _OutputCtf:
    def __init__(self, mic):
        self._mic = mic

    def getObjId(self):
        return self._mic.getObjId()

    def getMicrograph(self):
        return self._mic


class _LogicalCtfOutput(LogicalSetFake):
    def __init__(self):
        super().__init__([], streamClosed=False)

    def appendMic(self, mic):
        self._items.append(_OutputCtf(mic))


class _GctfOutputCompletionHarness(GctfStreamingBase):
    def __init__(self):
        mic = _Mic(1, "mic_001")
        self.micDict = OrderedDict([("mic_001", mic)])
        self._pendingMics = OrderedDict()
        self._steps = [
            _CtfStep(
                "estimateCtfStep",
                '["mic_001"]',
                True,
            ),
        ]
        self.streamClosed = False
        self.finished = False
        self.outputCTF = _LogicalCtfOutput()
        self.publishedBatches = []
        self.sleepCalls = 0

    def _getFinishedCtfMicNames(self):
        return ProtGctf._getFinishedCtfMicNames(self)

    def _updateOutputCTFSet(self, micList, streamMode):
        micList = list(micList)
        self.publishedBatches.append(
            [mic.getMicName() for mic in micList]
        )
        for mic in micList:
            self.outputCTF.appendMic(mic)
        return micList

    def _streamingSleepOnWait(self):
        self.sleepCalls += 1

    def _readDoneList(self):
        raise AssertionError("Gctf completion must not read DONE/all.TXT.")

    def _isMicDone(self, mic):
        raise AssertionError("Gctf completion must not read per-micrograph DONE files.")

    def _writeDoneList(self, micList):
        raise AssertionError("Gctf completion must not write DONE/all.TXT.")


class TestGctfStreamingCompletion(unittest.TestCase):
    def test_GctfPublishesFinishedStepWithoutDoneSidecars(self):
        protocol = _GctfOutputCompletionHarness()

        ProtGctf._checkNewOutput(protocol)
        ProtGctf._checkNewOutput(protocol)

        self.assertEqual(
            [["mic_001"]],
            protocol.publishedBatches,
        )
        self.assertFalse(protocol.finished)


class _GctfProcessingHarness:
    def __init__(self):
        self.micDict = {
            "mic_001": _Mic(1, "mic_001"),
            "mic_002": _Mic(2, "mic_002"),
        }
        self.singleCalls = []
        self.batchCalls = []
        self.infoMessages = []

    def _estimateCTF(self, mic, *args):
        self.singleCalls.append(mic.getMicName())

    def _estimateCtfList(self, micList, *args):
        self.batchCalls.append(
            [mic.getMicName() for mic in micList]
        )

    def _getMicrographDone(self, mic):
        raise AssertionError(
            "Gctf processing must not access per-micrograph DONE files."
        )

    def _isMicDone(self, mic):
        raise AssertionError(
            "Gctf processing must not inspect per-micrograph DONE files."
        )

    def isContinued(self):
        return True

    def info(self, message):
        self.infoMessages.append(message)


class TestGctfStreamingProcessing(unittest.TestCase):
    def test_GctfSingleProcessingDoesNotUseDoneSidecar(self):
        protocol = _GctfProcessingHarness()

        ProtGctf.estimateCtfStep(
            protocol,
            "mic_001",
        )

        self.assertEqual(
            ["mic_001"],
            protocol.singleCalls,
        )

    def test_GctfBatchProcessingDoesNotUseDoneSidecars(self):
        protocol = _GctfProcessingHarness()

        ProtGctf.estimateCtfListStep(
            protocol,
            ["mic_001", "mic_002"],
        )

        self.assertEqual(
            [["mic_001", "mic_002"]],
            protocol.batchCalls,
        )


class _InsertedGeneratorStep:
    def __init__(self, funcName):
        self.funcName = funcName


class _GctfGeneratorHarness:
    def __init__(self):
        self.insertedSteps = []
        self.donePathAccesses = 0
        self.inputLoads = 0

    def resumableStepGeneratorStep(self, timestamp):
        pass

    def _insertFunctionStep(self, func, *args, **kwargs):
        funcName = func.__name__ if callable(func) else func
        self.insertedSteps.append(
            _InsertedGeneratorStep(funcName)
        )
        return len(self.insertedSteps)

    def _getAllDone(self):
        self.donePathAccesses += 1
        raise AssertionError(
            "Generator orchestration must not create or access DONE/all.TXT."
        )

    def _loadInputList(self):
        self.inputLoads += 1
        raise AssertionError(
            "Generator orchestration must not snapshot input during _insertAllSteps()."
        )


class TestGctfGeneratorOrchestration(unittest.TestCase):
    def test_GctfStreamingUsesSingleResumableGenerator(self):
        protocol = _GctfGeneratorHarness()

        ProtGctf._insertAllSteps(protocol)

        self.assertEqual(
            ["resumableStepGeneratorStep"],
            [step.funcName for step in protocol.insertedSteps],
        )
        self.assertEqual(0, protocol.donePathAccesses)
        self.assertEqual(0, protocol.inputLoads)


class _GctfResumeHarness(GctfStreamingBase):
    def __init__(self):
        self._inputMics = _LogicalMicrographSet(
            [
                _Mic(1, "mic_001"),
                _Mic(2, "mic_002"),
                _Mic(3, "mic_003"),
            ],
            streamClosed=False,
        )
        self.outputCTF = _LogicalCtfOutput()
        self.outputCTF.appendMic(_Mic(1, "mic_001"))
        self._steps = [
            _CtfStep(
                "estimateCtfStep",
                '["mic_002"]',
                False,
            ),
        ]
        self.micDict = OrderedDict()
        self._pendingMics = OrderedDict()
        self._knownMicIds = set()
        self._lastInputId = 0
        self.streamClosed = False
        self.scheduledMicNames = []
        self.updateStepsCalls = 0

    def getInputMicrographs(self):
        return self._inputMics

    def _getFinishedCtfMicNames(self):
        return ProtGctf._getFinishedCtfMicNames(self)

    def _getScheduledCtfMicNames(self):
        return ProtGctf._getScheduledCtfMicNames(self)

    def _loadSet(self, inputSet, SetClass, getKeyFunc):
        return ProtGctf._loadSet(
            self,
            inputSet,
            SetClass,
            getKeyFunc,
        )

    def _loadInputList(self):
        return ProtGctf._loadInputList(self)

    def _insertNewMicsSteps(self, newMics):
        newMics = list(newMics)
        self.scheduledMicNames.extend(
            mic.getMicName()
            for mic in newMics
        )
        for mic in newMics:
            self.micDict[mic.getMicName()] = mic
        return [
            201 + index
            for index, _ in enumerate(newMics)
        ]

    def updateSteps(self):
        self.updateStepsCalls += 1

    def debug(self, message):
        pass


class TestGctfGeneratorResumeSafety(unittest.TestCase):
    def test_GctfResumeDoesNotReschedulePublishedOrScheduledMicrographs(self):
        protocol = _GctfResumeHarness()

        ProtGctf._restoreProcessedMicsFromPersistentState(protocol)
        ProtGctf._checkNewInput(protocol)

        self.assertEqual(
            ["mic_003"],
            protocol.scheduledMicNames,
        )
        # micDict now holds only what still has to be published, so the
        # micrograph a previous run already published stays out of it.
        self.assertEqual(
            {"mic_002", "mic_003"},
            set(protocol.micDict),
        )
        self.assertEqual(1, protocol.updateStepsCalls)


class _ThreadParam:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _StreamingInput:
    def __init__(self, streamOpen):
        self._streamOpen = streamOpen

    def isStreamOpen(self):
        return self._streamOpen


class _GctfThreadValidationHarness:
    def __init__(self, threads, streamOpen=True):
        self.numberOfThreads = _ThreadParam(threads)
        self._inputMics = _StreamingInput(streamOpen)

    def getInputMicrographs(self):
        return self._inputMics


class TestGctfStreamingThreadValidation(unittest.TestCase):
    def test_GctfStreamingRejectsTwoThreads(self):
        protocol = _GctfThreadValidationHarness(
            threads=2,
            streamOpen=True,
        )

        errors = ProtGctf._validateStreamingThreads(protocol)

        self.assertTrue(errors)

    def test_GctfStreamingAcceptsThreeThreads(self):
        protocol = _GctfThreadValidationHarness(
            threads=3,
            streamOpen=True,
        )

        errors = ProtGctf._validateStreamingThreads(protocol)

        self.assertEqual([], errors)

    def test_GctfNonStreamingStillAcceptsOneThread(self):
        protocol = _GctfThreadValidationHarness(
            threads=1,
            streamOpen=False,
        )

        errors = ProtGctf._validateStreamingThreads(protocol)

        self.assertEqual([], errors)


class _Value:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _GctfArgsProtocol:
    def __init__(self):
        self.minDefocus = _Value(5000.0)
        self.maxDefocus = _Value(90000.0)
        self.stepDefocus = _Value(500.0)
        self.astigmatism = 1000.0
        self.doEPA = True
        self.plotResRing = True
        self.bfactor = 150
        self.overlap = 0.5
        self.convsize = 85
        self.doHighRes = False
        self.smoothResL = 1000
        self.EPAsmp = 4
        self.doPhShEst = False
        self.doValidate = False

    def getCtfParamsDict(self):
        return {
            "samplingRate": 1.2,
            "voltage": 300.0,
            "sphericalAberration": 2.7,
            "ampContrast": 0.1,
            "scannedPixelSize": None,
            "lowRes": 50.0,
            "highRes": 4.0,
            "windowSize": 1024,
        }


class TestGctfOptionalDetectorPixelSize(unittest.TestCase):
    def test_GctfOmitsDstepWhenScannedPixelSizeIsUnavailable(self):
        protocol = _GctfArgsProtocol()

        program = ProgramGctf.__new__(ProgramGctf)
        args, params = program._getArgs(protocol)
        params["GPU"] = "0"

        self.assertNotIn("--dstep", args)
        formatted = args % params
        self.assertIn("--apix 1.200000", formatted)


class _ProcessingValue:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _ProcessingMic:
    def __init__(self, objId, fileName):
        self._objId = objId
        self._fileName = fileName

    def getObjId(self):
        return self._objId

    def getFileName(self):
        return self._fileName


class _FailingCommandGctfProgram:
    def getCommand(self, **kwargs):
        raise RuntimeError("synthetic Gctf command failure")


class _MissingOutputGctfProgram:
    def getCommand(self, **kwargs):
        return "gctf", "--synthetic"

    def getExt(self):
        return ".epa"


class _GctfFailurePropagationHarness:
    def __init__(self, root, program):
        self._root = root
        self._gctfProgram = program
        self.ctfDownFactor = _ProcessingValue(1)
        self._params = {"scannedPixelSize": 1.0}
        self.errors = []

    def _getMicrographDir(self, mic):
        return os.path.join(self._root, "mic_tmp")

    def runJob(self, program, params, env=None):
        pass

    def error(self, message):
        self.errors.append(message)

    def _getPsdPath(self, micFn):
        return os.path.join(self._root, "out_ctf.mrc")

    def _getCtfOutPath(self, micFn):
        return os.path.join(self._root, "out_ctf.log")

    def _getCtfFitOutPath(self, micFn):
        return os.path.join(self._root, "out_EPA.log")


class _PartialOutputHarness(_GctfFailurePropagationHarness):
    """Only the second micrograph of the batch produces results."""

    def __init__(self, root):
        super().__init__(root, _MissingOutputGctfProgram())

    def runJob(self, program, params, env=None):
        # A multi-micrograph batch gets its own suffixed directory, so take
        # the real one from the command the protocol just built.
        micPath = params.rsplit(" ", 1)[-1][:-len("/*.mrc")]

        for suffix in (".epa", "_gctf.log", "_EPA.log"):
            with open(os.path.join(micPath, "mic_002" + suffix), "w") as fh:
                fh.write("result\n")

    def collectedFor(self, micBase):
        return os.path.join(self._root, micBase + "_out_ctf.mrc")

    def _getPsdPath(self, micFn):
        return self.collectedFor(pwutils.removeBaseExt(micFn))

    def _getCtfOutPath(self, micFn):
        return os.path.join(
            self._root, pwutils.removeBaseExt(micFn) + "_out_ctf.log")

    def _getCtfFitOutPath(self, micFn):
        return os.path.join(
            self._root, pwutils.removeBaseExt(micFn) + "_out_EPA.log")


class TestGctfProcessingFailurePropagation(unittest.TestCase):
    def test_GctfCommandFailureFailsProcessingStep(self):
        with tempfile.TemporaryDirectory() as root:
            micFn = os.path.join(root, "mic_001.mrc")
            open(micFn, "wb").close()
            protocol = _GctfFailurePropagationHarness(
                root,
                _FailingCommandGctfProgram(),
            )

            with self.assertRaisesRegex(
                RuntimeError,
                "synthetic Gctf command failure",
            ):
                ProtGctf._estimateCtfList(
                    protocol,
                    [_ProcessingMic(1, micFn)],
                )

    def test_GctfMissingOutputFailsProcessingStep(self):
        # Nothing came out of the batch at all, so the step has to fail
        # rather than report the micrograph as estimated.
        with tempfile.TemporaryDirectory() as root:
            micFn = os.path.join(root, "mic_001.mrc")
            open(micFn, "wb").close()
            protocol = _GctfFailurePropagationHarness(
                root,
                _MissingOutputGctfProgram(),
            )

            with self.assertRaises(Exception):
                ProtGctf._estimateCtfList(
                    protocol,
                    [_ProcessingMic(1, micFn)],
                )

    def test_GctfBatchWithOneGoodResultDoesNotFailTheWholeStep(self):
        # The counterpart of the test above, and the reason the step cannot
        # simply re-raise: one bad micrograph must not discard the others.
        with tempfile.TemporaryDirectory() as root:
            micFns = []

            for index in (1, 2):
                micFn = os.path.join(root, "mic_%03d.mrc" % index)
                open(micFn, "wb").close()
                micFns.append(micFn)

            protocol = _PartialOutputHarness(root)

            ProtGctf._estimateCtfList(
                protocol,
                [_ProcessingMic(index + 1, micFn)
                 for index, micFn in enumerate(micFns)],
            )

            self.assertEqual(1, len(protocol.errors))
            self.assertTrue(os.path.exists(protocol.collectedFor("mic_002")))



class _GctfCostHarness(_GctfLogicalInputHarness):
    """Discovery over a Set that reports what it hydrated."""

    def __init__(self, inputSet):
        super().__init__(inputSet)
        self.streamClosed = False
        self.insertedMicNames = []
        self.updateStepsCalls = 0

    def _getScheduledCtfMicNames(self):
        return ProtGctf._getScheduledCtfMicNames(self)

    def _getFinishedCtfMicNames(self):
        return ProtGctf._getFinishedCtfMicNames(self)

    def _loadInputList(self):
        return ProtGctf._loadSet(
            self, self._inputSet, None, lambda mic: mic.getMicName())

    def _insertNewMicsSteps(self, newMics):
        newMics = list(newMics)
        self.insertedMicNames.extend(mic.getMicName() for mic in newMics)

        for mic in newMics:
            self.micDict[mic.getMicName()] = mic

        return []

    def updateSteps(self):
        self.updateStepsCalls += 1


class TestGctfStreamingDiscoveryCost(unittest.TestCase):
    """A poll must cost what just arrived, not everything seen so far."""

    @staticmethod
    def _mics(firstId, count):
        return [_Mic(micId, "mic_%03d" % micId)
                for micId in range(firstId, firstId + count)]

    def test_GctfPollOnlyHydratesTheMicrographsThatJustArrived(self):
        inputSet = _LogicalMicrographSet(self._mics(1, 500))
        protocol = _GctfCostHarness(inputSet)

        ProtGctf._checkNewInput(protocol)

        self.assertEqual(500, inputSet.hydratedItems)

        inputSet.addItems(self._mics(501, 3))
        hydratedBefore = inputSet.hydratedItems

        ProtGctf._checkNewInput(protocol)

        self.assertEqual(3, inputSet.hydratedItems - hydratedBefore)
        self.assertEqual(0, inputSet.fullScans)

    def test_GctfIdlePollHydratesNothingAtAll(self):
        inputSet = _LogicalMicrographSet(self._mics(1, 500))
        protocol = _GctfCostHarness(inputSet)

        ProtGctf._checkNewInput(protocol)
        hydratedBefore = inputSet.hydratedItems

        ProtGctf._checkNewInput(protocol)

        self.assertEqual(hydratedBefore, inputSet.hydratedItems)

    def test_GctfStepGraphScanDoesNotReparseStepsItAlreadyRead(self):
        protocol = _GctfCostHarness(_LogicalMicrographSet([]))

        parsed = []
        original = GctfStreamingBase._parseStepArgKeys

        def countingParse(step, dictField, keyType):
            parsed.append(step)
            return original(step, dictField, keyType)

        # A plain function, not staticmethod(): an instance attribute is
        # never bound, and a staticmethod object is only callable itself
        # from Python 3.10 on.
        protocol._parseStepArgKeys = countingParse

        protocol._steps = [
            _CtfStep("estimateCtfStep", '["mic_%03d"]' % micId, True)
            for micId in range(1, 201)
        ]

        self.assertEqual(200, len(protocol._getFinishedCtfMicNames()))
        self.assertEqual(200, len(parsed))

        protocol._steps.append(_CtfStep("estimateCtfStep", '["mic_201"]', True))

        self.assertEqual(201, len(protocol._getFinishedCtfMicNames()))
        # Only the new step is parsed again, not the 200 already read.
        self.assertEqual(201, len(parsed))

    def test_GctfStepGraphScanReadsStepsRestoredOnResume(self):
        # Work finished before a Continue lives in _prevSteps; reading only
        # _steps would make it invisible and schedule it all over again.
        protocol = _GctfCostHarness(_LogicalMicrographSet([]))
        protocol._steps = [_CtfStep("estimateCtfStep", '["mic_002"]', True)]
        protocol._prevSteps = [_CtfStep("estimateCtfStep", '["mic_001"]', True)]

        self.assertEqual({"mic_001", "mic_002"},
                         protocol._getFinishedCtfMicNames())

    def test_GctfStepGraphScanSeesAPendingStepFinishAfterTheListIsRebuilt(self):
        protocol = _GctfCostHarness(_LogicalMicrographSet([]))
        protocol._steps = [_CtfStep("estimateCtfStep", '["mic_001"]', False)]

        self.assertEqual(set(), protocol._getFinishedCtfMicNames())

        # pyworkflow can rebuild the step list from the database, handing
        # back new objects for the same positions.
        protocol._steps = [_CtfStep("estimateCtfStep", '["mic_001"]', True)]

        self.assertEqual({"mic_001"}, protocol._getFinishedCtfMicNames())

    def test_GctfResumeStartsAboveWhatIsPublishedButKeepsTheGaps(self):
        inputSet = _LogicalMicrographSet(self._mics(1, 100))
        protocol = _GctfCostHarness(inputSet)

        # A previous run published everything except mic_042.
        publishedIds = set(range(1, 101)) - {42}

        watermark, gapIds = protocol._resumeWatermarkWithGaps(
            inputSet, publishedIds)

        self.assertEqual(100, watermark)
        self.assertEqual({42}, gapIds)
        self.assertEqual(0, inputSet.hydratedItems)
