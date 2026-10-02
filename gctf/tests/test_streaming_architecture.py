# **************************************************************************
# *
# * Focused regression tests for Gctf streaming architecture.
# *
# **************************************************************************

import unittest

from gctf.protocols import ProtGctf


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


class _LogicalMicrographSet:
    def __init__(self, items, streamClosed=False):
        self._items = list(items)
        self._streamClosed = streamClosed
        self.loadCalls = 0

    def getFileName(self):
        raise AssertionError(
            "Gctf streaming discovery must not depend on a SQLite/storage filename."
        )

    def loadAllProperties(self):
        self.loadCalls += 1

    def iterItems(self):
        return iter(self._items)

    def isStreamClosed(self):
        return self._streamClosed


class _GctfLogicalInputHarness:
    def __init__(self, inputSet):
        self._inputSet = inputSet
        self.micDict = {}
        self.debugMessages = []

    def getInputMicrographs(self):
        return self._inputSet

    def debug(self, message):
        self.debugMessages.append(message)

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


class _JoinStep:
    def __init__(self):
        self.prerequisites = []

    def addPrerequisites(self, *deps):
        self.prerequisites.extend(deps)


class _GctfInputCheckHarness(_GctfLogicalInputHarness):
    def __init__(self, inputSet):
        super().__init__(inputSet)
        self.streamClosed = False
        self.joinStep = _JoinStep()
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
        return [
            101 + index
            for index, _ in enumerate(newMics)
        ]

    def _getFirstJoinStep(self):
        return self.joinStep

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
        self.assertEqual(
            [101, 102],
            protocol.joinStep.prerequisites,
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


class _GctfFinishedStepsHarness:
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

    def getMicrograph(self):
        return self._mic


class _LogicalCtfOutput:
    def __init__(self):
        self._items = []

    def iterItems(self):
        return iter(self._items)

    def appendMic(self, mic):
        self._items.append(_OutputCtf(mic))


class _GctfOutputCompletionHarness:
    def __init__(self):
        mic = _Mic(1, "mic_001")
        self.micDict = {"mic_001": mic}
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

    def _getPublishedCtfMicNames(self):
        return ProtGctf._getPublishedCtfMicNames(self)

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


class _GctfResumeHarness:
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
        self.micDict = {}
        self.streamClosed = False
        self.scheduledMicNames = []
        self.updateStepsCalls = 0

    def getInputMicrographs(self):
        return self._inputMics

    def _getPublishedCtfMicNames(self):
        return ProtGctf._getPublishedCtfMicNames(self)

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

    def _getFirstJoinStep(self):
        return None

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
        self.assertEqual(
            {"mic_001", "mic_002", "mic_003"},
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

