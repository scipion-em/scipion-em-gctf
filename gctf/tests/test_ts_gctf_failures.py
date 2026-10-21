# **************************************************************************
# *
# * Gctf on tilt series: a tilt series with no usable CTF is left out and
# * estimated again on Continue, the run fails only when every one failed,
# * and the gctf logs outlive the run's tmp folder.
# *
# **************************************************************************

import os
import tempfile
import unittest
from unittest import mock

from pyworkflow.utils import weakImport

with weakImport('tomo'):
    from gctf.protocols.protocol_ts_gctf import ProtTsGctf, TS_STEPS


    class _Ts:
        def __init__(self, tsId):
            self._tsId = tsId

        def getTsId(self):
            return self._tsId


    class _Harness:
        _stepResultPath = ProtTsGctf._stepResultPath
        _stepResultFile = ProtTsGctf._stepResultFile
        _prepareFailedForRetry = ProtTsGctf._prepareFailedForRetry
        _failedTsIds = ProtTsGctf._failedTsIds
        _createOutput = ProtTsGctf._createOutput
        closeOutputSetsStep = ProtTsGctf.closeOutputSetsStep

        def __init__(self, root, published=()):
            self._root = root
            self.tsDict = {'TS_01': _Ts('TS_01'), 'TS_02': _Ts('TS_02')}
            self._published = set(published)
            self.errors = []
            self._closeOutputSet = mock.Mock()
            self.getOutputCtfTomoSet = mock.Mock(
                side_effect=AssertionError("Nothing must be published for it."))

        def _getExtraPath(self, *paths):
            return os.path.join(self._root, 'extra', *paths)

        def _publishedTsIds(self):
            return set(self._published)

        def error(self, message):
            self.errors.append(message)


    class TestATiltSeriesWithoutCtfIsLeftOut(unittest.TestCase):
        def test_ItIsNotPublished(self):
            with tempfile.TemporaryDirectory() as tmp:
                prot = _Harness(tmp)
                prot._readTsCtfs = lambda ts: None

                prot._createOutput('TS_01')

            self.assertEqual({'TS_01'}, prot._failedTsIds())
            self.assertTrue(prot.errors)

        def test_EveryTiltSeriesFailedFailsTheRunButClosesTheOutput(self):
            with tempfile.TemporaryDirectory() as tmp:
                prot = _Harness(tmp)
                prot._failedTsIds().update({'TS_01', 'TS_02'})

                with self.assertRaisesRegex(Exception, "failed for every one of the 2"):
                    prot.closeOutputSetsStep()

            prot._closeOutputSet.assert_called_once()

        def test_SomePublishedTiltSeriesAreNoFailure(self):
            with tempfile.TemporaryDirectory() as tmp:
                prot = _Harness(tmp, published=['TS_01'])
                prot._failedTsIds().add('TS_02')

                prot.closeOutputSetsStep()


    class TestContinueEstimatesAgainWhatFailed(unittest.TestCase):
        def test_OnlyTheStepsOfTheFailedTiltSeriesRunAgain(self):
            with tempfile.TemporaryDirectory() as tmp:
                prot = _Harness(tmp, published=['TS_01'])
                markers = {}
                for tsId in prot.tsDict:
                    for kind in TS_STEPS:
                        name = '%s_%s' % (kind, tsId)
                        markers[name] = prot._stepResultFile(name)
                markers['close'] = prot._stepResultFile('close')

                prot._prepareFailedForRetry()

                remaining = {name for name, path in markers.items() if os.path.exists(path)}

            self.assertEqual({'%s_TS_01' % kind for kind in TS_STEPS}, remaining,
                             "TS_02 is not published: its steps, and the closing one, run again.")


    class TestTheGctfLogOutlivesTmp(unittest.TestCase):
        def test_TheLogIsReadFromExtra(self):
            read = []

            class _Program:
                def parseOutputAsCtf(self, fileName, psdFile=None):
                    read.append(fileName)
                    return mock.Mock()

            class _Ti:
                def getTsId(self):
                    return 'TS_01'

                def getObjId(self):
                    return 4

            with tempfile.TemporaryDirectory() as tmp:
                prot = mock.Mock()
                prot.getTiPrefix = ProtTsGctf.getTiPrefix
                prot._getExtraPath = lambda *paths: os.path.join(tmp, 'extra', *paths)
                prot._getTmpPath = lambda *paths: os.path.join(tmp, 'tmp', *paths)
                prot._gctfProgram = _Program()
                os.makedirs(prot._getExtraPath())
                open(prot._getExtraPath('TS_01_004_gctf.log'), 'w').close()

                with mock.patch('gctf.protocols.protocol_ts_gctf.CTFTomo'):
                    ProtTsGctf.getCtf(prot, _Ti())

            self.assertEqual([os.path.join(tmp, 'extra', 'TS_01_004_gctf.log')], read)


if __name__ == '__main__':
    unittest.main()
