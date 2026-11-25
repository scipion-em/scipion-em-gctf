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
"""Per-micrograph artefacts are named after the micrograph's basename.

A Set can perfectly well hold /data/sessionA/mic001.mrc and
/data/sessionB/mic001.mrc at once - different micrographs, different
science, one basename. Everything this protocol writes for a micrograph
is keyed on that basename alone, so the second one silently overwrites
the first and the CTF read back afterwards belongs to the wrong image.
"""
import os
import shlex
import shutil
import tempfile
import unittest

from gctf.protocols.protocol_gctf import ProtGctf


class _Mic:
    def __init__(self, objId, fileName):
        self._objId = objId
        self._fileName = fileName

    def getObjId(self):
        return self._objId

    def getFileName(self):
        return self._fileName


class _PathHarness:
    """Exercises the real per-micrograph path helpers."""

    _getPsdPath = ProtGctf._getPsdPath
    _getCtfOutPath = ProtGctf._getCtfOutPath
    _getCtfFitOutPath = ProtGctf._getCtfFitOutPath
    _itemScopedPath = ProtGctf._itemScopedPath
    _getMicArtefactBase = ProtGctf._getMicArtefactBase

    def __init__(self, extraDir):
        self._extraDir = extraDir

    def _getExtraPath(self, *paths):
        return os.path.join(self._extraDir, *paths)


SAME_BASENAME_A = '/data/sessionA/mic001.mrc'
SAME_BASENAME_B = '/data/sessionB/mic001.mrc'


class TestTwoMicrographsNeverShareAnArtefact(unittest.TestCase):

    def setUp(self):
        self.extraDir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.extraDir, True)
        self.harness = _PathHarness(self.extraDir)

    def _paths(self, mic):
        return (
            self.harness._getPsdPath(mic),
            self.harness._getCtfOutPath(mic),
            self.harness._getCtfFitOutPath(mic),
        )

    def testThePsdOfOneDoesNotOverwriteTheOther(self):
        first = self.harness._getPsdPath(_Mic(1, SAME_BASENAME_A))
        second = self.harness._getPsdPath(_Mic(2, SAME_BASENAME_B))

        self.assertNotEqual(
            first,
            second,
            "Both micrographs write their power spectrum to the same "
            "file: whichever finishes last wins and the other's CTF is "
            "read back from an image it does not belong to.",
        )

    def testTheCtfLogOfOneDoesNotOverwriteTheOther(self):
        first = self.harness._getCtfOutPath(_Mic(1, SAME_BASENAME_A))
        second = self.harness._getCtfOutPath(_Mic(2, SAME_BASENAME_B))

        self.assertNotEqual(first, second)

    def testTheFitLogOfOneDoesNotOverwriteTheOther(self):
        first = self.harness._getCtfFitOutPath(_Mic(1, SAME_BASENAME_A))
        second = self.harness._getCtfFitOutPath(_Mic(2, SAME_BASENAME_B))

        self.assertNotEqual(first, second)

    def testEveryArtefactOfOneMicrographIsDistinct(self):
        psd, log, fit = self._paths(_Mic(1, SAME_BASENAME_A))

        self.assertEqual(len({psd, log, fit}), 3)

    def testThePathIsStableForTheSameMicrograph(self):
        """Resume re-derives these paths and must land on the same files."""
        mic = _Mic(7, SAME_BASENAME_A)

        self.assertEqual(self._paths(mic), self._paths(_Mic(7, SAME_BASENAME_A)))

    def testAnArtefactFromAnOlderRunIsStillFound(self):
        """Continuing a project written before the paths were scoped.

        Those runs wrote the unscoped name, and the files are the user's
        results: a Continue has to keep finding them rather than treat
        the micrograph as never estimated.
        """
        legacy = os.path.join(self.extraDir, 'mic001_ctf.log')
        with open(legacy, 'w') as handle:
            handle.write('legacy')

        self.assertEqual(
            self.harness._getCtfOutPath(_Mic(1, SAME_BASENAME_A)),
            legacy,
            "A result file an earlier run already wrote must still be "
            "the one this run reads.",
        )

    def testAScopedArtefactWinsOverALegacyOne(self):
        legacy = os.path.join(self.extraDir, 'mic001_ctf.log')
        with open(legacy, 'w') as handle:
            handle.write('legacy')

        mic = _Mic(1, SAME_BASENAME_A)
        scoped = self.harness._getCtfOutPath(mic)
        self.assertEqual(scoped, legacy)

        # Once this run has written its own scoped file, that is the one.
        with open(os.path.join(self.extraDir, '000001__mic001_ctf.log'),
                  'w') as handle:
            handle.write('scoped')

        self.assertTrue(
            self.harness._getCtfOutPath(mic).endswith(
                '000001__mic001_ctf.log'),
            "Once the scoped artefact exists it is the authoritative one.",
        )


if __name__ == '__main__':
    unittest.main()


class _BatchWorkspaceHarness:
    """Exercises the name a micrograph gets inside the batch folder."""

    _getMicArtefactBase = ProtGctf._getMicArtefactBase
    _getBatchMicBase = ProtGctf._getBatchMicBase


class TestOneBatchNeverConvertsTwoMicrographsOntoOneFile(unittest.TestCase):
    """Within a batch every micrograph is converted into the same folder.

    gctf is then pointed at that folder with a wildcard, so two
    micrographs sharing a basename become a single input file: one is
    converted on top of the other before gctf even runs, and both end up
    reading the surviving one's results.
    """

    def setUp(self):
        self.harness = _BatchWorkspaceHarness()

    def testTwoMicrographsInOneBatchGetDifferentFiles(self):
        first = self.harness._getBatchMicBase(_Mic(1, SAME_BASENAME_A))
        second = self.harness._getBatchMicBase(_Mic(2, SAME_BASENAME_B))

        self.assertNotEqual(
            first,
            second,
            "Both micrographs are converted to the same file inside the "
            "batch folder, so one overwrites the other before gctf runs.",
        )

    def testTheNameIsStableForTheSameMicrograph(self):
        """The collector re-derives it to pick gctf's output back up."""
        mic = _Mic(3, SAME_BASENAME_A)

        self.assertEqual(
            self.harness._getBatchMicBase(mic),
            self.harness._getBatchMicBase(_Mic(3, SAME_BASENAME_A)),
        )

    def testTheNameStillCarriesTheOriginalBasename(self):
        """Keep it recognisable in logs and in gctf's own output."""
        self.assertIn(
            'mic001',
            self.harness._getBatchMicBase(_Mic(1, SAME_BASENAME_A)),
        )


class _CommandHarness:
    """Captures the command line the protocol hands to the shell."""

    _getMicArtefactBase = ProtGctf._getMicArtefactBase
    _getBatchMicBase = ProtGctf._getBatchMicBase

    def __init__(self, micPath):
        self._micPath = micPath

    def _batchInputArgument(self):
        return ProtGctf._batchInputArgument(self, self._micPath)


class TestTheBatchArgumentSurvivesASpace(unittest.TestCase):
    """The batch folder is handed to a shell as part of the command.

    Scipion projects live wherever the user put them, and a directory
    with a space in its name is ordinary. Unquoted, the shell reads it as
    two arguments and gctf is pointed at something else entirely.
    """

    def testAPathWithASpaceStaysOneArgument(self):
        harness = _CommandHarness('/home/me/my project/tmp/mic_0001')

        argument = harness._batchInputArgument()

        # Split it the way the shell would.
        tokens = shlex.split(argument)

        self.assertEqual(
            tokens,
            ['/home/me/my project/tmp/mic_0001/*.mrc'],
            "The shell reads this as more than one argument, so gctf is "
            "pointed at something other than the batch folder: %r"
            % argument,
        )

    def testTheWildcardStaysOutsideTheQuotes(self):
        harness = _CommandHarness('/tmp/mic_0001')

        argument = harness._batchInputArgument()

        self.assertTrue(
            argument.endswith('/*.mrc'),
            "The shell still has to expand the wildcard, so it must not "
            "be quoted along with the folder: %r" % argument,
        )

    def testAnOrdinaryPathIsUnchangedInSubstance(self):
        harness = _CommandHarness('/tmp/mic_0001')

        self.assertIn('/tmp/mic_0001', harness._batchInputArgument())
