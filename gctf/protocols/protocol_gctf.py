# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import shlex
import time
from collections import OrderedDict
from datetime import datetime

import pyworkflow.utils as pwutils
from pyworkflow.object import Boolean, Set
from pyworkflow.constants import PROD
from pwem import emlib
from pwem.objects import CTFModel
from pwem.protocols import ProtCTFMicrographs


from gctf import Plugin
from gctf.protocols.program_gctf import ProgramGctf
from gctf.protocols.protocol_streaming_base import GctfStreamingBase

# Module level: these helpers are called unbound on light test harnesses.
CTF_STEP_NAMES = ('estimateCtfStep', 'estimateCtfListStep')


class ProtGctf(GctfStreamingBase, ProtCTFMicrographs):
    """ Estimates CTF on a set of micrographs using Gctf.

    To find more information about Gctf go to:
    https://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/zhang-software/#gctf
    """
    _label = 'ctf estimation'
    _devStatus = PROD
    recalculate = Boolean(False, objDoStore=False)  # Legacy Sep 2024: to fake old recalculate param
    # that is still used in the ProtCTFMicrographs (to be removed)

    def _defineCtfParamsDict(self):
        ProtCTFMicrographs._defineCtfParamsDict(self)
        self._gctfProgram = ProgramGctf(self)

    def _defineParams(self, form):
        ProgramGctf.defineInputParams(form)
        ProgramGctf.defineProcessParams(form)
        form.getParam('numberOfThreads').config(default=3)
        self._defineStreamingParams(form)

    # -------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        """Insert only the resumable streaming generator."""
        self._insertFunctionStep(
            self.resumableStepGeneratorStep,
            str(datetime.now()),
            needsGPU=False,
        )

    def resumableStepGeneratorStep(self, timestamp):
        """Run the generator as a unique step on every resume."""
        self.stepsGeneratorStep()

    def stepsGeneratorStep(self):
        """Discover, process and publish CTFs incrementally."""
        self._defineCtfParamsDict()
        self.micDict = OrderedDict()
        self._pendingMics = OrderedDict()
        self._knownMicIds = set()
        self.streamClosed = False
        self.finished = False
        self.initialIds = self._insertInitialSteps()
        self._restoreProcessedMicsFromPersistentState()

        while not self.finished:
            # A failed step makes the executor stop and then join every
            # thread, this generator included: keep polling and the run
            # hangs for good with nothing left to do.
            if self._streamingMustStop():
                break

            self._checkNewInput()
            self._checkNewOutput()

            if self.finished:
                break

            sleepOnWait = self._getStreamingSleepOnWait()
            if sleepOnWait > 0:
                self._streamingSleepOnWait()
            else:
                time.sleep(1)

    def _stepsCheck(self):
        """Persist steps created by the generator without legacy polling."""
        if getattr(self, '_newSteps', False):
            self.updateSteps()

    def _insertInitialSteps(self):
        """No filesystem initialization is required for streaming state."""
        return []

    def _loadSet(self, inputSet, SetClass, getKeyFunc):
        """Discover the items added since the last poll.

        Only ids above the watermark are queried and only those items are
        hydrated, so a poll costs what just arrived instead of everything
        the stream has produced. Items no batch has taken yet stay in
        _pendingMics and are offered again next time, since the watermark
        will never look back at them.
        """
        self.debug("Discovering new items from the logical input set.")

        newItems, producerClosed, terminalConsistent = (
            self._discoverNewInputItems(inputSet, '_lastInputId',
                                        self._knownMicIds))

        gapIds = getattr(self, '_resumeGapIds', None)

        if gapIds:
            newItems = (self._loadLogicalSetItemsByIds(inputSet, gapIds)
                        + newItems)
            self._resumeGapIds = set()

        scheduledMicNames = self._getScheduledCtfMicNames()

        for item in newItems:
            itemId = item.getObjId()

            if itemId in self._knownMicIds:
                continue

            self._knownMicIds.add(itemId)

            micKey = getKeyFunc(item)

            if micKey in self.micDict:
                continue

            if micKey in scheduledMicNames:
                # Already has a step from an earlier run: it only needs
                # publishing, so it must not be scheduled a second time.
                self.micDict[micKey] = item
                continue

            self._pendingMics[micKey] = item

        return OrderedDict(self._pendingMics), producerClosed and terminalConsistent

    def _getScheduledCtfMicNames(self):
        """Return micrograph names represented by persisted CTF steps."""
        return self._collectStepArgKeys(CTF_STEP_NAMES, onlyFinished=False,
                                        keyType=str)

    def _getFinishedCtfMicNames(self):
        """Return micrograph names represented by finished CTF steps."""
        return self._collectStepArgKeys(CTF_STEP_NAMES, keyType=str)

    def _restoreProcessedMicsFromPersistentState(self):
        """Place the watermark past what a previous run already published.

        A published CTF carries the object id of its micrograph, so the
        output Set alone says where discovery can restart - no need to walk
        the input looking for names. Micrographs that were estimated but
        never published come back as gaps and are published on the first
        output check instead of being kept out of the stream.
        """
        self._lastInputId = getattr(self, '_lastInputId', 0)

        publishedIds = self._getOutputIdSet(getattr(self, 'outputCTF', None))

        if not publishedIds:
            return

        watermark, gapIds = self._resumeWatermarkWithGaps(
            self.getInputMicrographs(), publishedIds)

        self._lastInputId = max(self._lastInputId, watermark)
        self._knownMicIds.update(publishedIds)
        self._resumeGapIds = gapIds

    def _checkNewOutput(self):
        """Publish finished CTF steps without filesystem sidecars."""
        if getattr(self, 'finished', False):
            return

        finishedNames = self._getFinishedCtfMicNames()

        # micDict holds what has been scheduled and not published yet, so
        # only that has to be looked at - never every micrograph seen.
        newDone = [
            mic for micName, mic in self.micDict.items()
            if micName in finishedNames
        ]

        allDone = (len(newDone) == len(self.micDict)
                   and not self._pendingMics)

        self.finished = self.streamClosed and allDone
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            self._updateOutputCTFSet(newDone, streamMode)

            for mic in newDone:
                self.micDict.pop(mic.getMicName(), None)
        elif not self.finished:
            if allDone:
                self._streamingSleepOnWait()
            return

        if self.finished:
            self._updateStreamState(streamMode)

    def _checkNewInput(self):
        """Discover new micrographs directly from the logical input Set."""
        micDict, self.streamClosed = self._loadInputList()
        newMics = list(micDict.values())

        if newMics:
            self._insertNewMicsSteps(newMics)

            # pwem only takes whole batches; whatever it left out has to be
            # offered again, because discovery will not find it twice.
            for mic in newMics:
                if mic.getMicName() in self.micDict:
                    self._pendingMics.pop(mic.getMicName(), None)

            self.updateSteps()

    def estimateCtfStep(self, micName, *args):
        """Estimate one CTF without filesystem completion markers."""
        mic = self.micDict[micName]
        self.info("Estimating CTF of micrograph: %s " % mic.getObjId())
        self._estimateCTF(mic, *args)

    def estimateCtfListStep(self, micNameList, *args):
        """Estimate a CTF batch without filesystem completion markers."""
        micList = [self.micDict[micName] for micName in micNameList]
        self.info("Estimating CTF for micrographs: %s"
                  % [mic.getObjId() for mic in micList])
        self._estimateCtfList(micList, *args)

    def _estimateCTF(self, mic, *args):
        self._estimateCtfList([mic], *args)

    def _estimateCtfList(self, micList, *args, **kwargs):
        """ Estimate several micrographs at once, probably a bit more
        efficient. """
        try:
            micPath = self._getMicrographDir(micList[0])
            if len(micList) > 1:
                micPath += ('-%04d' % micList[-1].getObjId())

            pwutils.makePath(micPath)
            ih = emlib.image.ImageHandler()

            for mic in micList:
                micFn = mic.getFileName()
                # We convert the input micrograph on demand if not in .mrc

                if not os.path.exists(micFn):
                    self.error("Missing input micrograph: %s. Skipping..." % micFn)
                    continue

                downFactor = self.ctfDownFactor.get()
                micFnMrc = os.path.join(micPath,
                                        self._getBatchMicBase(mic) + '.mrc')

                if downFactor != 1:
                    # Replace extension by 'mrc' cause there are some formats
                    # that cannot be written (such as dm3)
                    ih.scaleFourier(micFn, micFnMrc, downFactor)
                    sps = self._params['scannedPixelSize'] * downFactor
                    kwargs['scannedPixelSize'] = sps
                else:
                    if micFn.endswith('.mrc'):
                        pwutils.createAbsLink(os.path.abspath(micFn), micFnMrc)
                    else:
                        ih.convert(micFn, micFnMrc, emlib.DT_FLOAT)

            program, params = self._gctfProgram.getCommand(**kwargs)
            params += self._batchInputArgument(micPath)
            self.runJob(program, params, env=Plugin.getEnviron())

            def _getFile(micBase, suffix):
                return os.path.join(micPath, micBase + suffix)

            collected = 0
            attempted = 0

            for mic in micList:
                micFn = mic.getFileName()

                if not os.path.exists(micFn):
                    continue

                attempted += 1
                micBase = self._getBatchMicBase(mic)
                micFnMrc = _getFile(micBase, '.mrc')
                # Let's clean the temporary mrc micrograph
                pwutils.cleanPath(micFnMrc)

                # move output from tmp to extra
                micFnCtf = _getFile(micBase, self._gctfProgram.getExt())
                micFnCtfLog = _getFile(micBase, '_gctf.log')
                micFnCtfFit = _getFile(micBase, '_EPA.log')

                micFnCtfOut = self._getPsdPath(mic)
                micFnCtfLogOut = self._getCtfOutPath(mic)
                micFnCtfFitOut = self._getCtfFitOutPath(mic)

                try:
                    pwutils.moveFile(micFnCtf, micFnCtfOut)
                    pwutils.moveFile(micFnCtfLog, micFnCtfLogOut)
                    pwutils.moveFile(micFnCtfFit, micFnCtfFitOut)
                    collected += 1
                except Exception:
                    # One micrograph without results must not discard the
                    # ones later in the batch that do have them, and
                    # publication already skips a micrograph whose CTF
                    # cannot be built.
                    self.error("ERROR: Gctf has failed for %s" % micFn)
                    import traceback
                    traceback.print_exc()

            pwutils.cleanPath(micPath)

            if attempted and not collected:
                # Nothing at all came out of this batch: the step really did
                # fail and has to say so, rather than reporting micrographs
                # as estimated when none of them are.
                raise Exception("Gctf produced no output for any micrograph "
                                "in %s" % micPath)

        except Exception:
            self.error("ERROR: Gctf has failed for %s/*.mrc" % micPath)
            import traceback
            traceback.print_exc()
            raise

    def _createCtfModel(self, mic, updateSampling=True):
        #  When downsample option is used, we need to update the
        # sampling rate of the micrograph associated with the CTF
        # since it could be downsampled
        if updateSampling:
            newSampling = mic.getSamplingRate() * self.ctfDownFactor.get()
            mic.setSamplingRate(newSampling)

        ctf = self._gctfProgram.parseOutputAsCtf(
            self._getCtfOutPath(mic), psdFile=self._getPsdPath(mic))
        ctf.setMicrograph(mic)

        return ctf

    def _createOutputStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validateStreamingThreads(self):
        inputMics = self.getInputMicrographs()

        if inputMics is not None and inputMics.isStreamOpen() and self.numberOfThreads.get() < 3:
            return ['Gctf streaming requires at least 3 threads.']

        return []

    def _validate(self):
        errors = []
        nprocs = max(self.numberOfMpi.get(), self.numberOfThreads.get())

        if nprocs < len(self.getGpuList()):
            errors.append("Multiple GPUs can not be used by a single process. "
                          "Make sure you specify more processors than GPUs. ")

        errors.extend(self._validateStreamingThreads())
        return errors

    def _methods(self):
        if self.inputMicrographs.get() is None:
            return ['Input micrographs not available yet.']
        methods = "We calculated the CTF of "
        methods += self.getObjectTag('inputMicrographs')
        methods += " using Gctf [Zhang2016]. "
        methods += self.methodsVar.get('')

        if self.hasAttribute('outputCTF'):
            methods += 'Output CTFs: %s' % self.getObjectTag('outputCTF')

        return [methods]

    # -------------------------- UTILS functions ------------------------------
    def _getPsdPath(self, mic):
        return self._itemScopedPath(mic, self._getMicArtefactBase(mic)
                                    + '_ctf.mrc')

    def _getCtfOutPath(self, mic):
        return self._itemScopedPath(mic, self._getMicArtefactBase(mic)
                                    + '_ctf.log')

    def _getCtfFitOutPath(self, mic):
        return self._itemScopedPath(mic, self._getMicArtefactBase(mic)
                                    + '_ctf_EPA.log')

    def _getMicArtefactBase(self, mic):
        """The basename every artefact of this micrograph is built on."""
        return pwutils.removeBaseExt(mic.getFileName())

    def _batchInputArgument(self, micPath):
        """The batch folder as gctf's input argument.

        This goes to a shell, and a Scipion project can perfectly well
        live under a directory with a space in its name - unquoted, the
        shell would read that as two arguments and point gctf somewhere
        else. The folder is quoted; the wildcard stays outside the quotes
        so the shell still expands it.
        """
        return " %s/*.mrc" % shlex.quote(micPath)

    def _getBatchMicBase(self, mic):
        """Name this micrograph takes inside its batch folder.

        Every micrograph of a batch is converted into one folder and gctf
        is pointed at it with a wildcard, so two micrographs whose files
        share a basename would become a single input: one converted on
        top of the other before gctf even runs, and both reading back the
        survivor's results. The id keeps them apart; the original
        basename stays in the name so logs remain readable.
        """
        return '%06d__%s' % (mic.getObjId(), self._getMicArtefactBase(mic))

    def _parseOutput(self, filename):
        """ Try to find the output estimation parameters
        from filename. It search for lines containing: Final Values
        and Resolution limit.
        """
        return self._gctfProgram.parseOutput(filename)

    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)

        return ctf
