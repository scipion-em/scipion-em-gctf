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
import time
import json
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


class ProtGctf(ProtCTFMicrographs):
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
        self.streamClosed = False
        self.finished = False
        self.initialIds = self._insertInitialSteps()
        self._restoreProcessedMicsFromPersistentState()

        while not self.finished:
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
        """Load new items directly from the logical input Set."""
        self.debug("Loading logical input set.")
        inputSet.loadAllProperties()

        newItemDict = {}
        for item in inputSet.iterItems():
            micKey = getKeyFunc(item)
            if micKey not in self.micDict:
                newItemDict[micKey] = item.clone()

        streamClosed = inputSet.isStreamClosed()
        return newItemDict, streamClosed

    def _getPublishedCtfMicNames(self):
        """Return micrographs already present in the logical CTF output."""
        outputCtf = getattr(self, 'outputCTF', None)

        if outputCtf is None:
            return set()

        iterator = (
            outputCtf.iterItems()
            if hasattr(outputCtf, 'iterItems')
            else iter(outputCtf)
        )

        publishedNames = set()

        for ctf in iterator:
            mic = ctf.getMicrograph()

            if mic is not None:
                publishedNames.add(mic.getMicName())

        return publishedNames

    def _getScheduledCtfMicNames(self):
        """Return micrograph names represented by persisted CTF steps."""
        scheduledNames = set()

        for step in getattr(self, '_steps', []):
            funcName = getattr(step, 'funcName', None)
            if hasattr(funcName, 'get'):
                funcName = funcName.get()

            if funcName not in ('estimateCtfStep', 'estimateCtfListStep'):
                continue

            argsStr = getattr(step, 'argsStr', None)
            if hasattr(argsStr, 'get'):
                argsStr = argsStr.get('[]')

            try:
                args = json.loads(argsStr or '[]')
            except (TypeError, ValueError):
                continue

            if not args:
                continue

            if funcName == 'estimateCtfStep':
                if isinstance(args[0], str):
                    scheduledNames.add(args[0])
            elif isinstance(args[0], list):
                scheduledNames.update(
                    micName for micName in args[0]
                    if isinstance(micName, str)
                )

        return scheduledNames

    def _getFinishedCtfMicNames(self):
        """Return micrograph names represented by finished CTF steps."""
        finishedNames = set()

        for step in getattr(self, '_steps', []):
            if not step.isFinished():
                continue

            funcName = getattr(step, 'funcName', None)
            if hasattr(funcName, 'get'):
                funcName = funcName.get()

            if funcName not in ('estimateCtfStep', 'estimateCtfListStep'):
                continue

            argsStr = getattr(step, 'argsStr', None)
            if hasattr(argsStr, 'get'):
                argsStr = argsStr.get('[]')

            try:
                args = json.loads(argsStr or '[]')
            except (TypeError, ValueError):
                continue

            if not args:
                continue

            if funcName == 'estimateCtfStep':
                if isinstance(args[0], str):
                    finishedNames.add(args[0])
            elif isinstance(args[0], list):
                finishedNames.update(
                    micName for micName in args[0]
                    if isinstance(micName, str)
                )

        return finishedNames

    def _restoreProcessedMicsFromPersistentState(self):
        """Restore published and already scheduled inputs before discovery."""
        processedNames = (
            self._getPublishedCtfMicNames()
            | self._getScheduledCtfMicNames()
        )

        if not processedNames:
            return

        inputMics = self.getInputMicrographs()
        inputMics.loadAllProperties()

        for mic in inputMics.iterItems():
            micName = mic.getMicName()
            if micName in processedNames:
                self.micDict[micName] = mic.clone()

    def _checkNewOutput(self):
        """Publish finished CTF steps without filesystem sidecars."""
        if getattr(self, 'finished', False):
            return

        publishedNames = self._getPublishedCtfMicNames()
        finishedNames = self._getFinishedCtfMicNames()
        micList = list(self.micDict.values())

        newDone = [
            mic for mic in micList
            if (
                mic.getMicName() in finishedNames
                and mic.getMicName() not in publishedNames
            )
        ]

        completedNames = publishedNames | finishedNames
        allDone = all(
            mic.getMicName() in completedNames
            for mic in micList
        )

        self.finished = self.streamClosed and allDone
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            self._updateOutputCTFSet(newDone, streamMode)
        elif not self.finished:
            if allDone:
                self._streamingSleepOnWait()
            return

        if self.finished:
            self._updateStreamState(streamMode)
            outputStep = self._getFirstJoinStep()

            if outputStep and outputStep.isWaiting():
                from pyworkflow.protocol.constants import STATUS_NEW
                outputStep.setStatus(STATUS_NEW)

    def _checkNewInput(self):
        """Discover new micrographs directly from the logical input Set."""
        micDict, self.streamClosed = self._loadInputList()
        newMics = micDict.values()
        outputStep = self._getFirstJoinStep()

        if newMics:
            dependencies = self._insertNewMicsSteps(newMics)
            if outputStep is not None:
                outputStep.addPrerequisites(*dependencies)
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
                micFnMrc = os.path.join(micPath, pwutils.replaceBaseExt(micFn, 'mrc'))

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
            params += ' %s/*.mrc' % micPath
            self.runJob(program, params, env=Plugin.getEnviron())

            def _getFile(micBase, suffix):
                return os.path.join(micPath, micBase + suffix)

            for mic in micList:
                micFn = mic.getFileName()

                if not os.path.exists(micFn):
                    continue

                micBase = pwutils.removeBaseExt(micFn)
                micFnMrc = _getFile(micBase, '.mrc')
                # Let's clean the temporary mrc micrograph
                pwutils.cleanPath(micFnMrc)

                # move output from tmp to extra
                micFnCtf = _getFile(micBase, self._gctfProgram.getExt())
                micFnCtfLog = _getFile(micBase, '_gctf.log')
                micFnCtfFit = _getFile(micBase, '_EPA.log')

                micFnCtfOut = self._getPsdPath(micFn)
                micFnCtfLogOut = self._getCtfOutPath(micFn)
                micFnCtfFitOut = self._getCtfFitOutPath(micFn)

                try:
                    pwutils.moveFile(micFnCtf, micFnCtfOut)
                    pwutils.moveFile(micFnCtfLog, micFnCtfLogOut)
                    pwutils.moveFile(micFnCtfFit, micFnCtfFitOut)
                except Exception:
                    self.error("ERROR: Gctf has failed for %s" % micFn)
                    import traceback
                    traceback.print_exc()

            pwutils.cleanPath(micPath)

        except:
            self.error("ERROR: Gctf has failed for %s/*.mrc" % micPath)
            import traceback
            traceback.print_exc()

    def _createCtfModel(self, mic, updateSampling=True):
        #  When downsample option is used, we need to update the
        # sampling rate of the micrograph associated with the CTF
        # since it could be downsampled
        if updateSampling:
            newSampling = mic.getSamplingRate() * self.ctfDownFactor.get()
            mic.setSamplingRate(newSampling)

        micFn = mic.getFileName()
        ctf = self._gctfProgram.parseOutputAsCtf(
            self._getCtfOutPath(micFn), psdFile=self._getPsdPath(micFn))
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
    def _getPsdPath(self, micFn):
        micFnBase = pwutils.removeBaseExt(micFn)
        return self._getExtraPath(micFnBase + '_ctf.mrc')

    def _getCtfOutPath(self, micFn):
        micFnBase = pwutils.removeBaseExt(micFn)
        return self._getExtraPath(micFnBase + '_ctf.log')

    def _getCtfFitOutPath(self, micFn):
        micFnBase = pwutils.removeBaseExt(micFn)
        return self._getExtraPath(micFnBase + '_ctf_EPA.log')

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
