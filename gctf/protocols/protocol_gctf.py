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


def isUsableCtf(ctfModel):
    """False for the -999 placeholder of a failed or missing result."""
    defocusU, defocusV = ctfModel.getDefocusU(), ctfModel.getDefocusV()

    return (defocusU is not None and defocusV is not None
            and defocusU >= 0 and defocusV >= 0)


class ProtGctf(GctfStreamingBase, ProtCTFMicrographs):
    """ Estimates CTF on a set of micrographs using Gctf.

    To find more information about Gctf go to:
    https://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/zhang-software/#gctf

    AI Generated:

    Gctf CTF Estimation (ProtGctf) — User Manual
        Overview

        The Gctf protocol performs Contrast Transfer Function estimation on
        cryo-EM micrographs using the GPU-accelerated Gctf software package.
        Its primary purpose is to characterize the optical effects introduced
        by the electron microscope during image acquisition so that these
        effects can later be corrected or compensated during reconstruction
        and refinement procedures.

        In cryo-EM workflows, accurate CTF estimation is one of the earliest
        and most biologically important preprocessing steps. The protocol
        determines parameters such as defocus values, astigmatism, and power
        spectrum characteristics from each micrograph. Reliable estimation is
        essential because downstream operations including particle alignment,
        three-dimensional reconstruction, and high-resolution refinement all
        depend strongly on the quality of the CTF model.

        Inputs and General Workflow

        The protocol accepts a set of micrographs as input. These micrographs
        may originate from single-particle cryo-EM acquisitions collected
        under a wide variety of imaging conditions. The protocol processes
        each image independently and produces a corresponding CTF estimation
        together with diagnostic outputs useful for visual inspection and
        validation.

        During execution, the protocol prepares the micrographs for analysis,
        applies optional downsampling if requested, and estimates the CTF
        parameters using Gctf. The resulting outputs are organized so they
        can be directly consumed by later processing stages inside Scipion
        workflows.

        Downsampling and Frequency Considerations

        One of the key configurable options is the CTF downsampling factor.
        Downsampling reduces the effective image size used during estimation,
        which can significantly accelerate processing and improve visibility
        of the Thon rings under certain conditions.

        From a practical biological perspective, moderate downsampling is
        often beneficial when working with very large micrographs or when the
        Thon rings are overly concentrated near the origin of the power
        spectrum. However, excessive downsampling may reduce the visibility
        of high-frequency information and can introduce aliasing effects that
        negatively affect estimation accuracy.

        In most routine workflows, users begin with no downsampling and only
        introduce moderate factors when computational efficiency or ring
        visibility becomes problematic.

        GPU Acceleration and High-Throughput Processing

        Gctf is optimized for GPU execution, making this protocol especially
        suitable for modern cryo-EM facilities and automated processing
        pipelines. The protocol can process multiple micrographs efficiently
        and is compatible with streaming-oriented workflows where data are
        analyzed during acquisition.

        This capability is particularly important in high-throughput
        environments where rapid feedback about microscope performance and
        dataset quality is required. Early detection of imaging problems can
        prevent unnecessary collection of unusable data and improve overall
        experimental efficiency.

        Interpretation of the Outputs

        The protocol produces CTF models associated with each input
        micrograph together with diagnostic power spectrum files and fitting
        information. These outputs allow users to evaluate whether the
        estimated defocus and astigmatism values are physically reasonable
        and whether the Thon ring fitting quality is sufficient for
        high-resolution analysis.

        In biological practice, users typically inspect the visibility and
        continuity of the fitted Thon rings as indicators of image quality.
        Strong and extended rings usually suggest good preservation of
        high-resolution signal, whereas weak or poorly fitted rings may
        indicate problems such as drift, charging, contamination, excessive
        ice thickness, or low signal-to-noise ratio.

        The estimated parameters also provide important information about the
        acquisition conditions across the dataset. Large defocus variations
        or systematic astigmatism trends may reveal microscope alignment
        issues or acquisition inconsistencies that should be considered
        before further analysis.

        Practical Recommendations

        For most datasets, the default Gctf settings provide a good starting
        point. Users should initially verify the estimated defocus values and
        visually inspect a representative subset of power spectra before
        processing the entire dataset.

        Moderate downsampling may improve performance for large datasets, but
        users interested in the highest possible resolution should verify
        that high-frequency information remains well fitted. GPU resources
        should also be configured carefully to ensure stable execution and
        efficient workload distribution.

        When processing heterogeneous datasets acquired over multiple
        sessions, it is advisable to monitor trends in the estimated CTF
        parameters because abrupt changes may indicate acquisition problems
        or microscope instability.

        Final Perspective

        CTF estimation is not simply a technical preprocessing step but a
        fundamental component of reliable cryo-EM analysis. Accurate optical
        characterization directly influences the interpretability and
        attainable resolution of the final reconstruction. Careful validation
        of the estimated parameters, thoughtful use of downsampling, and
        consistent quality monitoring are essential practices for obtaining
        biologically meaningful results in cryo-EM workflows.
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
                if self._ctfResultIsUsable(item):
                    # Already estimated by an earlier run: it only needs
                    # publishing, so it must not be scheduled a second time.
                    self.micDict[micKey] = item
                    continue

                # Its earlier estimation failed: estimate it again.
                self._retriedMicNames = getattr(
                    self, '_retriedMicNames', set())
                self._retriedMicNames.add(micKey)

            self._pendingMics[micKey] = item

        return OrderedDict(self._pendingMics), producerClosed and terminalConsistent

    def _getScheduledCtfMicNames(self):
        """Return micrograph names represented by persisted CTF steps."""
        return self._collectStepArgKeys(CTF_STEP_NAMES, onlyFinished=False,
                                        keyType=str)

    def _hasActiveStreamingWork(self):
        """Micrographs scheduled or waiting for a batch, not yet published."""
        return bool(getattr(self, 'micDict', None)
                    or getattr(self, '_pendingMics', None))

    def _getFinishedCtfMicNames(self):
        """Return micrograph names represented by finished CTF steps.

        A micrograph estimated again in this run (its earlier result was
        unusable) only counts once its new step has finished: its old step
        finished too, and publishing on that would drop it once more.
        """
        finished = self._collectStepArgKeys(CTF_STEP_NAMES, keyType=str)
        retried = getattr(self, '_retriedMicNames', None)

        if retried:
            finishedNow = self._collectStepArgKeys(CTF_STEP_NAMES, keyType=str,
                                                   includePrevious=False)
            finished = {name for name in finished
                        if name not in retried or name in finishedNow}

        return finished

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
            self._failIfNoCtfWasEstimated()

    def _failIfNoCtfWasEstimated(self):
        """A failed micrograph is left out, not published as -999; the run
        fails only when that happened to every micrograph."""
        if not getattr(self, '_knownMicIds', None):
            return

        output = getattr(self, 'outputCTF', None)

        if output is None or output.getSize() == 0:
            raise Exception("No CTF could be estimated for any of the %d "
                            "micrographs; see the errors above. Continue "
                            "estimates them again."
                            % len(self._knownMicIds))

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
                self.error("ERROR: Gctf produced no output for any "
                           "micrograph in %s" % micPath)

        except Exception:
            # One failed batch must not stop the stream: its micrographs
            # have no usable result, so they are left out of the output
            # and estimated again on Continue.
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

        ctf = self._gctfProgram.parseOutputAsCtf(
            self._getCtfOutPath(mic), psdFile=self._getPsdPath(mic))

        # A missing or unusable result raises instead of becoming a -999
        # CTF: pwem then leaves the micrograph out and lists it as failed.
        if not isUsableCtf(ctf):
            raise ValueError("No usable CTF estimation for micrograph %s"
                             % mic.getObjId())

        ctf.setMicrograph(mic)

        return ctf

    def _ctfResultIsUsable(self, mic):
        """Whether an earlier run left a usable estimation for ``mic``."""
        try:
            self._createCtfModel(mic.clone(), updateSampling=False)
        except Exception:
            return False

        return True

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
