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

import pyworkflow.utils as pwutils
from pyworkflow.object import Boolean
from pyworkflow.constants import PROD
from pwem import emlib
from pwem.objects import CTFModel
from pwem.protocols import ProtCTFMicrographs


from gctf import Plugin
from gctf.protocols.program_gctf import ProgramGctf


class ProtGctf(ProtCTFMicrographs):
    """
    Estimates the Contrast Transfer Function (CTF) of cryo-EM micrographs
    using the Gctf package. The protocol is designed to determine optical
    parameters such as defocus and astigmatism from electron microscopy
    images, enabling downstream processing steps including particle picking,
    refinement, reconstruction, and quality assessment. More info:
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
        self._defineStreamingParams(form)

    # -------------------------- STEPS functions ------------------------------
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

                pwutils.moveFile(micFnCtf, micFnCtfOut)
                pwutils.moveFile(micFnCtfLog, micFnCtfLogOut)
                pwutils.moveFile(micFnCtfFit, micFnCtfFitOut)

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
    def _validate(self):
        errors = []
        nprocs = max(self.numberOfMpi.get(), self.numberOfThreads.get())

        if nprocs < len(self.getGpuList()):
            errors.append("Multiple GPUs can not be used by a single process. "
                          "Make sure you specify more processors than GPUs. ")

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
