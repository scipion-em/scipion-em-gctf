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
from pwem import emlib
from pwem.objects import CTFModel
from pwem.protocols import ProtCTFMicrographs


from .program_gctf import ProgramGctf


class ProtGctf(ProtCTFMicrographs):
    """ Estimates CTF on a set of micrographs using Gctf.

    To find more information about Gctf go to:
    https://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/zhang-software/#gctf
    """
    _label = 'ctf estimation'

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
                downFactor = self.ctfDownFactor.get()
                micFnMrc = pwutils.join(micPath, pwutils.replaceBaseExt(micFn, 'mrc'))

                if downFactor != 1:
                    # Replace extension by 'mrc' cause there are some formats
                    # that cannot be written (such as dm3)
                    ih.scaleFourier(micFn, micFnMrc, downFactor)
                    sps = self._params['scannedPixelSize'] * downFactor
                    kwargs['scannedPixelSize'] = sps
                else:
                    ih.convert(micFn, micFnMrc, emlib.DT_FLOAT)

            program, args = self._gctfProgram.getCommand(**kwargs)
            args += ' %s/*.mrc' % micPath
            self.runJob(program, args)  # , env=gctf.Plugin.getEnviron())

            def _getFile(micBase, suffix):
                return os.path.join(micPath, micBase + suffix)

            for mic in micList:
                micFn = mic.getFileName()
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
            print("ERROR: Gctf has failed on %s/*.mrc" % micPath)
            import traceback
            traceback.print_exc()

    def _reEstimateCTF(self, mic, ctf):
        """ Re-run gctf with required parameters """
        self._estimateCtfList([mic], **self._getRecalCtfParamsDict(ctf))

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
    def _getRecalCtfParamsDict(self, ctfModel):
        values = [float(x) for x in ctfModel.getObjComment().split()]
        sampling = ctfModel.getMicrograph().getSamplingRate()
        return {
            'step_focus': 500.0,
            'lowRes': sampling / values[3],
            'highRes': sampling / values[4],
            'minDefocus': min([values[0], values[1]]),
            'maxDefocus': max([values[0], values[1]])
        }

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
