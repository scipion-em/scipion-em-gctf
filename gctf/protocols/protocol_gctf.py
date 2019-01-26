# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
import pyworkflow.em as em
import pyworkflow.protocol.params as params
from pyworkflow.utils.properties import Message
from pyworkflow.protocol import STEPS_PARALLEL

import gctf
from gctf.convert import readCtfModel, parseGctfOutput
from gctf.constants import CCC, MAXRES


class ProtGctf(em.ProtCTFMicrographs):
    """ Estimates CTF on a set of micrographs using Gctf.

    To find more information about Gctf go to:
    http://www.mrc-lmb.cam.ac.uk/kzhang
    """
    _label = 'ctf estimation'

    def __init__(self, **kwargs):
        em.ProtCTFMicrographs.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_CTF_ESTI)
        form.addParam('recalculate', params.BooleanParam, default=False,
                      condition='recalculate',
                      label="Do recalculate ctf?")

        form.addParam('continueRun', params.PointerParam, allowsNull=True,
                      condition='recalculate', label="Input previous run",
                      pointerClass=self.getClassName())
        form.addHidden('sqliteFile', params.FileParam, condition='recalculate',
                       allowsNull=True)

        form.addParam('inputMicrographs', params.PointerParam, important=True,
                      condition='not recalculate', label=Message.LABEL_INPUT_MIC,
                      pointerClass='SetOfMicrographs')
        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      condition='not recalculate',
                      help='Set to 1 for no downsampling. Non-integer '
                           'downsample factors are possible. This downsampling '
                           'is only used for estimating the CTF and it does not '
                           'affect any further calculation. Ideally the estimation '
                           'of the CTF is optimal when the Thon rings are not too '
                           'concentrated at the origin (too small to be seen) and '
                           'not occupying the whole power spectrum (since this '
                           'downsampling might entail aliasing).')

        line = form.addLine('Resolution', condition='not recalculate',
                            help='Give a value in digital frequency (i.e. between '
                                 '0.0 and 0.5). These cut-offs prevent the typical '
                                 'peak at the center of the PSD and high-resolution '
                                 'terms where only noise exists, to interfere with '
                                 'CTF estimation. The default lowest value is 0.05 '
                                 'but for micrographs with a very fine sampling this '
                                 'may be lowered towards 0. The default highest '
                                 'value is 0.35, but it should be increased for '
                                 'micrographs with signals extending beyond this '
                                 'value. However, if your micrographs extend further '
                                 'than 0.35, you should consider sampling them at a '
                                 'finer rate.')
        line.addParam('lowRes', params.FloatParam, default=0.05,
                      label='Lowest')
        line.addParam('highRes', params.FloatParam, default=0.35,
                      label='Highest')

        line = form.addLine('Defocus search range (microns)',
                            expertLevel=params.LEVEL_ADVANCED,
                            condition='not recalculate',
                            help='Select _minimum_ and _maximum_ values for '
                                 'defocus search range (in microns). '
                                 'Underfocus is represented by a positive '
                                 'number.')
        line.addParam('minDefocus', params.FloatParam, default=0.25,
                      label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=4.,
                      label='Max')

        form.addParam('astigmatism', params.FloatParam, default=100.0,
                      label='Expected (tolerated) astigmatism',
                      help='Estimated astigmatism in Angstroms',
                      expertLevel=params.LEVEL_ADVANCED)
        form.addParam('windowSize', params.IntParam, default=512,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Window size', condition='not recalculate',
                      help='The PSD is estimated from small patches of this '
                           'size. Bigger patches allow identifying more '
                           'details. However, since there are fewer windows, '
                           'estimations are noisier.')
        form.addParam('plotResRing', params.BooleanParam, default=True,
                      label='Plot a resolution ring on a PSD file',
                      help='Whether to plot an estimated resolution ring '
                           'on the power spectrum',
                      expertLevel=params.LEVEL_ADVANCED)

        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Choose GPU IDs",
                      help="GPU may have several cores. Set it to zero"
                           " if you do not know what we are talking about."
                           " First core index is 0, second 1 and so on."
                           " You can use multiple GPUs - in that case"
                           " set to i.e. *0 1 2*.")

        form.addSection(label='Advanced')
        group = form.addGroup('EPA')
        group.addParam('doEPA', params.BooleanParam, default=False,
                       label="Do EPA",
                       help='Do Equiphase average used for output CTF file. '
                            'Only for nice output, will NOT be used for CTF '
                            'determination.')

        group.addParam('EPAsmp', params.IntParam, default=4,
                       condition='doEPA',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Over-sampling factor for EPA")
        group.addParam('doBasicRotave', params.BooleanParam, default=False,
                       condition='doEPA',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Do rotational average",
                       help='Do rotational average used for output CTF file. '
                            'Only for nice output, will NOT be used for CTF '
                            'determination.')

        group.addParam('overlap', params.FloatParam, default=0.5,
                       condition='doEPA',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Overlap factor",
                       help='Overlapping factor for grid boxes sampling, '
                            'for windowsize=512, 0.5 means 256 pixels overlapping.')
        group.addParam('convsize', params.IntParam, default=85,
                       condition='doEPA',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Boxsize for smoothing",
                       help='Boxsize to be used for smoothing, '
                            'suggested 1/5 ~ 1/20 of window size in pixel, '
                            'e.g. 99 for 512 window')

        if self._isVersion118():
            group.addParam('smoothResL', params.IntParam, default=1000,
                           expertLevel=params.LEVEL_ADVANCED,
                           condition='doEPA',
                           label='Resolution for smoothing',
                           help='Provide a reasonable resolution for low '
                                'frequency background smoothing; 20 '
                                'angstrom suggested, 10-50 is proper range')

        form.addParam('bfactor', params.IntParam, default=150,
                      label="B-factor",
                      help='B-factors used to decrease high resolution '
                           'amplitude, A^2; suggested range 50~300 except '
                           'using REBS method  (see the paper for the details).')

        group = form.addGroup('High-res refinement')
        group.addParam('doHighRes', params.BooleanParam, default=False,
                       label="Do high-resolution refinement",
                       help='Whether to do High-resolution refinement or not, '
                            'very useful for selecting high quality micrographs. '
                            'Especially useful when your data has strong '
                            'low-resolution bias')
        group.addParam('HighResL', params.FloatParam, default=15.0,
                       condition='doHighRes',
                       label="Lowest resolution",
                       help='Lowest resolution  to be used for High-resolution '
                            'refinement, in Angstroms')
        group.addParam('HighResH', params.FloatParam, default=4.0,
                       condition='doHighRes',
                       label="Highest resolution",
                       help='Highest resolution  to be used for High-resolution '
                            'refinement, in Angstroms')
        group.addParam('HighResBf', params.IntParam, default=50,
                       condition='doHighRes',
                       label="B-factor",
                       help='B-factor to be used for High-resolution '
                            'refinement, in Angstroms')

        form.addParam('doValidate', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Do validation",
                      help='Whether to validate the CTF determination.')

        form.addSection(label='Phase shift')
        form.addParam('doPhShEst', params.BooleanParam, default=False,
                      label="Estimate phase shift?",
                      help='For micrographs collected with phase-plate. '
                           'It is suggested to import such micrographs with '
                           'amplitude contrast = 0. Also, using smaller '
                           '_lowest resolution_ (e.g. 15A) and smaller '
                           '_boxsize for smoothing_ (e.g. 50 for 1024 '
                           'window size) might be better.')

        line = form.addLine('Phase shift range range (deg)',
                            condition='doPhShEst',
                            help='Select _lowest_ and _highest_ phase shift '
                                 '(in degrees).')
        line.addParam('phaseShiftL', params.FloatParam, default=0.0,
                      condition='doPhShEst',
                      label="Min")
        line.addParam('phaseShiftH', params.FloatParam, default=180.0,
                      condition='doPhShEst',
                      label="Max")

        form.addParam('phaseShiftS', params.FloatParam, default=10.0,
                      condition='doPhShEst',
                      label="Step",
                      help='Phase shift search step. Do not worry about '
                           'the accuracy; this is just the search step, '
                           'Gctf will refine the phase shift anyway.')
        form.addParam('phaseShiftT', params.EnumParam, default=CCC,
                      condition='doPhShEst',
                      label='Target',
                      choices=['CCC', 'Resolution limit'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Phase shift target in the search: CCC or '
                            'resolution limit. Second option might generate '
                            'more accurate estimation if results are '
                            'essentially correct, but it tends to overfit high '
                            'resolution noise and might have the potential '
                            'possibility to generate completely wrong results. '
                            'The accuracy of CCC method might not be as '
                            'good, but it is more stable in general cases.')

        if self._isVersion118():
             form.addParam('coSearchRefine', params.BooleanParam,
                           default=False, condition='doPhShEst',
                           label='Search and refine simultaneously?',
                           help='Specify this option to do refinement during '
                                'phase shift search. Default approach is to do '
                                'refinement after search.')
             form.addParam('refine2DT', params.IntParam,
                           validators=[params.Range(1, 3, "value should be "
                                                          "1, 2 or 3. ")],
                           default=1, condition='doPhShEst',
                           label='Refinement type',
                           help='Refinement type: 1, 2, 3 allowed.\n NOTE:  '
                                'This parameter is different from Target and is'
                                'optional for different types of refinement algorithm, '
                                'in general cases they work similar. In challenging '
                                'case, they might converge to different results, '
                                'try to see which works best in your case. '
                                'My suggestion is running as default first, and '
                                'then try new refinement on the micrographs '
                                'which failed.')

        self._defineStreamingParams(form)

        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- STEPS functions ------------------------------
    def _estimateCTF(self, mic, *args):
        self._estimateCtfList([mic], *args)

    def _getExt(self):
        if self._isVersion118():
            ext = '.pow' if not self.doEPA else '.epa'
        else:
            ext = '.ctf'
        return ext

    def _estimateCtfList(self, micList, *args):
        """ Estimate several micrographs at once, probably a bit more
        efficient. """
        micPath = self._getMicrographDir(micList[0])
        if len(micList) > 1:
            micPath += ('-%04d' % micList[-1].getObjId())

        pwutils.makePath(micPath)
        ih = em.ImageHandler()

        def _getFile(micBase, suffix):
            return os.path.join(micPath, micBase + suffix)

        for mic in micList:
            micFn = mic.getFileName()
            # We convert the input micrograph on demand if not in .mrc
            downFactor = self.ctfDownFactor.get()
            micFnMrc = pwutils.join(micPath, pwutils.replaceBaseExt(micFn, 'mrc'))

            if downFactor != 1:
                # Replace extension by 'mrc' cause there are some formats
                # that cannot be written (such as dm3)
                ih.scaleFourier(micFn, micFnMrc, downFactor)
                sps = self.inputMicrographs.get().getScannedPixelSize() * downFactor
                self._params['scannedPixelSize'] = sps
            else:
                ih.convert(micFn, micFnMrc, em.DT_FLOAT)

        try:
            args = self._args % self._params
            args += ' %s/*.mrc' % micPath
            self.runJob(gctf.Plugin.getProgram(), args,
                        env=gctf.Plugin.getEnviron())

            for mic in micList:
                micFn = mic.getFileName()
                micBase = pwutils.removeBaseExt(micFn)
                micFnMrc = _getFile(micBase, '.mrc')
                # Let's clean the temporary mrc micrograph
                pwutils.cleanPath(micFnMrc)

                # move output from tmp to extra
                micFnCtf = _getFile(micBase, self._getExt())
                micFnCtfLog = _getFile(micBase, '_gctf.log')
                micFnCtfFit = _getFile(micBase, '_EPA.log')

                micFnCtfOut = self._getPsdPath(micFn)
                micFnCtfLogOut = self._getCtfOutPath(micFn)
                micFnCtfFitOut = self._getCtfFitOutPath(micFn)

                pwutils.moveFile(micFnCtf, micFnCtfOut)
                pwutils.moveFile(micFnCtfLog, micFnCtfLogOut)
                pwutils.moveFile(micFnCtfFit, micFnCtfFitOut)

        except:
            print("ERROR: Gctf has failed on %s/*.mrc" % micPath)
            import traceback
            traceback.print_exc()

    def _restimateCTF(self, ctfId):
        ih = em.ImageHandler()
        ctfModel = self.recalculateSet[ctfId]
        mic = ctfModel.getMicrograph()
        micFn = mic.getFileName()
        micBase = pwutils.removeBaseExt(micFn)
        micFnMrc = self._getTmpPath(micBase + '.mrc')

        # We convert the input micrograph on demand if not in .mrc
        downFactor = self.ctfDownFactor.get()

        if downFactor != 1:
            # Replace extension by 'mrc' cause there are some formats
            # that cannot be written (such as dm3)
            ih.scaleFourier(micFn, micFnMrc, downFactor)
            sps = self.inputMicrographs.get().getScannedPixelSize() * downFactor
            self._params['scannedPixelSize'] = sps
        else:
            ih.convert(micFn, micFnMrc, em.DT_FLOAT)

        # Update _params dictionary
        self._prepareRecalCommand(ctfModel)

        try:
            args = self._args % self._params
            args += ' %s' % micFnMrc
            self.runJob(gctf.Plugin.getProgram(), args,
                        env=gctf.Plugin.getEnviron())
        except:
            print("ERROR: Gctf has failed for micrograph %s" % micFnMrc)
            import traceback
            traceback.print_exc()

        # Let's clean the temporary mrc micrograph
        pwutils.cleanPath(micFnMrc)

        # move output from tmp to extra

        micFnCtf = self._getTmpPath(micBase + self._getExt())
        micFnCtfLog = self._getTmpPath(micBase + '_gctf.log')
        micFnCtfFit = self._getTmpPath(micBase + '_EPA.log')

        micFnCtfOut = self._getPsdPath(micFn)
        micFnCtfLogOut = self._getCtfOutPath(micFn)
        micFnCtfFitOut = self._getCtfFitOutPath(micFn)

        pwutils.moveFile(micFnCtf, micFnCtfOut)
        pwutils.moveFile(micFnCtfLog, micFnCtfLogOut)
        pwutils.moveFile(micFnCtfFit, micFnCtfFitOut)

    def _createCtfModel(self, mic, updateSampling=True):
        #  When downsample option is used, we need to update the
        # sampling rate of the micrograph associated with the CTF
        # since it could be downsampled
        if updateSampling:
            newSampling = mic.getSamplingRate() * self.ctfDownFactor.get()
            mic.setSamplingRate(newSampling)

        micFn = mic.getFileName()
        out = self._getCtfOutPath(micFn)
        psdFile = self._getPsdPath(micFn)

        ctfModel2 = em.CTFModel()
        readCtfModel(ctfModel2, out)
        ctfModel2.setPsdFile(psdFile)
        ctfModel2.setMicrograph(mic)

        return ctfModel2

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
    def _prepareCommand(self):
        sampling = self.inputMics.getSamplingRate() * self.ctfDownFactor.get()
        # Convert digital frequencies to spatial frequencies
        self._params['sampling'] = sampling
        self._params['lowRes'] = sampling / self._params['lowRes']
        if self._params['lowRes'] > 50:
            self._params['lowRes'] = 50
        self._params['highRes'] = sampling / self._params['highRes']
        self._params['step_focus'] = 500.0
        self._argsGctf()

    def _prepareRecalCommand(self, ctfModel):
        line = ctfModel.getObjComment().split()
        self._defineRecalValues(ctfModel)

        # get the size and the image of psd
        imgPsd = ctfModel.getPsdFile()
        imgh = em.ImageHandler()
        size, _, _, _ = imgh.getDimensions(imgPsd)

        mic = ctfModel.getMicrograph()

        # Convert digital frequencies to spatial frequencies
        sampling = mic.getSamplingRate()
        self._params['step_focus'] = 1000.0
        self._params['sampling'] = sampling
        self._params['lowRes'] = sampling / float(line[3])
        self._params['highRes'] = sampling / float(line[4])
        self._params['minDefocus'] = min([float(line[0]), float(line[1])])
        self._params['maxDefocus'] = max([float(line[0]), float(line[1])])
        self._params['windowSize'] = size

        self._argsGctf()

    def _argsGctf(self):
        self._args = " --apix %f " % self._params['sampling']
        self._args += "--kV %f " % self._params['voltage']
        self._args += "--cs %f " % self._params['sphericalAberration']
        self._args += "--ac %f " % self._params['ampContrast']
        self._args += "--dstep %f " % self._params['scannedPixelSize']
        self._args += "--defL %f " % self._params['minDefocus']
        self._args += "--defH %f " % self._params['maxDefocus']
        self._args += "--defS %f " % self._params['step_focus']
        self._args += "--astm %f " % self.astigmatism.get()
        self._args += "--resL %f " % self._params['lowRes']
        self._args += "--resH %f " % self._params['highRes']
        self._args += "--do_EPA %d " % (1 if self.doEPA else 0)
        self._args += "--boxsize %d " % self._params['windowSize']
        self._args += "--plot_res_ring %d " % (1 if self.plotResRing else 0)
        self._args += "--gid %%(GPU)s "  # Use %% to escape when formatting
        self._args += "--bfac %d " % self.bfactor.get()
        self._args += "--B_resH %f " % (2 * self._params['sampling'])
        self._args += "--overlap %f " % self.overlap.get()
        self._args += "--convsize %d " % self.convsize.get()
        self._args += "--do_Hres_ref %d " % (1 if self.doHighRes else 0)

        if self._isVersion118():
            self._args += "--smooth_resL %d " % self.smoothResL.get()

        self._args += "--EPA_oversmp %d " % self.EPAsmp.get()

        if self.doPhShEst:
            self._args += "--phase_shift_L %f " % self.phaseShiftL.get()
            self._args += "--phase_shift_H %f " % self.phaseShiftH.get()
            self._args += "--phase_shift_S %f " % self.phaseShiftS.get()
            self._args += "--phase_shift_T %d " % (1 + self.phaseShiftT.get())

            if self._isVersion118():
                self._args += "--cosearch_refine_ps %d " % (1 if self.coSearchRefine else 0)
                self._args += "--refine_2d_T %d " % self.refine2DT.get()

        if self.doHighRes:
            self._args += "--Href_resL %d " % self.HighResL.get()
            self._args += "--Href_resH %d " % self.HighResH.get()
            self._args += "--Href_bfac %d " % self.HighResBf.get()

        self._args += "--ctfstar NONE --do_validation %d " % (1 if self.doValidate else 0)

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
        return parseGctfOutput(filename)

    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = em.CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)

        return ctf

    def _isVersion118(self):
        # specific case for v1.18
        return gctf.Plugin.getActiveVersion() in ['1.18']
