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

from collections import OrderedDict

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.constants import RELATION_CTF
from pwem import emlib
import pwem.emlib.metadata as md
from pwem.protocols import EMProtocol, ProtParticles

from .. import Plugin
from ..convert import CoordinatesWriter, rowToCtfModel, getShifts
from ..constants import *


class ProtGctfRefine(ProtParticles):
    """
    Refines local CTF of a set of particles using Gctf.

    To find more information about Gctf go to:
    https://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/zhang-software/#gctf
    """

    _label = 'ctf refinement'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self._params = {}
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      important=True,
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles for local CTF refinement.')
        form.addParam('applyShifts', params.BooleanParam, default=False,
                      label='Apply particle shifts?',
                      help='Apply particle shifts from 2D alignment to '
                           'recalculate new coordinates. This can be useful '
                           'for re-centering particle coordinates.')
        form.addParam('inputMicrographs', params.PointerParam, important=True,
                      label='Input micrographs',
                      pointerClass='SetOfMicrographs',
                      help='Select the SetOfMicrographs related to input particles.')
        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      help='Set to 1 for no downsampling. Non-integer '
                           'downsample factors are possible. This downsampling '
                           'is only used for estimating the CTF and it does not '
                           'affect any further calculation. Ideally the estimation '
                           'of the CTF is optimal when the Thon rings are not too '
                           'concentrated at the origin (too small to be seen) and '
                           'not occupying the whole power spectrum (since this '
                           'downsampling might entail aliasing).')

        form.addParam('windowSize', params.IntParam, default=1024,
                      label='Box size (px)',
                      help='Boxsize in pixels to be used for FFT, 512 or '
                           '1024 highly recommended')

        group = form.addGroup('Search limits')
        line = group.addLine('Resolution (A)',
                             help='The CTF model will be fit to regions '
                                  'of the amplitude spectrum corresponding '
                                  'to this range of resolution.')
        line.addParam('lowRes', params.FloatParam, default=50., label='Min')
        line.addParam('highRes', params.FloatParam, default=4., label='Max')

        line = group.addLine('Defocus search range (A)',
                             help='Select _minimum_ and _maximum_ values for '
                                  'defocus search range (in A). Underfocus'
                                  ' is represented by a positive number.')
        line.addParam('minDefocus', params.FloatParam, default=5000.,
                      label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=90000.,
                      label='Max')
        group.addParam('stepDefocus', params.FloatParam, default=500.,
                       label='Defocus step (A)',
                       help='Step size for the defocus search.')

        form.addParam('astigmatism', params.FloatParam, default=1000.0,
                      label='Expected (tolerated) astigmatism',
                      help='Estimated astigmatism in Angstroms',
                      expertLevel=params.LEVEL_ADVANCED)

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

        form.addSection(label='Local refinement')
        line = form.addLine('Local resolution (A)',
                            help='Select _lowest_ and _highest_ resolution '
                                 'to be used for local CTF (in Angstrom).')
        line.addParam('locResL', params.IntParam, default=15, label='Low')
        line.addParam('locResH', params.IntParam, default=5, label='High')

        form.addParam('locRad', params.IntParam, default=1024,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Radius for local refinement (px)',
                      help='Radius for local refinement, no weighting '
                           'if the distance is larger than that')
        form.addParam('locAveType', params.EnumParam,
                      default=WEIGHT_BOTH,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Local average type',
                      choices=['Equal weights', 'Distance', 'Distance and freq'],
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='_Equal weights_: equal weights for all local '
                           'areas, neither distance nor frequency is '
                           'weighted\n_Distance_: single weight for each '
                           'local area, only distance is weighted\n'
                           '_Distance and freq_: Guassian weighting for '
                           'both distance and frequency')
        form.addParam('locBoxSize', params.IntParam, default=512,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Boxsize (px)',
                      help='Boxsize for local refinement (in pixels)')
        form.addParam('locOverlap', params.FloatParam, default=0.5,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Overlap',
                      help='Overlapping factor for grid boxes sampling')
        form.addParam('locAstm', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Refine astigmatism?',
                      help='By default (False) only refine Z-height '
                           'changes in local area (suggested). If True, '
                           'refine local astigmatism (not suggested unless '
                           'SNR is very good).')

        form.addSection(label='CTF refinement')
        form.addParam('useInputCtf', params.BooleanParam, default=False,
                      label="Refine input CTFs",
                      help='Input CTF will be taken from input micrographs. '
                           'By default Gctf wil NOT refine user-provided '
                           'CTF parameters but do ab initial determination.')
        form.addParam('ctfRelations', params.RelationParam, allowsNull=True,
                      condition='useInputCtf',
                      relationName=RELATION_CTF,
                      attributeName='_getMicrographs',
                      label='Input CTF estimation',
                      help='Choose some CTF estimation related to input '
                           'micrographs.')
        form.addParam('defUerr', params.FloatParam, default=500.0,
                      condition='useInputCtf',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='DefocusU error (nm)',
                      help='Estimated error of input initial defocus_U.')
        form.addParam('defVerr', params.FloatParam, default=500.0,
                      condition='useInputCtf',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='DefocusV error (nm)',
                      help='Estimated error of input initial defocus_V.')
        form.addParam('defAerr', params.FloatParam, default=15.0,
                      condition='useInputCtf',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Defocus angle error',
                      help='Estimated error of input initial defocus angle.')
        form.addParam('Berr', params.FloatParam, default=50.0,
                      condition='useInputCtf',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='B-factor error',
                      help='Estimated error of input initial B-factor.')

        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- STEPS functions -------------------------------
    def _createMicDict(self):
        """ Create a dictionary with all micrographs that
         are both in the input micrographs set and there
         are particles belonging to it.
         micName will be the key to that dict.
         """
        inputParticles = self.inputParticles.get()
        firstCoord = inputParticles.getFirstItem().getCoordinate()
        self.hasMicName = firstCoord.getMicName() is not None
        inputMicDict = {mic.getMicName(): mic.clone()
                        for mic in self._getMicrographs()}
        # Check now which if these mics have particles belonging
        self.micDict = OrderedDict()
        # match the mic from coord with micDict
        lastMicId = None
        # TODO: If this loop is too expensive for very large input datasets,
        # we could consider using the aggregate functions in the mapper
        for particle in inputParticles.iterItems(orderBy='_micId'):
            micId = particle.getMicId()
            if micId != lastMicId:  # Do no repeat check when this is the same mic
                micName = particle.getCoordinate().getMicName()
                if micName in inputMicDict:
                    self.micDict[micName] = inputMicDict[micName]
                lastMicId = micId

    def _insertAllSteps(self):
        self._createMicDict()
        self._defineValues()
        self._prepareCommand()

        convIdDeps = [self._insertFunctionStep('convertInputStep')]
        refineDeps = []

        for micName, mic in self.micDict.items():
            stepId = self._insertFunctionStep('refineCtfStep', mic.getFileName(), micName,
                                              prerequisites=convIdDeps)
            refineDeps.append(stepId)

        self._insertFunctionStep('createOutputStep', prerequisites=refineDeps)

    def _iterParticlesMic(self, newMicCallback):
        """ Iterate through particles sorting by micId and only for
        those that are present in the input set of micrographs. """
        inputParts = self.inputParticles.get()
        lastMicId = None

        for particle in inputParts.iterItems(orderBy=['_micId', 'id']):
            coord = particle.getCoordinate()
            micId = particle.getMicId()
            micName = coord.getMicName()

            if micId != lastMicId:  # Do no repeat check when this is the same mic
                mic = self.micDict.get(micName, None)
                if mic is None:
                    print("Skipping all particles from micrograph, "
                          "key %s not found" % micName)
                else:
                    newMicCallback(mic)  # Notify about a new micrograph found
                lastMicId = micId

            if mic is not None:
                yield particle

    def convertInputStep(self):
        inputParts = self.inputParticles.get()
        alignType = inputParts.getAlignment()
        inputMics = self._getMicrographs()

        scale = inputParts.getSamplingRate() / inputMics.getSamplingRate()
        doScale = abs(scale - 1.0 > 0.00001)
        if doScale:
            print("Scaling coordinates by a factor *%0.2f*" % scale)

        self._lastWriter = None
        coordDir = self._getTmpPath()

        def _newMic(mic):
            if self._lastWriter:
                self._lastWriter.close()
            micBase = pwutils.removeBaseExt(mic.getFileName())
            posFn = pwutils.join(coordDir, micBase, micBase + '_coords.star')
            self._lastWriter = CoordinatesWriter(posFn)

        for particle in self._iterParticlesMic(newMicCallback=_newMic):
            coord = particle.getCoordinate()
            x, y = coord.getPosition()
            if self.applyShifts:
                shifts = getShifts(particle.getTransform(), alignType)
                x, y = x - int(shifts[0]), y - int(shifts[1])
            if doScale:
                x, y = x * scale, y * scale
            self._lastWriter.writeCoord(x, y)

        if self._lastWriter:
            self._lastWriter.close()  # Close file writing for last mic

    def refineCtfStep(self, micFn, micKey):
        micPath = self._getTmpPath(pwutils.removeBaseExt(micFn))
        # We convert the input micrograph on demand if not in .mrc

        downFactor = self.ctfDownFactor.get()
        ih = emlib.image.ImageHandler()
        micFnMrc = pwutils.join(micPath, pwutils.replaceBaseExt(micFn, 'mrc'))

        if downFactor != 1:
            # Replace extension by 'mrc' cause there are some formats
            # that cannot be written (such as dm3)
            ih.scaleFourier(micFn, micFnMrc, downFactor)
            sps = self.inputMicrographs.get().getScannedPixelSize() * downFactor
            self._params['scannedPixelSize'] = sps
        else:
            ih.convert(micFn, micFnMrc, emlib.DT_FLOAT)

        # Refine input CTFs, match ctf by micName
        if self.useInputCtf and self.ctfRelations.hasValue():
            ctfs = self._getCtfs()

            for ctf in ctfs:
                ctfMicName = ctf.getMicrograph().getMicName()
                ctfMicId = ctf.getMicrograph().getObjId()
                if micKey == ctfMicName or micKey == ctfMicId:
                    # add CTF refine options
                    self._params.update({'refine_input_ctf': 1,
                                         'defU_init': ctf.getDefocusU(),
                                         'defV_init': ctf.getDefocusV(),
                                         'defA_init': ctf.getDefocusAngle(),
                                         'B_init': self.bfactor.get()
                                         })
                    self._args += "--refine_input_ctf %d " % self._params['refine_input_ctf']
                    self._args += "--defU_init %f " % self._params['defU_init']
                    self._args += "--defV_init %f " % self._params['defV_init']
                    self._args += "--defA_init %f " % self._params['defA_init']
                    self._args += "--B_init %f " % self._params['B_init']
                    self._args += "--defU_err %f " % self.defUerr.get()
                    self._args += "--defV_err %f " % self.defVerr.get()
                    self._args += "--defA_err %f " % self.defAerr.get()
                    self._args += "--B_err %f " % self.Berr.get()
                    break

        # Run Gctf refine
        try:
            args = self._args % self._params
            args += ' %s' % micFnMrc
            self.runJob(Plugin.getProgram(), args,
                        env=Plugin.getEnviron())

            # Let's clean the temporary mrc micrograph
            pwutils.cleanPath(micFnMrc)
    
            # move output from tmp to extra
            micFnCtf = pwutils.join(micPath, pwutils.replaceBaseExt(micFn, 'ctf'))
            micFnCtfLog = pwutils.join(micPath, pwutils.removeBaseExt(micFn) + '_gctf.log')
            micFnCtfFit = pwutils.join(micPath, pwutils.removeBaseExt(micFn) + '_EPA.log')
            micFnCtfLocal = pwutils.join(micPath, pwutils.removeBaseExt(micFn) + '_local.star')
    
            micFnCtfOut = self._getPsdPath(micFn)
            micFnCtfLogOut = self._getCtfOutPath(micFn)
            micFnCtfFitOut = self._getCtfFitOutPath(micFn)
            micFnCtfLocalOut = self._getCtfLocalOutPath(micFn)
    
            pwutils.moveFile(micFnCtf, micFnCtfOut)
            pwutils.moveFile(micFnCtfLog, micFnCtfLogOut)
            pwutils.moveFile(micFnCtfFit, micFnCtfFitOut)
            pwutils.moveFile(micFnCtfLocal, micFnCtfLocalOut)
        except:
            print("ERROR: Gctf has failed on %s" % micFnMrc)
            import traceback
            traceback.print_exc()

    def createOutputStep(self):
        inputParts = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        partSet.copyInfo(inputParts)
        self._rowList = None
        self._rowCounter = 0

        def _newMic(mic):
            micBase = pwutils.removeBaseExt(mic.getFileName())
            ctfFn = self._getCtfLocalOutPath(micBase)
            self._rowCounter = 0
            if pwutils.exists(ctfFn):
                self._rowList = [row.clone() for row in md.iterRows(ctfFn)]
            else:
                self._rowList = None

        for particle in self._iterParticlesMic(newMicCallback=_newMic):
            if self._rowList is None:  # Ignore particles if not CTF
                continue
            newPart = particle.clone()
            row = self._rowList[self._rowCounter]
            self._rowCounter += 1
            rowToCtfModel(row, newPart.getCTF())
            partSet.append(newPart)

        self._defineOutputs(outputParticles=partSet)
        self._defineTransformRelation(self.inputParticles, partSet)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        if Plugin.getActiveVersion() in ['1.18']:
            errors.append('Gctf version 1.18 does not support local refinement.'
                          ' Please use version 1.06.')

        if self.useInputCtf and not self._getCtfs:
            errors.append("Please provide input CTFs for refinement.")

        return errors

    def _summary(self):
        summary = []

        if not hasattr(self, 'outputParticles'):
            summary.append("Output is not ready yet.")
        else:
            summary.append("CTF refinement of %d particles."
                           % self.inputParticles.get().getSize())

        return summary

    def _methods(self):
        if self.inputParticles.get() is None:
            return ['Input particles not available yet.']
        methods = "We refined the CTF of "
        methods += self.getObjectTag('inputParticles')
        methods += " using Gctf [Zhang2016]. "
        methods += self.methodsVar.get('')

        if self.hasAttribute('outputParticles'):
            methods += 'Output particles: %s' % self.getObjectTag('outputParticles')

        return [methods]

    # -------------------------- UTILS functions -------------------------------
    def _defineValues(self):
        """ This function get some parameters of the micrographs"""
        self.inputMics = self._getMicrographs()
        acq = self.inputMics.getAcquisition()

        self._params = {'voltage': acq.getVoltage(),
                        'sphericalAberration': acq.getSphericalAberration(),
                        'magnification': acq.getMagnification(),
                        'ampContrast': acq.getAmplitudeContrast(),
                        'samplingRate': self.inputMics.getSamplingRate(),
                        'scannedPixelSize': self.inputMics.getScannedPixelSize(),
                        'windowSize': self.windowSize.get(),
                        'lowRes': self.lowRes.get(),
                        'highRes': self.highRes.get(),
                        'minDefocus': self.minDefocus.get(),
                        'maxDefocus': self.maxDefocus.get()
                        }

    def _prepareCommand(self):
        sampling = self._getMicrographs().getSamplingRate() * self.ctfDownFactor.get()
        self._params['sampling'] = sampling
        if self._params['lowRes'] > 50:
            self._params['lowRes'] = 50
        self._params['step_focus'] = self.stepDefocus.get()

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

        # local refine options
        self._args += "--do_local_refine 1 --boxsuffix _coords.star "
        self._args += "--local_radius %d " % self.locRad.get()
        self._args += "--local_avetype %d " % self.locAveType.get()
        self._args += "--local_boxsize %d " % self.locBoxSize.get()
        self._args += "--local_overlap % 0.2f " % self.locOverlap.get()
        self._args += "--local_resL %d " % self.locResL.get()
        self._args += "--local_resH %d " % self.locResH.get()
        self._args += "--refine_local_astm %d " % (1 if self.locAstm else 0)

        self._args += "--EPA_oversmp %d " % self.EPAsmp.get()

        if self.doPhShEst:
            self._args += "--phase_shift_L %f " % self.phaseShiftL.get()
            self._args += "--phase_shift_H %f " % self.phaseShiftH.get()
            self._args += "--phase_shift_S %f " % self.phaseShiftS.get()
            self._args += "--phase_shift_T %d " % (1 + self.phaseShiftT.get())

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

    def _getCtfLocalOutPath(self, micFn):
        micFnBase = pwutils.removeBaseExt(micFn)
        return self._getExtraPath(micFnBase + '_local.star')

    def _getMicrographs(self):
        return self.inputMicrographs.get()

    def _getCtfs(self):
        return self.ctfRelations.get() if self.ctfRelations.hasValue() else None
