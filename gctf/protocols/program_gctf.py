# **************************************************************************
# *
# * Authors:     Josue Gomez BLanco (josue.gomez-blanco@mcgill.ca) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Department of Anatomy and Cell Biology, McGill University
# * [2] SciLifeLab, Stockholm University
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

from pwem.objects import CTFModel
import pyworkflow.protocol.params as params

from .. import Plugin
from ..convert import readCtfModel, parseGctfOutput
from ..constants import CCC


class ProgramGctf:
    """
    Wrapper of Gctf program that will handle parameters definition
    and also execution of the program with the proper arguments.
    This class is not a Protocol, but it is related, since it can be used from
    protocols that perform CTF estimation.
    """
    def __init__(self, protocol):
        self._args, self._params = self._getArgs(protocol)  # Load general arguments
        if self.isVersion118():
            self._ext = '.pow' if not protocol.doEPA else '.epa'
        else:
            self._ext = '.ctf'

    @classmethod
    def defineInputParams(cls, form):
        """ Define input parameters from this program into the given form. """
        form.addSection(label='Input')
        form.addParam('recalculate', params.BooleanParam, default=False,
                      condition='recalculate',
                      label="Do recalculate ctf?")
        form.addParam('continueRun', params.PointerParam, allowsNull=True,
                      condition='recalculate', label="Input previous run",
                      pointerClass='ProtGctf')
        form.addHidden('sqliteFile', params.FileParam, condition='recalculate',
                       allowsNull=True)
        form.addParam('inputMicrographs', params.PointerParam, important=True,
                      condition='not recalculate',
                      label='Input micrographs',
                      pointerClass='SetOfMicrographs')
        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      help='Set to 1 for no downsampling. Non-integer downsample '
                           'factors are possible. This downsampling is only used '
                           'for estimating the CTF and it does not affect any '
                           'further calculation. Ideally the estimation of the '
                           'CTF is optimal when the Thon rings are not too '
                           'concentrated at the origin (too small to be seen) '
                           'and not occupying the whole power spectrum (since '
                           'this downsampling might entail aliasing).')

    @classmethod
    def defineProcessParams(cls, form):
        form.addParam('windowSize', params.IntParam, default=1024,
                      label='Box size (px)', condition='not recalculate',
                      help='Boxsize in pixels to be used for FFT, 512 or '
                           '1024 highly recommended')

        group = form.addGroup('Search limits')
        line = group.addLine('Resolution (A)', condition='not recalculate',
                             help='The CTF model will be fit to regions '
                                  'of the amplitude spectrum corresponding '
                                  'to this range of resolution.')
        line.addParam('lowRes', params.FloatParam, default=50., label='Min')
        line.addParam('highRes', params.FloatParam, default=4., label='Max')

        line = group.addLine('Defocus search range (A)',
                             condition='not recalculate',
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
        group.addParam('doEPA', params.BooleanParam, default=True,
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

        if cls.isVersion118():
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

        if cls.isVersion118():
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

        form.addParallelSection(threads=1, mpi=1)

    @classmethod
    def getVersion(cls):
        return Plugin.getActiveVersion()

    @classmethod
    def isVersion118(cls):
        return cls.getVersion() in ['1.18']

    def getExt(self):
        return self._ext

    @staticmethod
    def _getProgram():
        """ Return the program to be used. """
        return Plugin.getProgram()

    def getCommand(self, **kwargs):
        """ Return the program and arguments to be run.
        The input keywords argument should contain key-values for
        one micrograph or group of micrographs.
        """
        params = dict(self._params)
        params.update(kwargs)
        return self._getProgram(), self._args % params

    def parseOutput(self, filename):
        """ Retrieve defocus U, V and angle from the
        output file of the program execution.
        """
        return parseGctfOutput(filename)

    def parseOutputAsCtf(self, ctfFile, psdFile=None):
        ctf = CTFModel()
        readCtfModel(ctf, ctfFile)
        if psdFile:
            ctf.setPsdFile(psdFile)

        return ctf

    def _getArgs(self, protocol):
        # Update first the _params dict
        params = protocol.getCtfParamsDict()

        if params['lowRes'] > 50:
            params['lowRes'] = 50

        # defocus is in Angstroms now
        params['minDefocus'] = protocol.minDefocus.get()
        params['maxDefocus'] = protocol.maxDefocus.get()
        params['step_focus'] = protocol.stepDefocus.get()

        args = " --apix %(samplingRate)f "
        args += "--kV %(voltage)f "
        args += "--cs %(sphericalAberration)f "
        args += "--ac %(ampContrast)f "
        args += "--dstep %(scannedPixelSize)f "
        args += "--defL %(minDefocus)f "
        args += "--defH %(maxDefocus)f "
        args += "--defS %(step_focus)f "
        args += "--astm %f " % protocol.astigmatism
        args += "--resL %(lowRes)f "
        args += "--resH %(highRes)f "
        args += "--do_EPA %d " % (1 if protocol.doEPA else 0)
        args += "--boxsize %(windowSize)d "
        args += "--plot_res_ring %d " % (1 if protocol.plotResRing else 0)
        args += "--gid %%(GPU)s "  # Use %% to escape when formatting
        args += "--bfac %d " % protocol.bfactor
        args += "--B_resH %f " % (2 * params['samplingRate'])
        args += "--overlap %f " % protocol.overlap
        args += "--convsize %d " % protocol.convsize
        args += "--do_Hres_ref %d " % (1 if protocol.doHighRes else 0)

        if self.isVersion118():
            args += "--smooth_resL %d " % protocol.smoothResL

        args += "--EPA_oversmp %d " % protocol.EPAsmp

        if protocol.doPhShEst:
            args += "--phase_shift_L %f " % protocol.phaseShiftL
            args += "--phase_shift_H %f " % protocol.phaseShiftH
            args += "--phase_shift_S %f " % protocol.phaseShiftS
            args += "--phase_shift_T %d " % (1 + protocol.phaseShiftT.get())

            if self.isVersion118():
                args += ("--cosearch_refine_ps %d "
                         % (1 if protocol.coSearchRefine else 0))
                args += "--refine_2d_T %d " % protocol.refine2DT

        if protocol.doHighRes:
            args += "--Href_resL %0.3f " % protocol.HighResL
            args += "--Href_resH %0.3f " % protocol.HighResH
            args += "--Href_bfac %d " % protocol.HighResBf

        args += ("--ctfstar NONE --do_validation %d "
                 % (1 if protocol.doValidate else 0))

        return args, params
