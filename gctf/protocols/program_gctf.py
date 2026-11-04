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

from gctf import Plugin
from gctf.convert.convert import readCtfModel, parseGctfOutput
from gctf.constants import CCC


class ProgramGctf:
    """
    Provides an interface for configuring and executing Gctf-based Contrast Transfer Function estimation workflows in cryo-EM image processing. The class centralizes the definition of acquisition, optical, and refinement parameters required for reliable CTF determination, allowing different protocols to perform consistent and reproducible analyses of electron micrographs. Its main objective is to simplify the interaction with the Gctf engine while exposing biologically relevant controls that influence the quality and interpretability of CTF estimation results. More info: https://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/zhang-software/#gctf

    AI Generated:

    Gctf Program Wrapper (ProgramGctf) - User Manual
        Overview

        The Gctf Program Wrapper provides a unified environment for preparing,
        configuring, and executing CTF estimation tasks using the Gctf software
        package. In cryo-EM workflows, accurate CTF estimation is essential
        because it determines how microscope optics have altered the recorded
        image and directly influences the quality of downstream reconstruction,
        particle alignment, and refinement procedures.

        This wrapper is intended to support protocols that require automated
        CTF determination while maintaining flexibility for advanced users who
        need precise control over acquisition parameters, search ranges, and
        refinement strategies. The class acts as an integration layer between
        biological image-processing workflows and the external Gctf executable,
        ensuring that microscope metadata and processing settings are applied
        consistently across datasets.

        Biological Context and Importance

        In transmission electron microscopy, the CTF modifies image contrast in
        a frequency-dependent manner. Correct estimation of defocus and
        astigmatism is therefore critical for recovering high-resolution
        structural information. Poor CTF estimation can reduce map quality,
        introduce reconstruction artifacts, or prevent accurate interpretation
        of flexible or heterogeneous biological assemblies.

        This wrapper is particularly relevant for high-throughput cryo-EM
        facilities and automated processing pipelines where large collections
        of micrographs must be analyzed reproducibly. It supports both routine
        datasets and more challenging acquisitions involving phase plates,
        strong astigmatism, or high-resolution refinement strategies.

        Input Parameters and Acquisition Information

        The wrapper manages all acquisition-related parameters required for CTF
        estimation, including sampling rate, acceleration voltage, spherical
        aberration, amplitude contrast, and detector pixel size. These
        parameters define the physical microscope model used during estimation
        and should accurately reflect the experimental setup.

        For biological users, ensuring correct acquisition metadata is one of
        the most important prerequisites for obtaining meaningful CTF values.
        Incorrect voltage or pixel size values can propagate systematic errors
        throughout the entire reconstruction workflow.

        The wrapper also supports optional downsampling during CTF estimation.
        Downsampling is often beneficial when Thon rings are concentrated at
        very low spatial frequencies or when computational efficiency is
        important. In practical workflows, moderate downsampling can accelerate
        processing substantially without compromising estimation quality for
        medium-resolution datasets.

        Defocus and Resolution Search Strategy

        The wrapper exposes biologically meaningful search ranges for defocus
        and resolution refinement. These parameters determine the frequency
        interval over which the CTF model is optimized and strongly influence
        estimation stability.

        Wide defocus ranges are useful when imaging conditions vary strongly
        across a dataset or when acquisition metadata is uncertain. Narrower
        ranges improve efficiency and may stabilize estimation in homogeneous
        datasets collected under controlled conditions.

        Resolution limits are equally important. Including excessively high
        frequencies in noisy datasets may lead to unstable fitting, whereas
        restricting the analysis to lower frequencies can improve robustness
        for difficult samples or low-dose acquisitions.

        Astigmatism and Optical Refinement

        The wrapper supports astigmatism estimation and refinement, which are
        critical for accurately modeling imperfections in microscope optics.
        Biological datasets acquired under suboptimal alignment conditions or
        from older instruments may exhibit substantial astigmatism, making
        refinement especially important.

        Additional refinement controls are available for high-resolution
        analysis. These options become particularly useful when processing
        datasets intended for near-atomic reconstructions, where small errors
        in CTF estimation can noticeably affect map interpretability and model
        building accuracy.

        EPA Visualization and Spectral Interpretation

        The Equiphase Averaging functionality improves the visual appearance
        and interpretability of power spectra generated during CTF analysis.
        Although this enhancement does not directly alter the numerical CTF
        determination, it can significantly help users visually inspect the
        quality of Thon ring fitting and identify problematic micrographs.

        In biological practice, these visual diagnostics are frequently used
        during data collection sessions to assess ice quality, contamination,
        beam-induced motion, or microscope stability. Cleaner spectra often
        correlate with higher-quality downstream reconstructions.

        Phase Plate Data Processing

        The wrapper includes dedicated support for phase-shift estimation in
        datasets acquired using phase plates. Such experiments require special
        treatment because the phase plate introduces additional optical shifts
        that must be estimated accurately for proper reconstruction.

        Phase-plate cryo-EM is particularly valuable for weakly scattering
        biological specimens, small proteins, and low-contrast complexes.
        However, the additional optical complexity increases the importance of
        robust parameter selection and refinement stability.

        Validation and Refinement Options

        Validation tools are available to help assess the reliability of the
        estimated parameters. These checks are especially useful in automated
        pipelines where large numbers of micrographs are processed without
        continuous manual inspection.

        The wrapper also supports different refinement strategies for difficult
        datasets. In challenging biological cases, such as thick ice, strong
        preferred orientation, or low particle density, testing alternative
        refinement approaches may substantially improve estimation quality.

        GPU Acceleration and High-Throughput Processing

        The class supports GPU-oriented execution strategies intended for
        modern cryo-EM processing environments. High-throughput facilities
        frequently rely on GPU acceleration to process large datasets in near
        real time during microscope acquisition sessions.

        Efficient GPU usage becomes especially important for large-scale
        single-particle studies, tomography workflows, and screening campaigns
        where thousands of micrographs must be analyzed rapidly.

        Outputs and Downstream Usage

        The resulting CTF information includes defocus estimation, astigmatism,
        spectral diagnostics, and associated metadata required for subsequent
        cryo-EM processing stages. These outputs are typically consumed by
        particle picking, motion correction validation, particle refinement,
        and 3D reconstruction protocols.

        From a biological perspective, accurate CTF estimation contributes
        directly to improved map resolution, better density interpretability,
        and more reliable structural conclusions.

        Practical Recommendations

        For routine cryo-EM workflows, it is generally advisable to begin with
        conservative search ranges and standard refinement settings. Visual
        inspection of the resulting spectra remains important, particularly for
        heterogeneous datasets or difficult imaging conditions.

        Moderate downsampling often improves efficiency without major loss of
        information. High-resolution refinement should typically be reserved
        for datasets that already demonstrate strong high-frequency signal.

        For phase-plate datasets, users should carefully tune phase-shift
        search ranges and monitor fitting stability. In difficult datasets,
        comparing multiple refinement strategies may help identify the most
        reliable solution.

        Final Perspective

        Reliable CTF estimation is one of the foundational steps in cryo-EM
        image processing. The ProgramGctf wrapper provides a structured and
        flexible interface for controlling this process while maintaining
        compatibility with automated and high-throughput workflows. Careful
        selection of optical parameters, refinement settings, and validation
        strategies is essential for achieving biologically meaningful and
        high-resolution structural results.
    """
    def __init__(self, protocol):
        self._args, self._params = self._getArgs(protocol)  # Load general arguments
        self._ext = '.pow' if not protocol.doEPA else '.epa'

    @classmethod
    def defineInputParams(cls, form):
        """ Define input parameters from this program into the given form. """
        form.addSection(label='Input')
        form.addParam('inputMicrographs', params.PointerParam, important=True,
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
        args += "--smooth_resL %d " % protocol.smoothResL
        args += "--EPA_oversmp %d " % protocol.EPAsmp

        if protocol.doPhShEst:
            args += "--phase_shift_L %f " % protocol.phaseShiftL
            args += "--phase_shift_H %f " % protocol.phaseShiftH
            args += "--phase_shift_S %f " % protocol.phaseShiftS
            args += "--phase_shift_T %d " % (1 + protocol.phaseShiftT.get())

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
