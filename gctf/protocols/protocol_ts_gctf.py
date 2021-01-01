# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

from pwem.protocols import EMProtocol, pwutils
from pyworkflow.protocol import STEPS_PARALLEL
import pyworkflow.protocol.params as params

try:
    from tomo.protocols import ProtTsEstimateCTF
except ImportError:
    raise ImportError(
        'To use a Tomography protocol scipion-em-tomo plugin is required.'
        ' See https://github.com/scipion-em/scipion-em-tomo for further details')

from .program_gctf import ProgramGctf
import gctf


class ProtTsGctf(ProtTsEstimateCTF):
    """
    CTF estimation on Tilt-Series using GCTF.
    """
    _label = 'tiltseries gctf'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _initialize(self):
        ProtTsEstimateCTF._initialize(self)
        self._gctfProgram = ProgramGctf(self)

    def _defineProcessParams(self, form):
        form.addParam('recalculate', params.BooleanParam, default=False,
                      condition='recalculate',
                      label="Do recalculate ctf?")
        form.addParam('continueRun', params.PointerParam, allowsNull=True,
                      condition='recalculate', label="Input previous run",
                      pointerClass='ProtTsCtffind')
        form.addHidden('sqliteFile', params.FileParam,
                       condition='recalculate',
                       allowsNull=True)
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

        ProgramGctf.defineProcessParams(form)

    # --------------------------- STEPS functions ----------------------------
    def _estimateCtf(self, workingDir, tiFn, ti):
        try:
            program, args = self._gctfProgram.getCommand(
                scannedPixelSize=self._params['scannedPixelSize'])
            args += ' %s/*.mrc' % workingDir
            self.runJob(program, args, env=gctf.Plugin.getEnviron())

            ext = self._gctfProgram.getExt()

            # Move files we want to keep
            micBase = pwutils.removeBaseExt(tiFn)

            def _getFile(suffix):
                print("File: %s" % os.path.join(workingDir, micBase + suffix))
                return os.path.join(workingDir, micBase + suffix)

            # move output from tmp to extra
            pwutils.moveFile(_getFile(ext),
                             self._getExtraPath(micBase + '_ctf.mrc'))
            pwutils.moveFile(_getFile('_gctf.log'), self._getTmpPath())
            pwutils.moveFile(_getFile('_EPA.log'), self._getTmpPath())

        except Exception as ex:
            print("ERROR: Gctf has failed for %s" % tiFn)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        return [self.summaryVar.get('')]

    # --------------------------- UTILS functions ----------------------------
    def _getArgs(self):
        """ Return a list with parameters that will be passed to the process
        TiltSeries step. It can be redefined by subclasses.
        """
        return []

    def getCtf(self, ti):
        """ Parse the CTF object estimated for this Tilt-Image
        """
        prefix = self.getTiPrefix(ti)
        psd = self._getExtraPath(prefix + '_ctf.mrc')
        outCtf = self._getTmpPath(prefix + '_gctf.log')
        return self._gctfProgram.parseOutputAsCtf(outCtf, psdFile=psd)
