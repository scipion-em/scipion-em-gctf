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

import pwem
import pyworkflow.utils as pwutils

from .constants import *


__version__ = '3.0.12'
_logo = "gctf_logo.png"
_references = ['Zhang2016']


class Plugin(pwem.Plugin):
    _homeVar = GCTF_HOME
    _pathVars = [GCTF_HOME]
    _supportedVersions = ['1.06', '1.18']
    _url = "https://github.com/scipion-em/scipion-em-gctf"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(GCTF_HOME, 'gctf-1.18')
        cls._defineVar(GCTF, 'Gctf_v1.18_sm30-75_cu10.1')
        cls._defineVar(GCTF_CUDA_LIB, pwem.Config.CUDA_LIB)

    @classmethod
    def getEnviron(cls):
        """ Return the environ settings to run Gctf program. """
        environ = pwutils.Environ(os.environ)
        # Get Gctf CUDA library path if defined
        cudaLib = cls.getVar(GCTF_CUDA_LIB, pwem.Config.CUDA_LIB)
        environ.addLibrary(cudaLib)
        return environ

    @classmethod
    def defineBinaries(cls, env):
        env.addPackage('gctf', version='1.06',
                       tar='Gctf_v1.06.tgz')

        env.addPackage('gctf', version='1.18',
                       tar='Gctf_v1.18.tgz',
                       default=True)

    @classmethod
    def getProgram(cls):
        """ Return the program binary that will be used. """
        return os.path.join(cls.getHome('bin'),
                            os.path.basename(cls.getVar(GCTF)))
