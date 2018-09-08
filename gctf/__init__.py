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

import pyworkflow.em
import pyworkflow.utils as pwutils

from .constants import *


_logo = "gctf_logo.png"
_references = ['Zhang2016']


class Plugin(pyworkflow.em.Plugin):
    _homeVar = GCTF_HOME
    _pathVars = [GCTF_HOME]
    _supportedVersions = ['0.50', '1.06', '1.18']

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(GCTF_HOME, 'gctf-1.06')
        cls._defineVar(GCTF, 'Gctf-v1.06_sm_20_cu8.0_x86_64')

    @classmethod
    def getEnviron(cls):
        """ Return the environ settings to run Gctf program. """
        environ = pwutils.Environ(os.environ)

        # Take Scipion CUDA library path
        cudaLib = environ.getFirst((GCTF_CUDA_LIB, CUDA_LIB))
        environ.addLibrary(cudaLib)

        return environ

    @classmethod
    def getProgram(cls):
        """ Return the program binary that will be used. """
        if (GCTF not in os.environ or
            GCTF_HOME not in os.environ):
            return None

        return os.path.join(cls.getHome('bin'),
                            os.path.basename(os.environ[GCTF]))

    @classmethod
    def isNewVersion(cls):
        # exclude the oldest 0.50 version
        return not cls.getActiveVersion().startswith("0.50")


pyworkflow.em.Domain.registerPlugin(__name__)
