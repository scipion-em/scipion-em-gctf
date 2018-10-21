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
    def defineBinaries(cls, env):
        env.addPackage('gctf', version='0.50',
                       tar='Gctf_v0.50.tgz')

        env.addPackage('gctf', version='1.06',
                       tar='Gctf_v1.06.tgz',
                       default=True)

        env.addPackage('gctf', version='1.18',
                       tar='Gctf_v1.18.tgz')

    @classmethod
    def getProgram(cls):
        """ Return the program binary that will be used. """
        return os.path.join(cls.getHome('bin'), cls.getVar(GCTF))

    @classmethod
    def isNewVersion(cls):
        return not  pyworkflow.em.Plugin.getActiveVersion().startswith("0.50")


pyworkflow.em.Domain.registerPlugin(__name__)
