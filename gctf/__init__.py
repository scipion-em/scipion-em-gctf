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
from pyworkflow.config import VarTypes
import pyworkflow.utils as pwutils

from gctf.constants import *


__version__ = '3.2.1'
_logo = "gctf_logo.png"
_references = ['Zhang2016']


class Plugin(pwem.Plugin):
    _homeVar = GCTF_HOME
    _pathVars = [GCTF_HOME]
    _supportedVersions = [V1_18]
    _url = "https://github.com/scipion-em/scipion-em-gctf"

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(GCTF_HOME, f'gctf-{V1_18}',
                         description='Path to the Gctf installation folder',
                         var_type=VarTypes.STRING)
        cls._defineVar(GCTF, f'Gctf_v{V1_18}_sm30-75_cu10.1',
                       description='Gctf binary filename',
                       var_type=VarTypes.STRING)
        cls._defineVar(GCTF_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD,
                       description='GCTF environment activation command',
                       var_type=VarTypes.STRING)

    @classmethod
    def getEnviron(cls):
        """ Return the environ settings to run Gctf program. """
        environ = pwutils.Environ(os.environ)
        environ.update({'PATH': Plugin.getHome('bin')},
                       position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def getGctfEnvActivation(cls):
        return cls.getVar(GCTF_ENV_ACTIVATION)

    @classmethod
    def getDependencies(cls):
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    @classmethod
    def defineBinaries(cls, env):
        from scipion.install.funcs import CondaCommandDef
        installCmd = CondaCommandDef("gctf", cls.getCondaActivationCmd())
        installCmd.create(extraCmds='-y cudatoolkit=10.1')

        env.addPackage('gctf', version=V1_18,
                       tar=f'Gctf_v{V1_18}.tgz',
                       commands=installCmd.getCommands(),
                       neededProgs=cls.getDependencies(),
                       default=True)

    @classmethod
    def getProgram(cls):
        """ Return the program binary that will be used. """
        return " ".join([
            cls.getCondaActivationCmd(),
            cls.getGctfEnvActivation(),
            "&& LD_LIBRARY_PATH=$CONDA_PREFIX/lib &&",
            os.path.basename(cls.getVar(GCTF))
        ])
