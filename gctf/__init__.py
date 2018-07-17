# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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
"""
This sub-package contains data and protocol classes
wrapping Kai Zhang's GCTF program
"""
import os

import pyworkflow.em
import pyworkflow.utils as pwutils


_logo = "gctf_logo.png"

GCTF_HOME = 'GCTF_HOME'


# The following class is required for Scipion to detect this Python module
# as a Scipion Plugin. It needs to specify the PluginMeta __metaclass__
# Some function related to the underlying package binaries need to be
# implemented
class Plugin:
    __metaclass__ = pyworkflow.em.PluginMeta

    @classmethod
    def getEnviron(cls):
        """ Return the environ settings to run Gctf program. """
        environ = pwutils.Environ(os.environ)

        # Take Scipion CUDA library path
        cudaLib = environ.getFirst(('GCTF_CUDA_LIB', 'CUDA_LIB'))
        environ.addLibrary(cudaLib)

        return environ

    @classmethod
    def getActiveVersion(cls):
        """ Return the version of the Gctf binary that is currently active. """
        path = os.environ[GCTF_HOME]
        for v in cls.getSupportedVersions():
            if v in path or v in os.path.realpath(path):
                return v
        return ''

    @classmethod
    def getSupportedVersions(cls):
        """ Return the list of supported binary versions. """
        return ['0.50', '1.06']

    @classmethod
    def validateInstallation(cls):
        """ This function will be used to check if package is
        properly installed."""
        environ = cls.getEnviron()

        missingPaths = ["%s: %s" % (var, environ[var])
                        for var in [GCTF_HOME]
                        if not os.path.exists(environ[var])]

        return (["Missing variables:"] + missingPaths) if missingPaths else []

