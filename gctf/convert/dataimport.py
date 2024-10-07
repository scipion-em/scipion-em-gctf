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

import pyworkflow.utils as pwutils
from pwem.objects import CTFModel

from gctf.convert.convert import readCtfModel


class GctfImportCTF:
    """ Import CTF estimated with GCTF. """
    def __init__(self, protocol):
        self.protocol = protocol
        self.copyOrLink = self.protocol.getCopyOrLink()

    def importCTF(self, mic, fileName):
        ctf = CTFModel()
        ctf.setMicrograph(mic)
        readCtfModel(ctf, fileName)
        
        fnBase = pwutils.removeExt(fileName)
        psdFile = self._findPsdFile(fnBase)
        ctf.setPsdFile(psdFile)

        return ctf

    @staticmethod
    def _findPsdFile(fnBase):
        """ Try to find the given PSD file associated with the cttfind log file
        We handle special cases of .ctf extension and _ctffindX prefix for Relion runs
        """
        for suffix in ['_psd.mrc', '.mrc', '_ctf.mrcs',
                       '.mrcs', '.ctf']:
            psdPrefixes = [fnBase,
                           fnBase.replace('_ctffind3', ''),
                           fnBase.replace('_gctf', '')]
            for prefix in psdPrefixes:
                psdFile = prefix + suffix
                if os.path.exists(psdFile):
                    if psdFile.endswith('.ctf'):
                        psdFile += ':mrc'
                    return psdFile
        return None
