# **************************************************************************
# *
# * Authors:    Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
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
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr

from tomo.protocols import ProtImportTs
from ..protocols import ProtTsGctf


class TestBase(BaseTest):
    @classmethod
    def runImportTiltSeries(cls, **kwargs):
        cls.protImportTS = cls.newProtocol(ProtImportTs, **kwargs)
        cls.launchProtocol(cls.protImportTS)
        return cls.protImportTS


class TestGctfTs(TestBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('tutorialDataImodCTF')
        cls.inputSoTS = cls.inputDataSet.getFile('tsCtf1')

        print(magentaStr("\n==> Importing data - tilt series:"))
        cls.protImportTS = cls.runImportTiltSeries(filesPath=os.path.dirname(cls.inputSoTS),
                                                   filesPattern="WTI042413_1series4.mdoc",
                                                   voltage=300,
                                                   sphericalAberration=2.7,
                                                   amplitudeContrast=0.07,
                                                   anglesFrom=2)

    def testGctfTs(self):
        print(magentaStr("\n==> Testing gctf:"))
        protCTF = ProtTsGctf()
        protCTF.inputTiltSeries.set(self.protImportTS.outputTiltSeries)
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputSetOfCTFTomoSeries, "SetOfCTFTomoSeries has not been produced.")
