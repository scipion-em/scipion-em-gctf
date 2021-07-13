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

try:
    from tomo.protocols import ProtImportTs
except ImportError as e:
    if "'tomo'" not in str(e):
        raise e

from ..protocols import ProtTsGctf


class TestBase(BaseTest):
    @classmethod
    def runImportTiltSeries(cls, filesPath, pattern, voltage, magnification,
                            sphericalAberration, amplitudeContrast,
                            samplingRate, stepAngle, anglesFrom=0,
                            minAngle=-60.0, maxAngle=60.0):
        cls.protImportTS = cls.newProtocol(ProtImportTs,
                                           filesPath=filesPath,
                                           filesPattern=pattern,
                                           voltage=voltage,
                                           anglesFrom=anglesFrom,
                                           magnification=magnification,
                                           sphericalAberration=sphericalAberration,
                                           amplitudeContrast=amplitudeContrast,
                                           samplingRate=samplingRate,
                                           minAngle=minAngle,
                                           maxAngle=maxAngle,
                                           stepAngle=stepAngle)
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
                                                   pattern="WTI042413_1series4.st",
                                                   voltage=300,
                                                   magnification=33000,
                                                   sphericalAberration=2.7,
                                                   amplitudeContrast=0.07,
                                                   samplingRate=6.73981,
                                                   stepAngle=2.0)

    def testGctfTs(self):
        print(magentaStr("\n==> Testing gctf:"))
        protCTF = ProtTsGctf()
        protCTF.inputTiltSeries.set(self.protImportTS.outputTiltSeries)
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputSetOfCTFTomoSeries, "SetOfCTFTomoSeries has not been produced.")
