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
from pyworkflow.utils import weakImport

with weakImport('tomo'):
    from pyworkflow.tests import BaseTest, DataSet, setupTestProject
    from pyworkflow.utils import magentaStr, cyanStr

    from tomo.protocols import ProtImportTs
    from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
    from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer

    from gctf.protocols import ProtTsGctf


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

            self.assertIsNotNone(protCTF.CTFs,
                                 "SetOfCTFTomoSeries has not been produced.")

    class TestGctfTsTCL(TestBaseCentralizedLayer):
        importedTs = None
        unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value

        @classmethod
        def setUpClass(cls):
            setupTestProject(cls)
            cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
            cls.runPrevProtocols()

        @classmethod
        def runPrevProtocols(cls):
            print(cyanStr('--------------------------------- RUNNING PREVIOUS PROTOCOLS ---------------------------------'))
            cls._runPreviousProtocols()
            print(cyanStr('\n-------------------------------- PREVIOUS PROTOCOLS FINISHED ---------------------------------'))

        @classmethod
        def _runPreviousProtocols(cls):
            cls.importedTs = cls._runImportTs()

        @classmethod
        def _runImportTs(cls, filesPattern=DataSetRe4STATuto.tsPattern.value,
                         exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value):
            print(magentaStr("\n==> Importing the tilt series:"))
            protImportTs = cls.newProtocol(ProtImportTs,
                                           filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                           filesPattern=filesPattern,
                                           exclusionWords=exclusionWords,
                                           anglesFrom=2,  # From tlt file
                                           voltage=DataSetRe4STATuto.voltage.value,
                                           magnification=DataSetRe4STATuto.magnification.value,
                                           sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                           amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                           samplingRate=cls.unbinnedSRate,
                                           doseInitial=DataSetRe4STATuto.initialDose.value,
                                           dosePerFrame=DataSetRe4STATuto.dosePerTiltImgWithTltFile.value,
                                           tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

            cls.launchProtocol(protImportTs)
            tsImported = getattr(protImportTs, 'outputTiltSeries', None)
            return tsImported

        @classmethod
        def _runEstimateCtf(cls, inTsSet, objLabel=None):
            print(magentaStr("\n==> Running the CTF estimation:"))
            protEstimateCtf = cls.newProtocol(ProtTsGctf,
                                              inputTiltSeries=inTsSet,
                                              lowRes=50,
                                              highRes=4,
                                              minDefocus=15000,
                                              maxDefocus=50000)
            if objLabel:
                protEstimateCtf.setObjLabel(objLabel)
            cls.launchProtocol(protEstimateCtf)
            outTsSet = getattr(protEstimateCtf, protEstimateCtf._possibleOutputs.CTFs.name, None)
            return outTsSet

        def _checkCtfs(self, inCtfSet):
            expectedSetSize = 2  # TS_03 and TS_54
            self.checkCTFs(inCtfSet, expectedSetSize=expectedSetSize)

        def testEstimateCtf01(self):
            ctfs = self._runEstimateCtf(self.importedTs)
            self._checkCtfs(ctfs)
