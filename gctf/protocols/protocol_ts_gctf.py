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
from enum import Enum

from pyworkflow.object import Set, Boolean
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.convert.headers import MRC, getFileFormat
from pwem.emlib.image import DT_FLOAT, ImageHandler
from pwem.protocols import EMProtocol

from tomo.objects import CTFTomo, SetOfCTFTomoSeries, TiltImage, CTFTomoSeries
from tomo.protocols.protocol_ts_estimate_ctf import createCtfParams

from gctf.protocols.program_gctf import ProgramGctf
from gctf import Plugin


class TsGctfOutputs(Enum):
    CTFs = SetOfCTFTomoSeries


class ProtTsGctf(EMProtocol):
    """ CTF estimation on a set of tilt series using GCTF. """
    _label = 'tilt-series gctf'
    _devStatus = PROD
    _possibleOutputs = TsGctfOutputs
    recalculate = Boolean(False, objDoStore=False)  # Legacy Sep 2024: to fake old recalculate param
    # that is still used in the ProtCTFMicrographs (to be removed)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        self._gctfProgram = None
        self.inTsSet = None
        self.tsDict = None
        self._params = None
        self.ih = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=pwutils.Message.LABEL_INPUT)
        form.addParam('inputTiltSeries', params.PointerParam,
                      important=True,
                      pointerClass='SetOfTiltSeries, SetOfCTFTomoSeries',
                      label='Tilt series')
        form.addParam('ctfDownFactor', params.FloatParam,
                      default=1.,
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
    def _insertAllSteps(self):
        self._initialize()
        pIdList = []
        for tsId in self.tsDict.keys():
            pidProcess = self._insertFunctionStep(self.processTiltSeriesStep,
                                                  tsId, prerequisites=[])
            pidCreateOutput = self._insertFunctionStep(self.createOutputStep,
                                                       tsId, prerequisites=pidProcess)
            pIdList.append(pidCreateOutput)
        self._insertFunctionStep(self.closeOutputSetsStep, prerequisites=pIdList)

    def _initialize(self):
        self.ih = ImageHandler()
        self.inTsSet = self._getInputTs()
        self._params = createCtfParams(self.inTsSet, self.windowSize.get(),
                                       self.lowRes.get(), self.highRes.get(),
                                       self.minDefocus.get(), self.maxDefocus.get(),
                                       downFactor=self.getAttributeValue('ctfDownFactor', 1.0))
        self._gctfProgram = ProgramGctf(self)
        self.tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.inTsSet.iterItems()}

    def processTiltSeriesStep(self, tsId):
        ts = self.tsDict[tsId]
        for ti in ts.iterItems():
            workingDir = self._getTiWorkingDir(ti)
            tiFnMrc = os.path.join(workingDir, self.getTiPrefix(ti) + '.mrc')
            pwutils.makePath(workingDir)
            self._convertInputTi(ti, tiFnMrc)
            self._estimateCtf(workingDir, tiFnMrc)

    def createOutputStep(self, tsId):
        with self._lock:
            outCtfSet = self.getOutputCtfTomoSet()
            # Generate the current CTF tomo series item
            ts = self.tsDict[tsId]
            newCTFTomoSeries = CTFTomoSeries()
            newCTFTomoSeries.copyInfo(ts)
            newCTFTomoSeries.setTiltSeries(ts)
            newCTFTomoSeries.setObjId(ts.getObjId())
            newCTFTomoSeries.setTsId(ts.getTsId())
            outCtfSet.append(newCTFTomoSeries)

            # Generate the ti CTF and populate the corresponding CTF tomo series
            for i, tiltImage in enumerate(ts.iterItems()):
                ctfTomo = self.getCtf(tiltImage)
                ctfTomo.setIndex(tiltImage.getIndex())
                ctfTomo.setAcquisitionOrder(tiltImage.getAcquisitionOrder())
                newCTFTomoSeries.append(ctfTomo)

            outCtfSet.update(newCTFTomoSeries)
            self._store()

    def closeOutputSetsStep(self):
        self._closeOutputSet()

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions ----------------------------
    def getCtf(self, ti: TiltImage) -> CTFTomo:
        """ Parse the CTF object estimated for this Tilt-Image. """
        prefix = self.getTiPrefix(ti)
        psd = self._getExtraPath(prefix + '_ctf.mrc')
        outCtf = self._getTmpPath(prefix + '_gctf.log')
        ctfModel = self._gctfProgram.parseOutputAsCtf(outCtf, psdFile=psd)
        ctfTomo = CTFTomo.ctfModelToCtfTomo(ctfModel)

        return ctfTomo

    def _estimateCtf(self, workingDir, tiFn):
        try:
            program, args = self._gctfProgram.getCommand(
                scannedPixelSize=self._params['scannedPixelSize'])
            args += ' %s/*.mrc' % workingDir
            self.runJob(program, args, env=Plugin.getEnviron())

            ext = self._gctfProgram.getExt()

            # Move files we want to keep
            micBase = pwutils.removeBaseExt(tiFn)

            def _getFile(suffix):
                return os.path.join(workingDir, micBase + suffix)

            # move output from tmp to extra
            pwutils.moveFile(_getFile(ext),
                             self._getExtraPath(micBase + '_ctf.mrc'))
            pwutils.moveFile(_getFile('_gctf.log'), self._getTmpPath())
            pwutils.moveFile(_getFile('_EPA.log'),
                             self._getExtraPath())

        except Exception:
            self.error("ERROR: Gctf has failed for %s" % tiFn)

    def _getInputTs(self, pointer=False):
        if isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries):
            return self.inputTiltSeries.get().getSetOfTiltSeries(pointer=pointer)
        return self.inputTiltSeries.get() if not pointer else self.inputTiltSeries

    def _convertInputTi(self, ti, tiFn):
        """ This function will convert the input tilt-image
        taking into account the downFactor.
        It can be overwritten in subclasses if another behaviour is required.
        """
        downFactor = self.ctfDownFactor.get()

        if not self.ih.existsLocation(ti):
            raise Exception("Missing input file: %s" % ti)

        tiFName = ti.getFileName()
        # Make xmipp considers the input object as TS to work as expected
        if getFileFormat(tiFName) == MRC:
            tiFName = tiFName.split(':')[0] + ':mrcs'
        tiFName = str(ti.getIndex()) + '@' + tiFName

        if downFactor != 1:
            self.ih.scaleFourier(tiFName, tiFn, downFactor)
        else:
            self.ih.convert(tiFName, tiFn, DT_FLOAT)

    @staticmethod
    def getTiRoot(ti: TiltImage):
        return '%s_%02d' % (ti.getTsId(), ti.getObjId())

    def _getTiWorkingDir(self, ti: TiltImage):
        return self._getTmpPath(self.getTiRoot(ti))

    @staticmethod
    def getTiPrefix(ti: TiltImage):
        return '%s_%03d' % (ti.getTsId(), ti.getObjId())

    def getOutputCtfTomoSet(self) -> SetOfCTFTomoSeries:
        outCtfSet = getattr(self, TsGctfOutputs.CTFs.name, None)
        if outCtfSet:
            outCtfSet.enableAppend()
        else:
            inTsPointer = self._getInputTs(pointer=True)
            outCtfSet = SetOfCTFTomoSeries.create(self._getPath(), template='ctfTomoSeries%s.sqlite')
            outCtfSet.setSetOfTiltSeries(inTsPointer)
            outCtfSet.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{self._possibleOutputs.CTFs.name: outCtfSet})
            self._defineSourceRelation(inTsPointer, outCtfSet)
        return outCtfSet

    def getCtfParamsDict(self):
        """ Return a copy of the global params dict,
        to avoid overwriting values. """
        return self._params
