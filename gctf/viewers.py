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
from matplotlib.figure import Figure

from pwem.emlib.image import ImageHandler
from pyworkflow.gui.project import ProjectWindow
import pyworkflow.utils as pwutils
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pwem.viewers import EmPlotter, CtfView, showj
from tomo.viewers.viewers_data import CtfEstimationTomoViewer

from . import Plugin
from .protocols import ProtGctf, ProtTsGctf


def createCtfPlot(ctfSet, ctfId):
    ctfModel = ctfSet[ctfId]
    psdFn = ctfModel.getPsdFile()
    fn = pwutils.removeExt(psdFn) + "_EPA.log"
    xplotter = EmPlotter(windowTitle='CTF Fitting')
    plot_title = getPlotSubtitle(ctfModel)
    a = xplotter.createSubPlot(plot_title, 'Resolution (Angstroms)', 'CTF')
    a.invert_xaxis()
    version = Plugin.getActiveVersion()
    curves = [1, 4, 5] if version == '1.18' else [1, 3, 4]

    for i in curves:
        _plotCurve(a, i, fn)
    xplotter.showLegend(['simulated CTF',
                         # 'equiphase avg.',
                         # 'bg', #  only for v1.18
                         'equiphase avg. - bg',
                         'cross correlation'])
    a.grid(True)
    xplotter.show()


OBJCMD_GCTF = "GCTF plot results"

ProjectWindow.registerObjectCommand(OBJCMD_GCTF, createCtfPlot)


def _plotCurve(a, i, fn):
    freqs = _getValues(fn, 0)
    curv = _getValues(fn, i)
    a.plot(freqs, curv)


def _getValues(fn, col):
    values = list()
    with open(fn) as f:
        for line in f:
            if not line.startswith('Resolution', 2, 12):
                values.append(float(line.split()[col]))
    return values


def getPlotSubtitle(ctf):
    """ Create plot subtitle using CTF values. """
    ang = u"\u212B"
    deg = u"\u00b0"
    def1, def2, angle = ctf.getDefocus()
    phSh = ctf.getPhaseShift()
    score = ctf.getFitQuality()
    res = ctf.getResolution()

    title = "Def1: %d %s | Def2: %d %s | Angle: %0.1f%s | " % (
        def1, ang, def2, ang, angle, deg)

    if phSh is not None:
        title += "Phase shift: %0.2f %s | " % (phSh, deg)

    title += "Fit: %0.1f %s | Score: %0.3f" % (res, ang, score)

    return title


class GctfViewer(Viewer):
    """ Visualization of Gctf results. """
    _environments = [DESKTOP_TKINTER]
    _targets = [ProtGctf]

    def _visualize(self, prot, **kwargs):
        outputCTF = getattr(prot, 'outputCTF', None)

        if outputCTF is not None:
            ctfView = CtfView(self._project, outputCTF)
            viewParams = ctfView.getViewParams()
            viewParams[showj.OBJCMDS] = "'%s'" % OBJCMD_GCTF
            return [ctfView]
        else:
            return [self.infoMessage("The output SetOfCTFs has not been "
                                     "produced", "Missing output")]


class CtfEstimationTomoViewerGctf(CtfEstimationTomoViewer):
    """ This class implements a view using Tkinter CtfEstimationListDialog
    and the CtfEstimationTreeProvider.
    """
    _targets = [ProtTsGctf]

    def plot1D(self, ctfSet, ctfId):
        ctfModel = ctfSet[ctfId]
        psdFn = ctfModel.getPsdFile()
        fn = os.path.join(pwutils.removeExt(psdFn).replace("_ctf", "") + '_EPA.log')

        xplotter = EmPlotter(windowTitle='GCTF results')
        plot_title = '%s # %d\n' % (ctfSet.getTsId(), ctfId) + getPlotSubtitle(ctfModel)
        a = xplotter.createSubPlot(plot_title, 'Resolution (Angstroms)', 'CTF')
        a.invert_xaxis()
        version = Plugin.getActiveVersion()
        curves = [1, 4, 5] if version == '1.18' else [1, 3, 4]

        for i in curves:
            _plotCurve(a, i, fn)
        xplotter.showLegend(['simulated CTF',
                             # 'equiphase avg.',
                             # 'bg', #  only for v1.18
                             'equiphase avg. - bg',
                             'cross correlation'])
        a.grid(True)

        return xplotter

    def plot2D(self, ctfSet, ctfId):
        ctfModel = ctfSet[ctfId]
        psdFn = ctfModel.getPsdFile()
        img = ImageHandler().read(psdFn)
        fig = Figure(figsize=(7, 7), dpi=100)
        psdPlot = fig.add_subplot(111)
        psdPlot.get_xaxis().set_visible(False)
        psdPlot.get_yaxis().set_visible(False)
        psdPlot.set_title('%s # %d\n' % (ctfSet.getTsId(), ctfId) + getPlotSubtitle(ctfModel))
        psdPlot.imshow(img.getData(), cmap='gray')

        return fig
