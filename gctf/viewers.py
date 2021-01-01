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

from pyworkflow.gui.project import ProjectWindow
import pyworkflow.utils as pwutils
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pwem.viewers import EmPlotter, CtfView, showj

from . import Plugin
from .protocols import ProtGctf


def createCtfPlot(ctfSet, ctfId):
    ctfModel = ctfSet[ctfId]
    psdFn = ctfModel.getPsdFile()
    fn = pwutils.removeExt(psdFn) + "_EPA.log"
    xplotter = EmPlotter(windowTitle='CTF Fitting')
    plot_title = "CTF Fitting"
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


OBJCMD_GCTF = "Display Ctf Analysis"

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
