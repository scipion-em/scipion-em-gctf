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
import re
import logging
logger = logging.getLogger(__name__)


def parseGctfOutput(filename):
    """ Retrieve defocus U, V, angle, crossCorrelation
    and ctfResolution from the output file of the Gctf execution.
    """

    if os.path.exists(filename):
        # Create an empty list with: defU, defV, angle, CC and resolution
        result = [0.] * 6
        ansi_escape = re.compile(r'\x1b[^m]*m')
        f = open(filename)
        for line in f:
            if 'Final Values' in line:
                parts = line.strip().split()
                # line = DefocusU, DefocusV, Angle, crossCorrelation, Final, Values
                # OR
                # line = DefocusU, DefocusV, Angle, ctfPhaseShift, crossCorrelation, Final, Values
                # Always map defU, defV and angle
                result[0:3] = map(float, parts[0:3])

                if parts[4] == 'Final':  # no ctfPhaseShift
                    result[3] = float(parts[3])
                else:
                    result[3] = float(parts[4])  # CC is now in position 4
                    result[4] = float(parts[3])  # get ctfPhaseShift
            if 'Resolution limit estimated by EPA' in line:
                # Take ctfResolution as a tuple
                # that is the last value in the line
                # but remove escape characters first
                resol = ansi_escape.sub('', line)
                result[5] = float(resol.strip().split()[-1])
                break
        f.close()
    else:
        result = None
        logger.warning(f"Warning: Missing file: {filename}")

    return result


def setWrongDefocus(ctfModel):
    ctfModel.setDefocusU(-999)
    ctfModel.setDefocusV(-1)
    ctfModel.setDefocusAngle(-999)


def readCtfModel(ctfModel, filename):
    result = parseGctfOutput(filename)
    if result is None:
        setWrongDefocus(ctfModel)
        ctfFit, ctfResolution, ctfPhaseShift = -999, -999, 0
    else:
        defocusU, defocusV, defocusAngle, ctfFit, ctfPhaseShift, ctfResolution = result
        ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)
    ctfModel.setFitQuality(ctfFit)
    ctfModel.setResolution(ctfResolution)

    # Avoid creation of phaseShift
    if ctfPhaseShift != 0:
        ctfModel.setPhaseShift(ctfPhaseShift)
