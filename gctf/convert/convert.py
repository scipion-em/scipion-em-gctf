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
import numpy
import logging
logger = logging.getLogger(__name__)

import pyworkflow.utils as pwutils
from pwem.constants import ALIGN_2D, ALIGN_3D, ALIGN_PROJ, ALIGN_NONE
from pwem.convert.transformations import translation_from_matrix


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


class CoordinatesWriter:
    """ Simple class to write coordinates (in star file as in Relion). """
    HEADER = """
data_

loop_
_rlnCoordinateX #1
_rlnCoordinateY #2
"""

    def __init__(self, filename):
        """ Filename where to write the coordinates. """
        pwutils.makePath(os.path.dirname(filename))  # Ensure path exists
        self._f = open(filename, 'w')
        self._f.write(self.HEADER)

    def writeCoord(self, pos):
        self._f.write("%d %d\n" % (pos[0], pos[1]))

    def close(self):
        self._f.close()


def rowToCtfModel(row, ctf):
    """ Create a CTFModel from a row of a meta """
    ctf.setDefocusU(row.rlnDefocusU)
    ctf.setDefocusV(row.rlnDefocusV)
    ctf.setDefocusAngle(row.rlnDefocusAngle)
    ctf.setResolution(row.get('rlnCtfMaxResolution', 0))
    ctf.setFitQuality(row.get('rlnCtfFigureOfMerit', 0))

    if hasattr(row, 'rlnPhaseShift'):
        ctf.setPhaseShift(row.rlnPhaseShift)
    ctf.standardize()

    return ctf


def getShifts(transform, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
                          -> for xmipp implies alignment
    """
    if alignType == ALIGN_NONE:
        return None

    inverseTransform = alignType == ALIGN_PROJ
    # only flip is meaningful if 2D case
    # in that case the 2x2 determinant is negative
    matrix = transform.getMatrix()
    if alignType == ALIGN_2D:
        # get 2x2 matrix and check if negative
        flip = bool(numpy.linalg.det(matrix[0:2, 0:2]) < 0)
        if flip:
            matrix[0, :2] *= -1.  # invert only the first two columns keep x
            matrix[2, 2] = 1.  # set 3D rot

    elif alignType == ALIGN_3D:
        flip = bool(numpy.linalg.det(matrix[0:3, 0:3]) < 0)
        if flip:
            matrix[0, :4] *= -1.  # now, invert first line including x
            matrix[3, 3] = 1.  # set 3D rot

    shifts = geometryFromMatrix(matrix, inverseTransform)

    return shifts


def geometryFromMatrix(matrix, inverseTransform):

    if inverseTransform:
        matrix = numpy.linalg.inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    return shifts
