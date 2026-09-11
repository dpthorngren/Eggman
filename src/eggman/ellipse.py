# flake8: noqa: E302
from typing import Optional

import cython
import numpy as np
from cython.cimports.eggman import cy_eggman as cye  # type: ignore


@cython.cclass
class Ellipse:
    '''An ellipse centered at the origin, defined by two conjugate diameters e1 and e1.'''

    @property
    def e1(self):
        '''The first of the two conjugate diameters that define the ellipse.'''
        return np.array([self.cell.e1.x, self.cell.e1.y, self.cell.e1.z])

    @property
    def e2(self):
        '''The second of the two conjugate diameters that define the ellipse.'''
        return np.array([self.cell.e2.x, self.cell.e2.y, self.cell.e2.z])

    @property
    def det(self):
        '''The determinant of the matrix of conjugate diameters (with the third dimension dropped).'''
        return self.cell.det

    @property
    def x_size(self):
        '''The on-sky half-width of the ellipse.'''
        return self.cell.x_size

    @property
    def y_size(self):
        '''The on-sky half-height of the ellipse.'''
        return self.cell.y_size

    def __init__(self, e1, e2):
        '''Define an ellipse in terms of its conjugate diameters.

        Args:
            e1: The first conjugate diameter (a 1-d array of length 3).
            e2: The second conjugate diameter (a 1-d array of length 3).
        '''
        cython.declare(ve1=cye.Vec3, ve2=cye.Vec3)
        ve1 = cye.Vec3(e1[0], e1[1], e1[2])
        ve2 = cye.Vec3(e2[0], e2[1], e2[2])
        self.cell = cye.CEllipse(ve1, ve2)

    @classmethod
    def create_from_rot_radii(cls, a: float, b: float, rot: np.ndarray):
        '''Define an ellipse in terms of a rotation matrix and its semimajor and semiminor axes.

        Args:
            a: The semimajor axis of the ellipse.
            b: The semiminor axis of the ellipse.
            rot: A 3x3 rotation matrix to apply to the Ellipse.
        '''
        assert a > 0 and b > 0
        assert type(rot) is np.ndarray and rot.shape == (3, 3)
        assert np.all(rot.T @ rot - np.eye(3) < 1e-9), "Invalid rotation matrix."
        e1 = (rot @ np.vstack([a, 0, 0])).flat
        e2 = (rot @ np.vstack([0, b, 0])).flat
        return Ellipse(e1, e2)

    def get_ybounds(self, x: float):
        '''The upper and lower bounding points of the intersection of the ellipse and the plane
        defined by the given x.

        Args:
            x: The x value to find the y-bounds for.

        Returns:
            A 2x3 array representing the min and max bounds of the intersection.  If there is
            no intsersection, the array is filled with nans.'''
        out_min: cye.Vec3 = cye.Vec3()
        out_max: cye.Vec3 = cye.Vec3()
        self.cell.get_ybounds(x, out_min, out_max)
        return np.array([
            [out_min.x, out_min.y, out_min.z],
            [out_max.x, out_max.y, out_max.z],
        ])

    def line_intersects(self, x: float, y: float) -> tuple[cython.bint, np.ndarray]:
        '''Checks whether a line across z at (x, y) passes through the ellipse.

        Args:
            x: the x position of the line to check for an intersection.
            y: the y position of the line to check for an intersection.

        Returns: A tuple with two components, a bool indicating whether the line
            intersects the ellipsoid and a length-3 Numpy array representing
            the point of intersection (nans if there wasn't one).
        '''
        result: cye.Vec3 = cye.Vec3()
        hit = self.cell.line_intersects(x, y, cython.address(result))
        return hit, np.array([result.x, result.y, result.z])

    def nearest_to_line(self, xt: float, yt: float):
        '''Locate the point on the ellipse nearest to the line xt, yt.'''
        result: cye.Vec3 = self.cell.nearest_to_line(xt, yt)
        return np.array([result.x, result.y, result.z])

    def outline(self, res=200, dir: Optional[np.ndarray] = None, method: str = "angles"):
        '''Generates an array of points along the ellipse.'''
        if method == "angles":
            # TODO: Set start/end points at break plane
            tmin = 0
            tmax = 2 * np.pi
            t = np.linspace(tmin, tmax, res)[:, None]
            xyz = self.e1 * np.cos(t) + self.e2 * np.sin(t)
        elif method == "ybounds":
            x = np.linspace(-self.x_size, self.x_size, res)
            ylow = []
            yhigh = []
            for xi in x:
                y = self.get_ybounds(xi)
                ylow.append(y[0])
                yhigh.append(y[1])
            xyz = np.vstack([yhigh, ylow[-1::-1]])
        else:
            raise ValueError(f"Unrecognized method {method}, must be 'angles' or 'ybounds'")
        if dir is not None:
            return xyz[xyz.dot(dir) >= 0]
        return xyz

    def plot_outline(self, res=200, origin=(0., 0.), dir=None, method="angles", **args):
        '''Plots the ellipse in the x-y plane (ignores z).'''
        from matplotlib import pyplot as plt

        xyz = self.outline(res, dir, method=method)
        plt.plot(xyz[:, 0] + origin[0], xyz[:, 1] + origin[1], **args)

    def plot_vectors(self, origin=(0., 0.), **args):
        '''Plots the conjugate diameters of the ellipse in the x-y plane.'''
        from matplotlib import pyplot as plt

        plt.plot([origin[0], origin[0] + self.e1[0]], [origin[1], origin[1] + self.e1[1]], **args)
        plt.plot([origin[0], origin[0] + self.e2[0]], [origin[1], origin[1] + self.e2[1]], **args)
