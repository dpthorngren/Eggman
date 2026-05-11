import numpy as np

cdef class EllipseWrap:
    @property
    def e1(self):
        return np.array([self.ell.e1.x, self.ell.e1.y, self.ell.e1.z])

    @property
    def e2(self):
        return np.array([self.ell.e2.x, self.ell.e2.y, self.ell.e2.z])

    @property
    def x_size(self):
        return self.ell.x_size

    @property
    def y_size(self):
        return self.ell.y_size

    def __init__(self, e1, e2):
        cdef Vec3 ve1 = Vec3(e1[0], e1[1], e1[2])
        cdef Vec3 ve2 = Vec3(e2[0], e2[1], e2[2])
        self.ell = Ellipse(ve1, ve2)

    @classmethod
    def create_from_rot_radii(cls, a, b, rot):
        assert a > 0 and b > 0
        assert type(rot) is np.ndarray and rot.shape == (3, 3)
        assert np.all(rot.T @ rot - np.eye(3) < 1e-9), "Invalid rotation matrix."
        e1 = (rot @ np.row_stack([a, 0, 0])).flat
        e2 = (rot @ np.row_stack([0, b, 0])).flat
        return EllipseWrap(e1, e2)

    def get_ybounds(self, double x):
        cdef Vec3 out_min
        cdef Vec3 out_max
        self.ell.get_ybounds(x, &out_min, &out_max)
        return np.array([
            [out_min.x, out_min.y, out_min.z],
            [out_max.x, out_max.y, out_max.z],
        ])

    def line_intersects(self, x, y):
        return self.ell.line_intersects(x, y)

    def nearest_to_line(self, xt, yt):
        cdef Vec3 result = self.ell.nearest_to_line(xt, yt)
        return np.array([result.x, result.y, result.z])

    def outline(self, res=200, dir=None, method="angles"):
        # TODO: Set start/end points at break plane
        if method == "angles":
            tmin = 0
            tmax = 2*np.pi
            t = np.linspace(tmin, tmax, res)[:, None]
            xyz = self.e1 * np.cos(t) + self.e2*np.sin(t)
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
        from matplotlib import pyplot as plt

        xyz = self.outline(res, dir, method=method)
        plt.plot(xyz[:, 0] + origin[0], xyz[:, 1] + origin[1], **args)

    def plot_vectors(self, origin=(0., 0.), **args):
        from matplotlib import pyplot as plt

        plt.plot([origin[0], origin[0] + self.e1[0]], [origin[1], origin[1] + self.e1[1]], **args)
        plt.plot([origin[0], origin[0] + self.e2[0]], [origin[1], origin[1] + self.e2[1]], **args)
