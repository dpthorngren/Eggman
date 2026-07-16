cdef class BiellipsoidWrap:
    '''A wrapper class for the C++ class Biellipsoid, with some added convenience functions.'''

    @property
    def r_forward(self):
        return self.bell.r_forward

    @property
    def r_back(self):
        return self.bell.r_back

    @property
    def r_up(self):
        return self.bell.r_up

    @property
    def r_side(self):
        return self.bell.r_side

    @property
    def position(self):
        return np.array([self.bell.position.x, self.bell.position.y, self.bell.position.z])

    @property
    def radii(self):
        '''The forward, backward, up, and side radii of the bielipsoid.'''
        return np.array([self.bell.r_forward, self.bell.r_back, self.bell.r_up, self.bell.r_side])

    @property
    def rot(self):
        '''The rotation matrix of the biellipsoid, transforming from aligned to view space.'''
        return np.array([
            [self.bell.rot.xx, self.bell.rot.xy, self.bell.rot.xz],
            [self.bell.rot.yx, self.bell.rot.yy, self.bell.rot.yz],
            [self.bell.rot.zx, self.bell.rot.zy, self.bell.rot.zz],
        ])

    @property
    def forward(self):
        '''The forward vector in view space.'''
        return np.array([self.bell.rot.xx, self.bell.rot.yx, self.bell.rot.zx])

    @property
    def up(self):
        '''The biellipsoid's up (polar) vector in view space.'''
        return np.array([self.bell.rot.xy, self.bell.rot.yy, self.bell.rot.zy])

    @property
    def side(self):
        '''The biellipsoid's side (forward cross up) vector in view space.'''
        return np.array([self.bell.rot.xz, self.bell.rot.yz, self.bell.rot.zz])

    @property
    def limb_forward(self):
        cdef Ellipse ell = self.bell.f_limb
        return EllipseWrap(
            np.array([ell.e1.x, ell.e1.y, ell.e1.z]),
            np.array([ell.e2.x, ell.e2.y, ell.e2.z]),
        )

    @property
    def limb_back(self):
        cdef Ellipse ell = self.bell.b_limb
        return EllipseWrap(
            np.array([ell.e1.x, ell.e1.y, ell.e1.z]),
            np.array([ell.e2.x, ell.e2.y, ell.e2.z]),
        )

    @property
    def joint(self):
        cdef Ellipse split = self.bell.joint
        return EllipseWrap(
            np.array([split.e1.x, split.e1.y, split.e1.z]),
            np.array([split.e2.x, split.e2.y, split.e2.z]),
        )

    def __init__(self, double x, double y, double z, double theta, double phi, double gamma,
                 double r_forward, double r_back, double r_up, double r_side, double ci=-2):
        self.bell = Biellipsoid(r_forward, r_back, r_up, r_side)
        self.set_position(x, y, z)
        self.set_rotation(theta, phi, gamma, ci)

    def set_rotation(self, double theta, double phi, double gamma, double ci):
        self.bell.set_rotation(theta, phi, gamma, ci)

    def set_position(self, double x, double y, double z):
        cdef Vec3 position = Vec3(x, y, z)
        self.bell.set_position(position)

    def is_visible(self, double x, double y, double z):
        cdef Vec3 loc = Vec3(x, y, z)
        return self.bell.is_visible(loc)

    def is_forward(self, double x, double y, double z):
        cdef Vec3 loc = Vec3(x, y, z)
        return self.bell.is_forward(loc)

    def is_forward_2d(self, double x, double y, bint local=0):
        return self.bell.is_forward_2d(x, y, local)

    def x_bounds(self):
        cdef Bounds b = self.bell.x_bounds()
        return np.array([b.min, b.max])

    def y_bounds(self):
        cdef Bounds b = self.bell.y_bounds()
        return np.array([b.min, b.max])

    def slice_ylimits(self, double x):
        cdef Bounds b = self.bell.slice_ylimits(x)
        return np.array([b.min, b.max])

    def line_intersects(self, double x, double y):
        return self.bell.line_intersects(x, y)

    def raycast(self, double x, double y):
        cdef Vec3 loc
        cdef double mu
        cdef bint hit = self.bell.raycast(x, y, &mu, &loc)
        return hit, np.array([loc.x, loc.y, loc.z]), mu

    def world_to_aligned(self, loc):
        cdef Vec3 v = Vec3(loc[0], loc[1], loc[2])
        v = self.bell.world_to_aligned(v)
        return np.array([v.x, v.y, v.z])

    def world_to_sphere(self, loc):
        cdef Vec3 v = Vec3(loc[0], loc[1], loc[2])
        v = self.bell.world_to_sphere(v)
        return np.array([v.x, v.y, v.z])

    def aligned_to_world(self, loc):
        cdef Vec3 v = Vec3(loc[0], loc[1], loc[2])
        v = self.bell.aligned_to_world(v)
        return np.array([v.x, v.y, v.z])

    def aligned_to_sphere(self, loc):
        cdef Vec3 v = Vec3(loc[0], loc[1], loc[2])
        v = self.bell.aligned_to_sphere(v)
        return np.array([v.x, v.y, v.z])

    def sphere_to_world(self, loc):
        cdef Vec3 v = Vec3(loc[0], loc[1], loc[2])
        v = self.bell.sphere_to_world(v)
        return np.array([v.x, v.y, v.z])

    def sphere_to_aligned(self, loc):
        cdef Vec3 v = Vec3(loc[0], loc[1], loc[2])
        v = self.bell.sphere_to_aligned(v)
        return np.array([v.x, v.y, v.z])

    def principal_axes(self):
        '''Get the principle axes from a biellipsoid, for visualization and testing purposes.'''
        forward = self.rot @ np.diag([self.r_forward, self.r_up, self.r_side])
        back = self.rot @ np.diag([-self.r_back, self.r_up, self.r_side])
        return np.hstack([forward[:, :1], back])

    def meridians(self, axis=0, res=200):
        '''Get the positions of the equator and 90-degree meridians on the biellipsoid.'''
        axes = self.principal_axes()
        forward = axes[:, 0]
        t = np.linspace(0, 2 * np.pi, res)[:, None]
        if axis == 0:
            u = np.cos(t) * axes[:, 2] + np.sin(t) * axes[:, 3]
            v = np.cos(t) * axes[:, 2] + np.sin(t) * axes[:, 3]
        elif axis == 1:
            u = np.cos(t) * axes[:, 0] + np.sin(t) * axes[:, 3]
            v = np.cos(t) * axes[:, 1] + np.sin(t) * axes[:, 3]
        elif axis == 2:
            u = np.cos(t) * axes[:, 0] + np.sin(t) * axes[:, 2]
            v = np.cos(t) * axes[:, 1] + np.sin(t) * axes[:, 2]
        else:
            raise ValueError(f"Axis must be 0, 1, or 2, not {axis}.")
        u = u[u.dot(forward) > 0]
        v = v[v.dot(forward) <= 0]
        result = np.vstack([u, v])
        result = np.vstack([result, result[0]])
        result = np.array(sorted(result, key=lambda x: np.arctan2(x[1], x[0])))
        return result.T + self.position[:, None]

    def outline(self, res=200, yrange=False):
        if yrange:
            xmin, xmax = self.x_bounds()
            x = np.linspace(xmin, xmax, res)
            ymin = []
            ymax = []
            for xi in x:
                bounds = self.slice_ylimits(xi)
                ymin.append(bounds[0])
                ymax.append(bounds[1])
            return x, np.array(ymin), np.array(ymax)
        else:
            xyz1 = self.limb_forward.outline(res, dir=self.forward, method="angles")
            xyz2 = self.limb_back.outline(res, dir=-self.forward, method="angles")
            result = np.vstack([xyz1, xyz2])
            result = sorted(result, key=lambda x: np.arctan2(x[1], x[0]))
            result = np.vstack([result, result[0]])
            return np.array(result).T + self.position[:, None]

    def plot_ellipses(self, res=200, method="angles", f_args=dict(), b_args=dict()):
        '''Similar to plot_outline, but plots each ellipse that makes up the bielipsoid's limb
        separately (including on the other ellipse's side), for debugging visualization purposes.'''
        origin = (self.bell.position.x, self.bell.position.y)
        self.limb_forward.plot_outline(res, origin, dir=self.forward, method=method, **f_args)
        self.limb_back.plot_outline(res, origin, dir=-self.forward, method=method, **b_args)

    def plot_axes(self, colors=None, **args):
        from matplotlib import pyplot as plt
        x0, y0, _ = self.position
        axes = self.principal_axes()
        if colors is None:
            colors = ['red', 'green', 'blue', "orange"]
        for i in range(4):
            plt.plot([x0, x0 + axes[0, i]], [y0, y0 + axes[1, i]], zorder=1, color=colors[i], alpha=.3, **args)

    def plot_outline(self, res=200, **args):
        from matplotlib import pyplot as plt
        x, y, _ = self.outline(res)
        plt.plot(x, y, **args)

    def plot_area(self, res=200, **args):
        from matplotlib import pyplot as plt
        x, ymin, ymax = self.outline(res, True)
        plt.fill_between(x, ymin, ymax, **args)

    def plot_meridians(self, res=200, eq_args=dict(), prime_args=dict(), front_args=dict(), visible_only=True):
        from matplotlib import pyplot as plt

        for i, args in enumerate([prime_args, eq_args, front_args]):
            if args is not None:
                merid = self.meridians(i, res)
                if visible_only:
                    mask = np.array([not self.is_visible(i[0], i[1], i[2]+1e-6) for i in merid.T])
                    mask = np.vstack([mask]*3)
                    merid = np.ma.masked_where(mask, merid)
                plt.plot(merid[0], merid[1], **args)
