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
        return np.array([self.bell.r_forward, self.bell.r_back, self.bell.r_up, self.bell.r_side]),

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
        return np.row_stack([self.bell.f_limb.e1, self.bell.f_limb.e2])

    @property
    def limb_back(self):
        return np.row_stack([self.bell.b_limb.e1, self.bell.b_limb.e2])

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

    def x_bounds(self):
        cdef Bounds b = self.bell.x_bounds()
        return np.array([b.min, b.max])

    def y_bounds(self):
        cdef Bounds b = self.bell.y_bounds()
        return np.array([b.min, b.max])

    def slice_ylimits(self, double x):
        cdef Bounds b = self.bell.slice_ylimits(x)
        return np.array([b.min, b.max])

    def principal_axes(self):
        '''Get the principle axes from a biellipsoid, for visualization and testing purposes.'''
        forward = self.rot @ np.diag([self.r_forward, self.r_up, self.r_side])
        back = self.rot @ np.diag([-self.r_back, self.r_up, self.r_side])
        return np.row_stack([forward[:, 0], back])

    def meridians(self, axis=0, res=200):
        '''Get the positions of the equator and 90-degree meridians on the biellipsoid.'''
        axes = self.principal_axes()
        forward = axes[:, 0]
        t = np.linspace(0, 2 * np.pi, res)[:, None]
        if axis == 0:
            v = np.cos(t) * axes[:, 2] + np.sin(t) * axes[:, 3]
            return v.T + self.position[:, None]
        elif axis == 1:
            u = np.cos(t) * axes[:, 0] + np.sin(t) * axes[:, 3]
            v = np.cos(t) * axes[:, 1] + np.sin(t) * axes[:, 3]
        elif axis == 2:
            u = np.cos(t) * axes[:, 0] + np.sin(t) * axes[:, 2]
            v = np.cos(t) * axes[:, 1] + np.sin(t) * axes[:, 2]
        else:
            raise ValueError(f"Axis must be 0, 1, or 2, not {axis}.")
        u = u[np.dot(u, forward) > 0]
        v = v[np.dot(v, forward) < 0]
        result = np.vstack([u, v])
        result = sorted(result, key=lambda x: np.arctan2(x[1], x[0]))
        result = np.vstack([result, result[0]])
        return np.array(result).T + self.position[:, None]

    def outline(self, res=200):
        '''Get an outline of the biellipsoid at the requested resolution.'''
        assert res > 4 and res % 2 == 0, "Resolution must be even and at least 4."
        xmin, xmax = self.x_bounds()
        x = np.linspace(xmin, xmax, 1 + (res//2))
        ymin = []
        ymax = []
        for xi in x:
            bounds = self.slice_ylimits(xi)
            ymin.append(bounds[0])
            ymax.append(bounds[1])
        # TODO: fix?
        # x = np.concatenate([x, x[-1:1:-1]]).T
        # y = np.concatenate([ymax, ymin[-1:1:-1]]).T
        return x, ymin, ymax

    def plot_axes(self):
        from matplotlib import pyplot as plt
        x0, y0, _ = self.position
        axes = self.principal_axes()
        for i in range(4):
            color = ['red', 'green', 'blue', "orange"][i]
            plt.plot([x0, x0 + axes[i, 0]], [y0, y0 + axes[i, 1]], zorder=1, color=color, alpha=.3)

    def plot_outline(self, res=200):
        from matplotlib import pyplot as plt
        x, y = self.outline(res)

        plt.plot(x, y)
