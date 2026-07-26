cdef class BiellipsoidWrap:
    '''A wrapper class for the C++ class Biellipsoid, with some added convenience functions.

    The biellipsoid consists of two tri-axial ellipsoids which share radii along the y and z axis
    but not necessarily the x axis, which are joined along the y-z plane, rotated arbitrarily, and then
    shifted to some position.

    For these calculations, three reference frames are used -- the conversion functions are exposed
    fmainly for debugging purposes.  The frames are:
        :View Frame: This is the native frame users will be operating in. It is easiest to think of
            +x as being "right", +y as "up", and +z as "towards the viewer".
        :Aligned Frame: This is the frame in which the forward vector points in the +x,
            the up vector points in the +y, and the side vector points in the +z directions.  The
            rotation matrix of the biellipsoid rotates points from this frame to the view frame.
        :Sphere Frame: In this frame, the biellipsoid is a unit sphere.  This is a non-linear
            transform from the aligned frame, scaling y and z by r_up and r_side respectively,
            and x by r_forward or r_backward depending on whether x > 0.  Many geometric calculations
            amount to transforming to this frame, finding some point, line, or plane, and then
            transforming back to the view space.
    '''

    @property
    def r_forward(self):
        '''The forward radius of the object -- along the +x direction for an identity rotation matrix.'''
        return self.bell.r_forward

    @property
    def r_back(self):
        '''The forward radius of the object -- along the -x direction for an identity rotation matrix.'''
        return self.bell.r_back

    @property
    def r_up(self):
        '''The upward radius of the object -- along the y axis for an identity rotation matrix.'''
        return self.bell.r_up

    @property
    def r_side(self):
        '''The upward radius of the object -- along the z axis for an identity rotation matrix.'''
        return self.bell.r_side

    @property
    def position(self):
        '''The position of the center of the biellipsoid, as a numpy array.'''
        return np.array([self.bell.position.x, self.bell.position.y, self.bell.position.z])

    @property
    def radii(self):
        '''The forward, backward, up, and side radii of the bielipsoid, as a numpy array.'''
        return np.array([self.bell.r_forward, self.bell.r_back, self.bell.r_up, self.bell.r_side])

    @property
    def rot(self):
        '''The rotation matrix of the biellipsoid, a 3x3 numpy array, which transforms points
            from aligned to view space.'''
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
        '''The biellipsoid's up (along the poles) vector in view space.'''
        return np.array([self.bell.rot.xy, self.bell.rot.yy, self.bell.rot.zy])

    @property
    def side(self):
        '''The biellipsoid's side (forward cross up) vector in view space.'''
        return np.array([self.bell.rot.xz, self.bell.rot.yz, self.bell.rot.zz])

    @property
    def limb_forward(self):
        '''The ellipse defining the limb (visible edge) of the Biellipsoid in the forward direction.'''
        cdef Ellipse ell = self.bell.f_limb
        return EllipseWrap(
            np.array([ell.e1.x, ell.e1.y, ell.e1.z]),
            np.array([ell.e2.x, ell.e2.y, ell.e2.z]),
        )

    @property
    def limb_back(self):
        '''The ellipse defining the limb (visible edge) of the Biellipsoid in the backwards direction.'''
        cdef Ellipse ell = self.bell.b_limb
        return EllipseWrap(
            np.array([ell.e1.x, ell.e1.y, ell.e1.z]),
            np.array([ell.e2.x, ell.e2.y, ell.e2.z]),
        )

    @property
    def joint(self):
        '''The ellipse around the surface of the planet within the plane that joins the two ellipsoids
        that make up the biellipsoid.'''
        cdef Ellipse split = self.bell.joint
        return EllipseWrap(
            np.array([split.e1.x, split.e1.y, split.e1.z]),
            np.array([split.e2.x, split.e2.y, split.e2.z]),
        )

    def __init__(self, double r_forward, double r_back, double r_up, double r_side, double theta=0, double phi=0,
                 double gamma=0, double ci=-2, double x=0, double y=0, double z=0):
        '''Define a biellipsoid in terms of it's origin position, rotation angles, and radii.

        Args:
            r_forward: The forward radius of the object -- along the +x direction for an identity rotation matrix.
            r_back: The forward radius of the object -- along the -x direction for an identity rotation matrix.
            r_up: The upward radius of the object -- along the y axis for an identity rotation matrix.
            r_side: The upward radius of the object -- along the z axis for an identity rotation matrix.
            theta: The counter-clockwise rotation of the biellipsoid around the z axis.
            phi: The counter-clockwise rotation of the biellipsoid around the y axis.
            gamma: The counter-clockwise rotation of the biellipsoid around the x axis.
            ci: The cosine of the inclination of the orbit.  This is used only to reorient the object
                along its orbit at the current position.  If the object should not be rotated with its
                current position, set as -2.
            x: The x position of the biellipsoid, in view space.
            y: The y position of the biellipsoid, in view space.
            z: The z position of the biellipsoid, in view space.
        '''
        self.bell = Biellipsoid(r_forward, r_back, r_up, r_side)
        cdef Vec3 loc = Vec3(x, y, z)
        self.bell.set_position(loc)
        self.bell.set_rotation(theta, phi, gamma, ci)

    def set_rotation(self, double theta, double phi, double gamma, double ci=-2):
        '''Sets the rotation of the biellipsoid to the new given values. The first three arguments are Euler angles
            around the z, y, and x axes respectively, applied in reverse order.  If the final argument is between 0
            and 1, the biellipsoid will further be rotated along its orbit.

        Args:
            theta: The counter-clockwise rotation of the biellipsoid around the z axis.
            phi: The counter-clockwise rotation of the biellipsoid around the y axis.
            gamma: The counter-clockwise rotation of the biellipsoid around the x axis.
            ci: The cosine of the inclination of the orbit.  This is used only to reorient the object
                along its orbit at the current position.  If the object should not be rotated with its
                current position, set as -2.'''
        self.bell.set_rotation(theta, phi, gamma, ci)

    def set_position(self, double x, double y, double z):
        '''Set the position of the biellipsoid.  Note that this will not update the rotation of the object with its
        orbit; if this is desired, users should call set_rotation next.'''
        cdef Vec3 position = Vec3(x, y, z)
        self.bell.set_position(position)

    def is_visible(self, double x, double y, double z):
        '''Determines whether a the given point is 'visible' from the view plane; that is, whether a straight ray from
        the given point extending in the -z direction does NOT pass through the biellipsoid.'''
        cdef Vec3 loc = Vec3(x, y, z)
        return self.bell.is_visible(loc)

    def is_forward(self, double x, double y, double z):
        '''Determines whether a point is on the forward half of the biellipsoid or the backwards half.  The point does
        not need to be on the surface of the biellipsoid, as it will be compared against the joining plane.'''
        cdef Vec3 loc = Vec3(x, y, z)
        return self.bell.is_forward(loc)

    def is_forward_2d(self, double x, double y, bint local=0):
        '''Determines whether the point closest to the viewer (lowest z) on the surface of the biellipsoid at [x, y]
            is on the forward half of the biellipsoid or the backwards half. This is similar to `is_forward` but is
            faster than using raycast to find the intersection point and subsequently is_forward.

            Returns false if there are no points on the Biellipsoid at [x, y].'''
        return self.bell.is_forward_2d(x, y, local)

    def x_bounds(self):
        '''Returns the maximum and minimum extent of the biellipsoid in the x direction.'''
        cdef Bounds b = self.bell.x_bounds()
        return np.array([b.min, b.max])

    def y_bounds(self):
        '''Returns the maximum and minimum extent of the biellipsoid in the y direction.'''
        cdef Bounds b = self.bell.y_bounds()
        return np.array([b.min, b.max])

    def slice_ylimits(self, double x):
        '''For a given x, returns the largest and smallest values of y for which there is some z such that [x, y, z]
        is on the biellipsoid.  This is useful for computing integrals in view-space.'''
        cdef Bounds b = self.bell.slice_ylimits(x)
        return np.array([b.min, b.max])

    def line_intersects(self, double x, double y):
        '''Determines whether the line passing through [x, y, 0] in the +/-z direction passes through the biellipsoid.'''
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
