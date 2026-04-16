import numpy as np


cdef double darkening(double rSq, int limbType, double limb0, double limb1, double limb2, double limb3) noexcept:
    cdef double x
    if limbType == 0:
        # Quadratic; x = 1 - mu
        x = 1. - math.sqrt(max(1. - rSq, 0.))
        return (1 - limb0*x - limb1*x*x)
    elif limbType == 1:
        # Nonlinear; x = math.sqrt(mu)
        x = math.sqrt(math.sqrt(1 - rSq))
        return (1. - limb0*(1. - x) - limb1*(1. - x**2) - limb2*(1. - x**3) - limb3*(1. - x**4))
    # Invalid limbType
    return nan


cdef double darkening_normalization(int limbType, double limb0, double limb1, double limb2, double limb3) noexcept:
    if limbType == 0:
        return (1. - limb0/3. - limb1/6.)*pi
    elif limbType == 1:
        return (-limb0/5. - limb1/3. - 3.*limb2/7. - limb3/2. + 1.)*pi
    # Invalid limbType
    return nan


cpdef object prolate_transit(double rPrograde, double rSubstellar, double rPolar, double[:] t, double t0, double period, double semimajor, double inclination, str limbType, object limb):
    cdef double xp, yp
    cdef int i = 0
    cdef int nTimes = t.shape[0]
    output = np.empty_like(t)
    cdef double[:] outputView = output
    cdef double result
    cdef double theta

    # Check inputs
    cdef int limbCode
    cdef double[:] limbParams = np.zeros(5)
    limbType = limbType.lower().strip()
    if limbType == "quadratic":
        assert len(limb) == 2, "Error: Quadratic limb darkening takes exactly two parameters."
        limbCode = 0
        for i in range(2):
            limbParams[i] = limb[i]
    elif limbType == "nonlinear":
        assert len(limb) == 4, "Error: Nonlinear limb darkening takes exactly two parameters."
        limbCode = 1
        for i in range(4):
            limbParams[i] = limb[i]
    else:
        raise ValueError("limbType not recognized, must be one of quadratic, nonlinear.")

    # Do not crash the program due to lack of precision
    gsl_set_error_handler_off()

    for i in range(nTimes):
        theta = 2*pi * (t[i]-t0) / period

        if abs(theta) > pi/2.:
            # No transit if the planet is behind the star
            outputView[i] = 1.
            continue

        # Calculate the planet's position in the sky (aligned with orbit) frame
        majorLen, minorLen, xp, yp = orbit_geometry(rPrograde, rSubstellar, rPolar, semimajor, theta, inclination)
        # xe, ye, ze = orbit_to_position(t[i], semimajor, period, eccen, inclination, lonPeriapse)

        # Axis-aligned bounding box check to rule out trivial non-transits
        if (xp + majorLen < -1.) or (xp - majorLen > 1) or (yp + minorLen < -1) or (yp - minorLen > 1):
            outputView[i] = 1.
            continue

        result = brute_integrate(majorLen, minorLen, xp, yp, limbParams, limbType=limbCode, limitsMode=0)
        outputView[i] = 1. - result / darkening_normalization(limbCode, limbParams[0], limbParams[1], limbParams[2], limbParams[3])
    return output


cpdef object asymmetric_transit(double rMorning, double rEvening, double rPole, double[:] t, double t0, double period, double semimajor, double inclination, str limbType, object limb, double eccen=0, double lonPeriapse=90.):
    '''Calculates the transit of a piecewise-elliptical planet.  Assumes the same projected shape regardless of its position
    in the orbit.  Uses the same model as catwoman (two spheres split down the middle) if rPole is negative.

    Parameters:
        rMorning        The radius of the planet at the morning (right) side of the planet at the equator relative to
                            the stellar radius.
        rEvening        The radius of the planet at the evening (left) side of the planet at the equator relative top
                            the stellar radius.
        rPole           The radius of the planet at the poles (top and bottom) relative to the stellar radius.  If -1,
                            rPole is set to rMorning one morning side and rEvening on the evening side.
        t               A 1-d Numpy array of observation times to be simulated (must be the same unit as t0 and period).
        t0              The time of a mid-transit (needn't be the observed one).
        period          The orbital period of the planet.
        semimajor       The semimajor axis of the planet relative to the stellar radius.
        inclination     The inclination of the orbit in degrees (near 90 for transiting planets).
        limbType        The type of limb darkening to use: 'quadratic' or 'nonlinear'.
        limb            An iterable of limb darkening parameters, length 2 for quadratic and 4 for nonlinear.
        eccen           The orbital eccentricity of the planet.  Must be 0 <= e < 1, defaults to zero.
        lonPeriapse     The longitude of the periapse of the planet's orbit, in degrees.  The default of 90 means the pariapse occurs at mid-transit.

    Returns:
        The relative flux from the star at the time given, so 1 if the planet is out of transit.
    '''
    cdef double xe, ye, ze
    cdef int i = 0
    cdef int nTimes = t.shape[0]
    output = np.empty_like(t)
    cdef double[:] outputView = output
    cdef double result

    # Check inputs
    cdef int limbCode
    cdef double[:] limbParams = np.zeros(5)
    limbType = limbType.lower().strip()
    if limbType == "quadratic":
        assert len(limb) == 2, "Error: Quadratic limb darkening takes exactly two parameters."
        limbCode = 0
        for i in range(2):
            limbParams[i] = limb[i]
    elif limbType == "nonlinear":
        assert len(limb) == 4, "Error: Nonlinear limb darkening takes exactly two parameters."
        limbCode = 1
        for i in range(4):
            limbParams[i] = limb[i]
    else:
        raise ValueError("limbType not recognized, must be one of quadratic, nonlinear.")
    assert eccen >= 0. and eccen < 1.

    # Handle negative rPole (asymmetric pole radius locked to circles)
    cdef double rPoleMorning = rPole
    cdef double rPoleEvening = rPole
    if rPole < 0:
        rPoleMorning = rMorning
        rPoleEvening = rEvening

    # Do not crash the program due to lack of precision
    gsl_set_error_handler_off()

    for i in range(nTimes):
        # Calculate the planet's position in the sky (aligned with orbit) frame
        xe, ye, ze = orbit_to_position(t[i], semimajor, period, eccen, inclination, lonPeriapse)

        # No transit if the planet is behind the star
        if ze < 0:
            outputView[i] = 1.
            continue

        # Axis-aligned bounding box check to rule out trivial non-transits
        if (xe + rMorning < -1.) or (xe - rEvening > 1) or \
                (ye + max(rPoleMorning, rPoleEvening) < -1) or (ye - max(rPoleMorning, rPoleEvening) > 1):
            outputView[i] = 1.
            continue

        result = 0.
        if xe > -1.:
            result += brute_integrate(rEvening, rPoleEvening, xe, ye, limbParams, limbType=limbCode, limitsMode=1)
        if xe < 1.:
            result += brute_integrate(rMorning, rPoleMorning, xe, ye, limbParams, limbType=limbCode, limitsMode=2)
        outputView[i] = 1. - result / darkening_normalization(limbCode, limbParams[0], limbParams[1], limbParams[2], limbParams[3])
    return output


cdef double brute_integrand(double y, void* params) noexcept:
    cdef BruteIntegralParams* g = <BruteIntegralParams*>params
    cdef double rsq = g.x*g.x + y*y
    return darkening(rsq, g.limbType, g.limb0, g.limb1, g.limb2, g.limb3)


cdef double brute_integrateY(double x, void* params) noexcept:
    cdef double result, err
    cdef BruteIntegralParams* g = <BruteIntegralParams*>params
    g.x = x

    # Calculate the bounds of integration in y
    cdef double yMin = (x-g.xe)/g.a
    yMin = g.b * math.sqrt(1.-yMin*yMin)
    cdef double yMax = g.ye + yMin
    yMin = g.ye - yMin
    # Clip to within the bounds of the star
    cdef double yStarEdge = math.sqrt(1 - x*x)
    yMin = max(yMin, -yStarEdge)
    yMax = min(yMax, yStarEdge)
    if yMin >= yMax:
        return 0.

    cdef int code = gsl_integration_qag(g.integrand, yMin, yMax, 1e-7, 1e-7, 100, 1, g.work, &result, &err)
    if code != 0:
        return nan
    return result


cpdef double brute_integrate(double a, double b, double xe, double ye, double[:] limb, int limbType=0, int limitsMode=0) except -999.:
    cdef double result, err
    cdef int code = 0
    cdef BruteIntegralParams g = BruteIntegralParams(a, b, xe, ye, 0., limbType, limb[0], limb[1], limb[2], limb[3], NULL, NULL)

    # Prepare the inner (y) integral variables
    cdef gsl_integration_workspace* workspaceInner = gsl_integration_workspace_alloc(100)
    cdef gsl_function integInner
    integInner.function = &brute_integrand
    integInner.params = &g
    g.work = workspaceInner
    g.integrand = &integInner

    # Now prepare the outer (x) integral variables
    cdef gsl_integration_workspace* workspaceOuter = gsl_integration_workspace_alloc(100)
    cdef gsl_function integOuter
    integOuter.function = &brute_integrateY
    integOuter.params = &g

    # Compute limits of integration on the x axis
    cdef double x0, x1
    if limitsMode == 0:     # Integrate both sides of the planet
        x0 = max(xe-a, -1)
        x1 = min(xe+a, 1)
    elif limitsMode == 1:   # Integrate only the evening (left) side of the planet
        x0 = max(xe-a, -1)
        x1 = min(xe, 1)
    elif limitsMode == 2:   # Integrate only the morning (right) side of the planet
        x0 = max(xe, -1)
        x1 = min(xe+a, 1)
    if x0 >= x1:
        result = 0.
    else:
        code = gsl_integration_qag(&integOuter, x0, x1, 1e-7, 1e-7, 100, 1, workspaceOuter, &result, &err)

    gsl_integration_workspace_free(workspaceInner)
    gsl_integration_workspace_free(workspaceOuter)

    if code != 0:
        raise RuntimeError("Integration failed in brute_integrate().  Code: " + str(code) + ", message: " + bytes(gsl_strerror(code)).decode("ascii"),
                           (a, b, xe, ye, np.asarray(limb), limbType, limitsMode))

    return result


cpdef (double, double, double, double) orbit_geometry(double a, double b, double c, double semimajor, double theta, double phi):
    ''' a, b, and c are the planet radii relative to the stellar radius
        semimajor is the semimajor axis / stellar radius
        Theta is 2 * pi * time / orbital period
        phi is inclination'''
    # Calculate the sine and cosine of the angles: we'll use them repeatedly
    cdef SpiceDouble ct = math.cos(theta)
    cdef SpiceDouble st = math.sin(theta)
    phi = (90 - phi)*pi/180
    cdef SpiceDouble cp = math.cos(phi)
    cdef SpiceDouble sp = math.sin(phi)

    # Calculate the view angle in the planet's frame
    cdef SpiceDouble[3] view = [ct*cp, -st*cp, sp]

    # Generate the plane of the planet limb.
    cdef SpiceDouble[3] limbNormal = [view[0]/a/a, view[1]/b/b, view[2]/c/c]
    cdef SpicePlane limbPlane
    nvc2pl_c(limbNormal, 0., &limbPlane)

    # Get the ellipse of the limb from the ellipsoid and plane.
    cdef SpiceEllipse limbEllipse
    cdef SpiceBoolean found
    inedpl_c(a, b, c, &limbPlane, &limbEllipse, &found)
    if not found:
        print("ERROR: Failed to find intersection between the ellipse and limb plane.")

    # Project the limb ellipse onto the plane perpendicular to the view vector.
    cdef SpicePlane viewPlane
    cdef SpiceEllipse projEllipse
    nvc2pl_c(view, 0., &viewPlane)
    pjelpl_c(&limbEllipse, &viewPlane, &projEllipse)

    # Get the major and minor axis vectors of the ellipse
    cdef SpiceDouble[3] center
    cdef SpiceDouble[3] gv1
    cdef SpiceDouble[3] gv2
    el2cgv_c(&projEllipse, center, gv1, gv2)

    cdef SpiceDouble[3] majorAx
    cdef SpiceDouble[3] minorAx
    saelgv_c(gv1, gv2, majorAx, minorAx)

    cdef double majorLen = math.sqrt(majorAx[0]*majorAx[0] + majorAx[1]*majorAx[1] + majorAx[2]*majorAx[2])
    cdef double minorLen = math.sqrt(minorAx[0]*minorAx[0] + minorAx[1]*minorAx[1] + minorAx[2]*minorAx[2])

    # print("=== Planet Frame ===")
    if abs(majorAx[1]*view[1]+majorAx[2]*view[2] + majorAx[0]*view[0]) > 1e-14 or \
       abs(minorAx[1]*view[1]+minorAx[2]*view[2] + minorAx[0]*view[0]) > 1e-14:
        print("ERROR: projected axes were not perpendicular to the view vector!")
    # print("View         ", view)
    # print("MajorAx      ", majorAx)
    # print("MinorAx      ", minorAx)

    # Rotate the major axis into the unaligned view frame
    majorAx = [
        -st*majorAx[0]      - ct*majorAx[1],                        # no-cython-lint
        -ct*sp*majorAx[0]   + st*sp*majorAx[1]  + cp*majorAx[2],    # no-cython-lint
        0.  # Zero by construction: -ct*cp*majorAx[0]    + st*cp*majorAx[1]  - sp*majorAx[2]
    ]
    if majorAx[0] < 0:
        majorAx[0] = -majorAx[0]
        majorAx[1] = -majorAx[1]

    # Calculate the planet's position in the unaligned view frame
    cdef double xpu = semimajor * st
    cdef double ypu = semimajor * ct * sp

    # print("=== Unaligned View Frame ===")
    # print("MajorAx      ", majorAx)
    # print("Planet       ", xpu, ypu)

    # Rotate the planet's position into the aligned frame
    cdef double xp = (xpu*majorAx[0] + ypu*majorAx[1])/majorLen
    cdef double yp = (-xpu*majorAx[1] + ypu*majorAx[0])/majorLen

    # Write out the results and return
    return majorLen, minorLen, xp, yp


cdef class Biellipsoid:
    def __init__(self, double x, double y, double z, double theta, double phi, double gamma,
                 double r_forward, double r_back, double r_up, double r_side):
        self.position = [x, y, z]
        self.scale_f = [r_forward, r_up, r_side]
        self.scale_b = [r_back, r_up, r_side]
        cdef int i
        theta = theta * pi / 180
        phi = phi * pi / 180
        gamma = gamma * pi / 180
        cdef double ct = math.cos(theta)
        cdef double st = math.sin(theta)
        cdef double cp = math.cos(phi)
        cdef double sp = math.sin(phi)
        cdef double cg = math.cos(gamma)
        cdef double sg = math.sin(gamma)
        self.rot = [
            cp*ct,  -cp*st*cg + sp*sg, cp*st*sg + sp*cg,
            st,     ct*cg,             -ct*sg,
            -sp*ct, sp*st*cg + cp*sg,  -sp*st*sg + cp*cg
        ]

        # Calculate limb vectors
        cdef double[3] plane = [self.rot[2] / r_forward, self.rot[5] / r_up, self.rot[8] / r_side]
        normalize3(plane)
        self.f1 = [1., 0., 0.]
        self.f2 = [0., 1., 0.]
        # TODO: Fix wrong condition?  Also below
        if (plane[0]*plane[0] + plane[1]*plane[1] > 1e-12):
            self.f1 = [plane[1], -plane[0], 0.]
            normalize3(self.f1)
            cross(plane, self.f1, out=self.f2)
        # Calculate limb vectors
        plane = [self.rot[2] / r_back, self.rot[5] / r_up, self.rot[8] / r_side]
        normalize3(plane)
        self.b1 = [1., 0., 0.]
        self.b2 = [0., 1., 0.]
        if (plane[0]*plane[0] + plane[1]*plane[1] > 1e-12):
            self.b1 = [plane[1], -plane[0], 0.]
            normalize3(self.b1)
            cross(plane, self.b1, out=self.b2)

        # Transform back to view space
        for i in range(3):
            self.f1[i] *= self.scale_f[i]
            self.f2[i] *= self.scale_f[i]
            self.b1[i] *= self.scale_b[i]
            self.b2[i] *= self.scale_b[i]
        matmul3x3(self.rot, self.f1, self.f1)
        matmul3x3(self.rot, self.f2, self.f2)
        matmul3x3(self.rot, self.b1, self.b1)
        matmul3x3(self.rot, self.b2, self.b2)

    cpdef double get_half_width(self, int axis=0, int side=0) noexcept:
        if side == 0:
            return math.sqrt(self.f1[axis]*self.f1[axis] + self.f2[axis]*self.f2[axis])
        else:
            return math.sqrt(self.b1[axis]*self.b1[axis] + self.b2[axis]*self.b2[axis])

    cpdef (double, double) get_ybounds(self, double x) noexcept:
        cdef int i
        cdef double[3] xf1, xf2, xb1, xb2

        cdef double xwidth2 = self.f1[0]*self.f1[0] + self.f2[0]*self.f2[0]
        cdef double disc = self.f1[0]**4 - self.f1[0]**2 * x**2 + self.f1[0]**2 * self.f2[0]**2
        cdef double B = (x * self.f2[0] + math.sqrt(disc)) / xwidth2
        cdef double A = (x - self.f2[0] * B) / self.f1[0]
        for i in range(3):
            xf1[i] = A * self.f1[i] + B * self.f2[i]
        B = (x * self.f2[0] - math.sqrt(disc)) / xwidth2
        A = (x - self.f2[0] * B) / self.f1[0]
        for i in range(3):
            xf2[i] = A * self.f1[i] + B * self.f2[i]
        B = (x * self.f2[0] - math.sqrt(disc)) / xwidth2
        A = (x - self.f2[0] * B) / self.f1[0]
        for i in range(3):
            xb1[i] = A * self.f1[i] + B * self.f2[i]
        B = (x * self.f2[0] + math.sqrt(disc)) / xwidth2
        A = (x - self.f2[0] * B) / self.f1[0]
        for i in range(3):
            xb2[i] = A * self.f1[i] + B * self.f2[i]

        # Select bounds based on their location
        cdef double ymin=0.
        cdef double ymax=0.
        cdef double[3] forward = [self.rot[0], self.rot[3], self.rot[5]]
        if dot3(xf1, forward) >= 0:
            if xf1[1] < ymin:
                ymin = xf1[1]
            elif xf1[1] > ymax:
                ymax = xf1[1]
        if dot3(xf2, forward) >= 0:
            if xf2[1] < ymin:
                ymin = xf2[1]
            elif xf2[1] > ymax:
                ymax = xf2[1]
        if dot3(xb1, forward) < 0:
            if xb1[1] < ymin:
                ymin = xb1[1]
            elif xb1[1] > ymax:
                ymax = xb1[1]
        if dot3(xb2, forward) < 0:
            if xb2[1] < ymin:
                ymin = xb2[1]
            elif xb2[1] > ymax:
                ymax = xb2[1]
        return ymin, ymax


cpdef double solve_kepler(double mean_anomaly, double eccen):
    cdef double e_anom = mean_anomaly
    cdef double e_new = mean_anomaly
    cdef int i = 0
    if eccen > 0.95:
        # Direct iteration
        for i in range(50):
            e_new = mean_anomaly + eccen * math.sin(e_anom)
            if abs(e_new - e_anom) < 1e-12:
                break
            e_anom = e_new
    if eccen > 0.:
        # Newton-Raphson
        for i in range(50):
            e_new -= (e_anom - eccen*math.sin(e_anom) - mean_anomaly) / (1 - eccen * math.cos(e_anom))
            if abs(e_new - e_anom) < 1e-12:
                break
            e_anom = e_new
    return e_new


cpdef (double, double, double) orbit_to_position(double t, double semimajor, double period, double eccen, double inclination, double lon_periapse):
    # Notation: t = time, ta = true anomaly, ea = eccentric anomaly, ma = mean anomaly, x,y,z = position
    cdef double ta_c = (lon_periapse-90)*pi/180.
    cdef double sta_c = math.sin(ta_c)
    cdef double cta_c = math.cos(ta_c)
    # Get time of conjunction given longitude of periapse
    cdef double ea_c = math.atan2(math.sqrt(1 - eccen*eccen) * sta_c, eccen + cta_c)
    cdef double ma_c = ea_c - eccen * math.sin(ea_c)
    cdef double t_c = ma_c * period / (2*pi)
    # Solve Kepler's equation in the rotated and inclined frame, relative to inferior conjunction
    cdef double ma = 2*pi * (t-t_c) / period
    cdef double ea = solve_kepler(ma, eccen)
    # Position in the rotated frame where lon_periapse = 0, inclination = 0
    cdef double x_p = semimajor * math.sqrt(1 - eccen*eccen) * math.sin(ea)
    cdef double y_p = semimajor * (math.cos(ea) - eccen)
    # Transform back to frame where lon_periapse != 0 (still inclination = 0)
    cdef double x_inc = x_p * cta_c + y_p * sta_c
    cdef double y_inc = x_p * -sta_c + y_p * cta_c
    # Transform back to the view frame (x = right, y = up, z = towards observer) centered on star
    # Longitude of ascending node is 270 degrees by construction, so x = x_inc
    inclination = inclination*pi/180
    cdef double y = y_inc * math.cos(inclination)
    cdef double z = y_inc * math.sin(inclination)
    return x_inc, y, z


# ========== Helpers ==========
cdef inline double normalize3(double[3] vec) noexcept:
    cdef double length = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]
    length = math.sqrt(length)
    vec[0] /= length
    vec[1] /= length
    vec[2] /= length
    return length


cdef inline double dot3(double[3] a, double[3] b) noexcept:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


cdef inline void cross(double[3] a, double[3] b, double[3] out) noexcept:
    out[0] = a[1]*b[2] - b[1]*a[2]
    out[1] = a[2]*b[0] - b[2]*a[0]
    out[2] = a[0]*b[1] - b[0]*a[1]
    return


cdef inline void matmul3x3(double[9] m, double[3] v, double[3] out) noexcept:
    # Separated in case v is out
    cdef double[3] result
    result[0] = m[0]*v[0] + m[1]*v[1] + m[2]*v[2]
    result[1] = m[3]*v[0] + m[4]*v[1] + m[5]*v[2]
    result[2] = m[6]*v[0] + m[7]*v[1] + m[8]*v[2]
    out[0] = result[0]
    out[1] = result[1]
    out[2] = result[2]
