cdef double darkening(double r, void* params) noexcept:
    cdef IntegralParams* g = <IntegralParams*>params
    cdef double mu = 1. - sqrt(max(1. - r*r, 0.))
    return g.limb[0] - g.limb[1]*mu - g.limb[2]*mu*mu
 

cdef double originDist(double t, void* params) noexcept:
    cdef IntegralParams* g = <IntegralParams*>params
    return (g.a*cos(t)+g.xe)**2 + (g.b*sin(t)+g.ye)**2 - g.rsq


cdef double originDistDiff(double theta, void* params) noexcept:
    # Note: I left out a constant factor of 2 because it doesn't matter
    cdef IntegralParams* g = <IntegralParams*>params
    cdef double ct = cos(theta)
    cdef double st = sin(theta)
    return -g.a*st*(g.xe + g.a*ct) + g.b*ct*(g.ye + g.b*st)


cdef double integrand(double rad, void* params) noexcept:
    cdef IntegralParams* g = <IntegralParams*>params
    g.rsq = rad*rad

    # Find right-side intersection
    gsl_root_fsolver_set(g.solver, g.func, g.tFar, g.tNear+2*pi)
    while (gsl_root_fsolver_x_upper(g.solver) - gsl_root_fsolver_x_lower(g.solver)) > 1e-9:
        gsl_root_fsolver_iterate(g.solver)
    cdef double tLeft = gsl_root_fsolver_root(g.solver)

    # Find right-side intersection
    gsl_root_fsolver_set(g.solver, g.func, g.tNear, g.tFar)
    while (gsl_root_fsolver_x_upper(g.solver) - gsl_root_fsolver_x_lower(g.solver)) > 1e-9:
        gsl_root_fsolver_iterate(g.solver)
    cdef double tRight = gsl_root_fsolver_root(g.solver)

    # Calculate the actual integrand, which is the integral of darkening across the arc
    cdef double thetaLeft = atan2(g.b*sin(tLeft)+g.ye, g.a*cos(tLeft)+g.xe)
    cdef double thetaRight = atan2(g.b*sin(tRight)+g.ye, g.a*cos(tRight)+g.xe)
   
    # Ensure that we're getting the length of the correct arc (the one inside the ellipse)
    if thetaLeft < thetaRight:
        thetaLeft += 2*pi

    return (thetaLeft-thetaRight) * rad * darkening(rad, g)


cpdef double transitDepth(double a, double b, double c, double semimajor, double theta, double phi, double[:] limb):
    theta = (theta+pi)%(2*pi) - pi
    if abs(theta) > pi/2.:
        return 0.
    # Get the sky-projected coordinates
    a, b, xe, ye  = orbitGeometry(a, b, c, semimajor, theta, phi)
    # Feed them into the transit depth integral
    return transitIntegral(a, b, xe, ye, limb) / pi


cpdef double transitIntegral(double a, double b, double xe, double ye, double[:] limb, int preferBrute=1):
    # Ensure the planet is in the first quadrant for consistency
    ye = abs(ye)
    xe = abs(xe)

    # Simple, blazing fast axis-aligned bounding box check to rule out many non-transits
    if xe - a > 1. or ye-b > 1:
        return 0.

    # If the ellipse contains the center, jump straight to the brute integral
    if (((xe/a)**2 + (ye/b)**2) < 1) or ():
        return bruteIntegrate(a, b, xe, ye, limb)

    # Init root-finding requirements
    gsl_set_error_handler_off()
    cdef IntegralParams g = IntegralParams(a, b, xe, ye, -1., -1., 0, [limb[i] for i in range(5)], NULL, NULL)
    cdef gsl_function dist
    dist.function = &originDistDiff
    dist.params = &g
    g.solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent)
    g.func = &dist

    # Find the pseudo-angle t that minimizes the distance to the origin
    if gsl_root_fsolver_set(g.solver, &dist, pi-.1, .1+3*pi/2):
        print("Root solver error while finding min distance!")
        gsl_root_fsolver_free(g.solver)
        return nan
    while (gsl_root_fsolver_x_upper(g.solver) - gsl_root_fsolver_x_lower(g.solver)) > 1e-9:
        gsl_root_fsolver_iterate(g.solver)
    g.tNear = gsl_root_fsolver_root(g.solver)
    cdef double rNear = sqrt(abs(originDist(g.tNear, &g)))
    if rNear > 1.:
        # Planet and star do not overlap -> zero transit depth.
        gsl_root_fsolver_free(g.solver)
        return 0.

    # Now find the t of max distance
    if gsl_root_fsolver_set(g.solver, &dist, 0., pi/2.):
        # Try again but with a narrow focus around pi/2.
        if gsl_root_fsolver_set(g.solver, &dist, pi/2-.1, pi/2+.1):
            if gsl_root_fsolver_set(g.solver, &dist, -.1, .1):
            # print("Root solver error while finding max distance!")
                print("Max error", xe, ye, originDistDiff(-.1, <void*>&g), originDistDiff(.1,<void*>&g), flush=True)
                gsl_root_fsolver_free(g.solver)
                return nan
    while (gsl_root_fsolver_x_upper(g.solver) - gsl_root_fsolver_x_lower(g.solver)) > 1e-9:
        gsl_root_fsolver_iterate(g.solver)
    g.tFar = gsl_root_fsolver_root(g.solver)
    cdef double rFar = sqrt(abs(originDist(g.tFar, &g)))

    # If the planet is fully in-transit, we can use either integrator
    # Which one is faster depends on how well-behaved the limb-darkening is
    if (rFar < 1.0) and preferBrute:
        gsl_root_fsolver_free(g.solver)
        return bruteIntegrate(a, b, xe, ye, limb)

    # Ensure that tNear and tFar are sorted
    if g.tNear > g.tFar:
        g.tNear -= 2*pi

    # Prepare the integrator
    cdef double result, err
    dist.function = &originDist
    cdef gsl_integration_workspace* workspace
    workspace = gsl_integration_workspace_alloc(100)
    cdef gsl_function integ
    integ.function = &integrand
    integ.params = &g

    gsl_integration_qag(&integ, rNear, min(1., rFar), 0., 1e-7, 100, 1, workspace, &result, &err)

    gsl_root_fsolver_free(g.solver)
    gsl_integration_workspace_free(workspace)
    return result


cdef double bruteIntegrand(double y, void* params) noexcept:
    cdef BruteIntegralParams* g = <BruteIntegralParams*>params
    cdef double mu = 1. - g.x*g.x - y*y
    mu = 1 - sqrt(max(mu, 0.))
    return g.limb[0] - g.limb[1]*mu - g.limb[2]*mu*mu


cdef double bruteIntegrateY(double x, void* params) noexcept:
    cdef double result, err
    cdef BruteIntegralParams* g = <BruteIntegralParams*>params

    # Calculate the bounds of integration in y
    cdef double yMin = (x-g.xe)/g.a
    yMin = g.b * sqrt(1.-yMin*yMin)
    cdef double yMax = g.ye + yMin
    yMin = g.ye - yMin
    g.x = x

    gsl_integration_qag(g.integrand, yMin, yMax, 0., 1e-7, 100, 1, g.work, &result, &err)
    return result


cpdef double bruteIntegrate(double a, double b, double xe, double ye, double[:] limb):
    cdef double result, err
    cdef BruteIntegralParams g = BruteIntegralParams(a, b, xe, ye, 0., [limb[i] for i in range(5)], NULL, NULL)
    
    # Prepare the inner (y) integral variables
    cdef gsl_integration_workspace* workspaceInner = gsl_integration_workspace_alloc(100)
    cdef gsl_function integInner
    integInner.function = &bruteIntegrand
    integInner.params = &g
    g.work = workspaceInner
    g.integrand = &integInner

    # Now prepare the outer (x) integral variables
    cdef gsl_integration_workspace* workspaceOuter = gsl_integration_workspace_alloc(100)
    cdef gsl_function integOuter
    integOuter.function = &bruteIntegrateY
    integOuter.params = &g

    gsl_integration_qag(&integOuter, xe-a, xe+a, 0., 1e-7, 100, 1, workspaceOuter, &result, &err)

    gsl_integration_workspace_free(workspaceInner)
    gsl_integration_workspace_free(workspaceOuter)
    return result


cpdef (double, double, double, double) orbitGeometry(double a, double b, double c, double semimajor, double theta, double phi):
    ''' a, b, and c are the planet radii relative to the stellar radius
        semimajor is the semimajor axis / stellar radius
        Theta is 2 * pi * time / orbital period
        phi is inclination - 90 degrees'''
    # Calculate the sine and cosine of the angles: we'll use them repeatedly
    cdef SpiceDouble ct = cos(theta)
    cdef SpiceDouble st = sin(theta)
    cdef SpiceDouble cp = cos(phi)
    cdef SpiceDouble sp = sin(phi)

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

    cdef double majorLen = sqrt(majorAx[0]*majorAx[0] + majorAx[1]*majorAx[1] + majorAx[2]*majorAx[2])
    cdef double minorLen = sqrt(minorAx[0]*minorAx[0] + minorAx[1]*minorAx[1] + minorAx[2]*minorAx[2])

    # print("=== Planet Frame ===")
    if abs(majorAx[1]*view[1]+majorAx[2]*view[2] + majorAx[0]*view[0]) > 1e-14 or \
          abs(minorAx[1]*view[1]+minorAx[2]*view[2] + minorAx[0]*view[0]) > 1e-14:
        print("ERROR: projected axes were not perpendicular to the view vector!")
    # print("View         ", view)
    # print("MajorAx      ", majorAx)
    # print("MinorAx      ", minorAx)

    # Rotate the major axis into the unaligned view frame
    majorAx = [
        -st*majorAx[0]      - ct*majorAx[1],
        -ct*sp*majorAx[0]   + st*sp*majorAx[1]  + cp*majorAx[2],
        0. # Zero by contsruction: -ct*cp*majorAx[0]    + st*cp*majorAx[1]  - sp*majorAx[2]
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
