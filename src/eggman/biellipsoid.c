#include "eggman.h"

Biellipsoid create_biellipsoid(
    Vec3 position, double theta, double phi, double gamma, double r_forward, double r_back,
    double r_up, double r_side
) {
    double l;
    double ct, st, cp, sp, cg, sg;

    // Setup the rotation matrix (from ellipsoid-aligned to view frame)
    theta = theta * M_PI / 180;
    phi = phi * M_PI / 180;
    gamma = gamma * M_PI / 180;
    ct = cos(theta);
    st = sin(theta);
    cp = cos(phi);
    sp = sin(phi);
    cg = cos(gamma);
    sg = sin(gamma);
    Mat3 rot = {
        cp * ct,  -cp * st * cg + sp * sg, cp * st * sg + sp * cg,  st, ct * cg, -ct * sg,
        -sp * ct, sp * st * cg + cp * sg,  -sp * st * sg + cp * cg,
    };

    // Calculate limb planes and vectors
    Vec3 f_limb = (Vec3){-rot.zx / r_forward, -rot.zy / r_up, -rot.zz / r_side};
    l = LENGTH(f_limb);
    if (f_limb.z < 0.0) {
        l *= -1.0;
    }
    RESCALE(f_limb, l);

    Vec3 f1 = {1.0, 0.0, 0.0};
    Vec3 f2 = {0.0, 1.0, 0.0};
    if (f_limb.x * f_limb.x + f_limb.y * f_limb.y > 1e-12) {
        f1 = (Vec3){f_limb.y, -f_limb.x, 0.};
        l = LENGTH(f1);
        RESCALE(f1, l);
        CROSS(f1, f_limb, f2);
    }

    Vec3 b_limb = (Vec3){-rot.zx / r_back, -rot.zy / r_up, -rot.zz / r_side};
    l = LENGTH(b_limb);
    if (b_limb.z < 0.0) {
        l *= -1.0;
    }
    RESCALE(b_limb, l);

    Vec3 b1 = {1.0, 0.0, 0.0};
    Vec3 b2 = {0.0, 1.0, 0.0};
    if (b_limb.x * b_limb.x + b_limb.y * b_limb.y > 1e-12) {
        b1 = (Vec3){b_limb.y, -b_limb.x, 0.};
        l = LENGTH(b1);
        RESCALE(b1, l);
        CROSS(b1, b_limb, b2);
    }

    // Transform the limb planes into view space
    Vec3 temp = {f_limb.x / r_forward, f_limb.y / r_up, f_limb.z / r_side};
    MATMUL(rot, temp, f_limb);

    temp = (Vec3){b_limb.x / r_back, b_limb.y / r_up, b_limb.z / r_side};
    MATMUL(rot, temp, b_limb);

    // Transform the limb vectors into view space
    temp = (Vec3){f1.x * r_forward, f1.y * r_up, f1.z * r_side};
    MATMUL(rot, temp, f1);

    temp = (Vec3){f2.x * r_forward, f2.y * r_up, f2.z * r_side};
    MATMUL(rot, temp, f2);

    temp = (Vec3){b1.x * r_back, b1.y * r_up, b1.z * r_side};
    MATMUL(rot, temp, b1);

    temp = (Vec3){b2.x * r_back, b2.y * r_up, b2.z * r_side};
    MATMUL(rot, temp, b2);

    // Get the limb breakpoints
    double d1 = rot.xx * f1.x + rot.yx * f1.y + rot.zx * f1.z;
    double d2 = rot.xx * f2.x + rot.yx * f2.y + rot.zx * f2.z;
    double cbreak = -1. / sqrt(1 + d1 * d1 / (d2 * d2));
    double sbreak = 1. / sqrt(1 + d2 * d2 / (d1 * d1));
    Vec3 break_offset;
    WEIGHTED_SUM(cbreak, sbreak, f1, f2, break_offset);
    if (break_offset.x < 0) {
        RESCALE(break_offset, -1.);
    }
    Biellipsoid result = {0.,     0.,     0.,           r_forward, r_back, r_up, r_side, rot,
                          f_limb, b_limb, break_offset, f1,        f2,     b1,   b2};
    set_position(&result, position);
    return result;
}


void set_position(Biellipsoid *bell, Vec3 position) {
    // Sets the position and updates the bounding box for the biellipsoid
    bell->position = position;
    double f_width = sqrt(bell->f1.x * bell->f1.x + bell->f2.x * bell->f2.x);
    double b_width = sqrt(bell->b1.x * bell->b1.x + bell->b2.x * bell->b2.x);
    if (bell->rot.xx > 0) {
        bell->xbounds = (Bounds){position.x - b_width, position.x + f_width};
    } else {
        bell->xbounds = (Bounds){position.x - f_width, position.x + b_width};
    }
    f_width = sqrt(bell->f1.y * bell->f1.y + bell->f2.y * bell->f2.y);
    b_width = sqrt(bell->b1.y * bell->b1.y + bell->b2.y * bell->b2.y);
    if (bell->rot.yy > 0) {
        bell->ybounds = (Bounds){position.y - b_width, position.y + f_width};
    } else {
        bell->ybounds = (Bounds){position.y - f_width, position.y + b_width};
    }
    return;
}


Bounds get_ylimits(Biellipsoid *bell, double x) {
    double st, ct;
    Vec3 e1, e2;
    Vec3 forward = {bell->rot.xx, bell->rot.yx, bell->rot.zx};
    Bounds ybounds = {0., 0.};

    if ((x <= bell->xbounds.min) || (x >= bell->xbounds.max)) {
        return ybounds;
    }
    x -= bell->position.x;

    double xwidth2 = bell->f1.x * bell->f1.x + bell->f2.x * bell->f2.x;
    int aligned = bell->f1.x * bell->f1.x < xwidth2 * 1e-12;
    double disc = (x * x) * (bell->f2.x * bell->f2.x) - xwidth2 * (x * x - bell->f1.x * bell->f1.x);
    disc = sqrt(disc);
    // First and second forward ellipse intersections
    st = aligned ? x / bell->f2.x : (x * bell->f2.x + disc) / xwidth2;
    ct = aligned ? sqrt(1 - st * st) : (x - bell->f2.x * st) / bell->f1.x;
    WEIGHTED_SUM(ct, st, bell->f1, bell->f2, e1);
    st = aligned ? x / bell->f2.x : (x * bell->f2.x - disc) / xwidth2;
    ct = aligned ? -sqrt(1 - st * st) : (x - bell->f2.x * st) / bell->f1.x;
    WEIGHTED_SUM(ct, st, bell->f1, bell->f2, e2);

    // Select bounds ased on their location
    if (DOT3(e1, forward) >= 0) {
        if (e1.y < ybounds.min)
            ybounds.min = e1.y;
        else if (e1.y > ybounds.max)
            ybounds.max = e1.y;
    }
    if (DOT3(e2, forward) >= 0) {
        if (e2.y < ybounds.min)
            ybounds.min = e2.y;
        else if (e2.y > ybounds.max)
            ybounds.max = e2.y;
    }

    // First and second backward ellipse intersections
    xwidth2 = bell->b1.x * bell->b1.x + bell->b2.x * bell->b2.x;
    aligned = bell->b1.x * bell->b1.x < xwidth2 * 1e-12;
    disc = (x * x) * (bell->b2.x * bell->b2.x) - xwidth2 * (x * x - bell->b1.x * bell->b1.x);
    disc = sqrt(disc);
    st = aligned ? x / bell->b2.x : (x * bell->b2.x + disc) / xwidth2;
    ct = aligned ? sqrt(1 - st * st) : (x - bell->b2.x * st) / bell->b1.x;
    WEIGHTED_SUM(ct, st, bell->b1, bell->b2, e1);
    st = aligned ? x / bell->b2.x : (x * bell->b2.x - disc) / xwidth2;
    ct = aligned ? -sqrt(1 - st * st) : (x - bell->b2.x * st) / bell->b1.x;
    WEIGHTED_SUM(ct, st, bell->b1, bell->b2, e2);

    // Select bounds based on their location
    if (DOT3(e1, forward) < 0) {
        if (e1.y < ybounds.min)
            ybounds.min = e1.y;
        else if (e1.y > ybounds.max)
            ybounds.max = e1.y;
    }
    if (DOT3(e2, forward) < 0) {
        if (e2.y < ybounds.min)
            ybounds.min = e2.y;
        else if (e2.y > ybounds.max)
            ybounds.max = e2.y;
    }
    ybounds.min += bell->position.y;
    ybounds.max += bell->position.y;
    return ybounds;
}


Vec3 nearest_point_3d(Biellipsoid *bell) {
    // Returns the x, y, and squared dist of the nearest point to the origin on the
    // biellipsoid, unless it contains the origin, in which case returns [0, 0, -1]
    Vec3 e1, e2;

    // Which ellipsoid is relevant (forward or back?)
    Vec3 zero = {0., 0., 0.};
    if (BIELLIPSE_DIR(bell, zero) >= 0.) {
        e1 = bell->f1;
        e2 = bell->f2;
    } else {
        e1 = bell->b1;
        e2 = bell->b2;
    }

    // Inverse matrix to transform to circle space, check if point is inside
    double det = e1.x*e2.y - e2.x*e1.y;
    double u = (-e2.y*bell->position.x + e2.x*bell->position.y)/det;
    double v = (e1.y*bell->position.x - e1.x*bell->position.y)/det;
    if (u * u + v * v < 1.) {
        // Origin is inside ellipse
        return (Vec3){0., 0., -1};
    }

    // Is e2 on the near or far side of the ellipse from the origin?
    double sign = -DOT3(e2, bell->position);
    sign = SIGN(sign);

    // Search in cos(theta), sin(theta) space,
    double ct, st, x, y, d;
    // Bound with 1,0 to -1,0 pick branch from whichever of 0,1 0,-1 is larger.
    double ct0 = 1.0;
    double st0 = 0.;
    x = bell->position.x + e1.x * ct0;
    y = bell->position.y + e1.y * ct0;
    double d0 = x * x + y * y;

    double ct1 = -1.0;
    double st1 = 0.;
    x = bell->position.x + e1.x * ct1;
    y = bell->position.y + e1.y * ct1;
    double d1 = x * x + y * y;
    d = d1;

    // Midpoint method
    int i;
    for (i = 0; i < 64; i++) {
        ct = (ct0 + ct1) / 2.;
        st = sign * sqrt(1 - ct * ct);
        x = bell->position.x + e1.x * ct + e2.x * st;
        y = bell->position.y + e1.y * ct + e2.y * st;
        d = x * x + y * y;
        if (d1 < d0) {
            ct0 = ct;
            st0 = st;
            d0 = d;
        } else {
            ct1 = ct;
            st1 = st;
            d1 = d;
        }
    }
    return (Vec3){x, y, d};
}
