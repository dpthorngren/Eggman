#include "biellipsoid.hpp"
#include "math_utils.hpp"
#include <cmath>

Biellipsoid::Biellipsoid() {
    position = {0., 0., 0.};
    r_forward = 1.0;
    r_back = 1.0;
    r_up = 1.0;
    r_side = 1.0;
    theta = 0.;
    phi = 0.;
    gamma = 0.;
    rot = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    update_derived();
}

Biellipsoid::Biellipsoid(double r_forward, double r_back, double r_up, double r_side) {
    position = {0., 0., 0.};
    this->r_forward = r_forward;
    this->r_back = r_back;
    this->r_up = r_up;
    this->r_side = r_side;
    // Valid defaults, but expect user to usually call set_rotation next.
    theta = 0.;
    phi = 0.;
    gamma = 0.;
    rot = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    update_derived();
}

void Biellipsoid::set_radii(double r_forward, double r_back, double r_up, double r_side) {
    this->r_forward = r_forward;
    this->r_back = r_back;
    this->r_up = r_up;
    this->r_side = r_side;
    update_derived();
}

void Biellipsoid::set_rotation(double theta, double phi, double gamma, double ci) {
    // Setup the rotation matrix (from ellipsoid-aligned to view frame)
    this->theta = theta * M_PI / 180;
    this->phi = phi * M_PI / 180;
    this->gamma = gamma * M_PI / 180;
    double ct = cos(this->theta);
    double st = sin(this->theta);
    double cp = cos(this->phi);
    double sp = sin(this->phi);
    double cg = cos(this->gamma);
    double sg = sin(this->gamma);
    rot = {
        cp * ct,  -cp * st * cg + sp * sg, cp * st * sg + sp * cg,  st, ct * cg, -ct * sg,
        -sp * ct, sp * st * cg + cp * sg,  -sp * st * sg + cp * cg,
    };
    if (fabs(ci) <= 1.) {
        double rad = LENGTH(position);
        double si;
        Mat3 orbit_rot;
        Mat3 result;
        if (rad >= 0) {
            // TODO: Wrong sign if not 0 < i < pi/2
            st = position.x / rad;
            ct = SIGN(position.z) * sqrt(1 - st * st);
            si = sqrt(1 - ci * ci);
            orbit_rot = (Mat3){ct, 0, st, -st * ci, si, ct * ci, -st * si, -ci, ct * si};
            matmul3x3(rot, orbit_rot, result);
            rot = result;
        }
    }
    update_derived();
}

void Biellipsoid::position_from_orbit(double t, const Orbit &orb, bool rotate_with_orbit) {
    set_position(orb.get_position(t));
    if (rotate_with_orbit) {
        set_rotation(
            theta * 180. / M_PI, phi * 180. / M_PI, gamma * 180. / M_PI, orb.get_cos_inc()
        );
    }
}

void Biellipsoid::set_position(Vec3 new_position) { this->position = new_position; }

void Biellipsoid::update_derived() {
    // See if we can bypass these calculations
    is_sphere =
        isclose(r_forward, r_back) && isclose(r_forward, r_up) && isclose(r_forward, r_side);
    if (is_sphere) {
        f_limb = Ellipse((Vec3){r_forward, 0., 0.}, (Vec3){0., r_forward, 0.});
        b_limb = Ellipse({r_forward, 0., 0.}, {0., r_forward, 0.});
        return;
    }
    // Calculate limb planes and construct the limb ellipses
    Vec3 temp, e1, e2;
    Vec3 limb_plane = (Vec3){-rot.zx / r_forward, -rot.zy / r_up, -rot.zz / r_side};
    double l = LENGTH(limb_plane);
    if (limb_plane.z < 0.0) {
        l *= -1.0;
    }
    RESCALE(limb_plane, l);
    e1 = {1.0, 0.0, 0.0};
    e2 = {0.0, 1.0, 0.0};
    if (limb_plane.x * limb_plane.x + limb_plane.y * limb_plane.y > 1e-12) {
        e1 = normalized({limb_plane.y, -limb_plane.x, 0.});
        CROSS(e1, limb_plane, e2);
    }
    // Transform back to the view plane
    temp = (Vec3){e1.x * r_forward, e1.y * r_up, e1.z * r_side};
    MATMUL(rot, temp, e1);
    temp = (Vec3){e2.x * r_forward, e2.y * r_up, e2.z * r_side};
    MATMUL(rot, temp, e2);
    f_limb = Ellipse(e1, e2);

    // Again, but for the back ellipse
    limb_plane = (Vec3){-rot.zx / r_back, -rot.zy / r_up, -rot.zz / r_side};
    l = LENGTH(limb_plane);
    if (limb_plane.z < 0.0) {
        l *= -1.0;
    }
    RESCALE(limb_plane, l);
    // Transform back to the view plane
    e1 = {1.0, 0.0, 0.0};
    e2 = {0.0, 1.0, 0.0};
    if (limb_plane.x * limb_plane.x + limb_plane.y * limb_plane.y > 1e-12) {
        e1 = normalized({limb_plane.y, -limb_plane.x, 0.});
        CROSS(e1, limb_plane, e2);
    }

    // Transform the limb vectors into view space
    temp = (Vec3){e1.x * r_back, e1.y * r_up, e1.z * r_side};
    MATMUL(rot, temp, e1);
    temp = (Vec3){e2.x * r_back, e2.y * r_up, e2.z * r_side};
    MATMUL(rot, temp, e2);
    b_limb = Ellipse(e1, e2);

    // Construct the joint ellipse
    e1 = {r_up * rot.xy, r_up * rot.yy, r_up * rot.zy};
    e2 = {r_side * rot.xz, r_side * rot.yz, r_side * rot.zz};
    joint = Ellipse(e1, e2);
    if (fabs(joint.det) < 1e-12) {
        // Joint has no area (edge-on), can't rotate axes.
        return;
    }

    // Find the limb-joint intersection point
    double d1 = rot.xx * f_limb.e1.x + rot.yx * f_limb.e1.y + rot.zx * f_limb.e1.z;
    double d2 = rot.xx * f_limb.e2.x + rot.yx * f_limb.e2.y + rot.zx * f_limb.e2.z;
    double ct = -1. / sqrt(1 + d1 * d1 / (d2 * d2));
    double st = 1. / sqrt(1 + d2 * d2 / (d1 * d1));
    WEIGHTED_SUM(ct, st, f_limb.e1, f_limb.e2, e1);
    if (e1.x * rot.yx - e1.y * rot.xx < 0) {
        RESCALE(e1, -1.);
    }

    // Rotate the joint ellipse axes so that e1 points to the limb joint point
    double u = (joint.e2.y * e1.x - joint.e2.x * e1.y) / joint.det;
    double v = (-joint.e1.y * e1.x + joint.e1.x * e1.y) / joint.det;
    WEIGHTED_SUM(-v, u, joint.e1, joint.e2, e2);
    joint = Ellipse(e1, e2);
}

Bounds Biellipsoid::slice_ylimits(double x) const {
    Vec3 lower, upper;
    Bounds result = {INFINITY, -INFINITY};
    x -= position.x;

    if (is_sphere) {
        x /= r_forward;
        result.max = r_forward * sqrt(1 - x * x);
        result.min = position.y - result.max;
        result.max += position.y;
        return result;
    }

    f_limb.get_ybounds(x, lower, upper);
    if (is_forward_local(lower)) {
        result.min = fmin(result.min, lower.y);
        result.max = fmax(result.max, lower.y);
    }
    if (is_forward_local(upper)) {
        result.min = fmin(result.min, upper.y);
        result.max = fmax(result.max, upper.y);
    }

    b_limb.get_ybounds(x, lower, upper);
    if (!is_forward_local(lower)) {
        result.min = fmin(result.min, lower.y);
        result.max = fmax(result.max, lower.y);
    }
    if (!is_forward_local(upper)) {
        result.min = fmin(result.min, upper.y);
        result.max = fmax(result.max, upper.y);
    }

    result.min += position.y;
    result.max += position.y;
    return result;
}

bool Biellipsoid::line_intersects(double x, double y) const {
    bool hit;
    Vec3 loc;
    x -= position.x;
    y -= position.y;
    hit = f_limb.line_intersects(x, y, &loc);
    if (hit && is_forward_local(loc)) {
        return true;
    }
    hit = b_limb.line_intersects(x, y, &loc);
    if (hit && (!is_forward_local(loc))) {
        return true;
    }
    return false;
}

bool Biellipsoid::raycast(double x, double y, double *mu_out, Vec3 *hit_out) const {
    double mu, rsq, len_u_sq, det, offset, lusq;
    Vec3 hit, p0, u;
    x -= position.x;
    y -= position.y;
    if (is_sphere) {
        // Much simpler calculation for spheres
        hit.x = x / r_forward;
        hit.y = y / r_forward;
        rsq = hit.x * hit.x + hit.y * hit.y;
        if (rsq < 0) {
            return false; // No intersections
        }
        hit.z = sqrt(1 - rsq);
        mu = hit.z;
        if (hit_out != nullptr) {
            MATMUL_T(rot, hit, p0);
            hit = p0;
        }
        goto OUTPUT_RET;
    }
    // TODO: Determine forward or backwards ellipse directly, skip for symmetric case

    // Line origin in forward sphere-space: (x, y, 0) in world-space
    p0 = {
        (rot.xx * x + rot.yx * y) / r_forward, (rot.xy * x + rot.yy * y) / r_up,
        (rot.xz * x + rot.yz * y) / r_side
    };
    // Direction of the line in forward sphere-space
    u = {rot.zx / r_forward, rot.zy / r_up, rot.zz / r_side};
    // Saving a sqrt by not normalizing u, but it makes some of the algebra a little funky
    lusq = u.x * u.x + u.y * u.y + u.z * u.z;
    offset = u.x * p0.x + u.y * p0.y + u.z * p0.z;
    len_u_sq = p0.x * p0.x + p0.y * p0.y + p0.z * p0.z;
    det = offset * offset / lusq - len_u_sq + 1;
    if (det >= 0) {
        offset = (sqrt(lusq * det) - offset) / lusq;
        hit = {(p0.x + offset * u.x), (p0.y + offset * u.y), (p0.z + offset * u.z)};
        if (hit.x >= 0) {
            // TODO: Calculate mu
            goto OUTPUT_RET;
        }
    }

    // Find the backward near-point in backward sphere-space
    p0.x = p0.x * r_forward / r_back;
    u = {rot.zx / r_back, rot.zy / r_up, rot.zz / r_side};
    lusq = u.x * u.x + u.y * u.y + u.z * u.z;
    offset = u.x * p0.x + u.y * p0.y + u.z * p0.z;
    len_u_sq = p0.x * p0.x + p0.y * p0.y + p0.z * p0.z;
    det = offset * offset / lusq - len_u_sq + 1;
    if (det < 0) {
        return false; // No intersections
    }
    offset = (sqrt(lusq * det) - offset) / lusq;
    hit = {(p0.x + offset * u.x), (p0.y + offset * u.y), (p0.z + offset * u.z)};
    if (hit.x > 0) {
        return false; // No intersections
    }

OUTPUT_RET: // Yes, yes, goto's are bad.  Deal with it :)
    if (mu_out != nullptr) {
        *mu_out = mu;
    }
    if (hit_out != nullptr) {
        *hit_out = hit;
    }
    return true;
}

Vec3 Biellipsoid::nearest_to_line(double x, double y) const {
    Vec3 result;
    if (is_forward((Vec3){x, y, position.z})) {
        result = f_limb.nearest_to_line(x - position.x, y - position.y);
    } else {
        result = b_limb.nearest_to_line(x - position.x, y - position.y);
    }
    result.x += position.x;
    return result;
}

inline bool Biellipsoid::is_forward(Vec3 loc) const {
    return ((loc.x - position.x) * rot.xx + (loc.y - position.y) * rot.yx +
            (loc.z - position.z) * rot.zx) >= 0;
}

inline bool Biellipsoid::is_forward_local(Vec3 loc) const {
    return (loc.x * rot.xx + loc.y * rot.yx + loc.z * rot.zx) >= 0;
}

bool Biellipsoid::is_forward_2d(double x, double y, bool local) const {
    if (!local) {
        x -= position.x;
        y -= position.y;
    }
    if (fabs(joint.det) < 1e-12) {
        // Ellipse has no area, all points are outside it
        return x * rot.xx + y * rot.yx > 0.;
    }
    // Is the forward side closer to the viewer?
    bool forward_near = rot.zx >= 0;
    // Is the point on the forward-side of the joint line?
    bool side = joint.e1.x * y - joint.e1.y * x >= 0;
    if (side == forward_near) {
        // The point is on the near side, return that.
        return forward_near;
    }
    // The point is on the far side, making things harder
    // If the point overlaps the joint ellipse, then it hits the near side
    double u = (joint.e2.y * x - joint.e2.x * y) / joint.det;
    double v = (-joint.e1.y * x + joint.e1.x * y) / joint.det;
    bool in_ellipse = (u * u + v * v <= 1);
    // If in ellipse, return the forward side, otherwise return the back side
    // Forward side indicated by forward_near, so this is equivalent to an XNOR.
    return !(in_ellipse ^ forward_near);
}

bool Biellipsoid::is_visible(Vec3 loc) const {
    // In this niche world-aligned-but-scaled frame, planet is a unit sphere at the origin
    // but the view vector is (0, 0, -1)
    Vec3 loc_sph = world_to_sphere(loc);
    MATMUL(rot, loc_sph, loc);
    // Thus the visibilty test is simple: is the point behind a unit sphere at the origin?
    double xysq = loc.x * loc.x + loc.y * loc.y;
    if (xysq >= 1) {
        return true;
    }
    if (loc.z < 0) {
        return false;
    }
    return (xysq + loc.z * loc.z) > 1;
}

Bounds Biellipsoid::x_bounds() const {
    if (rot.xx > 0) {
        return {position.x - b_limb.x_size, position.x + f_limb.x_size};
    } else {
        return {position.x - f_limb.x_size, position.x + b_limb.x_size};
    }
}

Bounds Biellipsoid::y_bounds() const {
    if (rot.yy > 0) {
        return {position.y - b_limb.y_size, position.y + f_limb.y_size};
    } else {
        return {position.y - f_limb.y_size, position.y + b_limb.y_size};
    }
}

double Biellipsoid::get_area() const { return 0.5 * (f_limb.get_area() + b_limb.get_area()); }

inline Vec3 Biellipsoid::world_to_aligned(Vec3 loc) const {
    Vec3 result;
    loc.x -= position.x;
    loc.y -= position.y;
    loc.z -= position.z;
    MATMUL_T(rot, loc, result);
    return result;
}

inline Vec3 Biellipsoid::world_to_sphere(Vec3 loc) const {
    Vec3 result;
    loc.x -= position.x;
    loc.y -= position.y;
    loc.z -= position.z;
    MATMUL_T(rot, loc, result);
    result.x /= (result.x < 0 ? r_back : r_forward);
    result.y /= r_up;
    result.z /= r_side;
    return result;
}

inline Vec3 Biellipsoid::aligned_to_world(Vec3 loc) const {
    Vec3 result;
    MATMUL(rot, loc, result);
    result.x += position.x;
    result.y += position.y;
    result.z += position.z;
    return result;
}

inline Vec3 Biellipsoid::aligned_to_sphere(Vec3 loc) const {
    loc.x /= (loc.x < 0 ? r_back : r_forward);
    loc.y /= r_up;
    loc.z /= r_side;
    return loc;
}

inline Vec3 Biellipsoid::sphere_to_world(Vec3 loc) const {
    Vec3 result;
    loc.x *= (loc.x < 0 ? r_back : r_forward);
    loc.y *= r_up;
    loc.z *= r_side;
    MATMUL(rot, loc, result);
    result.x += position.x;
    result.y += position.y;
    result.z += position.z;
    return result;
}

inline Vec3 Biellipsoid::sphere_to_aligned(Vec3 loc) const {
    loc.x *= (loc.x < 0 ? r_back : r_forward);
    loc.y *= r_up;
    loc.z *= r_side;
    return loc;
}
