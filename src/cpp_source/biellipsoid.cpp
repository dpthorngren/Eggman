#include "biellipsoid.hpp"
#include "math_utils.hpp"
#include <cmath>

Biellipsoid::Biellipsoid() {
    r_forward = 1.0;
    r_back = 1.0;
    r_up = 1.0;
    r_side = 1.0;
    rot = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    f_limb = Ellipse();
    b_limb = Ellipse();
}

Biellipsoid::Biellipsoid(double r_forward, double r_back, double r_up, double r_side) {
    this->r_forward = r_forward;
    this->r_back = r_back;
    this->r_up = r_up;
    this->r_side = r_side;
    // Valid defaults, but expect user to usually call set_rotation next.
    rot = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    f_limb = Ellipse();
    b_limb = Ellipse();
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
    theta = theta * M_PI / 180;
    phi = phi * M_PI / 180;
    gamma = gamma * M_PI / 180;
    double ct = cos(theta);
    double st = sin(theta);
    double cp = cos(phi);
    double sp = sin(phi);
    double cg = cos(gamma);
    double sg = sin(gamma);
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

void Biellipsoid::set_position(Vec3 new_position) { this->position = new_position; }

void Biellipsoid::update_derived() {
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
        e1 = (Vec3){limb_plane.y, -limb_plane.x, 0.};
        l = LENGTH(e1);
        RESCALE(e1, l);
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
        e1 = (Vec3){limb_plane.y, -limb_plane.x, 0.};
        l = LENGTH(e1);
        RESCALE(e1, l);
        CROSS(e1, limb_plane, e2);
    }

    // Transform the limb vectors into view space
    temp = (Vec3){e1.x * r_back, e1.y * r_up, e1.z * r_side};
    MATMUL(rot, temp, e1);
    temp = (Vec3){e2.x * r_back, e2.y * r_up, e2.z * r_side};
    MATMUL(rot, temp, e2);
    b_limb = Ellipse(e1, e2);
}

Bounds Biellipsoid::slice_ylimits(double x) {
    Vec3 lower, upper;
    Bounds result = {INFINITY, -INFINITY};
    x -= position.x;

    f_limb.get_ybounds(x, &lower, &upper);
    if (is_forward_local(lower)) {
        result.min = fmin(result.min, lower.y);
        result.max = fmax(result.max, lower.y);
    }
    if (is_forward_local(upper)) {
        result.min = fmin(result.min, upper.y);
        result.max = fmax(result.max, upper.y);
    }

    b_limb.get_ybounds(x, &lower, &upper);
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

bool Biellipsoid::line_intersects(double x, double y) {
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

Vec3 Biellipsoid::line_project(double x, double y) {
    x -= position.x;
    y -= position.y;
    Vec3 hit;
    // Line origin in forward sphere-space: (x, y, 0) in world-space
    Vec3 p0 = {
        (rot.xx * x + rot.yx * y) / r_forward, (rot.xy * x + rot.yy * y) / r_up,
        (rot.xz * x + rot.yz * y) / r_side
    };
    // Find the forward near-point in forward sphere-space
    // Direction of the line in forward sphere-space
    Vec3 u = {rot.zx / r_forward, rot.zy / r_up, rot.zz / r_side};
    double len = LENGTH(u);
    RESCALE(u, len);
    double det, offset;
    offset = u.x * p0.x + u.y * p0.y + u.z * p0.z;
    len = LENGTH(p0);
    det = offset * offset - len * len + 1;
    if (det > 0) {
        offset = sqrt(det) - offset;
        hit = (Vec3){
            (p0.x + offset * u.x) * r_forward, (p0.y + offset * u.y) * r_up,
            (p0.z + offset * u.z) * r_side
        };
        if (hit.x > 0) {
            return aligned_to_world(hit);
        }
    }

    // Find the backward near-point in backward sphere-space
    p0.x = p0.x * r_forward / r_back;
    u = {rot.zx / r_back, rot.zy / r_up, rot.zz / r_side};
    len = LENGTH(u);
    RESCALE(u, len);
    offset = u.x * p0.x + u.y * p0.y + u.z * p0.z;
    len = LENGTH(p0);
    det = offset * offset - len * len + 1;
    if (det < 0) {
        // No intersections
        return (Vec3){NAN, NAN, NAN};
    }
    offset = sqrt(det) - offset;
    hit = {
        (p0.x + offset * u.x) * r_back, (p0.y + offset * u.y) * r_up, (p0.z + offset * u.z) * r_side
    };
    if (hit.x > 0) {
        // No intersections
        return (Vec3){NAN, NAN, NAN};
    }
    return aligned_to_world(hit);
}

Vec3 Biellipsoid::nearest_to_line(double x, double y) {
    Vec3 result;
    if (is_forward((Vec3){x, y, position.z})) {
        result = f_limb.nearest_to_line(x - position.x, y - position.y);
    } else {
        result = b_limb.nearest_to_line(x - position.x, y - position.y);
    }
    result.x += position.x;
    return result;
}

inline bool Biellipsoid::is_forward(Vec3 loc) {
    return ((loc.x - position.x) * rot.xx + (loc.y - position.y) * rot.yx +
            (loc.z - position.z) * rot.zx) >= 0;
}

inline bool Biellipsoid::is_forward_local(Vec3 loc) {
    return (loc.x * rot.xx + loc.y * rot.yx + loc.z * rot.zx) >= 0;
}

bool Biellipsoid::is_visible(Vec3 loc) {
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

Bounds Biellipsoid::x_bounds() {
    if (rot.xx > 0) {
        return {position.x - b_limb.x_size, position.x + f_limb.x_size};
    } else {
        return {position.x - f_limb.x_size, position.x + b_limb.x_size};
    }
}

Bounds Biellipsoid::y_bounds() {
    if (rot.yy > 0) {
        return {position.y - b_limb.y_size, position.y + f_limb.y_size};
    } else {
        return {position.y - f_limb.y_size, position.y + b_limb.y_size};
    }
}

inline Vec3 Biellipsoid::world_to_aligned(Vec3 loc) {
    Vec3 result;
    loc.x -= position.x;
    loc.y -= position.y;
    loc.z -= position.z;
    MATMUL_T(rot, loc, result);
    return result;
}

inline Vec3 Biellipsoid::world_to_sphere(Vec3 loc) {
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

inline Vec3 Biellipsoid::aligned_to_world(Vec3 loc) {
    Vec3 result;
    MATMUL(rot, loc, result);
    result.x += position.x;
    result.y += position.y;
    result.z += position.z;
    return result;
}

inline Vec3 Biellipsoid::aligned_to_sphere(Vec3 loc) {
    loc.x /= (loc.x < 0 ? r_back : r_forward);
    loc.y /= r_up;
    loc.z /= r_side;
    return loc;
}

inline Vec3 Biellipsoid::sphere_to_world(Vec3 loc) {
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

inline Vec3 Biellipsoid::sphere_to_aligned(Vec3 loc) {
    loc.x *= (loc.x < 0 ? r_back : r_forward);
    loc.y *= r_up;
    loc.z *= r_side;
    return loc;
}
