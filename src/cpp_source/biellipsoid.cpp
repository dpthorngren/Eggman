#include "biellipsoid.hpp"
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
            st = position.x / rad;
            ct = SIGN(position.z) * sqrt(1 - st * st);
            si = sqrt(1 - ci * ci);
            orbit_rot = (Mat3){ct, 0, st, -st * ci, si, ct * ci, -st * si, -ci, ct * si};
            matmul3x3(orbit_rot, rot, result);
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
    Bounds result = {0., 0.};
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
    if (is_forward((Vec3){x, y, position.z})) {
        return f_limb.line_intersects(x - position.x, y - position.y);
    } else {
        return b_limb.line_intersects(x - position.x, y - position.y);
    }
}

Vec3 Biellipsoid::nearest_to_line(double x, double y) {
    if (is_forward((Vec3){x, y, position.z})) {
        return f_limb.nearest_to_line(x - position.x, y - position.y);
    } else {
        return b_limb.nearest_to_line(x - position.x, y - position.y);
    }
}

inline bool Biellipsoid::is_forward(Vec3 loc) {
    return ((loc.x - position.x) * rot.xx + (loc.y - position.y) * rot.yx +
            (loc.z - position.z) * rot.zx) >= 0;
}

inline bool Biellipsoid::is_forward_local(Vec3 loc) {
    return (loc.x * rot.xx + loc.y * rot.yx + loc.z * rot.zx) >= 0;
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
