#include "ellipse.hpp"
#include "math_utils.hpp"
#include <cmath>

Ellipse::Ellipse() {
    e1 = {1., 0., 0.};
    e2 = {0., 1., 0.};
    x_size = 1.0;
    y_size = 1.0;
}

Ellipse::Ellipse(Vec3 e1, Vec3 e2) {
    this->e1 = e1;
    this->e2 = e2;
    x_size = sqrt(e1.x * e1.x + e2.x * e2.x);
    y_size = sqrt(e1.y * e1.y + e2.y * e2.y);
}

void Ellipse::get_ybounds(double x, Vec3 *out_min, Vec3 *out_max) {
    x = CLAMP(x, -x_size, x_size);
    double disc = (x * x) * (e2.x * e2.x) - (x_size * x_size) * (x * x - e1.x * e1.x);
    // If the discriminant is negative return the boundary instead
    disc = disc > 0. ? sqrt(disc) : 0.;

    // Check for the axis-aligned case where the normal formula would cause a divide by zero
    bool aligned = fabs(e1.x) < x_size * 1e-6;
    // First intersection
    double st, ct;
    Vec3 p1, p2;
    st = aligned ? x / e2.x : (x * e2.x + disc) / (x_size * x_size);
    ct = aligned ? sqrt(1 - st * st) : (x - e2.x * st) / e1.x;
    WEIGHTED_SUM(ct, st, e1, e2, p1);
    // Second intersection
    st = aligned ? x / e2.x : (x * e2.x - disc) / (x_size * x_size);
    ct = aligned ? -sqrt(1 - st * st) : (x - e2.x * st) / e1.x;
    WEIGHTED_SUM(ct, st, e1, e2, p2);

    // Write the results, in sorted order, to the output locations
    if (p1.y < p2.y) {
        *out_min = p1;
        *out_max = p2;
    } else {
        *out_min = p2;
        *out_max = p1;
    }
}

bool Ellipse::line_intersects(double x, double y) {
    // Inverse matrix to transform to circle space, check if point is inside
    double det = e1.x * e2.y - e2.x * e1.y;
    double u = (e2.y * x - e1.y * y) / det;
    double v = (-e2.x * x + e1.x * y) / det;
    return u * u + v * v < 1.;
}

Vec3 Ellipse::plane_intersection(Vec3 normal) {
    // Get the limb breakpoints
    double d1 = DOT3(normal, e1);
    double d2 = DOT3(normal, e2);
    double cbreak = -1. / sqrt(1 + d1 * d1 / (d2 * d2));
    double sbreak = 1. / sqrt(1 + d2 * d2 / (d1 * d1));
    Vec3 break_offset;
    WEIGHTED_SUM(cbreak, sbreak, e1, e2, break_offset);
    if (break_offset.x < 0) {
        RESCALE(break_offset, -1.);
    }
    return break_offset;
}


Vec3 Ellipse::nearest_to_line(double x0, double y0) {
    if (line_intersects(x0, y0)) {
        return {x0, y0, 0.0};
    }

    // Initial guess
    double det = e1.x * e2.y - e2.x * e1.y;
    double u = (e2.y * x0 - e1.y * y0) / det;
    double v = (-e2.x * x0 + e1.x * y0) / det;
    double t = atan2(v, u);

    // Newton's method in t to find the point where the normal vector points to (x0, y0)
    double ct, st, x, y, dxdt, dydt;
    double dldt, d2ldt2;

    for (int i = 0; i < 16; i++) {
        ct = cos(t);
        st = sin(t);
        x = e1.x * ct + e2.x * st;
        y = e1.y * ct + e2.y * st;
        dxdt = -e1.x * st + e2.x * ct;
        dydt = -e1.y * st + e2.y * ct;

        dldt = -x0 * y + y0 * x;
        d2ldt2 = -x0 * dydt + y0 * dxdt;

        t -= CLAMP(dldt / d2ldt2, -.2, .2);
    }
    return (Vec3){x, y, 0.};
}
