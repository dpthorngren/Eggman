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

bool Ellipse::line_intersects(double x, double y, Vec3 *out) {
    // Inverse matrix to transform to circle space, check if point is inside
    double det = e1.x * e2.y - e2.x * e1.y;
    double u = (e2.y * x - e2.x * y) / det;
    double v = (-e1.y * x + e1.x * y) / det;
    if (out != nullptr) {
        WEIGHTED_SUM(u, v, e1, e2, (*out));
    }
    return u * u + v * v <= 1.;
}

Vec3 Ellipse::nearest_to_line(double x0, double y0) {
    Vec3 result = {x0, y0, 0.0};
    if (line_intersects(x0, y0, &result)) {
        return result;
    }

    // Initial guess
    double det = e1.x * e2.y - e2.x * e1.y;
    double u = (e2.y * x0 - e2.x * y0) / det;
    double v = (-e1.y * x0 + e1.x * y0) / det;
    double t = atan2(v, u);

    // Newton's method to find zero of l = dot(dr/dt, (r0-r))
    double ct, st, x, y, dxdt, dydt;
    double l, dldt;

    for (int i = 0; i < 32; i++) {
        ct = cos(t);
        st = sin(t);
        x = e1.x * ct + e2.x * st;
        y = e1.y * ct + e2.y * st;
        dxdt = -e1.x * st + e2.x * ct;
        dydt = -e1.y * st + e2.y * ct;

        l = (x0 - x) * dxdt + (y0 - y) * dydt;
        dldt = x * x + y * y - x0 * x - y0 * y - dxdt * dxdt - dydt * dydt;
        t -= CLAMP(l / dldt, -.2, .2);
    }
    result = (Vec3){x, y, e1.z * ct + e2.z * st};
    return result;
}
