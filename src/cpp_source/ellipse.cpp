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
    double u = (e2.y * x + e2.x * y) / det;
    double v = (e1.y * x - e1.x * y) / det;
    return u * u + v * v < 1.;
}

Vec3 Ellipse::nearest_to_line(double xt, double yt) {
    if (line_intersects(xt, yt)) {
        return {xt, yt, 0.0};
    }
    // Is e2 on the near or far side of the ellipse from the point?
    Vec3 pos = {xt, yt, 0.};
    double sign = DOT3(e2, pos);
    sign = SIGN(sign);

    // Search in cos(theta), sin(theta) space,
    double ct, st, x, y, d, ct0, ct1, d0, d1;
    // Bound with 1,0 to -1,0 pick branch from whichever of 0,1 0,-1 is larger.
    ct0 = 1.0;
    x = e1.x * ct0 - xt;
    y = e1.y * ct0 - yt;
    d0 = x * x + y * y;

    ct1 = -1.0;
    x = e1.x * ct1 - xt;
    y = e1.y * ct1 - yt;
    d1 = x * x + y * y;

    // Midpoint method
    for (int i = 0; i < 64; i++) {
        ct = (ct0 + ct1) / 2.;
        st = sign * sqrt(1 - ct * ct);
        x = e1.x * ct + e2.x * st - xt;
        y = e1.y * ct + e2.y * st - yt;
        d = x * x + y * y;
        if (d1 < d0) {
            ct0 = ct;
            d0 = d;
        } else {
            ct1 = ct;
            d1 = d;
        }
    }
    return (Vec3){x, y, d};
}
