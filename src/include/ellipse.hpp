#ifndef ELLIPSE_HPP
#define ELLIPSE_HPP

#include "math_utils.hpp"

/* Describes an ellipse in 3d space in terms of the constructing vectors e1 and e2.  The ellipse
 * comprises all points along a*e1 + b*e2 where a^2 + b^2 = 1.  The axis-aligned sizes are
 * pre-computed and member functions provide y-bounds are a given x for integration purposes.  All
 * members are public for simplicity, but should not be written to else outputs may be incorrect --
 * instantiate a new ellipse instead.
 */
class Ellipse {
  public:
    Vec3 e1;       // The first constructing vector for the ellipse
    Vec3 e2;       // The second constructing vector for the ellipse
    double x_size; // Bounds the ellipse in x to [+x_size, -x_size]
    double y_size; // Bounds the ellipse to y to [+y_size, -y_size]

    Ellipse();
    Ellipse(Vec3 e1, Vec3 e2);
    void get_ybounds(double x, Vec3 *out_min, Vec3 *out_max);
    // Returns whether the line through x, y in view space intersects the ellipse
    bool line_intersects(double x, double y, Vec3* out=nullptr);
    // Finds the nearest point on/in the ellipse to the line through x0, y0
    Vec3 nearest_to_line(double x0, double y0);
};

#endif
