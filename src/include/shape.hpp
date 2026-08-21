#ifndef STAR_SYSTEM_HPP
#define STAR_SYSTEM_HPP

#include "ellipse.hpp"
#include "math_utils.hpp"
#include "orbit.hpp"

typedef enum {
    Biellipsoid, // General case: two ellipsoids smoothly and continuously joined along the yz plane
    Ellipsoid,   // Forward and backward radii are equal
    Sphere,      // All radii are equal
    Ring,        // Upwards radius is zero (changes interpretation of radii and limbs)
} ShapeType;

/* Describes a shape for use in transit and emission calculations, taking one of the ShapeType
 * types, rotated by rot and shifted by position. Member functions provide bounds and limits related
 * to integrating across it.  All members are public for simplicity, but mutation should be done
 * via the setters or outputs may be incorrect.
 *
 * If r_up == 0, then the Shape is interpreted as a ring (technically an annulus) and the member
 * variables are interpreted as follows:
 * - r_forward is outer radius
 * - r_back is inner radius.
 * - f_limb is the outer disk edge.
 * - b_limb is the inner disk edge.
 * - joint and r_side are unused.
 */
class Shape {
  public:
    Vec3 position;
    double r_forward;
    double r_back;
    double r_up;
    double r_side;
    double theta;
    double phi;
    double gamma;
    ShapeType shape_type;
    // Rotation matrix from ellipsoid space to view space.
    // Forward vector (view space) is the first column rot[::3].
    // View vector (ellipsoid space) is minus the last column -rot[2::3].
    Mat3 rot;
    Ellipse f_limb; // The ellipse defining the forward portion of the limb
    Ellipse b_limb; // The ellipse defining the backwards portion of the limb
    Ellipse joint;  // The ellipse connecting the two half-ellipsoids

    Shape();
    Shape(double r_forward, double r_backward, double r_up, double r_side);
    void set_radii(double r_forward, double r_back, double r_up, double r_side);
    void set_position(Vec3 new_position);
    // Rotate around z, y, and x; then if ci (cos(inclination)) is given, align with orbit
    void set_rotation(double theta, double phi, double gamma, double ci = -2.);
    void position_from_orbit(
        double t, const Orbit &orb, bool rotate_with_orbit = false, Vec3 origin = {0., 0., 0.}
    );
    void update_derived(); // Updates f_limb and b_limb when radii or rot change.

    // Derived info
    Bounds x_bounds() const;
    Bounds y_bounds() const;
    // Checks if loc is on the forward side of the biellipsoid or the back
    bool is_forward(Vec3 loc) const;
    bool is_forward_local(Vec3 loc) const;
    bool is_forward_2d(double x, double y, bool local = false) const;
    // Checks if loc is not behind / inside the biellipsoid
    bool is_visible(Vec3 loc) const;
    // Determines the y range occupied by the biellipsoid for this x value (min=max outside range)
    Bounds slice_ylimits(double x, Bounds *out2 = nullptr, int zcut = 0) const;
    // Returns whether the line through x, y intersects the biellipsoid
    bool line_intersects(double x, double y) const;
    // Finds the near intersection of a line through (x, y) with the biellipsoid
    // writing to the angle of incidence mu and hit position (Sphere space) if not null
    bool raycast(double x, double y, double *mu_out = nullptr, Vec3 *hit_out = nullptr) const;
    // Finds the nearest point on/in the biellipse to the line through x, y
    Vec3 nearest_to_line(double x, double y) const;
    // Gets the area of the biellipsoid in the x-y (view) plane
    double get_area() const;

    // Coordinate transforms
    Vec3 world_to_aligned(Vec3 loc) const;
    Vec3 world_to_sphere(Vec3 loc) const;
    Vec3 aligned_to_world(Vec3 loc) const;
    Vec3 aligned_to_sphere(Vec3 loc) const;
    Vec3 sphere_to_world(Vec3 loc) const;
    Vec3 sphere_to_aligned(Vec3 loc) const;
};

// Inlined functions

inline Vec3 Shape::world_to_aligned(Vec3 loc) const {
    Vec3 result;
    loc.x -= position.x;
    loc.y -= position.y;
    loc.z -= position.z;
    MATMUL_T(rot, loc, result);
    return result;
}

inline Vec3 Shape::world_to_sphere(Vec3 loc) const {
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

inline Vec3 Shape::aligned_to_world(Vec3 loc) const {
    Vec3 result;
    MATMUL(rot, loc, result);
    result.x += position.x;
    result.y += position.y;
    result.z += position.z;
    return result;
}

inline Vec3 Shape::aligned_to_sphere(Vec3 loc) const {
    loc.x /= (loc.x < 0 ? r_back : r_forward);
    loc.y /= r_up;
    loc.z /= r_side;
    return loc;
}

inline Vec3 Shape::sphere_to_world(Vec3 loc) const {
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

inline Vec3 Shape::sphere_to_aligned(Vec3 loc) const {
    loc.x *= (loc.x < 0 ? r_back : r_forward);
    loc.y *= r_up;
    loc.z *= r_side;
    return loc;
}

inline Bounds z_clamp(Vec3 lower, Vec3 upper, int zcut) {
    if ((zcut == 0) || (upper.z == lower.z)) {
        return {lower.y, upper.y};
    }
    double slope = (upper.y - lower.y) / (upper.z - lower.z);
    double ycut = upper.y - upper.z * slope;
    if ((zcut > 0) ^ (slope > 0)) {
        return {fmin(lower.y, ycut), fmin(upper.y, ycut)};
    }
    return {fmax(lower.y, ycut), fmax(upper.y, ycut)};
}

#endif
