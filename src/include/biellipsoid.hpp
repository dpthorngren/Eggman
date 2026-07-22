#ifndef BIELLIPSOID_HPP
#define BIELLIPSOID_HPP

#include "ellipse.hpp"
#include "math_utils.hpp"
#include "orbit.hpp"

/* Describes a biellipsoid -- two ellipsoids smoothly and continuously joined along the yz plane,
 * then rotated and shifted. Member functions provide bounds and limits related to integrating
 * across them.  All members are public for simplicity, but mutation should be done via the
 * setters or outputs may be incorrect.
 */
class Biellipsoid {
  public:
    Vec3 position;
    double r_forward;
    double r_back;
    double r_up;
    double r_side;
    bool is_sphere;
    double theta;
    double phi;
    double gamma;
    // Rotation matrix from ellipsoid space to view space.
    // Forward vector (view space) is the first column rot[::3].
    // View vector (ellipsoid space) is minus the last column -rot[2::3].
    Mat3 rot;
    Ellipse f_limb; // The ellipse defining the forward portion of the limb
    Ellipse b_limb; // The ellipse defining the backwards portion of the limb
    Ellipse joint;  // The ellipse connecting the two half-ellipsoids

    Biellipsoid();
    Biellipsoid(double r_forward, double r_backward, double r_up, double r_side);
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
    Bounds slice_ylimits(double x) const;
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

#endif
