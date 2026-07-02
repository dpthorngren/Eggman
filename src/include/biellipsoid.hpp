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

    Biellipsoid();
    Biellipsoid(double r_forward, double r_backward, double r_up, double r_side);
    void set_radii(double r_forward, double r_back, double r_up, double r_side);
    void set_position(Vec3 new_position);
    // Rotate around z, y, and x; then if ci (cos(inclination)) is given, align with orbit
    void set_rotation(double theta, double phi, double gamma, double ci = -2.);
    void position_from_orbit(double t, const Orbit &orb, bool rotate_with_orbit = false);
    void update_derived(); // Updates f_limb and b_limb when radii or rot change.

    // Derived info
    Bounds x_bounds();
    Bounds y_bounds();
    // Checks if loc is on the forward side of the biellipsoid or the back
    bool is_forward(Vec3 loc);
    bool is_forward_local(Vec3 loc);
    // Checks if loc is not behind / inside the biellipsoid
    bool is_visible(Vec3 loc);
    // Determines the y range occupied by the biellipsoid for this x value (min=max outside range)
    Bounds slice_ylimits(double x);
    // Returns whether the line through x, y intersects the biellipsoid
    bool line_intersects(double x, double y);
    // Finds the intersection of a line through (x, y) with the biellipsoid, choosing point with the
    // larger z and returning the location in world space or hit longitude, sin(latitude),
    // and cos(incidence angle) if mulatlon is set.
    Vec3 line_project(double x, double y, bool mulatlon = false);
    // Finds the nearest point on/in the biellipse to the line through x, y
    Vec3 nearest_to_line(double x, double y);
    // Gets the area of the biellipsoid in the x-y (view) plane
    double get_area();

    // Coordinate transforms
    Vec3 world_to_aligned(Vec3 loc);
    Vec3 world_to_sphere(Vec3 loc);
    Vec3 aligned_to_world(Vec3 loc);
    Vec3 aligned_to_sphere(Vec3 loc);
    Vec3 sphere_to_world(Vec3 loc);
    Vec3 sphere_to_aligned(Vec3 loc);
};

#endif
