#ifndef EGGMAN_H
#define EGGMAN_H
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

// 3-vector macros
#define LENGTH(arr) sqrt(arr.x*arr.x + arr.y*arr.y + arr.z*arr.z)
#define RESCALE(arr, coeff) arr.x /= coeff; arr.y /= coeff; arr.z /= coeff;
#define DOT3(a, b) a.x*b.x + a.y*b.y + a.z*b.z
#define MATMUL(mat, vec, out) \
out.x = mat.xx*vec.x + mat.xy*vec.y + mat.xz*vec.z;\
out.y = mat.yx*vec.x + mat.yy*vec.y + mat.yz*vec.z;\
out.z = mat.zx*vec.x + mat.zy*vec.y + mat.zz*vec.z;
#define CROSS(a, b, out) \
out.x = a.y*b.z - b.y*a.z; \
out.y = a.z*b.x - b.z*a.x; \
out.z = a.x*b.y - b.x*a.y;
#define WEIGHTED_SUM(w1, w2, v1, v2, out) \
out.x = w1*v1.x + w2*v2.x; \
out.y = w1*v1.y + w2*v2.y; \
out.z = w1*v1.z + w2*v2.z;
#define PRINT(vec) printf("%f, %f, %f\n", vec.x, vec.y, vec.z);


typedef struct{
    double x;
    double y;
    double z;
} Vec3;


typedef struct{
    double xx;
    double xy;
    double xz;
    double yx;
    double yy;
    double yz;
    double zx;
    double zy;
    double zz;
} Mat3;


typedef struct{
    double min;
    double max;
} Bounds;


typedef struct{
    Vec3 position;
    double radii[4];
    // Rotation matrix from ellipsoid space to view space
    // Forward vector (view space) is the first column rot[::3]
    // View vector (ellipsoid space) is minus the last column -rot[2::3]
    Mat3 rot;
    // Plane of the limb for the forward (backward) ellipsoid
    Vec3 f_limb;
    Vec3 b_limb;
    // x offset of the intersection of limb ellipse and break plane
    Vec3 break_offset;
    // Basis vectors 1 and 2 for the bounding ellipsoids f and b
    Vec3 f1;
    Vec3 f2;
    Vec3 b1;
    Vec3 b2;
    // View-axis-aligned bounding box relative to position
    Bounds xbounds;
    Bounds ybounds;
} Biellipsoid;


typedef struct{
    double period;
    double semimajor;
    double eccen;
    // Sin and cosine of the inclination
    double s_inc;
    double c_inc;
    // Sin and cosine of the true anomaly at periapse
    double sta_p;
    double cta_p;
    // Time of periapse passage
    double t_p;
} Orbit;


// Biellipsoid Functions
Biellipsoid create_biellipsoid(Vec3 position, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side);
Bounds get_ylimits(Biellipsoid* bell, double x);


// Orbit Functions
Orbit create_orbit(double period, double semimajor, double eccen, double inclination, double lon_periapse);
double solve_kepler(double mean_anomaly, double eccen);
Vec3 get_position(Orbit* orb, double t);
#endif
