#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

#include <cmath>
#include <gsl/gsl_errno.h>
#include <vector>

// Linear algebra structs and macros
#define SIGN(x) ((x < 0.) ? -1. : ((x > 0.) ? 1. : 0.))
#define CLAMP(x, xmin, xmax) (fmax(fmin(x, xmax), xmin))
#define LENGTH(arr) sqrt(arr.x *arr.x + arr.y * arr.y + arr.z * arr.z)
#define RESCALE(arr, coeff)                                                                        \
    arr.x /= coeff;                                                                                \
    arr.y /= coeff;                                                                                \
    arr.z /= coeff;
#define DOT3(a, b) (a.x * b.x + a.y * b.y + a.z * b.z)
#define MATMUL(mat, vec, out)                                                                      \
    out.x = mat.xx * vec.x + mat.xy * vec.y + mat.xz * vec.z;                                      \
    out.y = mat.yx * vec.x + mat.yy * vec.y + mat.yz * vec.z;                                      \
    out.z = mat.zx * vec.x + mat.zy * vec.y + mat.zz * vec.z;
#define MATMUL_T(mat, vec, out)                                                                    \
    out.x = mat.xx * vec.x + mat.yx * vec.y + mat.zx * vec.z;                                      \
    out.y = mat.xy * vec.x + mat.yy * vec.y + mat.zy * vec.z;                                      \
    out.z = mat.xz * vec.x + mat.yz * vec.y + mat.zz * vec.z;
#define CROSS(a, b, out)                                                                           \
    out.x = a.y * b.z - b.y * a.z;                                                                 \
    out.y = a.z * b.x - b.z * a.x;                                                                 \
    out.z = a.x * b.y - b.x * a.y;
#define WEIGHTED_SUM(w1, w2, v1, v2, out)                                                          \
    out.x = w1 * v1.x + w2 * v2.x;                                                                 \
    out.y = w1 * v1.y + w2 * v2.y;                                                                 \
    out.z = w1 * v1.z + w2 * v2.z;
#define BIELLIPSE_DIR(bell, pos)                                                                   \
    (pos.x * bell->rot.xx + pos.y * bell->rot.yx + pos.z * bell->rot.zx)
#define PRINT(vec) printf("%f, %f, %f\n", vec.x, vec.y, vec.z);
#define PRINTMAT(mat)                                                                              \
    printf(                                                                                        \
        "[(%6f, %6f, %6f)\n (%6f, %6f, %6f)\n (%6f, %6f, %6f)]\n", mat.xx, mat.xy, mat.xz, mat.yx, \
        mat.yy, mat.yz, mat.zx, mat.zy, mat.zz                                                     \
    );

typedef struct {
    double x;
    double y;
    double z;
} Vec3;

typedef struct {
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

typedef struct {
    double min;
    double max;
} Bounds;

inline bool isclose(double a, double b, double tol = 1e-9) {
    return fabs(a - b) < tol + tol * fabs(a);
}

inline void matmul3x3(Mat3 &a, Mat3 &b, Mat3 &out) {
    out.xx = a.xx * b.xx + a.xy * b.yx + a.xz * b.zx;
    out.xy = a.xx * b.xy + a.xy * b.yy + a.xz * b.zy;
    out.xz = a.xx * b.xz + a.xy * b.yz + a.xz * b.zz;

    out.yx = a.yx * b.xx + a.yy * b.yx + a.yz * b.zx;
    out.yy = a.yx * b.xy + a.yy * b.yy + a.yz * b.zy;
    out.yz = a.yx * b.xz + a.yy * b.yz + a.yz * b.zz;

    out.zx = a.zx * b.xx + a.zy * b.yx + a.zz * b.zx;
    out.zy = a.zx * b.xy + a.zy * b.yy + a.zz * b.zy;
    out.zz = a.zx * b.xz + a.zy * b.yz + a.zz * b.zz;
}

inline Vec3 normalized(Vec3 v) {
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return (Vec3){v.x / len, v.y / len, v.z / len};
}

inline bool integration_failed(
    int code, double result, double err, double atol, double rtol, bool print_error = true
) {
    if (code == 0) {
        // No errors, accept the result
        return false;
    }
    if ((code == GSL_EMAXITER) || (code == GSL_EROUND)) {
        // Integrator didn't meet its error targets...
        if (err > atol + rtol * fabs(result)) {
            // ...but it meets ours, so accept the result without error
            return false;
        }
    }
    if (print_error) {
        const char *msg = gsl_strerror(code);
        printf("INTEGRATION ERROR %i: %s. Result=%f, err=%f\n", code, msg, result, err);
    }
    return true;
}

inline void sphere_to_gridmap(Vec3 &loc, int n, int m) {
    loc.x = n * fmod(atan2(loc.x, loc.z) / (2 * M_PI) + 0.5, 1);
    loc.y = (m + 1) * (CLAMP(loc.y, -1., 1.) + 1.) / 2. - 1.;
    loc.z = 0.;
    return;
}

inline void gridmap_to_sphere(Vec3 &loc, int n, int m) {
    loc.y = 2 * ((loc.y + 1.) / (m + 1.)) - 1.;
    double dist = sqrt(1 - loc.y * loc.y);
    double theta = 2 * M_PI * (loc.x / n - .5);
    loc.z = dist * cos(theta);
    loc.x = dist * sin(theta);
    return;
}

inline double smootherstep(double x) {
    if (x < 0)
        return 0.;
    if (x > 1)
        return 1.;
    return x * x * x * (x * (6. * x - 15.) + 10.);
}

inline double interp_gridmap(Vec3 loc, int n, int m, std::vector<double> data) {
    if (data.size() != n * m + 2) {
        return NAN;
    }
    sphere_to_gridmap(loc, n, m);

    double low, high;
    double weight = smootherstep(fmod(loc.x, 1));
    int i0 = (int)loc.x;
    int i1 = (i0 + 1) % (n - 1);
    int j = (int)loc.y;
    if (loc.y < 0) {
        low = data[n * m];
        high = data[i0] * (1 - weight) + data[i1] * weight;
    } else if (loc.y > m - 1) {
        low = data[n * j + i0] * (1 - weight) + data[n * j + i1] * weight;
        high = data[n * m + 1];
    } else {
        low = data[n * j + i0] * (1 - weight) + data[n * j + i1] * weight;
        j += 1;
        high = data[n * j + i0] * (1 - weight) + data[n * j + i1] * weight;
    }
    weight = smootherstep(fmod(1. + loc.y, 1));
    return (1. - weight) * low + weight * high;
}

#endif
