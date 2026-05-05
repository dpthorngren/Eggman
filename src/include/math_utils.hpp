#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

// Linear algebra structs and macros
#define SIGN(x) ((x < 0.) ? -1. : ((x > 0.) ? 1. : 0.))
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

#endif
