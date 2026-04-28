#include <math.h>
#include <stdio.h>

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

Biellipsoid create_biellipsoid(double x, double y, double z, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side){
    double l;
    double ct, st, cp, sp, cg, sg;

    // Setup the rotation matrox (from ellipsoid-aligned to view frame)
    theta = theta * M_PI / 180;
    phi = phi*M_PI / 180;
    gamma = gamma * M_PI / 180;
    ct = cos(theta);
    st = sin(theta);
    cp = cos(phi);
    sp = sin(phi);
    cg = cos(gamma);
    sg = sin(gamma);
    Mat3 rot = {
        cp * ct,  -cp * st * cg + sp*sg, cp*st*sg + sp*cg,
        st,       ct * cg,               -ct * sg,
        -sp * ct, sp*st*cg + cp*sg,      -sp * st * sg + cp*cg,
    };

    // Calculate limb planes and vectors
    Vec3 f_limb = (Vec3){-rot.zx/r_forward, -rot.zy/r_up, -rot.zz/r_side};
    l = LENGTH(f_limb);
    if (f_limb.z < 0.0){
        l *= -1.0;
    }
    RESCALE(f_limb, l);

    Vec3 f1 = {1.0, 0.0, 0.0};
    Vec3 f2 = {0.0, 1.0, 0.0};
    if (f_limb.x*f_limb.x + f_limb.y*f_limb.y > 1e-12){
        f1 = (Vec3){f_limb.y, -f_limb.x, 0.};
        l = LENGTH(f1);
        RESCALE(f1, l);
        CROSS(f1, f_limb, f2);
    }

    Vec3 b_limb = (Vec3){-rot.zx/r_back, -rot.zy/r_up, -rot.zz / r_side};
    l = LENGTH(b_limb);
    if (b_limb.z < 0.0){
        l *= -1.0;
    }
    RESCALE(b_limb, l);

    Vec3 b1 = {1.0, 0.0, 0.0};
    Vec3 b2 = {0.0, 1.0, 0.0};
    if (b_limb.x*b_limb.x + b_limb.y*b_limb.y > 1e-12){
        b1 = (Vec3){b_limb.y, -b_limb.x, 0.};
        l = LENGTH(b1);
        RESCALE(b1, l);
        CROSS(b1, b_limb, b2);
    }

    // Transform the limb planes into view space
    Vec3 temp = {f_limb.x/r_forward, f_limb.y/r_up, f_limb.z/r_side};
    MATMUL(rot, temp, f_limb);

    temp = (Vec3){b_limb.x/r_back, b_limb.y/r_up, b_limb.z/r_side};
    MATMUL(rot, temp, b_limb);

    // Transform the limb vectors into view space
    temp = (Vec3){f1.x*r_forward, f1.y*r_up, f1.z*r_side};
    MATMUL(rot, temp, f1);

    temp = (Vec3){f2.x*r_forward, f2.y*r_up, f2.z*r_side};
    MATMUL(rot, temp, f2);

    temp = (Vec3){b1.x*r_back, b1.y*r_up, b1.z*r_side};
    MATMUL(rot, temp, b1);

    temp = (Vec3){b2.x*r_back, b2.y*r_up, b2.z*r_side};
    MATMUL(rot, temp, b2);

    // Get the limb breakpoints
    double d1 = rot.xx*f1.x + rot.yx*f1.y + rot.zx*f1.z;
    double d2 = rot.xx*f2.x + rot.yx*f2.y + rot.zx*f2.z;
    double cbreak = -1./sqrt(1 + d1*d1 / (d2*d2));
    double sbreak = 1./sqrt(1 + d2*d2 / (d1*d1));
    Vec3 break_offset;
    WEIGHTED_SUM(cbreak, sbreak, f1, f2, break_offset);
    if (break_offset.x < 0){
        RESCALE(break_offset, -1.);
    }

    // Get the bounds
    Bounds xbounds, ybounds;
    double f_width = sqrt(f1.x*f1.x + f2.x*f2.x);
    double b_width = sqrt(b1.x*b1.x + b2.x*b2.x);
    if (rot.xx > 0){
        xbounds = (Bounds){x - b_width, x + f_width};
    }
    else {
        xbounds = (Bounds){x - f_width, x + b_width};
    }
    f_width = sqrt(f1.y*f1.y + f2.y*f2.y);
    b_width = sqrt(b1.y*b1.y + b2.y*b2.y);
    if (rot.yy > 0){
        ybounds = (Bounds){y - b_width, y + f_width};
    }
    else {
        ybounds = (Bounds){y - f_width, y + b_width};
    }
    Biellipsoid result = {x, y, z, r_forward, r_back, r_up, r_side, rot, f_limb, b_limb, break_offset, f1, f2, b1, b2, xbounds, ybounds};
    return result;
}

Bounds get_ylimits(Biellipsoid* bell, double x){
    double st, ct;
    Vec3 e1, e2;
    Vec3 forward = {bell->rot.xx, bell->rot.yx, bell->rot.zx};
    Bounds ybounds = {0., 0.};

    if ((x <= bell->xbounds.min) || (x >= bell->xbounds.max)){
        return ybounds;
    }
    x -= bell->position.x;

    double xwidth2 = bell->f1.x*bell->f1.x + bell->f2.x*bell->f2.x;
    int aligned = bell->f1.x*bell->f1.x < xwidth2*1e-12;
    double disc = (x*x) * (bell->f2.x*bell->f2.x) - xwidth2 * (x*x - bell->f1.x * bell->f1.x);
    disc = sqrt(disc);
    // First and second forward ellipse intersections
    st = aligned? x / bell->f2.x: (x * bell->f2.x + disc) / xwidth2;
    ct = aligned? sqrt(1 - st*st): (x - bell->f2.x * st) / bell->f1.x;
    WEIGHTED_SUM(ct, st, bell->f1, bell->f2, e1);
    st = aligned? x / bell->f2.x: (x * bell->f2.x - disc) / xwidth2;
    ct = aligned? -sqrt(1 - st*st): (x - bell->f2.x * st) / bell->f1.x;
    WEIGHTED_SUM(ct, st, bell->f1, bell->f2, e2);

    // Select bounds ased on their location
    if (DOT3(e1, forward) >= 0){
        if (e1.y < ybounds.min)
            ybounds.min = e1.y;
        else if (e1.y > ybounds.max)
            ybounds.max = e1.y;
    }
    if (DOT3(e2, forward) >= 0){
        if (e2.y < ybounds.min)
            ybounds.min = e2.y;
        else if (e2.y > ybounds.max)
            ybounds.max = e2.y;
    }

    // First and second backward ellipse intersections
    xwidth2 = bell->b1.x*bell->b1.x + bell->b2.x*bell->b2.x;
    aligned = bell->b1.x*bell->b1.x < xwidth2*1e-12;
    disc = (x*x) * (bell->b2.x*bell->b2.x) - xwidth2 * (x*x - bell->b1.x * bell->b1.x);
    disc = sqrt(disc);
    st = aligned? x / bell->b2.x: (x * bell->b2.x + disc) / xwidth2;
    ct = aligned? sqrt(1 - st*st): (x - bell->b2.x * st) / bell->b1.x;
    WEIGHTED_SUM(ct, st, bell->b1, bell->b2, e1);
    st = aligned? x / bell->b2.x: (x * bell->b2.x - disc) / xwidth2;
    ct = aligned? -sqrt(1 - st*st): (x - bell->b2.x * st) / bell->b1.x;
    WEIGHTED_SUM(ct, st, bell->b1, bell->b2, e2);

    // Select bounds based on their location
    if (DOT3(e1, forward) < 0){
        if (e1.y < ybounds.min)
            ybounds.min = e1.y;
        else if (e1.y > ybounds.max)
            ybounds.max = e1.y;
    }
    if (DOT3(e2, forward) < 0){
        if (e2.y < ybounds.min)
            ybounds.min = e2.y;
        else if (e2.y > ybounds.max)
            ybounds.max = e2.y;
    }
    ybounds.min += bell->position.y;
    ybounds.max += bell->position.y;
    return ybounds;
}
