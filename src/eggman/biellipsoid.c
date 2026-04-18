#include <math.h>

// 3-vector macros
#define ASSIGN(arr, a, b, c) arr[0] = a; arr[1] = b; arr[2] = c;
#define LENGTH(arr) sqrt(arr[0]*arr[0] + arr[1]*arr[1] + arr[2]*arr[2])
#define SCALE(arr, coeff) arr[0] /= coeff; arr[1] /= coeff; arr[2] /= coeff;
#define DOT3(a, b) a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
#define MATMUL(mat, vec, out) \
out[0] = mat[0]*vec[0] + mat[1]*vec[1] + mat[2]*vec[2];\
out[1] = mat[3]*vec[0] + mat[4]*vec[1] + mat[5]*vec[2];\
out[2] = mat[6]*vec[0] + mat[7]*vec[1] + mat[8]*vec[2];
#define CROSS(a, b, out) \
out[0] = a[1]*b[2] - b[1]*a[2]; \
out[1] = a[2]*b[0] - b[2]*a[0]; \
out[2] = a[0]*b[1] - b[0]*a[1];
#define PRINT(vec) printf("%f, %f, %f\n", vec[0], vec[1], vec[2]);


typedef struct{
    double position[3];
    // Rotation matrix of the biellipsoid (forward vector is first column (rot[::3])
    double rot[9];
    // Basis vectors 1 and 2 for the bounding ellipsoids f and b
    double f1[3];
    double f2[3];
    double b1[3];
    double b2[3];
    // View-axis-aligned bounding box relative to position
    double bounds[4];
} Biellipsoid;


typedef struct{
    double min;
    double max;
} Bounds;


void init_biellipsoid(Biellipsoid* f, double x, double y, double z, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side){
    double plane[3];
    double l;
    double ct, st, cp, sp, cg, sg;
    double temp[3];

    ASSIGN(f->position, x, y, z);

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
    f->rot[0] = cp * ct;
    f->rot[1] = -cp * st * cg + sp*sg;
    f->rot[2] = cp*st*sg + sp*cg;
    f->rot[3] = st;
    f->rot[4] = ct * cg;
    f->rot[5] = -ct * sg;
    f->rot[6] = -sp * ct;
    f->rot[7] = sp*st*cg + cp*sg;
    f->rot[8] = -sp * st * sg + cp*cg;

    // Calculate limb vectors
    ASSIGN(plane, f->rot[2] / r_forward, f->rot[5] / r_up, f->rot[8] / r_side);
    l = LENGTH(plane);
    SCALE(plane, l);
    if (plane[0]*plane[0] + plane[1]*plane[1] < 1e-12){
        ASSIGN(f->f1, 1.0, 0.0, 0.0);
        ASSIGN(f->f2, 0.0, 1.0, 0.0);
    }
    else {
        ASSIGN(f->f1, plane[1], -plane[2], 0.)
        l = LENGTH(f->f1);
        SCALE(f->f1, l);
        CROSS(f->f1, plane, f->f2);
    }

    ASSIGN(plane, f->rot[2] / r_back, f->rot[5] / r_up, f->rot[8] / r_side);
    l = LENGTH(plane);
    SCALE(plane, l);
    if (plane[0]*plane[0] + plane[1]*plane[1] < 1e-12){
        ASSIGN(f->b1, 1.0, 0.0, 0.0);
        ASSIGN(f->b2, 0.0, 1.0, 0.0);
    }
    else {
        ASSIGN(f->b1, plane[1], -plane[2], 0.)
        l = LENGTH(f->b1);
        SCALE(f->b1, l);
        CROSS(f->b1, plane, f->b2);
    }

    // Transform the limb vectors into view space
    ASSIGN(temp, f->f1[0]/r_forward, f->f1[1]/r_up, f->f1[2]/r_back);
    MATMUL(f->rot, temp, f->f1);

    ASSIGN(temp, f->f2[0]/r_forward, f->f2[1]/r_up, f->f2[2]/r_back);
    MATMUL(f->rot, temp, f->f2);

    ASSIGN(temp, f->b1[0]/r_forward, f->b1[1]/r_up, f->b1[2]/r_back);
    MATMUL(f->rot, temp, f->b1);

    ASSIGN(temp, f->b2[0]/r_forward, f->b2[1]/r_up, f->b2[2]/r_back);
    MATMUL(f->rot, temp, f->b2);

    // Get the bounds
    f->bounds[0] = fabs(sqrt(f->f1[0]*f->f1[0] + f->f2[0]*f->f2[0]));
    f->bounds[1] = -f->bounds[0];
    f->bounds[2] = fabs(sqrt(f->b1[0]*f->f1[0] + f->b2[0]*f->b2[0]));
    f->bounds[3] = -f->bounds[2];
    return;
}

Bounds get_ybounds(Biellipsoid* bell, double x){
    int i = 0;
    double xf1[3];
    double xf2[3];
    double xb1[3];
    double xb2[3];

    double xwidth2 = bell->f1[0]*bell->f1[0] + bell->f2[0]*bell->f2[0];
    // Quadratic equation discriminant
    double disc = bell->f1[0] * bell->f1[0];
    disc = disc*disc - disc * x*x + disc * bell->f2[0] * bell->f2[0];
    double B = (x * bell->f2[0] + sqrt(disc)) / xwidth2;
    double A = (x - bell->f2[0] * B) / bell->f1[0];
    for(int i = 0; i < 3; i++){
        xf1[i] = A * bell->f1[i] + B * bell->f2[i];
    }
    B = (x * bell->f2[0] - sqrt(disc)) / xwidth2;
    A = (x - bell->f2[0] * B) / bell->f1[0];
    for(int i = 0; i < 3; i++){
        xf2[i] = A * bell->f1[i] + B * bell->f2[i];
    }
    B = (x * bell->f2[0] - sqrt(disc)) / xwidth2;
    A = (x - bell->f2[0] * B) / bell->f1[0];
    for(int i = 0; i < 3; i++){
        xb1[i] = A * bell->f1[i] + B * bell->f2[i];
    }
    B = (x * bell->f2[0] + sqrt(disc)) / xwidth2;
    A = (x - bell->f2[0] * B) / bell->f1[0];
    for(int i = 0; i < 3; i++){
        xb2[i] = A * bell->f1[i] + B * bell->f2[i];
    }

    // Select bounds based on their location
    double ymin=0.;
    double ymax=0;
    double forward[3] = {bell->rot[0], bell->rot[3], bell->rot[5]};
    if (DOT3(xf1, forward) >= 0){
        if (xf1[1] < ymin)
            ymin = xf1[1];
        else if (xf1[1] > ymax)
            ymax = xf1[1];
    }
    if (DOT3(xf2, forward) >= 0){
        if (xf2[1] < ymin)
            ymin = xf2[1];
        else if (xf2[1] > ymax)
            ymax = xf2[1];
    }
    if (DOT3(xb1, forward) >= 0){
        if (xb1[1] < ymin)
            ymin = xb1[1];
        else if (xb1[1] > ymax)
            ymax = xb1[1];
    }
    if (DOT3(xb2, forward) >= 0){
        if (xb2[1] < ymin)
            ymin = xb2[1];
        else if (xb2[1] > ymax)
            ymax = xb2[1];
    }
    Bounds b = {ymin, ymax};
    return b;
}
