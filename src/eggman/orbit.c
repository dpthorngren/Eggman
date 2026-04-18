#include <math.h>


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


void init_orbit(Orbit* orb, double period, double semimajor, double eccen, double inclination, double lon_periapse){
    /*Notation:
        t: time...
        ta: true anomaly...
        ea: eccentric anomaly...
        ma: mean anomaly...
        [no suffix]: ...at the requested time t
        _p: ...at periapse
        _t: ...at transit (x=0, z>0)
        x,y,z: position...
        _rot: ...in the inclined, rotated frame
        _inc: ...in the inclined frame
        [no suffix]: ...in the view frame
    */
    // Basic data
    orb->period = period;
    orb->semimajor = semimajor;
    orb->eccen = eccen;
    // Precompute inclination sin and cosine
    inclination *= M_PI/180.;
    orb->s_inc = sin(inclination);
    orb->c_inc = cos(inclination);
    // Precompute periapse data
    double ta_p = (lon_periapse-90.)*M_PI/180.;
    orb->sta_p = sin(ta_p);
    orb->cta_p = cos(ta_p);
    double ea_p = atan2(sqrt(1 - orb->eccen*orb->eccen) * orb->sta_p, eccen + orb->cta_p);
    double ma_p = ea_p - eccen * sin(ea_p);
    orb->t_p = ma_p * orb->period / (2*M_PI);
    return;
}


double solve_kepler(double mean_anomaly, double eccen){
    double e_anom = mean_anomaly;
    double e_new = mean_anomaly;
    int i = 0;
    if(eccen > 0.95){
        // Direct iteration
        for(i=0; i < 50; i++){
            e_new = mean_anomaly + eccen * sin(e_anom);
            if(fabs(e_new - e_anom) < 1e-12)
                break;
            e_anom = e_new;
        }
    }
    if (eccen > 0.){
        // Newton-Raphson
        for(i=0; i < 50; i++){
            e_new -= (e_anom - eccen*sin(e_anom) - mean_anomaly) / (1 - eccen * cos(e_anom));
            if (fabs(e_new - e_anom) < 1e-12)
                break;
            e_anom = e_new;
        }
    }
    return e_new;
}


void get_position(Orbit* orb, double t, double* x, double* y, double* z){
    // Solve Kepler's equation in the rotated and inclined frame, relative to inferior conjunction
    double ma = 2*M_PI * (t - orb->t_p) / orb->period;
    double ea = solve_kepler(ma, orb->eccen);
    // Position in the rotated frame where lon_periapse = 0, inclination = 0
    double x_rot = orb->semimajor * sqrt(1 - orb->eccen*orb->eccen) * sin(ea);
    double y_rot = orb->semimajor * (cos(ea) - orb->eccen);
    // Transform back to frame where lon_periapse != 0 (still inclination = 0)
    double x_inc = x_rot * orb->cta_p + y_rot * orb->sta_p;
    double y_inc = x_rot * -orb->sta_p + y_rot * orb->cta_p;
    // Transform back to the view frame (x = right, y = up, z = towards observer) centered on star
    *x = x_inc;
    *y = y_inc * orb->c_inc;
    *z = y_inc * orb->s_inc;
    return;
}
