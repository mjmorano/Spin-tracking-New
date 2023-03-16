#include "../include/dists.h"

double approx(const double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) / (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double normal01(const double u){
    if (u < 0.5)
        return -approx( sqrt(-2.0 * log(u) ) );
    else
        return approx( sqrt(-2.0 * log(1.0-u) ) );
}

double maxboltz(const double u, const double sqrtkT_m){
    return normal01(u) * sqrtkT_m;
}

double unif02pi(const double u){
    return u * 2.0 * M_PI;
}

double exponential(const double u, const double tc){
    return - tc * log(1.0 - u);
}
