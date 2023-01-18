#include <math.h>
#include "../include/desprng.h"

#pragma acc routine(approx) seq
double approx(const double t);

#pragma acc routine(normal01) seq
double normal01(const double u);

#pragma acc routine(maxboltz) seq
double maxboltz(const double u, const double T, const double m);

#pragma acc routine(unif02pi) seq
double unif02pi(const double u);

#pragma acc routine(exponential) seq
double exponential(const double u, const double tc);