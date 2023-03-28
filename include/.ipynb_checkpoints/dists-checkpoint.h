#include <math.h>
#include <cmath>
#include "../include/desprng.h"

#if(_OPENMP)
#pragma omp declare target
#elif(_OPENACC)
#pragma acc routine seq
#endif
double approx(const double t);

#if(_OPENMP)
#pragma omp declare target
#elif(_OPENACC)
#pragma acc routine seq
#endif
double normal01(const double u);

#if(_OPENMP)
#pragma omp declare target
#elif(_OPENACC)
#pragma acc routine seq
#endif
double maxboltz(const double u, const double sqrtkT_m);

#if(_OPENMP)
#pragma omp declare target
#elif(_OPENACC)
#pragma acc routine seq
#endif
double unif02pi(const double u);

#if(_OPENMP)
#pragma omp declare target
#elif(_OPENACC)
#pragma acc routine seq
#endif
double exponential(const double u, const double tc);