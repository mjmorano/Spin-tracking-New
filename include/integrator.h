#pragma once
#include <stdio.h>
#include <cmath>
#include "../include/double3.h"
#include "../include/options.h"
#include "../include/coeff.h"

#if __HIPCC__
#include <hip/hip_runtime.h>
#include <hiprand/hiprand.h>
#include <hiprand/hiprand_kernel.h>
#define __PREPROC__ __device__
#elif __NVCC__
#define __PREPROC__ __device__
#else
#include <random>
#define __PREPROC__
#endif

__PREPROC__ void obs(long nr, double xold, double x, double3 y, double3 pos, int* irtrn, 
	options OPT, double* lastOutput, unsigned int* lastIndex, outputDtype* outputArray);		

__PREPROC__ double3 pulse(const double t);

__PREPROC__ double3 findCrossTerm(const double t, const double3& y, const double B0, const double E, const double gamma, const double t0, const double tf ,const double3& p_old, const double3& p_new, const double3& v_old, const double3& v_new);

__PREPROC__ void Bloch(const double t, const double3& y, double3& f, const double B0, const double E, const double gamma, const double t0, const double tf ,const double3& p_old, const double3& p_new, const double3& v_old, const double3& v_new);

__PREPROC__ void interpolate(const double t, const double t0, const double tf, const double3& p_old, const double3& p_new, const double3& v_old, const double3& v_new, double3& p_out, double3& v_out);

__PREPROC__ double3 grad(double3&);

// double hinit(double, double*, double, double*, double*, double*, int, double, double, double);
__PREPROC__ int integrateDOP(double t0, double tf, double3& y, const double3& p_old, const double3& p_new, 
	const double3& v_old, const double3& v_new, options OPT,
	double& lastOutput, unsigned int& lastIndex, outputDtype* outputArray);

__PREPROC__ int integrateRK45Hybrid(double t0, double tf, double3& y, const double3& p_old, const double3& p_new, 
	const double3& v_old, const double3& v_new, options OPT, double& h,
	double& lastOutput, unsigned int& lastIndex, outputDtype* outputArray);

__PREPROC__ double sign(double, double);

__PREPROC__ double min_d(double, double);

__PREPROC__ double max_d(double, double);

