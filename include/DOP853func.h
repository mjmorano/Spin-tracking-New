#pragma once
#include <stdio.h>
#include <cmath>
#include "../include/double3.h"
#include "../include/options.h"

void obs(long nr, double xold, double x, double3 y, int* irtrn, 
	options OPT, double* lastOutput, unsigned int* lastIndex, float* outputArray);		
			
double pulse(const double t);

void Bloch(const double t, const double3& y, double3& f, const double B0, const double E, const double gamma, const double t0, const double tf ,const double3& p_old, const double3& p_new, const double3& v_old, const double3& v_new);

void interpolate(const double t, const double t0, const double tf, const double3& p_old, const double3& p_new, const double3& v_old, const double3& v_new, double3& p_out, double3& v_out);

void grad(double3&, double3&);

// double hinit(double, double*, double, double*, double*, double*, int, double, double, double);

int integrate(double t0, double tf, double3& y, const double3& p_old, const double3& p_new, 
	const double3& v_old, const double3& v_new, options OPT,
	double& lastOutput, unsigned int& lastIndex, float* outputArray);

double sign(double, double);

double min_d(double, double);

double max_d(double, double);

