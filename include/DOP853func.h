#pragma once
#include <stdio.h>
#include <cmath>
#include "../include/options.h"

void obs(long nr, double xold, double x, double* y, int* irtrn);

double pulse(const double t);

void Bloch(const double t, const double* y, double* f, const double B0, const double gamma, const double t0, const double tf ,const double* x_old, const double* x_new, const double* v_old, const double* v_new);

void interpolate(const double t, const double t0, const double tf, const double* x_old, const double* x_new, const double* v_old, const double* v_new, double* p_out, double* v_out);

void grad(double*, double*);

// double hinit(double, double*, double, double*, double*, double*, int, double, double, double);

int integrate(double t0, double tf, double* y0, const double* x_old, const double* x_new, const double* v_old, const double* v_new, options OPT);

double sign(double, double);

double min_d(double, double);

double max_d(double, double);