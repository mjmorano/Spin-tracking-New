#pragma once

#include "double3.h"

const double G_CONST = -9.8;
const double c2 = 299792458.0 * 299792458.0;

struct options{
	double3 L = {0.07, 0.1, 0.4};
	char dist = 'C';
	double m = 2.2*5e-27;
	bool gas_coll = true;
	double tc = 1e-4;
	double T = 4.2;
	bool diffuse = true;
	double gamma = -2.078e8;
	double V = 5.0;
	bool gravity = true;
	//these are x, y, z coordinates
	double3 B0 = {0.0, 0.0, 3e-6};
	double3 E = {0, 0, 75e5};
	double t0 = 0.0;
	double tf = 10.0;
	double rtol = 1e-12;
	double atol = 1e-12;
	double beta = 0.0;
	int iout = 2;
	double ioutInt = 0.05;
	double uround = 1e-16;
	double safe = 0.9;
	double fac1 = 0.333;
	double fac2 = 6.0;
	double hmax = 1.0;
	double h = 0.001;
	unsigned int nmax = 10000000;
	int integratorType = 0; //0 means DOP853, 1 means hybrid RK45 approach
	double swapStepSize = 1.0-4; //above this use rotations, below this use standard RK techniques
	int numParticles = 1000;
	int numPerGPUBlock = 128;
};
