#pragma once

const double G_CONST = -9.8;
const double c2 = 299792458.0 * 299792458.0;

struct options{
	const double Lx = 0.5;
	const double Ly = 0.5;
	const double Lz = 0.5;
	char dist = 'M';
	double m = 2.2*5.00823452e-27;
	bool gas_coll = false;
	const double T = 0.3;
	bool diffuse = true;
	double gamma = -2.037894569e8;
	double V = 5.0;
	const bool gravity = false;
	const double B0 = 3e-6;
	const double E = 0.0;
	const double t0 = 0.0;
	const double tf = 100.0;
	const double rtol = 1e-12;
	const double atol = 1e-12;
	const double beta = 0.0;
	const int iout = 2;
	const double ioutInt = 1.0;
	const double uround = 2.2e-16;
	const double safe = 0.9;
	const double fac1 = 0.333;
	const double fac2 = 6.0;
	const double hmax = 1.0;
	const double h = 0.001;
	const unsigned int nmax = 10000000;
	const double max_step = 0.001;
};
