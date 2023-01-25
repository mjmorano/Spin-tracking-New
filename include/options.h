#pragma once

const double G_CONST = -9.8;

struct options{
	const double Lx = 0.07;
	const double Ly = 0.1;
	const double Lz = 0.4;
	const char dist = 'C';
	const double m = 1.6e-27;
	const bool gas_coll = false;
	const double tc = 1e-3;
	const double T = 4.2;
	const bool diffuse = false;
	const double gamma = -2.078e8;
	const double V = 5.0;
	const bool gravity = true;
	const double B0 = 3e-6;
	const double E = 7500;
	const double t0 = 0.0;
	const double tf = 1.0;
	const double rtol = 1e-12;
	const double atol = 1e-12;
	const double beta = 0.0;
	const int iout = 2;
	const double ioutInt = 0.5;
	const double uround = 1e-16;
	const double safe = 0.9;
	const double fac1 = 0.333;
	const double fac2 = 6.0;
	const double hmax = 1.0;
	const double h = 0.001;
	const unsigned int nmax = 10000000;
	const double max_step = 0.001;
};
