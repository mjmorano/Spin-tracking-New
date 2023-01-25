#pragma once

const double G_CONST = -9.8;

struct options{
	double Lx = 0.07;
	double Ly = 0.1;
	double Lz = 0.4;
	char dist = 'C';
	double m = 1.6e-27;
	bool gas_coll = false;
	double tc = 1e-3;
	double T = 4.2;
	bool diffuse = false;
	double gamma = -2.078e8;
	double V = 5.0;
	bool gravity = true;
	double B0 = 3e-6;
	double E = 7500;
	double t0 = 0.0;
	double tf = 1.0;
	double rtol = 1e-12;
	double atol = 1e-12;
	double beta = 0.0;
	int iout = 2;
	double ioutInt = 1.0;
	double uround = 1e-16;
	double safe = 0.9;
	double fac1 = 0.333;
	double fac2 = 6.0;
	double hmax = 1.0;
	double h = 0.001;
	unsigned int nmax = 10000000;
	double max_step = 0.001;
};
