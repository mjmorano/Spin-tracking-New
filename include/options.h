#pragma once

struct options{
	double Lx = 0.5;
	double Ly = 0.5;
	double Lz = 0.5;
	char dist = 'C';
	double m = 1.6e-27;
	bool gas_coll = true;
	double tc = 1e-3;
	double T = 4.2;
    bool diffuse = true;
    double gamma = -2.078e8;
    double V = 5.0;
    bool gravity = false;
    double B0 = 3e-6;
    double E = 7500;
	double t0 = 0.0;
    double tf = 10.0;
    double rtol = 1e-12;
    double atol = 1e-12;
    double beta = 0.0;
    int iout = 1;
    double uround = 1e-16;
    double safe = 0.9;
    double fac1 = 0.333;
    double fac2 = 6.0;
    double hmax = 1.0;
    double h = 0.0;
    unsigned int nmax = 10000000;
    unsigned int n = 3;
};