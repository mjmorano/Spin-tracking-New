/* DOP853

This is my adaptation of the DOP853 code from J .Colinge (COLINGE@DIVSUN.UNIGE.CH) which in turn is an adaptation of E. Hairer's
original FORTRAN code. The main changes are the removal of pointer arrays for passing between functions
and using std::vector instead. The other main change is that the intergrator is now a class instead of a
function which relied on delcaring a bunch of static variables. I did this because it should make it cleaner
to use and expand in the future, and becuase the static variables will create huge problems when trying to
parallelize the spin simulation this is intended for. In addition to modernizing the C code I took the liberty
to remove certain parts that we don't really need for our purposes. Most of the error checking is still present
but more could probably be added.

Author: Matt Morano
Date: 1/4/2023

Credit for the original implementation goes to:

    Authors : E. Hairer & G. Wanner
      Universite de Geneve, dept. de Mathematiques
      CH-1211 GENEVE 4, SWITZERLAND
      e-mail : HAIRER@DIVSUN.UNIGE.CH, WANNER@DIVSUN.UNIGE.CH

The code is described in : E. Hairer, S.P. Norsett and G. Wanner, Solving
ordinary differential equations I, nonstiff problems, 2nd edition,
Springer Series in Computational Mathematics, Springer-Verlag (1993).

*/

#pragma once
#include "../include/coeff.h"
#include "../include/options.h"
#include <functional>
#include <vector>
#include <iostream>


// typedef void (*FCN)(const double x, const std::vector<double> y, std::vector<double>& f);
typedef void (*OBS)(long nr, double xold, double x, std::vector<double>& y, unsigned int n, int* irtrn);
// using std::function<void(double x, std::vector<double>& y, std::vector<double>& f)> FCN;

class DOP853
{
public:

    DOP853() {}

    DOP853(std::function<void(double x, std::vector<double>& y, std::vector<double>& f)> fcn, std::vector<double> y0, OBS obs, options OPT)
    : n(3), y(y0), fcn(fcn), rtol(OPT.rtol), atol(OPT.atol), obs(obs), iout(OPT.iout), uround(OPT.uround), safe(OPT.safe), 
      fac1(OPT.fac1), fac2(OPT.fac2), beta(OPT.beta), hmax(OPT.hmax), h(OPT.h), nmax(OPT.nmax), c(coeffs()), yy1(n), k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n), k8(n), 
      k9(n), k10(n), rcont1(n), rcont2(n), rcont3(n), rcont4(n), rcont5(n), rcont6(n), rcont7(n), rcont8(n) 
      {
        //printf("%f\t %f\t %f\n",y[0],y[1],y[2]);
      }

    int integrate(double x0, double xf);

    double sign(double a, double b);

    double min_d(double a, double b);

    double max_d(double a, double b);

    double hinit(std::function<void(double x, std::vector<double>& y, std::vector<double>& f)> fcn, double x, std::vector<double>& y,
      double posneg, std::vector<double>& f0, std::vector<double>& f1, std::vector<double>& yy1, int iord,
      double hmax, double atoler, double rtoler);

private:
    unsigned int n;
    std::function<void(double x, std::vector<double>& y, std::vector<double>& f)> fcn;
    std::vector<double> y;
    double x;
    double xend;
    double rtol;
    double atol;
    OBS obs;
    int iout;
    double uround;
    double safe;
    double fac1;
    double fac2;
    double beta;
    double hmax;
    double h;
    long nmax;
    double facold, expo1, fac, facc1, facc2, fac11, posneg, xph;
    double atoli, rtoli, hlamb, err, sk, hnew, yd0, ydiff, bspl;
    double stnum, stden, sqr, err2, erri, deno;
    int iasti, iord, irtrn, reject, last, nonsti;
    unsigned int i, j;
    unsigned int nfcn = 0, nstep = 0, naccpt = 0, nrejct = 0;
    unsigned int *indir;
    std::vector<double> yy1, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10;
    std::vector<double> rcont1, rcont2, rcont3, rcont4;
    std::vector<double> rcont5, rcont6, rcont7, rcont8;
    double hout, xold, xout;
    int arret, idid;
    coeffs c;
};