#include "../include/DOP853.h"
#include <cmath>
#include <cstring>
#include <iostream>

double DOP853::sign(double a, double b)
{
  return (b < 0.0)? -fabs(a) : fabs(a);
}

double DOP853::min_d(double a, double b)
{
  return (a < b)?a:b;
}

double DOP853::max_d(double a, double b)
{
  return (a > b)?a:b;
}

double DOP853::hinit(std::function<void(double x, std::vector<double>& y, std::vector<double>& f)> fcn, double x, std::vector<double>& y,
	    double posneg, std::vector<double>& f0, std::vector<double>& f1, std::vector<double>& yy1, int iord,
	    double hmax, double atoler, double rtoler)
{
    double   dnf, dny, atoli, rtoli, sk, h, h1, der2, der12, sqr;
    unsigned i;

    dnf = 0.0;
    dny = 0.0;
    atoli = atoler;
    rtoli = rtoler;

    for (i = 0; i < n; i++){
        sk = atoli + rtoli * fabs(y[i]);
        sqr = f0[i] / sk;
        dnf += sqr*sqr;
        sqr = y[i] / sk;
        dny += sqr*sqr;
    }

    if ((dnf <= 1.0E-10) || (dny <= 1.0E-10))
    h = 1.0E-6;
    else
    h = sqrt (dny/dnf) * 0.01;

    h = min_d(h, hmax);
    h = sign(h, posneg);

    /* perform an explicit Euler step */
    for (i = 0; i < n; i++)
        yy1[i] = y[i] + h * f0[i];
    fcn (x+h, yy1, f1);

    /* estimate the second derivative of the solution */
    der2 = 0.0;
    for (i = 0; i < n; i++){
        sk = atoli + rtoli * fabs(y[i]);
        sqr = (f1[i] - f0[i]) / sk;
        der2 += sqr*sqr;
    }
    der2 = sqrt (der2) / h;

    /* step size is computed such that h**iord * max_d(norm(f0),norm(der2)) = 0.01 */
    der12 = max_d(fabs(der2), sqrt(dnf));
    if (der12 <= 1.0E-15)
    h1 = max_d (1.0E-6, fabs(h)*1.0E-3);
    else
    h1 = pow (0.01/der12, 1.0/(double)iord);
    h = min_d (100.0 * fabs(h), min_d (h1, hmax));

    return sign (h, posneg);

}

int DOP853::integrate(double x0, double xf){
    x = x0;
    xend = xf;
    facold = 1.0E-4;
    expo1 = 1.0/8.0 - beta * 0.2;
    facc1 = 1.0 / fac1;
    facc2 = 1.0 / fac2;
    posneg = sign(1.0, xend-x);

    /* initial preparations */
    atoli = atol;
    rtoli = rtol;
    last  = 0;
    hlamb = 0.0;
    iasti = 0;
    fcn(x, y, k1);
    hmax = fabs(hmax);
    iord = 8;
    if (h == 0.0)
    h = hinit(fcn, x, y, posneg, k1, k2, k3, iord, hmax, atol, rtol);
    nfcn += 2;
    reject = 0;
    xold = x;

    if (iout){
        irtrn = 1;
        hout = 1.0;
        xout = x;
        obs(naccpt+1, xold, x, y, n, &irtrn);
    }

    while (1){

        if (nstep > nmax){
            xout = x;
            hout = h;
            return -2;
        }

        if (0.1 * fabs(h) <= fabs(x) * uround){
            xout = x;
            hout = h;
            return -3;
        }

        if ((x + 1.01*h - xend) * posneg > 0.0){
            h = xend - x;
            last = 1;
        }

        nstep++;

        /* the twelve stages */
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * c.a21 * k1[i];
        fcn (x+c.c2*h, yy1, k2);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (c.a31*k1[i] + c.a32*k2[i]);
        fcn (x+c.c3*h, yy1, k3);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (c.a41*k1[i] + c.a43*k3[i]);
        fcn (x+c.c4*h, yy1, k4);
        for (i = 0; i <n; i++)
            yy1[i] = y[i] + h * (c.a51*k1[i] + c.a53*k3[i] + c.a54*k4[i]);
        fcn (x+c.c5*h, yy1, k5);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (c.a61*k1[i] + c.a64*k4[i] + c.a65*k5[i]);
        fcn (x+c.c6*h, yy1, k6);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (c.a71*k1[i] + c.a74*k4[i] + c.a75*k5[i] + c.a76*k6[i]);
        fcn (x+c.c7*h, yy1, k7);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (c.a81*k1[i] + c.a84*k4[i] + c.a85*k5[i] + c.a86*k6[i] +
                    c.a87*k7[i]);
        fcn (x+c.c8*h, yy1, k8);
        for (i = 0; i <n; i++)
            yy1[i] = y[i] + h * (c.a91*k1[i] + c.a94*k4[i] + c.a95*k5[i] + c.a96*k6[i] +
                    c.a97*k7[i] + c.a98*k8[i]);
        fcn (x+c.c9*h, yy1, k9);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (c.a101*k1[i] + c.a104*k4[i] + c.a105*k5[i] + c.a106*k6[i] +
                    c.a107*k7[i] + c.a108*k8[i] + c.a109*k9[i]);
        fcn (x+c.c10*h, yy1, k10);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (c.a111*k1[i] + c.a114*k4[i] + c.a115*k5[i] + c.a116*k6[i] +
                    c.a117*k7[i] + c.a118*k8[i] + c.a119*k9[i] + c.a1110*k10[i]);
        fcn (x+c.c11*h, yy1, k2);
        xph = x + h;
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (c.a121*k1[i] + c.a124*k4[i] + c.a125*k5[i] + c.a126*k6[i] +
                    c.a127*k7[i] + c.a128*k8[i] + c.a129*k9[i] +
                    c.a1210*k10[i] + c.a1211*k2[i]);
        fcn (xph, yy1, k3);
        nfcn += 11;
        for (i = 0; i < n; i++){
            k4[i] = c.b1*k1[i] + c.b6*k6[i] + c.b7*k7[i] + c.b8*k8[i] + c.b9*k9[i] +
                c.b10*k10[i] + c.b11*k2[i] + c.b12*k3[i];
            k5[i] = y[i] + h * k4[i];
        }
        
        /* error estimation */
        err = 0.0;
        err2 = 0.0;

        for (i = 0; i < n; i++)
        {
            sk = atoli + rtoli * max_d (fabs(y[i]), fabs(k5[i]));
            erri = k4[i] - c.bhh1*k1[i] - c.bhh2*k9[i] - c.bhh3*k3[i];
            sqr = erri / sk;
            err2 += sqr*sqr;
            erri = c.er1*k1[i] + c.er6*k6[i] + c.er7*k7[i] + c.er8*k8[i] + c.er9*k9[i] +
                c.er10 * k10[i] + c.er11*k2[i] + c.er12*k3[i];
            sqr = erri / sk;
            err += sqr*sqr;
        }

        deno = err + 0.01 * err2;
        if (deno <= 0.0)
        deno = 1.0;
        err = fabs(h) * err * sqrt (1.0 / (deno*(double)n));

        /* computation of hnew */
        fac11 = pow (err, expo1);
        /* Lund-stabilization */
        fac = fac11 / pow(facold,beta);
        /* we require fac1 <= hnew/h <= fac2 */
        fac = max_d (facc2, min_d (facc1, fac/safe));
        hnew = h / fac;

        if (err <= 1.0)
        {
            /* step accepted */

            facold = max_d (err, 1.0E-4);
            naccpt++;
            fcn (xph, k5, k4);
            nfcn++;
            
            /* final preparation for dense output */
            if (iout == 2)
            {
            /* save the first function evaluations */
                for (i = 0; i < n; i++)
                {
                    rcont1[i] = y[i];
                    ydiff = k5[i] - y[i];
                    rcont2[i] = ydiff;
                    bspl = h * k1[i] - ydiff;
                    rcont3[i] = bspl;
                    rcont4[i] = ydiff - h*k4[i] - bspl;
                    rcont5[i] = c.d41*k1[i] + c.d46*k6[i] + c.d47*k7[i] + c.d48*k8[i] +
                        c.d49*k9[i] + c.d410*k10[i] + c.d411*k2[i] + c.d412*k3[i];
                    rcont6[i] = c.d51*k1[i] + c.d56*k6[i] + c.d57*k7[i] + c.d58*k8[i] +
                        c.d59*k9[i] + c.d510*k10[i] + c.d511*k2[i] + c.d512*k3[i];
                    rcont7[i] = c.d61*k1[i] + c.d66*k6[i] + c.d67*k7[i] + c.d68*k8[i] +
                        c.d69*k9[i] + c.d610*k10[i] + c.d611*k2[i] + c.d612*k3[i];
                    rcont8[i] = c.d71*k1[i] + c.d76*k6[i] + c.d77*k7[i] + c.d78*k8[i] +
                        c.d79*k9[i] + c.d710*k10[i] + c.d711*k2[i] + c.d712*k3[i];
                }

                /* the next three function evaluations */
                for (i = 0; i < n; i++)
                    yy1[i] = y[i] + h * (c.a141*k1[i] + c.a147*k7[i] + c.a148*k8[i] +
                                c.a149*k9[i] + c.a1410*k10[i] + c.a1411*k2[i] +
                                c.a1412*k3[i] + c.a1413*k4[i]);
                fcn (x+c.c14*h, yy1, k10);
                for (i = 0; i < n; i++)
                    yy1[i] = y[i] + h * (c.a151*k1[i] + c.a156*k6[i] + c.a157*k7[i] + c.a158*k8[i] +
                                c.a1511*k2[i] + c.a1512*k3[i] + c.a1513*k4[i] +
                                c.a1514*k10[i]);
                fcn (x+c.c15*h, yy1, k2);
                for (i = 0; i < n; i++)
                    yy1[i] = y[i] + h * (c.a161*k1[i] + c.a166*k6[i] + c.a167*k7[i] + c.a168*k8[i] +
                                c.a169*k9[i] + c.a1613*k4[i] + c.a1614*k10[i] +
                                c.a1615*k2[i]);
                fcn (x+c.c16*h, yy1, k3);
                nfcn += 3;

                /* final preparation */
                for (i = 0; i < n; i++){
                    rcont5[i] = h * (rcont5[i] + c.d413*k4[i] + c.d414*k10[i] +
                            c.d415*k2[i] + c.d416*k3[i]);
                    rcont6[i] = h * (rcont6[i] + c.d513*k4[i] + c.d514*k10[i] +
                            c.d515*k2[i] + c.d516*k3[i]);
                    rcont7[i] = h * (rcont7[i] + c.d613*k4[i] + c.d614*k10[i] +
                            c.d615*k2[i] + c.d616*k3[i]);
                    rcont8[i] = h * (rcont8[i] + c.d713*k4[i] + c.d714*k10[i] +
                            c.d715*k2[i] + c.d716*k3[i]);
                }
            }

            k1 = k4;
            y = k5;
            xold = x;
            x = xph;

            if (iout){
                hout = h;
                xout = x;
                obs(naccpt+1, xold, x, y, n, &irtrn);
                if (irtrn < 0)
                    return 2;
            }

            /* normal exit */
            if (last)
            {
                hout=hnew;
                xout = x;
                return 1;
            }

            if (fabs(hnew) > hmax)
                hnew = posneg * hmax;
            if (reject)
                hnew = posneg * min_d (fabs(hnew), fabs(h));

            reject = 0;
        }
        else{
            /* step rejected */
            hnew = h / min_d (facc1, fac11/safe);
            reject = 1;
            if (naccpt >= 1)
            nrejct=nrejct + 1;
            last = 0;
        }

        h = hnew;
    }

}