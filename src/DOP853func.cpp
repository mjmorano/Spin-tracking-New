#include "../include/DOP853func.h"
#include "../include/options.h"
#include "../include/coeff.h"

double sign(double a, double b)
{
  return (b < 0.0)? -fabs(a) : fabs(a);
}

double min_d(double a, double b)
{
  return (a < b)?a:b;
}

double max_d(double a, double b)
{
  return (a > b)?a:b;
}

double pulse(const double t){
    // return 38.7511230e-6*cos(6000*t);
    return 0.0;
}

void obs(long nr, double xold, double x, double* y, int* irtrn){
    printf("%f\t %f\n", x, y[0]);
}

void grad(double* pos, double* G){
    G[0] = 0.0;
    G[1] = 0.0;
    G[2] = 0.0;
}

void interpolate(const double t, const double t0, const double tf, const double* x_old, const double* x_new, const double* v_old, const double* v_new, double* p_out, double* v_out){
	p_out[0] = (x_old[0]*(tf - t) + x_new[0]*(t - t0))/(tf-t0);
	p_out[1] = (x_old[1]*(tf - t) + x_new[1]*(t - t0))/(tf-t0);
	p_out[2] = (x_old[2]*(tf - t) + x_new[2]*(t - t0))/(tf-t0);
	// printf("%f\t %f\t %f\n", x_old[0], p_out[0], x_new[0]);
	v_out[0] = v_old[0];
	v_out[1] = v_old[1];
	v_out[2] = v_old[2];
}

void Bloch(const double t, const double* y, double* f, const double B0, const double gamma, const double t0, const double tf , const double* x_old, const double* x_new, const double* v_old, const double* v_new){
    double p[3], v[3], G[3], Bx, By, Bz;
	interpolate(t,t0,tf,x_old,x_new,v_old,v_new,p,v);
	Bx = pulse(t);
	grad(p,G);
	Bx += G[0];
	By = G[1];
	Bz = B0 + G[2];
	// printf("%f\t %f\t %f\n",Bx,By,Bz);
    f[0] = gamma * (y[1]*Bz - y[2]*By);
	f[1] = gamma * (y[2]*Bx - y[0]*Bz);
	f[2] = gamma * (y[0]*By - y[1]*Bx);
}

// double hinit(double x, double* y, double posneg, double* f0, double* f1, double* yy1, int iord, options OPT)
// {
//     double dnf, dny, atoli, rtoli, sk, h, h1, der2, der12, sqr;
//     unsigned i;
//     int n = 3;

//     dnf = 0.0;
//     dny = 0.0;
//     atoli = OPT.atol;
//     rtoli = OPT.rtol;

//     for (i = 0; i < n; i++){
//         sk = atoli + rtoli * fabs(y[i]);
//         sqr = f0[i] / sk;
//         dnf += sqr*sqr;
//         sqr = y[i] / sk;
//         dny += sqr*sqr;
//     }

//     if ((dnf <= 1.0E-10) || (dny <= 1.0E-10))
//     h = 1.0E-6;
//     else
//     h = sqrt (dny/dnf) * 0.01;

//     h = min_d(h, OPT.hmax);
//     h = sign(h, posneg);

//     /* perform an explicit Euler step */
//     for (i = 0; i < n; i++)
//         yy1[i] = y[i] + h * f0[i];
//     Bloch (x+h, yy1, f1);

//     /* estimate the second derivative of the solution */
//     der2 = 0.0;
//     for (i = 0; i < n; i++){
//         sk = atoli + rtoli * fabs(y[i]);
//         sqr = (f1[i] - f0[i]) / sk;
//         der2 += sqr*sqr;
//     }
//     der2 = sqrt (der2) / h;

//     /* step size is computed such that h**iord * max_d(norm(f0),norm(der2)) = 0.01 */
//     der12 = max_d(fabs(der2), sqrt(dnf));
//     if (der12 <= 1.0E-15)
//     h1 = max_d (1.0E-6, fabs(h)*1.0E-3);
//     else
//     h1 = pow (0.01/der12, 1.0/(double)iord);
//     h = min_d (100.0 * fabs(h), min_d (h1, OPT.hmax));

//     return sign (h, posneg);

// }

int integrate(double t0, double tf, double* y0, const double* x_old, const double* x_new, const double* v_old, const double* v_new, options OPT){

    double* y = y0;
    double yy1[3], k1[3], k2[3], k3[3], k4[3], k5[3], k6[3], k7[3], k8[3], k9[3], k10[3];
    double rcont1[3], rcont2[3], rcont3[3], rcont4[3], rcont5[3], rcont6[3], rcont7[3], rcont8[3];
    int arret, idid;
    int iasti, iord, irtrn, reject, last, nonsti;
    double facold, expo1, fac, facc1, facc2, fac11, posneg, xph;
    double stnum, stden, sqr, err2, erri, deno;
    double atoli, rtoli, hlamb, err, sk, hnew, ydiff, bspl;
    unsigned int nfcn = 0, nstep = 0, naccpt = 0, nrejct = 0;
    double hout, xold, xout;
    double x = t0;
    double xf = tf;
    double h = OPT.h;
    unsigned int i;
    int n = 3;

    facold = 1.0E-4;
    expo1 = 1.0/8.0 - OPT.beta * 0.2;
    facc1 = 1.0 / OPT.fac1;
    facc2 = 1.0 / OPT.fac2;
    posneg = sign(1.0, tf-t0);

    /* initial preparations */
    atoli = OPT.atol;
    rtoli = OPT.rtol;
    last  = 0;
    hlamb = 0.0;
    iasti = 0;
    Bloch(x, y, k1, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
    double hmax = fabs(OPT.hmax);
    iord = 8;
    // if (OPT.h == 0.0)
    //     h = hinit(fcn, x0, y, posneg, k1, k2, k3, iord, hmax, OPT.atol, OPT.rtol);
    nfcn += 2;
    reject = 0;
    xold = x;

    if (OPT.iout){
        irtrn = 1;
        hout = 1.0;
        xout = t0;
        obs(naccpt+1, xold, x, y, &irtrn);
    }

    while (1){

        if (nstep > OPT.nmax){
            xout = x;
            hout = h;
            return -2;
        }

        if (0.1 * fabs(h) <= fabs(x) * OPT.uround){
            xout = x;
            hout = h;
            return -3;
        }

        if ((x + 1.01*h - xf) * posneg > 0.0){
            h = xf - x;
            last = 1;
        }

        nstep++;

        /* the twelve stages */
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * COEF::a21 * k1[i];
        Bloch(x+COEF::c2*h, yy1, k2, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (COEF::a31*k1[i] + COEF::a32*k2[i]);
        Bloch(x+COEF::c3*h, yy1, k3, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (COEF::a41*k1[i] + COEF::a43*k3[i]);
        Bloch(x+COEF::c4*h, yy1, k4, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        for (i = 0; i <n; i++)
            yy1[i] = y[i] + h * (COEF::a51*k1[i] + COEF::a53*k3[i] + COEF::a54*k4[i]);
        Bloch(x+COEF::c5*h, yy1, k5, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (COEF::a61*k1[i] + COEF::a64*k4[i] + COEF::a65*k5[i]);
        Bloch(x+COEF::c6*h, yy1, k6, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (COEF::a71*k1[i] + COEF::a74*k4[i] + COEF::a75*k5[i] + COEF::a76*k6[i]);
        Bloch(x+COEF::c7*h, yy1, k7, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (COEF::a81*k1[i] + COEF::a84*k4[i] + COEF::a85*k5[i] + COEF::a86*k6[i] +
                    COEF::a87*k7[i]);
        Bloch(x+COEF::c8*h, yy1, k8, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        for (i = 0; i <n; i++)
            yy1[i] = y[i] + h * (COEF::a91*k1[i] + COEF::a94*k4[i] + COEF::a95*k5[i] + COEF::a96*k6[i] +
                    COEF::a97*k7[i] + COEF::a98*k8[i]);
        Bloch(x+COEF::c9*h, yy1, k9, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (COEF::a101*k1[i] + COEF::a104*k4[i] + COEF::a105*k5[i] + COEF::a106*k6[i] +
                    COEF::a107*k7[i] + COEF::a108*k8[i] + COEF::a109*k9[i]);
        Bloch(x+COEF::c10*h, yy1, k10, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (COEF::a111*k1[i] + COEF::a114*k4[i] + COEF::a115*k5[i] + COEF::a116*k6[i] +
                    COEF::a117*k7[i] + COEF::a118*k8[i] + COEF::a119*k9[i] + COEF::a1110*k10[i]);
        Bloch(x+COEF::c11*h, yy1, k2, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        xph = x + h;
        for (i = 0; i < n; i++)
            yy1[i] = y[i] + h * (COEF::a121*k1[i] + COEF::a124*k4[i] + COEF::a125*k5[i] + COEF::a126*k6[i] +
                    COEF::a127*k7[i] + COEF::a128*k8[i] + COEF::a129*k9[i] +
                    COEF::a1210*k10[i] + COEF::a1211*k2[i]);
        Bloch(xph, yy1, k3, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
        nfcn += 11;
        for (i = 0; i < n; i++){
            k4[i] = COEF::b1*k1[i] + COEF::b6*k6[i] + COEF::b7*k7[i] + COEF::b8*k8[i] + COEF::b9*k9[i] +
                COEF::b10*k10[i] + COEF::b11*k2[i] + COEF::b12*k3[i];
            k5[i] = y[i] + h * k4[i];
        }
        
        /* error estimation */
        err = 0.0;
        err2 = 0.0;

        for (i = 0; i < n; i++)
        {
            sk = atoli + rtoli * max_d (fabs(y[i]), fabs(k5[i]));
            erri = k4[i] - COEF::bhh1*k1[i] - COEF::bhh2*k9[i] - COEF::bhh3*k3[i];
            sqr = erri / sk;
            err2 += sqr*sqr;
            erri = COEF::er1*k1[i] + COEF::er6*k6[i] + COEF::er7*k7[i] + COEF::er8*k8[i] + COEF::er9*k9[i] +
                COEF::er10 * k10[i] + COEF::er11*k2[i] + COEF::er12*k3[i];
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
        fac = fac11 / pow(facold,OPT.beta);
        /* we require fac1 <= hnew/h <= fac2 */
        fac = max_d (facc2, min_d (facc1, fac/OPT.safe));
        hnew = h / fac;

        if (err <= 1.0)
        {
            /* step accepted */

            facold = max_d (err, 1.0E-4);
            naccpt++;
            Bloch(xph, k5, k4, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
            nfcn++;
            
            /* final preparation for dense output */
            if (OPT.iout == 2)
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
                    rcont5[i] = COEF::d41*k1[i] + COEF::d46*k6[i] + COEF::d47*k7[i] + COEF::d48*k8[i] +
                        COEF::d49*k9[i] + COEF::d410*k10[i] + COEF::d411*k2[i] + COEF::d412*k3[i];
                    rcont6[i] = COEF::d51*k1[i] + COEF::d56*k6[i] + COEF::d57*k7[i] + COEF::d58*k8[i] +
                        COEF::d59*k9[i] + COEF::d510*k10[i] + COEF::d511*k2[i] + COEF::d512*k3[i];
                    rcont7[i] = COEF::d61*k1[i] + COEF::d66*k6[i] + COEF::d67*k7[i] + COEF::d68*k8[i] +
                        COEF::d69*k9[i] + COEF::d610*k10[i] + COEF::d611*k2[i] + COEF::d612*k3[i];
                    rcont8[i] = COEF::d71*k1[i] + COEF::d76*k6[i] + COEF::d77*k7[i] + COEF::d78*k8[i] +
                        COEF::d79*k9[i] + COEF::d710*k10[i] + COEF::d711*k2[i] + COEF::d712*k3[i];
                }

                /* the next three function evaluations */
                for (i = 0; i < n; i++)
                    yy1[i] = y[i] + h * (COEF::a141*k1[i] + COEF::a147*k7[i] + COEF::a148*k8[i] +
                                COEF::a149*k9[i] + COEF::a1410*k10[i] + COEF::a1411*k2[i] +
                                COEF::a1412*k3[i] + COEF::a1413*k4[i]);
                Bloch(x+COEF::c14*h, yy1, k10, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
                for (i = 0; i < n; i++)
                    yy1[i] = y[i] + h * (COEF::a151*k1[i] + COEF::a156*k6[i] + COEF::a157*k7[i] + COEF::a158*k8[i] +
                                COEF::a1511*k2[i] + COEF::a1512*k3[i] + COEF::a1513*k4[i] +
                                COEF::a1514*k10[i]);
                Bloch(x+COEF::c15*h, yy1, k2, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
                for (i = 0; i < n; i++)
                    yy1[i] = y[i] + h * (COEF::a161*k1[i] + COEF::a166*k6[i] + COEF::a167*k7[i] + COEF::a168*k8[i] +
                                COEF::a169*k9[i] + COEF::a1613*k4[i] + COEF::a1614*k10[i] +
                                COEF::a1615*k2[i]);
                Bloch(x+COEF::c16*h, yy1, k3, OPT.B0, OPT.gamma, t0, tf, x_old, x_new, v_old, v_new);
                nfcn += 3;

                /* final preparation */
                for (i = 0; i < n; i++){
                    rcont5[i] = h * (rcont5[i] + COEF::d413*k4[i] + COEF::d414*k10[i] +
                            COEF::d415*k2[i] + COEF::d416*k3[i]);
                    rcont6[i] = h * (rcont6[i] + COEF::d513*k4[i] + COEF::d514*k10[i] +
                            COEF::d515*k2[i] + COEF::d516*k3[i]);
                    rcont7[i] = h * (rcont7[i] + COEF::d613*k4[i] + COEF::d614*k10[i] +
                            COEF::d615*k2[i] + COEF::d616*k3[i]);
                    rcont8[i] = h * (rcont8[i] + COEF::d713*k4[i] + COEF::d714*k10[i] +
                            COEF::d715*k2[i] + COEF::d716*k3[i]);
                }
            }

            k1[0] = k4[0];
            k1[1] = k4[1];
            k1[2] = k4[2];
            y[0] = k5[0];
            y[1] = k5[1];
            y[2] = k5[2];
            xold = x;
            x = xph;

            if (OPT.iout){
                hout = h;
                xout = x;
                obs(naccpt+1, xold, x, y, &irtrn);
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
            hnew = h / min_d (facc1, fac11/OPT.safe);
            reject = 1;
            if (naccpt >= 1)
            nrejct=nrejct + 1;
            last = 0;
        }

        h = hnew;
    }

}