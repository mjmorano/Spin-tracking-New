#include <iostream>
#include <math.h>
#include <fstream>
#include <chrono>
#include <vector>
#include <iomanip>
#include <openacc.h>
#include "include/test.cuh"
#include "include/DOP853.h"
#include "include/particle.h"
#include "include/options.h"

std::ofstream outfile;

using namespace std;
using namespace std::chrono;

double pulse(double t){
    // return 64.7775066e-6*cos(10000*t);
    return 0.0;
}

void grad(const vector<double> pos, vector<double>&G){
    G[0] = 0.0;
    G[1] = 0.0;
    G[2] = 0.0;
}

void solout(long nr, double xold, double x, std::vector<double>& y, unsigned int n, int* irtrn){
    outfile << setprecision(16);
    outfile << x << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\n";
}

int main() {

    outfile.open("/mnt/c/Users/moran/Desktop/dressing.txt");
    std::vector<double> yi = {1.0, 0.0, 0.0};
	options opt;
    opt.iout = 1;
    opt.gas_coll = false;
    opt.rtol = 1e-14;
    opt.t0 = 0.0;
    opt.tf = 1.0;
    opt.B0 = 3e-6;
    particle p(pulse, grad, solout, yi, opt);
    p.run();
    
    outfile.close();
}