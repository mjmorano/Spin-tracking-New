#include <iostream>
#include <chrono>
#include "include/particle.h"
#include "include/options.h"
#include <openacc.h>

using namespace std;
using namespace std::chrono;

ofstream outfile;

int main() {

    auto start = high_resolution_clock::now();
    double yi[3] = {1.0, 0.0, 0.0};
	options opt;
    opt.iout = 0;
    opt.diffuse = false;
    opt.gas_coll = false;
    opt.t0 = 0.0;
    opt.tf = 1.0;
    opt.h = 0.0001;
    opt.B0 = 3e-6;

    particle p(yi, opt);
    p.run();

    auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - start);
	cout << "Execution Time: " << duration.count() << " ms\n";
}