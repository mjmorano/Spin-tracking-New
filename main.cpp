#include <iostream>
#include <chrono>
#include "include/particle.h"
#include "include/options.h"
#include <openacc.h>

using namespace std;
using namespace std::chrono;

int main() {

    auto start = high_resolution_clock::now();
    double3 yi = {1.0, 0.0, 0.0};
    options opt;
    opt.iout = 2;
    opt.diffuse = false;
    opt.gas_coll = false;
    opt.t0 = 0.0;
    opt.tf = 10.0;
    opt.ioutInt = 0.00001;
    opt.h = 0.0001;
    opt.B0 = 3e-6;
	
	int numOutput = (opt.tf - opt.t0)/opt.ioutInt;
	printf("numOutput = %d\n", numOutput);
	
	size_t outputSize = numOutput * sizeof(float)*4; //4 because timestamp + spin vector
	int numParticles = 1000;
	
	#pragma acc parallel loop
	for(int n = 0; n<numParticles;n++){
			
		float * outputArray = (float*)malloc(outputSize);
		
		particle p(yi, opt, outputArray);
		p.run();
		
		//#pragma serial
		int loc = 0;
		for(int i = 0; i < numOutput; i++){
			loc = i*4;
			//printf("%lf %lf %lf %lf\n", outputArray[loc], outputArray[loc+1], 
				//outputArray[loc+2], outputArray[loc+3]);
		}
	}
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - start);
		//cout << "Execution Time: " << duration.count() << " ms\n";
}
