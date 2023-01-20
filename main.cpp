#include <iostream>
#include <chrono>
#include "include/particle.h"
#include "include/options.h"
#include "include/dists.h"
#include <openacc.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[]) {
    	auto start = high_resolution_clock::now();
    	double3 yi = {1.0, 0.0, 0.0};
    	options opt;
    	opt.iout = 2;
    	opt.diffuse = false;
    	opt.gas_coll = false;
    	opt.t0 = 0.0;
    	opt.tf = 1.0;
    	opt.ioutInt = 0.01;
    	opt.h = 0.0001;
    	opt.B0 = 3e-6;
	int numParticles = 10000;

	char * outputFilename = "data.bin";
	
	int numOutput = (opt.tf - opt.t0)/opt.ioutInt*4; // 4 because timestamp + spin vector	
	size_t outputSize = numOutput * sizeof(float); //4 because timestamp + spin vector
	float * outputArray = (float*)malloc(outputSize*numParticles);
	std::cout<<numOutput<<" : "<<outputSize<<" : "<<outputSize*numParticles<<std::endl;

	unsigned long* nident = (unsigned long*)malloc(8 * numParticles);
	desprng_common_t *process_data;
    	desprng_individual_t *thread_data;
   	thread_data = (desprng_individual_t*)malloc(sizeof(desprng_individual_t) * numParticles);
    	process_data = (desprng_common_t*)malloc(sizeof(desprng_common_t));
	initialize_common(process_data);

	#pragma acc parallel loop
	for(int n = 0; n<numParticles;n++){
		nident[n] = n;	
		particle p(yi, opt, thread_data, process_data, n, &nident[n], &outputArray[n*numOutput]);
		p.run();
		//printf("reached on %d\n", n);
	}
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - start);
	cout << "Execution Time: " << duration.count() << " ms\n";
	cout << numOutput * 4 <<std::endl;
	FILE* f = fopen(outputFilename, "wb");
	fwrite(outputArray, outputSize * numParticles, 1, f);
	fclose(f);	

	free(outputArray);
	free(nident);
	free(process_data);
	free(thread_data);
	return 0;
}

