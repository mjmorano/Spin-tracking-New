#include <iostream>
#include <chrono>
#include <ctime>
#include "include/particle.h"
#include "include/options.h"
#include "include/dists.h"
#include <openacc.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[]) {

	double3 yi = {1.0, 0.0, 0.0};
	options opt;
	unsigned int timestamp = time(NULL);
	int numParticles = 10000;

	char * outputFilename = "data.bin";
	
	int numOutput = (opt.tf - opt.t0)/opt.ioutInt; 	
	size_t outputSize = numOutput * sizeof(outputDtype);
	
	outputDtype * outputArray = (outputDtype*)malloc(outputSize*numParticles);
	unsigned long* nident = (unsigned long*)malloc(8 * numParticles);
	desprng_common_t *process_data;
    desprng_individual_t *thread_data;
   	thread_data = (desprng_individual_t*)malloc(sizeof(desprng_individual_t) * numParticles);
    process_data = (desprng_common_t*)malloc(sizeof(desprng_common_t));
	initialize_common(process_data);

	#pragma acc parallel loop
	for(unsigned int n = 0; n<numParticles;n++){
		nident[n] = timestamp+n;	//this is assigning the seed to the RNG
		particle p(yi, opt, thread_data + n, process_data, n, nident, &outputArray[n*numOutput]);
		p.run();
	}
		
	FILE* f = fopen(outputFilename, "wb");
	fwrite(outputArray, outputSize * numParticles, 1, f);
	fclose(f);

	free(outputArray);
	free(nident);
	free(process_data);
	free(thread_data);
	return 0;
}

