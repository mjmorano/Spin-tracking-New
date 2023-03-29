#include <iostream>
#include <chrono>
#include <ctime>
#include "include/particle.h"
#include "include/options.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __HIPCC__
#include <hip/hip_runtime.h>
#include <hiprand/hiprand.h>
#include <hiprand/hiprand_kernel.h>
#endif

using namespace std;
using namespace std::chrono;


void mainAnalysis(options opt, int numParticles, int totalTime, char* outputName);

int main(int argc, char* argv[]) {

	options opt;
	opt.iout = 2;
	opt.gas_coll = false;
	opt.rtol = 1.0e-14;
	opt.t0 = 0.0;
	opt.tf = 100.0;
	opt.ioutInt = 0.01;
	opt.gravity=true;
	opt.B0 = {0, 0, 3.0e-6};
	opt.E = {0, 0, 75.0e5};
	opt.integratorType = 0;
	opt.h = opt.ioutInt/2.0; //default to half of the output interval spacing to address a bug
	
	//set either of these to -1 to ignore that option
	//if both are -1, then code breaks
	//if both are NOT -1, then it will stop when first limit is reached
	int numParticles = 100; //desired number of particles
	int totalTime = 3600; //total time allowed in seconds
	char * outputName = "../data"; //the start of the output name, will automatically end in .bin and follow outputName01.bin, outputName02.bin, etc
	
	mainAnalysis(opt, numParticles, totalTime, outputName);
	return 0;
}

//this functions does the actual analysis and integration
void mainAnalysis(options opt, int numParticles, int totalTime, char* outputName){
	#if defined __NVCC__ || defined __HIPCC__
	{
		//In this case we're going to use the GPU to do mostly everything
		printf("using the GPU for things\n");
		
		//start up the same way basically
		auto start = high_resolution_clock::now();//start the clock on the process
		//Do the single CPU version of the code that uses all cores/threads on a singular CPU
		unsigned int timestamp = time(NULL);

		int numOutput = (opt.tf - opt.t0)/opt.ioutInt;
		size_t outputSize = numOutput * sizeof(outputDtype);
		double3 yi = {1.0, 0.0, 0.0};
		
		outputDtype * outputArrayHost;
		outputDtype * outputArrayGPU;
		unsigned long* nident;
		desprng_common_t *process_data;
		desprng_individual_t *thread_data;
		
		//do the CPU allocation first
		outputArrayHost  = (outputDtype*)malloc(outputSize * numParticles);
		#if __NVCC__
		cudaMalloc(&outputArrayGPU, outputSize*numParticles); //output location
		curandState *rngStates;
		cudaMalloc(&rngStates, sizeof(curandState)*numParticles); //allocate one state per particle
		#else
		hipMalloc(&outputArrayGPU, outputSize*numParticles);
		hiprandStateXORWOW_t *rngStates;
		hipMalloc(&rngStates, sizeof(hiprandStateXORWOW_t)*numParticles);
		#endif
		
		//first we need to initialize the common data on the device
		
		
		
		#if __NVCC__
		cudaFree(outputArrayGPU);
		cudaFree(nident);
		cudaFree(process_data);
		cudaFree(thread_data);
		#else
		hipFree(outputArrayGPU);
		hipFree(nident);
		hipFree(process_data);
		hipFree(thread_data);
		#endif
		
		free(outputArrayHost);
		printf("finishing with the GPU\n");
	}
	#else
	{
		printf("using the CPU for everything\n");
		auto start = high_resolution_clock::now();//start the clock on the process
		//Do the single CPU version of the code that uses all cores/threads on a singular CPU
		unsigned int timestamp = time(NULL);

		int numOutput = (opt.tf - opt.t0)/opt.ioutInt;
		size_t outputSize = numOutput * sizeof(outputDtype);
		double3 yi = {1.0, 0.0, 0.0};
		outputDtype * outputArray = (outputDtype*)malloc(outputSize*numParticles); //allocate the total output size
		unsigned long seed = 0;
		#if defined(_OPENMP)
		#pragma omp parallel for
		#endif
		for(unsigned int n = 0; n<numParticles;n++){
			//printf("%d\n", omp_get_num_threads());
			particle p(yi, opt, seed, n, &outputArray[n*numOutput]);
			p.run();
		}
		auto end = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(end-start);
		//printf("%d\n", (int)duration);

		FILE* f = fopen(outputName, "wb");
		fwrite(outputArray, outputSize * numParticles, 1, f);
		fclose(f);

		free(outputArray);
	}
	#endif
	return;
}