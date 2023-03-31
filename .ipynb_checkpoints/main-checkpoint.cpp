#include <iostream>
#include <chrono>
#include <ctime>
#include "include/particle.h"
#include "include/options.h"
#include "include/optionsParser.h"

#if defined(_OPENMP)
#include <omp.h>
#endif
#if defined(__HIPCC__)
#include <hip/hip_runtime.h>
#include <hiprand/hiprand.h>
#include <hiprand/hiprand_kernel.h>
#endif

using namespace std;
using namespace std::chrono;


void mainAnalysis(options opt, int totalTime, char* outputName, unsigned int seed, double3 yi);
options parseUserInput(int argc, char * argv[]);

int main(int argc, char* argv[]){
	options opt = parseUserInput(argc, argv);
	int totalTime = 3600; //total time allowed in seconds
	char * outputName = "../data"; //the start of the output name, will automatically end in .bin and follow outputName01.bin, outputName02.bin, etc
	unsigned int seed = 0;
	double3 yi = {1.0, 0.0, 0.0};
	
	auto start = high_resolution_clock::now();
	mainAnalysis(opt, totalTime, outputName, seed, yi);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop-start).count();
	std::cout<<duration<<std::endl;
	return 0;
}

//this functions does the actual analysis and integration
void mainAnalysis(options opt, int totalTime, char* outputName, unsigned int seed, double3 yi){
	#if defined(__NVCC__) || defined(__HIPCC__)
	{
		//In this case we're going to use the GPU to do mostly everything
		//start up the same way basically
		//start the clock on the process
		//Do the single CPU version of the code that uses all cores/threads on a singular CPU
		unsigned int timestamp = time(NULL);
		int numParticles = opt.numParticles;
		int numOutput = (opt.tf - opt.t0)/opt.ioutInt;
		size_t outputSize = numOutput * sizeof(outputDtype);
		
		outputDtype * outputArrayHost;
		outputDtype * outputArrayGPU;
		unsigned long* nident;
		
		//do the CPU allocation first
		outputArrayHost  = (outputDtype*)malloc(outputSize * numParticles);
		#if defined(__NVCC__)
		cudaMalloc(&outputArrayGPU, outputSize*numParticles); //output location
		curandState *rngStates;
		cudaMalloc(&rngStates, sizeof(curandState)*numParticles); //allocate one state per particle
		cudaDeviceSynchronize();
		#else
		hipMalloc(&outputArrayGPU, outputSize*numParticles);
		hiprandStateXORWOW_t *rngStates;
		hipMalloc(&rngStates, sizeof(hiprandStateXORWOW_t)*numParticles);
		hipDeviceSynchronize();
		#endif
		
		//now actually do the kernel call
		int numPartsPerBlock = opt.numPerGPUBlock;
		int numBlocks = std::ceil((double)opt.numParticles/(double)numPartsPerBlock);
		runSimulation<<<numBlocks, numPartsPerBlock>>>(numParticles, outputArrayGPU, rngStates, opt, seed, yi, numOutput);
		
		#if defined(__NVCC__)
		cudaDeviceSynchronize();
		cudaMemcpy(outputArrayHost, outputArrayGPU, outputSize * numParticles, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		cudaFree(outputArrayGPU);
		#else
		hipDeviceSynchronize();
		hipMemcpy(outputArrayHost, outputArrayGPU, outputSize * numParticles, hipMemcpyDeviceToHost);
		hipDeviceSynchronize();
		hipFree(outputArrayGPU);
		#endif
		
		FILE* f = fopen(outputName, "wb");
		fwrite(outputArrayHost, outputSize * numParticles, 1, f);
		fclose(f);
		
		free(outputArrayHost);
	}
	#else
	{
		//Do the single CPU version of the code that uses all cores/threads on a singular CPU
		unsigned int timestamp = time(NULL);
		int numParticles = opt.numParticles;
		int numOutput = (opt.tf - opt.t0)/opt.ioutInt;
		size_t outputSize = numOutput * sizeof(outputDtype);
		outputDtype * outputArray = (outputDtype*)malloc(outputSize*numParticles); //allocate the total output size
		
		runSimulation(numParticles, outputArray, opt, seed, yi, numOutput);
		
		//printf("%d\n", (int)duration);
		
		FILE* f = fopen(outputName, "wb");
		fwrite(outputArray, outputSize * numParticles, 1, f);
		fclose(f);
		
		free(outputArray);
	}
	#endif
	return;
}

options parseUserInput(int argc, char *argv[]){
	//assume the first input is the file name
	if(argc < 2){
		std::cout<<"invalid input options"<<std::endl;
		exit(-1);
	}
	char *filename = argv[1];
	options opt = optionParser(filename);
	return opt;
	
}