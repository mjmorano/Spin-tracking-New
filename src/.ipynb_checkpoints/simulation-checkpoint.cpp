#include "../include/simulation.h"

//this functions does the actual analysis and integration
void mainAnalysis(options opt, int totalTime, char* outputName, unsigned int seed, double3 yi){
	#if defined(__NVCOMPILER) || defined(__HIPCC__)
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
		#if defined(__NVCOMPILER)
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
		
		#if defined(__NVCOMPILER)
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