#include <iostream>
#include <chrono>
#include <ctime>
#include "include/particle.h"
#include "include/options.h"
#include "include/optionsParser.h"
#include "include/simulation.h"

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