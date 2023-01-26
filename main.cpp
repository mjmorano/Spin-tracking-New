#include <ctime>
#include "include/particle.h"
#include "include/options.h"
#include "include/dists.h"
#include <openacc.h>

using namespace std;

int main(int argc, char* argv[]) {

	double3 yi = {1.0, 0.0, 0.0};
	options opt;
	unsigned int timestamp = time(NULL);
	int numParticles = 1;

	char * outputFilename = "data.bin";
	
	int numOutput = (opt.tf - opt.t0)/opt.ioutInt; 	
	size_t outputSize = numOutput * sizeof(outputDtype);
	
	outputDtype * outputArray = (outputDtype*)malloc(outputSize*numParticles);

	#pragma acc parallel loop
	for(unsigned int n = 0; n<numParticles;n++){
		particle p(yi, n+timestamp, opt, &outputArray[n*numOutput]);
		p.run();
	}
		
	FILE* f = fopen(outputFilename, "wb");
	fwrite(outputArray, outputSize * numParticles, 1, f);
	fclose(f);

	free(outputArray);
	return 0;
}