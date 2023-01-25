#include <iostream>
#include <chrono>
#include <time.h>
#include "include/particle.h"
#include "include/options.h"
#include "include/dists.h"
#include <openacc.h>
#include <mpi.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[]) {
	auto start = high_resolution_clock::now();
	
	MPI_Init(&argc, &argv);
	int num_devices, gpuid, local_rank, global_rank;
	MPI_Comm shmcomm;
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
	MPI_Comm_rank(shmcomm, &local_rank);
	MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
	acc_device_t device_type;
	device_type = acc_get_device_type();
	num_devices = acc_get_num_devices(device_type);
	gpuid = local_rank % num_devices;
	acc_set_device_num(gpuid, device_type);

	
	double3 yi = {1.0, 0.0, 0.0};
	options opt;
	int numParticles = 80 * 64 * 2;

	char outputFilename[20];
	sprintf(outputFilename, "data%d.bin", global_rank);
	
	int numOutput = (opt.tf - opt.t0)/opt.ioutInt; 	
	size_t outputSize = numOutput * sizeof(outputDtype);
	
	outputDtype * outputArray = (outputDtype*)malloc(outputSize*numParticles);
	unsigned long* nident = (unsigned long*)malloc(8 * numParticles);
	desprng_common_t *process_data;
	desprng_individual_t *thread_data;
   	thread_data = (desprng_individual_t*)malloc(sizeof(desprng_individual_t) * numParticles);
	process_data = (desprng_common_t*)malloc(sizeof(desprng_common_t));
	initialize_common(process_data);
	
	#pragma acc parallel loop create(thread_data[:numParticles], nident[:numParticles]) copyin(process_data[:1]), copyout(outputArray[:numParticles*numOutput])
	for(int n = 0; n<numParticles;n++){
		nident[n] = (unsigned long)n;	//this is assigning the seed to the RNG
		particle p(yi, opt, thread_data+n, process_data, n, nident, &outputArray[n*numOutput]);
		p.run();		
	}
		
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - start);
	cout << "Execution Time: " << duration.count() << " ms\n";
	FILE* f = fopen(outputFilename, "wb");
	fwrite(outputArray, outputSize * numParticles, 1, f);
	fclose(f);	

	free(outputArray);
	free(nident);
	free(process_data);
	free(thread_data);
	return 0;
}

