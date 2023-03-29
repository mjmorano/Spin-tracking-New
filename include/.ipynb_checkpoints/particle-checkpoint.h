#pragma once
#include <math.h>
#include <algorithm>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <random>
#include "integrator.h"
#include "options.h"

// typedef void (*GRAD)(const double* pos, double* G);

class particle
{
	const double k = 1.380649e-23;

public:

	size_t n_bounce = 0;
	size_t n_coll = 0;
	size_t n_steps = 0;
	bool finished = false;
	bool bad = true;
	double lastOutput = 0.0;
	unsigned int lastIndex = 0;
	outputDtype *outputArray;
	
	particle(double3 y0, options OPT, unsigned long seed, unsigned int ipart, outputDtype* storage) :
		L(OPT.L), m(OPT.m), tc(OPT.tc), 
		dist(OPT.dist), V_init(OPT.V), t0(OPT.t0), tf(OPT.tf), 
		diffuse(OPT.diffuse), gas_coll(OPT.gas_coll), 
		gravity(OPT.gravity), pos(), sqrtKT_m(sqrt(k*OPT.T/opt.m)), max_step(OPT.hmax),
		pos_old(), v(), v_old(), B0(OPT.B0), p_interp(), v_interp(), 
		gamma(OPT.gamma), G(), opt(OPT), ipart(ipart), seed(seed),
		integrationType(OPT.integratorType), h(OPT.h)
	{
		// First do the RNG initialization. This will depend heavily on the compiler being used
		#if defined(__HIPCC__)
		hiprand_init(seed, ipart, 0, rngState);
		#elif defined(__NVCC__)
		curand_init(seed, ipart, 0, rngState);
		#endif
		// printf("%u\n", &thread_data);
		S = y0;
		outputArray = storage;

		pos.x = uniform()*L.x-L.x/2.0;
		pos.y = uniform()*L.y-L.y/2.0;
		pos.z = uniform()*L.z-L.z/2.0;
		
		// printf("%f\n", pos.x);
		// printf("%f\n", pos.y);
		
		pos_old = pos;
		t = t0;
		if (gas_coll == true){
			next_gas_coll_time = exponential(tc);
		}
		else
			next_gas_coll_time = tf + 1.0;

		if (dist == 'C') {

			double3 vec;
			vec.x = normal01();
			vec.y = normal01();
			vec.z = normal01();

			double vec_norm = len(vec);
			v = V_init * vec/vec_norm;
		}
		else if (dist == 'M') {
			v.x = maxboltz(sqrtKT_m);
			v.y = maxboltz(sqrtKT_m);
			v.z = maxboltz(sqrtKT_m);
		}
		
		// printf("%f\t %f\t %f\n", v.x, v.y, v.z);
		Vel = len(v);
		v_old = v;
	}

	~particle() {};
	void calc_next_collision_time();
	template <typename T> double sgn(T val);
	void new_velocities();
	void move();
	void step();
	void run();
	double uniform();
	double normal01();
	double maxboltz(const double);
	double unif02pi();
	double exponential(const double);

private:
	options opt;
	const bool diffuse;
	const bool gas_coll;
	const bool gravity;
	double y = 0;
	double theta = 0;
	double phi = 0;
	double m;
	double tc;
	double3 S;
	double3 v;
	double3 v_old;
	double3 pos;
	double3 pos_old;
	double3 p_interp;
	double3 v_interp;
	double3 G;
	double sqrtKT_m;
	double V_init;
	double Vel = 0.0;
	const double3 L;
	char coll_type = 'W';
	char dist = 'C';
	const double t0;
	const double tf;
	double t;
	double t_old;
	double dt = 0.0;
	double next_gas_coll_time;
	double dx = 0.0;
	double dy = 0.0;
	double dz = 0.0;
	double dtx = 0.0;
	double dty = 0.0;
	double dtz = 0.0;
	double tbounce = 0.0;
	char wall_hit = 'x';
	double Temp = 4.2;
	double3 B0 = {0.0, 0.0, 0.0};
	double gamma;
	unsigned int ipart;
	unsigned int icount = 0;
	unsigned long iprn;
	double max_step = 0.001;
	int integrationType = 0;//default to DOP853
	double h;//keep track of hte step size in the particle;
	
	//all the various RNG related things
	unsigned long seed = 0;
	#if defined(__HIPCC__)
	hiprandState rngState;
	#elif defined(__NVCC__)
	curandState rngState;
	#else
	std::random_device dev;
	std::mt19937_64 gen64{dev()};
	std::normal_distribution<double> dist_normal{0.0, 1.0};
	std::uniform_real_distribution<double> dist_uniform{0.0, 1.0};
	#endif
};
