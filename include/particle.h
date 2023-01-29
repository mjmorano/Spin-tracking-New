#pragma once
#include <math.h>
#include <algorithm>
#include <openacc.h>
#include "openacc_curand.h"
#include "DOP853func.h"
#include "options.h"

class particle
{
	const double k = 1.380649e-23;

public:

	size_t n_bounce = 0;
	size_t n_coll = 0;
	size_t n_steps = 0;
	size_t n_int_steps = 0;
	bool finished = false;
	bool bad = true;
	double lastOutput = 0.0;
	unsigned int lastIndex = 0;
	outputDtype *outputArray;
	
	particle(double3 y0, unsigned int seed, options OPT, outputDtype* storage) :
		Lx(OPT.Lx), Ly(OPT.Ly), Lz(OPT.Lz), m(OPT.m), tc(OPT.tc), 
		dist(OPT.dist), V_init(OPT.V), t0(OPT.t0), tf(OPT.tf), 
		diffuse(OPT.diffuse), gas_coll(OPT.gas_coll), 
		gravity(OPT.gravity), pos(), sqrtKT_m(sqrt(k*OPT.T/OPT.m)), max_step(OPT.max_step),
		pos_old(), v(), v_old(), Bz(OPT.B0), B0(OPT.B0), p_interp(), v_interp(), 
		gamma(OPT.gamma), G(), opt(OPT)
	{	
		
		curand_init(seed, 0, 0, &state);
		S = y0;
		outputArray = storage;
		pos.x = curand_uniform_double(&state)*Lx-Lx/2.0;
		pos.y = curand_uniform_double(&state)*Ly-Ly/2.0;
		pos.z = curand_uniform_double(&state)*Lz-Lz/2.0;
		// printf("%f\n", pos.x);
		// printf("%f\n", pos.y);
		
		pos_old = pos;
		t = t0;

		if (gas_coll == true)
			next_gas_coll_time = - tc * log(1.0 - curand_uniform_double(&state));
		else
			next_gas_coll_time = tf + 1.0;

		if (dist == 'C') {

			double3 vec;
			vec.x = curand_normal_double(&state);
			vec.y = curand_normal_double(&state);
			vec.z = curand_normal_double(&state);

			double vec_norm = len(vec);
			v = V_init * vec/vec_norm;
		}
		else if (dist == 'M') {
			v.x = curand_normal_double(&state) * sqrtKT_m;
			v.y = curand_normal_double(&state) * sqrtKT_m;
			v.z = curand_normal_double(&state) * sqrtKT_m;
		}
		
		// printf("%f\t %f\t %f\n", v.x, v.y, v.z);
		Vel = len(v);
		v_old = v;
	}

	~particle() {
		// printf("%u\n",n_int_steps);
	};

	#pragma acc routine seq nohost
	void calc_next_collision_time();
	#pragma acc routine seq nohost
	template <typename T> double sgn(T val);
	#pragma acc routine seq nohost
	void new_velocities();
	#pragma acc routine seq nohost
	void move();
	#pragma acc routine seq nohost
	void step();
	#pragma acc routine seq nohost
	void run();

private:
	options opt;
	curandState_t state;
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
	const double Lx;
	const double Ly;
	const double Lz;
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
	double Bx = 0.0;
	double By = 0.0;
	double Bz;
	double B0;
	double gamma;
	double max_step = 100.0;
};