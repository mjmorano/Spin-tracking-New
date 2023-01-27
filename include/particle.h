#pragma once
#include <math.h>
#include <algorithm>
#include <openacc.h>
#include "DOP853func.h"
#include "options.h"
#include "dists.h"

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
	
	particle(double3 y0, options OPT, desprng_individual_t* thread_data, desprng_common_t* process_data, unsigned int ipart, unsigned long* nident, outputDtype* storage) :
		Lx(OPT.Lx), Ly(OPT.Ly), Lz(OPT.Lz), m(OPT.m), tc(OPT.tc), 
		dist(OPT.dist), V_init(OPT.V), t0(OPT.t0), tf(OPT.tf), 
		diffuse(OPT.diffuse), gas_coll(OPT.gas_coll), 
		gravity(OPT.gravity), pos(), sqrtKT_m(sqrt(k*OPT.T/opt.m)), max_step(OPT.max_step),
		pos_old(), v(), v_old(), Bz(OPT.B0), B0(OPT.B0), p_interp(), v_interp(), 
		gamma(OPT.gamma), G(), opt(OPT), ipart(ipart), thread_data(thread_data), process_data(process_data) 
	{	

		create_identifier(nident+ipart);
		initialize_individual(process_data, thread_data, nident[ipart]);
		// printf("%u\n", &thread_data);
		S = y0;
		outputArray = storage;
		pos.x = get_uniform_prn(process_data, thread_data, ++icount, &iprn)*Lx-Lx/2.0;
		pos.y = get_uniform_prn(process_data, thread_data, ++icount, &iprn)*Ly-Ly/2.0;
		pos.z = get_uniform_prn(process_data, thread_data, ++icount, &iprn)*Lz-Lz/2.0;
		// printf("%f\n", pos.x);
		// printf("%f\n", pos.y);
		
		pos_old = pos;
		t = t0;

		if (gas_coll == true)
			next_gas_coll_time = exponential(get_uniform_prn(process_data, thread_data, ++icount, &iprn),tc);
		else
			next_gas_coll_time = tf + 1.0;

		if (dist == 'C') {

			double3 vec;
			vec.x = normal01(get_uniform_prn(process_data, thread_data, ++icount, &iprn));
			vec.y = normal01(get_uniform_prn(process_data, thread_data, ++icount, &iprn));
			vec.z = normal01(get_uniform_prn(process_data, thread_data, ++icount, &iprn));

			double vec_norm = len(vec);
			v = V_init * vec/vec_norm;
		}
		else if (dist == 'M') {
			v.x = maxboltz(get_uniform_prn(process_data, thread_data, ++icount, &iprn), sqrtKT_m);
			v.y = maxboltz(get_uniform_prn(process_data, thread_data, ++icount, &iprn), sqrtKT_m);
			v.z = maxboltz(get_uniform_prn(process_data, thread_data, ++icount, &iprn), sqrtKT_m);
		}
		
		// printf("%f\t %f\t %f\n", v.x, v.y, v.z);
		Vel = len(v);
		v_old = v;
	}

	~particle() {};
	#pragma acc routine seq
	void calc_next_collision_time();
	#pragma acc routine seq
	template <typename T> double sgn(T val);
	#pragma acc routine seq
	void new_velocities();
	#pragma acc routine seq
	void move();
	#pragma acc routine seq
	void step();
	#pragma acc routine seq
	void run();

private:
	desprng_individual_t* thread_data; 
	desprng_common_t* process_data;
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
	unsigned int ipart;
	unsigned int icount = 0;
	unsigned long iprn;
	double max_step = 0.001;
};
