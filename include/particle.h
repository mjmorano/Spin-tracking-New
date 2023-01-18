#pragma once
#include <random>
#include <math.h>
#include <algorithm>
#include <openacc.h>
#include "DOP853func.h"
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

	particle(double* y0, options OPT) :
		Lx(OPT.Lx), Ly(OPT.Ly), Lz(OPT.Lz), m(OPT.m), tc(OPT.tc), dist(OPT.dist), V_init(OPT.V), t0(OPT.t0), tf(OPT.tf), gen(rd()), diffuse(OPT.diffuse), gas_coll(OPT.gas_coll), Temp(OPT.T), gravity(OPT.gravity), pos(),
		pos_old(), v(), v_old(), S(y0), Bz(OPT.B0), B0(OPT.B0), p_interp(), v_interp(), gamma(OPT.gamma), G(), opt(OPT),
		initx(-OPT.Lx/2,OPT.Lx/2), inity(-OPT.Ly/2,OPT.Ly/2), initz(-OPT.Lz/2,OPT.Lz/2), gen_coll_time(1/OPT.tc), gen_max_boltz(0.0,sqrt(k * OPT.T / OPT.m)), gen_norm(0.0,1.0), uni_0_1(0.0,1.0), uni_0_2pi(0.0,2*M_PI)
		// integrator(std::bind(&particle::Bloch,this,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3), y0, obs, OPT)
	{	
		pos[0] = initx(gen);
		pos[1] = inity(gen);
		pos[2] = initz(gen);
		pos_old[0] = pos[0];
		pos_old[1] = pos[1];
		pos_old[2] = pos[2];
        t = t0;

		if (gas_coll == true)
			next_gas_coll_time = gen_coll_time(gen);
		else
			next_gas_coll_time = tf + 1.0;

		if (dist == 'C') {

			double vec[3];
			for (int i = 0; i < 3; i++)
				vec[i] = gen_norm(gen);
			double vec_norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
			v[0] = V_init * vec[0] / vec_norm;
			v[1] = V_init * vec[1] / vec_norm;
			v[2] = V_init * vec[2] / vec_norm;
			Vel = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		}
		else if (dist == 'M') {
			v[0] = gen_max_boltz(gen);
			v[1] = gen_max_boltz(gen);
			v[2] = gen_max_boltz(gen);
			Vel = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		}

		v_old[0] = v[0];
		v_old[1] = v[1];
		v_old[2] = v[2];

	}

	~particle() {};

	void calc_next_collision_time();

	template <typename T> double sgn(T val);

	void new_velocities();

	void move();

	void step();

	void run();

private:
	options opt;
	const bool diffuse;
	const bool gas_coll;
	const bool gravity;
	std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> initx;
    std::uniform_real_distribution<double> inity;
    std::uniform_real_distribution<double> initz;
    std::exponential_distribution<double> gen_coll_time;
    std::normal_distribution<double> gen_norm;
    std::normal_distribution<double> gen_max_boltz;
    std::uniform_real_distribution<double> uni_0_1;
    std::uniform_real_distribution<double> uni_0_2pi;
	double y = 0;
	double theta = 0;
	double phi = 0;
	double m = 1.6e-27;
	double tc = 1.0;
	double* S;
	double v[3];
	double v_old[3];
	double pos[3];
	double pos_old[3];
	double p_interp[3];
	double v_interp[3];
	double G[3];
	double V_init;
	double Vel = 0.0;
	const double Lx;
	const double Ly;
	const double Lz;
	char coll_type = 'W';
	char dist = 'C';
	const double t0 = 0.0;
	const double tf = 0.0;
	const double t_total = tf - t0;
	double t = t0;
	double t_old = t0;
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
};