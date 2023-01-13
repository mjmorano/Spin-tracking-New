#pragma once
#include <vector>
#include <string>
#include <random>
#include <math.h>
#include <algorithm>
#include "DOP853.h"
#include "options.h"

typedef void (*GRAD)(const std::vector<double> pos, std::vector<double>& G);
typedef void (*OBS)(long nr, double xold, double x, std::vector<double>& y, unsigned int n, int* irtrn);
typedef std::function<void(double x, std::vector<double> y, std::vector<double>&f)> FCN;

class particle
{
	const double k = 1.380649e-23;

public:

	size_t n_bounce = 0;
	size_t n_coll = 0;
	size_t n_steps = 0;
	bool finished = false;
	bool bad = true;

	particle(std::function<double(double x)> pulse, GRAD grad, OBS obs, std::vector<double> y0, options OPT) :
		Lx(OPT.Lx), Ly(OPT.Ly), Lz(OPT.Lz), m(OPT.m), tc(OPT.tc), dist(OPT.dist), V_init(OPT.V), t0(OPT.t0), tf(OPT.tf), gen(rd()), diffuse(OPT.diffuse), gas_coll(OPT.gas_coll), Temp(OPT.T), gravity(OPT.gravity), pos(3),
		pos_old(3), v(3), v_old(3), S(y0), Bz(OPT.B0), p_interp(3), v_interp(3), gamma(OPT.gamma), pulse(pulse), grad(grad), G(3),
		initx(-OPT.Lx/2,OPT.Lx/2), inity(-OPT.Ly/2,OPT.Ly/2), initz(-OPT.Lz/2,OPT.Lz/2), gen_coll_time(1/OPT.tc), gen_max_boltz(0.0,sqrt(k * OPT.T / OPT.m)), gen_norm(0.0,1.0), uni_0_1(0.0,1.0), uni_0_2pi(0.0,2*M_PI),
		integrator(std::bind(&particle::Bloch,this,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3), y0, obs, OPT)
	{	
		pos = {initx(gen), inity(gen), initz(gen)};
		pos_old = pos;
        t = t0;

		if (gas_coll == true)
			next_gas_coll_time = gen_coll_time(gen);
		else
			next_gas_coll_time = tf + 1.0;

		if (dist == 'C') {

			std::vector<double> vec(3);
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

		v_old = v;

	}

	~particle() {};

	void calc_next_collision_time();

	template <typename T> double sgn(T val);

	void new_velocities();

	void move();

	void step();

	void run();

	void interpolate(const double ti, std::vector<double>& p_out, std::vector<double>& v_out);

	void integrate_step();

	void integrate_spin();

	void Bloch(const double t, const std::vector<double>& y, std::vector<double>& f);

private:

	const bool diffuse;
	const bool gas_coll;
	const bool gravity;
	DOP853 integrator;
	std::function<double(double x)> pulse;
	GRAD grad;
	std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> initx;
    std::uniform_real_distribution<double> inity;
    std::uniform_real_distribution<double> initz;
    std::exponential_distribution<double> gen_coll_time;
    std::uniform_real_distribution<double> gen_norm;
    std::uniform_real_distribution<double> gen_max_boltz;
    std::uniform_real_distribution<double> uni_0_1;
    std::uniform_real_distribution<double> uni_0_2pi;
	double y = 0;
	double theta = 0;
	double phi = 0;
	double m = 1.6e-27;
	double tc = 1.0;
	std::vector<double> S;
	std::vector<double> v;
	std::vector<double> v_old;
	std::vector<double> pos;
	std::vector<double> pos_old;
	std::vector<double> p_interp;
	std::vector<double> v_interp;
	std::vector<double> G;
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
	double gamma;
};