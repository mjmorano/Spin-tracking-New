#include "../include/particle.h"
/*
Outputs the sign of a number.
*/

template <typename T>
double particle::sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

/*
Calculates the timestep based on the next wall or gas collision.
*/

void particle::calc_next_collision_time() {
	dx = sgn(v.x) * Lx / 2.0 - pos.x;
	dy = sgn(v.y) * Ly / 2.0 - pos.y;
	dz = sgn(v.z) * Lz / 2.0 - pos.z;
	dtx = dx / v.x;
	dty = dy / v.y;
	dtz = dz / v.z;

	if (dtx < 1e-16)
		dtx = 1e6;
	else if (dty < 1e-16)
		dty = 1e6;
	else if (dtz < 1e-16)
		dtz = 1e6;

	tbounce = std::min( {dtx, dty, dtz} );

	if (t + tbounce < next_gas_coll_time && t + tbounce < tf) {
		dt = tbounce;
		n_bounce += 1;
		coll_type = 'W';
		if (tbounce == dtx)
			wall_hit = 'x';
		else if (tbounce == dty)
			wall_hit = 'y';
		else if (tbounce == dtz)
			wall_hit = 'z';
	}
	else if (t + tbounce > next_gas_coll_time && next_gas_coll_time < tf) {
		dt = next_gas_coll_time - t;
		next_gas_coll_time += exponential(get_uniform_prn(process_data, thread_data, ++icount, &iprn),tc);
		coll_type = 'G';
		n_coll += 1;
	}
	else {
		dt = tf - t;
		finished = true;
	}
}

/*
Calculates the new velocities after a wall or gas collision.
*/

void particle::new_velocities() {
	v_old = v;
	if (coll_type == 'W' && diffuse == false) {
		if (wall_hit == 'x')
			v.x *= -1.0;
		else if (wall_hit == 'y')
			v.y *= -1.0;
		else if (wall_hit == 'z')
			v.z *= -1.0;
	}
	else if (coll_type == 'W' && diffuse == true) {
		//V = sqrt(vx * vx + vy * vy + vz * vz);
		phi = acos(cbrt(1 - get_uniform_prn(process_data, thread_data, ++icount, &iprn)));
		theta = 2 * M_PI * get_uniform_prn(process_data, thread_data, ++icount, &iprn);
		if (wall_hit == 'x') {
			v.x = -1 * sgn(v.x) * Vel * cos(phi);
			v.y = -Vel * sin(phi) * cos(theta);
			v.z = Vel * sin(phi) * sin(theta);
		}
		else if (wall_hit == 'y') {
			v.x = Vel * sin(phi) * cos(theta);
			v.y = -1 * sgn(v.y) * Vel * cos(phi);
			v.z = Vel * sin(phi) * sin(theta);
		}
		else if (wall_hit == 'z') {
			v.x = Vel * sin(phi) * cos(theta);
			v.y = Vel * sin(phi) * sin(theta);
			v.z = -1 * sgn(v.z) * Vel * cos(phi);
		}
	}
	else if (coll_type == 'G' && dist == 'M') {
		v.x = maxboltz(get_uniform_prn(process_data, thread_data, ++icount, &iprn),temp,m);
		v.y = maxboltz(get_uniform_prn(process_data, thread_data, ++icount, &iprn),temp,m);
		v.z = maxboltz(get_uniform_prn(process_data, thread_data, ++icount, &iprn),temp,m);

		Vel = len(v);
	}
	else if (coll_type == 'G' && dist == 'C') {
		double3 vec;
		vec.x = normal01(get_uniform_prn(process_data, thread_data, ++icount, &iprn));
		vec.y = normal01(get_uniform_prn(process_data, thread_data, ++icount, &iprn));
		vec.z = normal01(get_uniform_prn(process_data, thread_data, ++icount, &iprn));
		double vec_norm = len(vec);
		v = Vel * vec/vec_norm;
		Vel = len(v);
	}
}

/*
Moves the particle based on the particle velocity and calcuated timestep.
*/

void particle::move() {
	// printf("%f\t %f\t %f\t %f\n",t,pos[0],pos[1],pos[2]);
	t_old = t;
	t += dt;
	pos_old = pos;
	pos = pos +  v * dt;
}

/*
Performs one particle and spin integration step.
*/

void particle::step() {
	calc_next_collision_time();
	move();
	new_velocities();
	integrate(t_old, t, S, pos_old, pos, v_old, v, opt, lastOutput, lastIndex, outputArray);
	// integrate_step();

	n_steps += 1;
}

/*
Convenience function that calls the step() function repeatedly until the end time is reached.
*/

void particle::run() {
	while (finished != true) {
		step();
	}
	// printf("%f\t %f\t %f\n", S[0], S[1], S[2]);
}

/*
void particle::interpolate(const double ti, double* p_out, double* v_out){
	p_out[0] = (pos_old[0]*(t - ti) + pos[0]*(ti - t_old))/(t-t_old);
	p_out[1] = (pos_old[1]*(t - ti) + pos[1]*(ti - t_old))/(t-t_old);
	p_out[2] = (pos_old[2]*(t - ti) + pos[2]*(ti - t_old))/(t-t_old);
	// printf("%f\t %f\t %f\n", pos_old[0], pos[0], p_out[0]);
	v_out[0] = v_old[0];
	v_out[1] = v_old[1];
	v_out[2] = v_old[2];
}
void particle::integrate_step(){
	integrator.integrate(t_old,t);
}
void particle::integrate_spin(){
	integrator.integrate(t0,tf);
}
void particle::Bloch(const double x, const double* y, double* f){
	interpolate(x,p_interp,v_interp);
	Bx = pulse(x);
	grad(p_interp,G);
	Bx += G[0];
	By = G[1];
	Bz = B0 + G[2];
	// printf("%f\t %f\t %f\n",Bx,By,Bz);
    f[0] = gamma * (y[1]*Bz - y[2]*By);
	f[1] = gamma * (y[2]*Bx - y[0]*Bz);
	f[2] = gamma * (y[0]*By - y[1]*Bx);
}

*/
