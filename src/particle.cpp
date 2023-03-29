#include "../include/particle.h"
/*
Outputs the sign of a number.
*/

template <typename T>
double particle::sgn(T val) {
	if(val<0)
		return -1.0;
	else
		return 1.0;
	//return (T(0) < val) - (val < T(0)); //we don't want the 0 case
}

/*
Calculates the timestep based on the next wall or gas collision.
*/

void particle::calc_next_collision_time() {
	if(gravity){
		//calculate distance to collision point
		dx = sgn(v.x) * Lx / 2.0 - pos.x;
		dz = sgn(v.z) * Lz / 2.0 - pos.z;
		
		//time to wall for x and z coordinate
		dtx = dx / v.x;
		dtz = dz / v.z;
		double y2 = v.y * v.y;
		if(sgn(v.y) <= 0.0){ //if the particle has negative y velocity
				dy = pos.y + Ly*0.5;
				double sqr = sqrt(-2.0*G_CONST*dy+y2);
				double temp1 = -(sqr+v.y)/G_CONST;
				double temp2 = (sqr-v.y)/G_CONST;
				//printf("temp1 = %f, temp2 = %f\n", temp1, temp2);
				dty = std::min({std::abs(temp1), std::abs(temp2)});
		}
		else{
				double maxHeight = -0.5 * y2/G_CONST + pos.y;
				if(maxHeight < 0.5 * Ly){ //in this case it can't hit the ceiling
						dy = pos.y+Ly*0.5;
						double sqr = sqrt(-2.0*G_CONST*dy+y2);
						double temp1 = -(sqr+v.y)/G_CONST;
						double temp2 = (sqr-v.y)/G_CONST;
						dty = std::max({temp1, temp2});
				}
				else{
						dy = Ly*0.5 - pos.y; //how far to ceiling
						double sqr = sqrt(-2.0*G_CONST*dy+y2);
						double temp1 = -(sqr+v.y)/G_CONST;
						double temp2 = (sqr-v.y)/G_CONST;
						dty = std::min({std::abs(temp1), std::abs(temp2)});
				}
		}
		if (dtx < 1e-16 || std::isnan(dtx))
			dtx = 1e6;
		else if (dty < 1e-16 || std::isnan(dty))
			dty = 1e6;
		else if (dtz < 1e-16 || std::isnan(dtz))
			dtz = 1e6;
		//printf("pos = %f, %f, %f, vel = %f %f %f, dt = %f %f %f\n", pos.x, pos.y, pos.z, v.x, v.y, v.z, dtx, dty, dtz);
	}
	else{	
		
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
	}
	int min_elm;
	if(dtx <= dty && dtx <= dtz){
		tbounce = dtx;
		min_elm = 0;
	}
	else if(dty <= dtx && dty <= dtz){
		tbounce = dty;
		min_elm = 1;
	}
	else if(dtz <= dtx && dtz <= dty){
		tbounce = dtz;
		min_elm = 2;
	}
	double timeToNextGas = next_gas_coll_time - t;
	//printf("t=%f v: %f %f %f, p = %f %f %f, tBounce = %f, minT = %f\n", t, v.x, v.y, v.z, pos.x, pos.y, pos.z, tbounce, max_step);
	if(max_step <= timeToNextGas && max_step <= tbounce && t + max_step < tf){ //check if the max step size is smaller than the next collision times
		//if so then just say we don't collide and keep going
		dt = max_step;
		coll_type = 'N';
	}
	else if(tbounce < timeToNextGas && t + tbounce < tf){ //is a wall bounce next
		dt = tbounce;
		n_bounce += 1;
		coll_type = 'W';
		if (min_elm == 0)
			wall_hit = 'x';
		else if (min_elm == 1)
			wall_hit = 'y';
		else if (min_elm == 2)
			wall_hit = 'z';
	}
	else if (t + tbounce > next_gas_coll_time && next_gas_coll_time < tf) { //is a gas collision next?
		dt = next_gas_coll_time - t;
		next_gas_coll_time += exponential(get_uniform_prn(process_data, thread_data, ++icount, &iprn), tc);
		coll_type = 'G';
		n_coll += 1;
	}
	else {
		//printf("huh?");
		coll_type = 'N';
		dt = tf - t;
		finished = true;
	}
}

/*
Calculates the new velocities after a wall or gas collision.
*/

void particle::new_velocities() {
	v_old = v;
	Vel = len(v);
	//printf("%c\n", coll_type);
	if (coll_type == 'N'){
		//in this case we don't have a wall collision and it's just iterating through space still
		//don't update the velocities they're fine
	}
	else if (coll_type == 'W' && diffuse == false) {
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
		v.x = maxboltz(get_uniform_prn(process_data, thread_data, ++icount, &iprn), sqrtKT_m);
		v.y = maxboltz(get_uniform_prn(process_data, thread_data, ++icount, &iprn), sqrtKT_m);
		v.z = maxboltz(get_uniform_prn(process_data, thread_data, ++icount, &iprn), sqrtKT_m);
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
	t_old = t; // update the time
	t += dt; //increment forward
	pos_old = pos; //update old position
	v_old = v; //update old velocity
	double3 a;
	if(gravity)
		a = (double3){0.0, G_CONST, 0.0}; //acceleration due to gravity
	else
		a = {0.0, 0.0, 0.0};
	pos = pos_old +  v * dt + 0.5 *a*dt*dt; //update position
	v = v + a * dt; //update velocity
	
}

/*
Performs one particle and spin integration step.
*/

void particle::step() {
	//printf("t = %f dt = %f, v = %f %f %f, pos = %f %f %f\n", t, dt, v.x, v.y, v.z, pos.x, pos.y, pos.z);
	calc_next_collision_time(); //when do we hit something next?
	move(); //move the particle forward that amount of time
	//printf("%f %f %f %f\n", t, pos.x, pos.y, pos.z);
	//printf("t = %f dt = %f, v = %f %f %f, pos = %f %f %f\n", t, dt, v.x, v.y, v.z, pos.x, pos.y, pos.z);
	new_velocities(); //update the velocity based on the collision type
	integrate(t_old, t, S, pos_old, pos, v_old, v, opt, lastOutput, lastIndex, max_index, outputArray); //integrate the spin along the path
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