#include "../include/particle.h"
/*
Outputs the sign of a number.
*/

#if defined(__HIPCC__)
#define __PREPROC__ __device__
#elif defined(__NVCOMPILER)
#define __PREPROC__ __device__
#else
#include <random>
#define __PREPROC__
#endif


#if defined(__HIPCC__)
__global__ void runSimulation(int numParticles, outputDtype* outputArray, hiprandStateXORWOW_t* states, options OPT, unsigned long seed, double3 yi, int numOutput){
	unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if(tid < numParticles){
		particle p(yi, OPT, seed, tid, &outputArray[tid*numOutput]);
		p.run();
	}
}
#elif defined(__NVCC__)
__global__ void runSimulation(int numParticles, outputDtype* outputArray, curandStateXORWOW_t* states, options OPT, unsigned long seed, double3 yi, int numOutput){
	unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if(tid < numParticles){
		particle p(yi, OPT, seed, tid, &outputArray[tid*numOutput]);
		p.run();
	}
}
#else
void runSimulation(int numParticles, outputDtype* outputArray, options OPT, unsigned long seed, double3 yi, int numOutput){
	#if defined(_OPENMP)
	#pragma omp parallel for
	#endif
	for(unsigned int n = 0; n<numParticles;n++){
		//printf("%d\n", omp_get_num_threads());
		particle p(yi, OPT, seed, n, &outputArray[n*numOutput]);
		p.run();
	}
}
#endif

template <typename T>
__PREPROC__ double particle::sgn(T val) {
	if(val<0)
		return -1.0;
	else
		return 1.0;
	//return (T(0) < val) - (val < T(0)); //we don't want the 0 case
}

/*
Calculates the timestep based on the next wall or gas collision.
*/

#if defined(__HIPCC__)
__PREPROC__ double particle::uniform(){
	return hiprand_uniform_double(&rngState);
}

__PREPROC__ double particle::normal01(){
	return hiprand_normal_double(&rngState);
}

__PREPROC__ double particle::maxboltz(const double sqrtkT_m){
    return normal01() * sqrtkT_m;
}

__PREPROC__ double particle::unif02pi(){
    return hiprand_uniform_double(&rngState) * 2.0 * M_PI;
}

__PREPROC__ double particle::exponential(const double tc){
    return - tc * log(1.0 - hiprand_uniform_double(&rngState));
}
#elif defined(__NVCC__)
__PREPROC__ double particle::uniform(){
	return curand_uniform_double(&rngState);
}

__PREPROC__ double particle::normal01(){
	return curand_normal_double(&rngState);
}

__PREPROC__ double particle::maxboltz(const double sqrtkT_m){
    return normal01() * sqrtkT_m;
}

__PREPROC__ double particle::unif02pi(){
    return curand_uniform_double(&rngState) * 2.0 * M_PI;
}

__PREPROC__ double particle::exponential(const double tc){
    return - tc * log(1.0 - curand_uniform_double(&rngState));
}
#else
__PREPROC__ double particle::uniform(){
	return dist_uniform(gen64);
}
__PREPROC__ double particle::normal01(){
	return dist_normal(gen64);
}

__PREPROC__ double particle::maxboltz(const double sqrtkT_m){
    return normal01() * sqrtkT_m;
}

__PREPROC__ double particle::unif02pi(){
    return dist_uniform(gen64) * 2.0 * M_PI;
}

__PREPROC__ double particle::exponential(const double tc){
    return - tc * log(1.0 - dist_uniform(gen64));
}
#endif

__PREPROC__ void particle::calc_next_collision_time() {
	if(gravity){
		//calculate distance to collision point
		dx = sgn(v.x) * L.x / 2.0 - pos.x;
		dz = sgn(v.z) * L.z / 2.0 - pos.z;
		
		//time to wall for x and z coordinate
		dtx = dx / v.x;
		dtz = dz / v.z;
		double y2 = v.y * v.y;
		if(sgn(v.y) <= 0.0){ //if the particle has negative y velocity
				dy = pos.y + L.y*0.5;
				double sqr = sqrt(-2.0*G_CONST*dy+y2);
				double temp1 = -(sqr+v.y)/G_CONST;
				double temp2 = (sqr-v.y)/G_CONST;
				//printf("temp1 = %f, temp2 = %f\n", temp1, temp2);
				dty = min(std::abs(temp1), std::abs(temp2));
		}
		else{
				double maxHeight = -0.5 * y2/G_CONST + pos.y;
				if(maxHeight < 0.5 * L.y){ //in this case it can't hit the ceiling
						dy = pos.y+L.y*0.5;
						double sqr = sqrt(-2.0*G_CONST*dy+y2);
						double temp1 = -(sqr+v.y)/G_CONST;
						double temp2 = (sqr-v.y)/G_CONST;
						dty = max(temp1, temp2);
				}
				else{
						dy = L.y*0.5 - pos.y; //how far to ceiling
						double sqr = sqrt(-2.0*G_CONST*dy+y2);
						double temp1 = -(sqr+v.y)/G_CONST;
						double temp2 = (sqr-v.y)/G_CONST;
						dty = min(std::abs(temp1), std::abs(temp2));
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
		
		dx = sgn(v.x) * L.x / 2.0 - pos.x;
		dy = sgn(v.y) * L.y / 2.0 - pos.y;
		dz = sgn(v.z) * L.z / 2.0 - pos.z;
		
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
		next_gas_coll_time += exponential(tc);
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

__PREPROC__ void particle::new_velocities() {
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
		phi = acos(cbrt(1.0 - uniform()));
		theta = unif02pi();
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
		v.x = maxboltz(sqrtKT_m);
		v.y = maxboltz(sqrtKT_m);
		v.z = maxboltz(sqrtKT_m);
		Vel = len(v);
	}
	else if (coll_type == 'G' && dist == 'C') {
		double3 vec;
		vec.x = normal01();
		vec.y = normal01();
		vec.z = normal01();
		double vec_norm = len(vec);
		v = Vel * vec/vec_norm;
		Vel = len(v);
	}
}

/*
Moves the particle based on the particle velocity and calcuated timestep.
*/
__PREPROC__ void particle::move() {
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

__PREPROC__ void particle::step() {
	//printf("t = %f dt = %f, v = %f %f %f, pos = %f %f %f\n", t, dt, v.x, v.y, v.z, pos.x, pos.y, pos.z);
	calc_next_collision_time(); //when do we hit something next?
	move(); //move the particle forward that amount of time
	//printf("%f %f %f %f\n", t, pos.x, pos.y, pos.z);
	//printf("t = %f dt = %f, v = %f %f %f, pos = %f %f %f\n", t, dt, v.x, v.y, v.z, pos.x, pos.y, pos.z);
	new_velocities(); //update the velocity based on the collision type
	//printf("h before = %lf, startTime = %lf, stopTime = %lf, lastOutput = %lf\n", h, t_old, t, lastOutput);
	if(integrationType == 0){
		//use the DOP853 algorithm for spin tracking
		integrateDOP(t_old, t, S, pos_old, pos, v_old, v, opt, lastOutput, lastIndex, outputArray);
	}
	else if(integrationType == 1){
		//use the hybrid RK45 method
		integrateRK45Hybrid(t_old, t, S, pos_old, pos, v_old, v, opt, h, lastOutput, lastIndex, outputArray);
	}
	else{
		//this is an unrecognized option so just don't integrate the spin in this case
	}
	//printf("h after = %lf, lastOutput = %lf\n", h, lastOutput);
	n_steps += 1;
}

/*
Convenience function that calls the step() function repeatedL.y until the end time is reached.
*/

__PREPROC__ void particle::run() {
	while (finished != true) {
		step();
	}
}
