#pragma once

#if(__NVCC__)
#define __PREPROC__ __host__ __device__
#elif(__HIPCC__)
#define __PREPROC__ __host__ __device__
#else
#define __PREPROC__
#endif

struct double3 {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

struct float3 {
	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;
};

struct outputDtype {
	double t = 0.0;
	double3 x;
	double3 s;
};

__PREPROC__ double3 operator+(const double3, const double3);
__PREPROC__ double3 operator+(const double3, const double);
__PREPROC__ double3 operator+(const double, const double3);
__PREPROC__ double3 operator-(const double3, const double3);
__PREPROC__ double3 operator*(const double3, const double);
__PREPROC__ double3 operator*(const double, const double3);
__PREPROC__ double3 operator*(const double3, const double3);
__PREPROC__ double3 operator/(const double3, const double);
__PREPROC__ double3 operator/(const double3, const double3);
__PREPROC__ double3 cross(const double3, const double3);
__PREPROC__ double sum(const double3);
__PREPROC__ double len(const double3);
__PREPROC__ double3 norm(const double3);
__PREPROC__ double3 fabs3(const double3);
__PREPROC__ double3 max_d3(const double3, const double3);
__PREPROC__ double max3(const double3);
__PREPROC__ double3 sgn(const double3);

struct quaternion{
	double w = 0.0;
	double x = 0.0; 
	double y = 0.0;
	double z = 0.0;
};

__PREPROC__ quaternion operator*(const quaternion, const quaternion);
__PREPROC__ quaternion conjugate(const quaternion);
__PREPROC__ quaternion qMult(const quaternion, const quaternion);
__PREPROC__ double3 qv_mult(const quaternion, const double3);
__PREPROC__ quaternion rodriguezQuat(const double3, const double);