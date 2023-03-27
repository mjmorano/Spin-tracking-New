#pragma once

struct double3 {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

struct float3 {
	float x;
	float y;
	float z;
};

struct outputDtype {
	double t;
	double3 x;
	double3 s;
};

#pragma acc routine seq
double3 operator+(const double3, const double3);
#pragma acc routine seq
double3 operator+(const double3, const double);
#pragma acc routine seq
double3 operator+(const double, const double3);
#pragma acc routine seq
double3 operator-(const double3, const double3);
#pragma acc routine seq
double3 operator*(const double3, const double);
#pragma acc routine seq
double3 operator*(const double, const double3);
#pragma acc routine seq
double3 operator*(const double3, const double3);
#pragma acc routine seq
double3 operator/(const double3, const double);
#pragma acc routine seq
double3 operator/(const double3, const double3);
#pragma acc routine seq
double3 cross(const double3, const double3);
#pragma acc routine seq
double sum(const double3);
#pragma acc routine seq
double len(const double3);
#pragma acc routine seq
double3 norm(const double3);
#pragma acc routine seq
double3 fabs3(const double3);
#pragma acc routine seq
double3 max_d3(const double3, const double3);
#pragma acc routine seq
double max3(const double3);
#pragma acc routine seq
double3 sgn(const double3);

struct quaternion{
	double w;
	double x; 
	double y;
	double z;
};

#pragma acc routine seq
quaternion operator*(const quaternion, const quaternion);
#pragma acc routine seq
quaternion conjugate(const quaternion);
quaternion qMult(const quaternion, const quaternion);
double3 qv_mult(const quaternion, const double3);
quaternion rodriguezQuat(const double3, const double);