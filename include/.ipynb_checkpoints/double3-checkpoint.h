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
	float t;
	float3 x;
	float3 s;
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
double3 sgn(const double3);

