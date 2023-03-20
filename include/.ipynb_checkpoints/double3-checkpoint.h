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
double3 sgn(const double3);

struct matrix{
	double xx = 0.0;
	double xy = 0.0;
	double xz = 0.0;
	double yx = 0.0;
	double yy = 0.0;
	double yz = 0.0;
	double zx = 0.0;
	double zy = 0.0;
	double zz = 0.0;
}

#pragma acc routine seq
matrix operator+(const matrix, const matrix);
#pragma acc routine seq
matrix operator+(const matrix, const double);
#pragma acc routine seq
matrix operator+(const double, const matrix);
#pragma acc routine seq
matrix operator-(const matrix, const matrix);
#pragma acc routine seq
matrix operator*(const matrix, const double);
#pragma acc routine seq
matrix operator*(const double, const matrix);
#pragma acc routine seq
matrix operator*(const matrix, const matrix);
#pragma acc routine seq
matrix operator*(const matrix, const double3);
#pragma acc routine seq
matrix operator/(const matrix, const double);
#pragma acc routine seq
matrix operator/(const matrix, const matrix);
#pragma acc routine seq
matrix determineMatrix(const double3);