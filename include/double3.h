#pragma once

struct double3 {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

double3 operator+(const double3, const double3);
double3 operator+(const double3, const double);
double3 operator+(const double, const double3);
double3 operator-(const double3, const double3);
double3 operator*(const double3, const double);
double3 operator*(const double, const double3);
double3 operator*(const double3, const double3);
double3 operator/(const double3, const double);
double3 operator/(const double3, const double3);
double3 cross(const double3, const double3);
double sum(const double3);
double len(const double3);
double3 norm(const double3);
double3 fabs3(const double3);
double3 max_d3(const double3, const double3);
double3 sgn(const double3);

