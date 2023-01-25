#include "../include/double3.h"
#include <cmath>

double asm_abs(double val){
	return std::abs(val);
}

double3 operator+(const double3 a, const double3 b){
	double3 out;
	out.x = a.x + b.x;
	out.y = a.y + b.y; 
	out.z = a.z + b.z;
	return out;
}

double3 operator+(const double3 a, const double b){
	double3 out;
	out.x = a.x + b;
	out.y = a.y + b;
	out.z = a.z + b;
	return out;
}

double3 operator+(const double a, const double3 b){
	return b+a;
}

double3 operator-(const double3 a, const double3 b){
	double3 out;
	out.x = a.x - b.x;
	out.y = a.y - b.y;
	out.z = a.z - b.z;
	return out;
}

double3 operator*(const double3 a, const double b){
	double3 out;
	out.x = a.x * b;
	out.y = a.y * b;
	out.z = a.z * b;
	return out;
}

double3 operator*(const double b, const double3 a){
	return a * b;
}

double3 operator*(const double3 a, const double3 b){
	double3 out;
	out.x = a.x * b.x;
	out.y = a.y * b.y;
	out.z = a.z * b.z;
	return out;
}

double3 operator/(const double3 a, const double b){
	double3 out;
	out.x = a.x / b;
	out.y = a.y / b;
	out.z = a.z / b;
	return out;
}
double3 operator/(const double3 a, const double3 b){
	double3 out;
	out.x = a.x/b.x;
	out.y = a.y/b.y;
	out.z = a.z/b.z;
	return out;
}

double3 cross(const double3 a, const double3 b){
	double3 out;
	out.x = a.y*b.z - a.z*b.y;
	out.y = a.z*b.x - a.x*b.z;
	out.z = a.x*b.y - a.y*b.x;
	return out;
}

double sum(const double3 a){
	return a.x+a.y+a.z;
}

double len(const double3 a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

double3 norm(const double3 a){	
	return a/len(a);
}

double3 fabs3(const double3 a){
	double3 out;
	out.x = asm_abs(a.x);
	out.y = asm_abs(a.y);
	out.z = asm_abs(a.z);
	return out;
}

double3 max_d3(const double3 a, const double3 b){
	double3 out;
	out.x = (a.x>b.x)?a.x:b.x;
	out.y = (a.y>b.y)?a.y:b.y;
	out.z = (a.z>b.z)?a.z:b.z;
	return out;
}

template <typename T>
int sgn (T val){
	return (T(0) < val) - (val < T(0));

  double3 sgn(const double3 a);
	double3 out;
	out.x = sgn(a.x);
	out.y = sgn(a.y);
	out.z = sgn(a.z);
	return out;
}
