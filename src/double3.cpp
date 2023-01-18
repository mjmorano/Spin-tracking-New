#include "../include/double3.h"
#include <cmath>

inline double3 operator+(const double3 a, const double3 b){
	double3 out;
	out.x = a.x + b.x;
	out.y = a.y + b.y; 
	out.z = a.z + b.z;
	return out;
}

inline double3 operator+(const double3 a, const double b){
	double3 out;
	out.x = a.x + b;
	out.y = a.y + b;
	out.z = a.z + b;
	return out;
}

inline double3 operator+(const double a, const double3 b){
	return b+a;
}

inline double3 operator-(const double3 a, const double3 b){
	double3 out;
	out.x = a.x - b.x;
	out.y = a.y - b.y;
	out.z = a.z - b.z;
	return out;
}

inline double3 operator*(const double3 a, const double b){
	double3 out;
	out.x = a.x * b;
	out.y = a.y * b;
	out.z = a.z * b;
	return out;
}

inline double3 operator*(const double b, const double3 a){
	return a * b;
}

inline double3 operator*(const double3 a, const double3 b){
	double3 out;
	out.x = a.x * b.x;
	out.y = a.y * b.y;
	out.z = a.z * b.z;
	return out;
}

inline double3 operator/(const double3 a, const double b){
	double3 out;
	out.x = a.x / b;
	out.y = a.y / b;
	out.z = a.z / b;
	return out;
}
inline double3 operator/(const double3 a, const double3 b){
	double3 out;
	out.x = a.x/b.x;
	out.y = a.y/b.y;
	out.z = a.z/b.z;
	return out;
}

inline double3 cross(const double3 a, const double3 b){
	double3 out;
	out.x = a.y*b.z - a.z*b.y;
	out.y = a.z*b.x - a.x*b.z;
	out.z = a.x*b.y - a.y*b.x;
	return out;
}

inline double sum(const double3 a){
	return a.x+a.y+a.z;
}

inline double len(const double3 a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

inline double3 norm(const double3 a){	
	return a/len(a);
}

inline double3 fabs(const double3 a){
	double3 out;
	out.x = fabs(a.x);
	out.y = fabs(a.y);
	out.z = fabs(a.z);
	return out;
}

inline double3 max_d3(const double3 a, const double3 b){
	double3 out;
	out.x = (a.x>b.x)?a.x:b.x;
	out.y = (a.y>b.y)?a.y:b.y;
	out.z = (a.z>b.z)?a.z:b.z;
	return out;
}

template <typename T>
int sgn (T val){
	return (T(0) < val) - (val < T(0));

inline double3 sgn(const double3 a);
	double3 out;
	out.x = sgn(a.x);
	out.y = sgn(a.y);
	out.z = sgn(a.z);
	return out;
}
