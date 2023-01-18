#ifndef __DOUBLE3_h__
#define __DOUBLE3_h__

struct double3 {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

inline double3 operator+(const double3, const double3);
inline double3 operator+(const double3, const double);
inline double3 operator+(const double, const double3);
inline double3 operator-(const double3, const double3);
inline double3 operator*(const double3, const double);
inline double3 operator*(const double, const double3);
inline double3 operator*(const double3, const double3);
inline double3 operator/(const double3, const double);
inline double3 operator/(const double3, const double3);
inline double3 cross(const double3, const double3);
inline double sum(const double3);
inline double len(const double3);
inline double3 norm(const double3);
inline double3 fabs(const double3);
inline double3 max_d3(const double3, const double3);
inline double3 sgn(const double3);


#endif
