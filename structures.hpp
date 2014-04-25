#include <iostream>
#include <iomanip>
#include <fstream>
#include <istream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>
#include <list>
#include <set>
#include <unordered_set>
#include <array>
#include <thread>
#include <iterator>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <gmpxx.h>
#include <mpfr.h>
#include <sys/time.h>
#include <functional>

const int verbose = 0;

unsigned long long int peano_hilbert_key(double x, double y, double z, int bits);

/*
	Point structure and its associated operators
*/

struct point {

	double x, y, z, q;

	unsigned long long int peanokey;

	point(void) : x(0), y(0), z(0), q(0), peanokey(0) { return; };
	point(double a, double b, double c) : x(a), y(b), z(c) {
		q = x*x + y*y + z*z;
		peanokey = peano_hilbert_key(x, y, z, 52);
		return; 
	};
};

inline std::ostream& operator<<(std::ostream &os, const point &p) {
	os << p.x << " " << p.y << " " << p.z; //<< " : pkey = " << p.peanokey;
	return os;
}

inline point operator+(const point &p1, const point &p2) {
	return point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}

inline point operator-(const point &p1, const point &p2) {
	return point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}

inline point operator*(const point &p1, const point &p2) {
	return point(p1.x * p2.x, p1.y * p2.y, p1.z * p2.z);
}

inline point operator+(const point &p1, double n) {
	return point(p1.x / n, p1.y / n, p1.z / n);
}

inline point& operator+=(point &p1, const point &p2) {

	p1.x += p2.x;
	p1.y += p2.y;
	p1.z += p2.z;

	return p1;
}

inline point operator*=(point &p, double t) {

	p.x *= t;
	p.y *= t;
	p.z *= t;

	return p;
}

struct tetra {

	std::array<point*, 4> p;
	std::array<tetra*, 4> ngb;
	std::array<unsigned int, 4> nid;

	std::array<std::array<point*, 3> , 4> f;

	/*
		Assume sane input is always given for construction
	*/

	tetra(point *a, point *b, point *c, point *d);
};

inline std::ostream& operator<<(std::ostream &os, const tetra &t) {
	os << &t << "\n" << "0 : " << *t.p[0] << "\n1 : " << *t.p[1] << "\n2 : " << *t.p[2] << "\n3 : " << *t.p[3] << "\n";

	for (int i = 0; i < 4; ++i)
		std::cout << "\tngb[" << i << "] : " << t.ngb[i] << std::endl;

	return os;
}

struct mesh {

	unsigned long int num_points;
	unsigned long int num_tetra;

	point *p_position;
	tetra *t_position;
	tetra *m_position; // mesh position (for walking)

	std::vector<point> p_buffer;
	std::vector<tetra> t_buffer;

	mesh(unsigned long int num_points) {
		this->num_points = num_points;
		num_tetra = 0;

		p_buffer.reserve(4 + num_points);
		t_buffer.reserve(16 * num_points);

		p_position = p_buffer.data();
		t_position = t_buffer.data();

		m_position = nullptr;

		return;
	};
};


/*
	Various global function declarations
*/

void write_sorted_set(std::vector<point> &buff, unsigned long int num_points, int box_length);

void walk(std::vector<std::pair<tetra*, unsigned int>> &leaves, const point *p, tetra *place, tetra *backup);

void fracture(mesh &m, std::vector<std::pair<tetra*, unsigned int>> &leaves, std::vector<tetra*> &pile, point *p);

int link_tetra(tetra *t1, tetra *t2);

unsigned int fast_triangle_ray_inclusion(const std::array<point*, 3> &f, const point *a, const point *b);

bool insphere(const tetra &t, const point &p);

void delaunay(mesh &m, std::vector<tetra*> &pile);

unsigned int halfspace_locator(const std::array<point*, 3> &f, const point *a, const point *b);
