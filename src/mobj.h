/*
 * mobj.h
 *
 *  Created on: 26 февр. 2016 г.
 *      Author: alex
 *
 *      Mathematical Objects
 */

#ifndef MOBJ_H_
#define MOBJ_H_

#include <vector>
#include <functional>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Grid/Grid.h"
#include <algorithm>

using std::abs;
using std::sqrt;

#ifdef __CUDACC__
#define CUDA_HOST_DEV_FUN __host__ __device__
#else
#define CUDA_HOST_DEV_FUN
#endif

template<int N, class T>
constexpr T pow(const T& x) {
	return N > 1 ? x*pow<(N-1)*(N>1)>(x)
			: N < 0 ? T(1)/pow<(-N)*(N<0)>(x)
					: N == 1 ? x
							: T(1);
}

struct limits {
	double lower;
	double upper;
	int n;
	double d() const { return (upper - lower) / (double)n; }
	double at(const int i) const { return lower + d()*i; }
	double atWh(const int i) const { return n>1 ? lower + dWh()*i : (upper + lower) / 2; }
	double dWh() const { return n>1 ? (upper - lower) / (double)(n - 1) : (upper - lower); }

	void scale(const double s) { lower *= s; upper *= s; }
	double width() const { return upper - lower; }

	CUDA_HOST_DEV_FUN int indd(const double v) const {
		return (int)floor(indR(v));
	}
	CUDA_HOST_DEV_FUN int indu(const double v) const {
		return (int)ceil(indR(v));
	}
private:
	CUDA_HOST_DEV_FUN double indR(const double v) const {
		if (v > upper) return n;
		if (v < lower) return 0;
		return double(n) * (v - lower) / (upper - lower);
	}
};

struct Field {
	limits Llim;
	limits Blim;
	std::vector<double> data;
	Field(const limits &Llim, const limits &Blim) : Llim(Llim), Blim(Blim) {
		data.resize(Llim.n * Blim.n);
	}
	Grid toGrid(void) {
		Grid g;
		g.data = data;
		g.nRow = Blim.n;
		g.nCol = Llim.n;
		g.xLL = Llim.lower;
		g.yLL = Blim.lower;
		g.xSize = Llim.dWh();
		g.ySize = Blim.dWh();
		return g;
	}
};

class Point {
public:
	double x;
	double y;
	double z;
	CUDA_HOST_DEV_FUN Point() : x(0), y(0), z(0) {}
	CUDA_HOST_DEV_FUN Point(const double x, const double y, const double z) : x(x), y(y), z(z) {}
	CUDA_HOST_DEV_FUN Point operator-(const Point& p) const {
		return {x-p.x, y-p.y, z-p.z};
	}
	CUDA_HOST_DEV_FUN Point operator+(const Point& p) const {
		return {x+p.x, y+p.y, z+p.z};
	}
	CUDA_HOST_DEV_FUN const Point& operator+=(const Point& p) {
		*this = *this + p;
		return *this;
	}
	CUDA_HOST_DEV_FUN Point operator*(const Point& p) const { //vector mul
		const double xt = y*p.z - z*p.y;
		const double yt = z*p.x - x*p.z;
		const double zt = x*p.y - y*p.x;
		return {xt,yt,zt};
	}
	CUDA_HOST_DEV_FUN Point operator*(const double a) const {
		return {x*a, y*a, z*a};
	}
	CUDA_HOST_DEV_FUN double operator^(const Point& p) const { //scalar mul
		return x*p.x + y*p.y + z*p.z;
	}
	CUDA_HOST_DEV_FUN Point operator/(const double a) const {
		return {x/a, y/a, z/a};
	}
	CUDA_HOST_DEV_FUN Point norm() const { //norm vector
		const double tmp = eqNorm();
		return {x/tmp, y/tmp, z/tmp};
	}
	CUDA_HOST_DEV_FUN double eqNorm() const {
		return sqrt(x*x + y*y + z*z);
	}
	CUDA_HOST_DEV_FUN bool operator==(const Point &p) const {
		return x == p.x && y == p.y && z == p.z;
	}
	CUDA_HOST_DEV_FUN bool operator!=(const Point &p) const {
		return !(*this == p);
	}
	friend std::ostream& operator<<(std::ostream& os, const Point& p);
};

class MassPoint : public Point {
public:
	double mass;
	CUDA_HOST_DEV_FUN MassPoint() : mass(0) {}
	CUDA_HOST_DEV_FUN MassPoint(const double x, const double y, const double z, const double mass) : Point(x, y, z), mass(mass) {}
	CUDA_HOST_DEV_FUN MassPoint(const Point &p, const double mass) : Point(p), mass(mass) {}
};

class Triangle {
public:
	Point p1;
	Point p2;
	Point p3;

	CUDA_HOST_DEV_FUN Triangle() {}
	CUDA_HOST_DEV_FUN Triangle(const Point a, const Point b, const Point c) : p1(a), p2(b), p3(c) {}

	CUDA_HOST_DEV_FUN Point convert(const double u, const double v) const {
		const double x = p1.x + (p3.x-p1.x)*u + (p2.x-p1.x)*v;
		const double y = p1.y + (p3.y-p1.y)*u + (p2.y-p1.y)*v;
		const double z = p1.z + (p3.z-p1.z)*u + (p2.z-p1.z)*v;
		return {x,y,z};
	}

	CUDA_HOST_DEV_FUN double intDeform() const {
		const double E = (p3.x-p1.x)*(p3.x-p1.x) + (p3.y-p1.y)*(p3.y-p1.y) + (p3.z-p1.z)*(p3.z-p1.z);
		const double F = (p3.x-p1.x)*(p2.x-p1.x) + (p3.y-p1.y)*(p2.y-p1.y) + (p3.z-p1.z)*(p2.z-p1.z);
		const double G = (p2.x-p1.x)*(p2.x-p1.x) + (p2.y-p1.y)*(p2.y-p1.y) + (p2.z-p1.z)*(p2.z-p1.z);
		return sqrt(E*G-F*F);
	}

	CUDA_HOST_DEV_FUN static bool check(const Triangle t) {
		return t.p1.eqNorm() + t.p2.eqNorm() >= t.p3.eqNorm() &&
				t.p2.eqNorm() + t.p3.eqNorm() >= t.p1.eqNorm() &&
				t.p1.eqNorm() + t.p3.eqNorm() >= t.p2.eqNorm();
	}

	CUDA_HOST_DEV_FUN Point normal() const {
		const Point vp = (p2-p1)*(p3-p1);
		return vp.norm();
	}

	friend std::ostream& operator<<(std::ostream& os, const Triangle& q);
};

class TriangleWid : public Triangle {
	bool isExtNorm;
	TriangleWid(const Triangle& t, bool isExtNorm) : Triangle(t), isExtNorm(isExtNorm) {}
};

class Quadrangle {
public:
	Point p1;
	Point p2;
	Point p3;
	Point p4;
	double val;
	Quadrangle() : Quadrangle(Point(), Point(), Point(), Point()) {}
	Quadrangle(const Point a, const Point b, const Point c, const Point d, const double val = 0)
		: p1(a), p2(b), p3(c), p4(d), val(val) {}
	CUDA_HOST_DEV_FUN Point normal() const {
		const Point vp = (p2-p1)*(p4-p1);
		return vp.norm();
	}
	friend std::ostream& operator<<(std::ostream& os, const Quadrangle& p);
	Triangle getT1() const {
		return { p1, p2, p3 };
	}
	Triangle getT2() const {
		return { p1, p3, p4 };
	}
};

class QuadrangleWid : public Quadrangle {
public:
	Point center;
	double mass = 0;
	QuadrangleWid() {}
	QuadrangleWid(const Quadrangle &q, const Point &center, const double mass) : Quadrangle(q), center(center), mass(mass) {}
};

class Tetrahedron {
public:
	Tetrahedron(const Point p, const Point b1, const Point b2, const Point b3) : p(p), base(b1, b2, b3) {}
	Point p;
	Triangle base;
	double volume() const {
		return abs((base.p2 - base.p1)*(base.p3 - base.p2)^(p - base.p1) / 6.);
	}
	Point massCenter() {
		return (p + base.p1 + base.p2 + base.p3) / 4.;
	}
};

class Hexahedron {
public:
	Point p[8];
	double dens;

	CUDA_HOST_DEV_FUN Hexahedron() : dens(0) {}
	Hexahedron(const std::vector<Point> &pin, const double dens = 0) : dens(dens) { set(pin); }
	void set(const std::vector<Point> &pin) {
		for (int i = 0; i < std::min(pin.size(), size_t(8)); ++i)
			p[i] = pin[i];
	}

	std::vector<Triangle> splitTr() const {
		std::vector<Triangle> trs = {
			Triangle(p[0], p[2], p[3]), //0
			Triangle(p[0], p[1], p[3]), //1
			Triangle(p[0], p[2], p[4]), //2
			Triangle(p[2], p[4], p[6]), //3
			Triangle(p[1], p[3], p[5]), //4
			Triangle(p[3], p[5], p[7]), //5
			Triangle(p[4], p[6], p[7]), //6
			Triangle(p[4], p[5], p[7]), //7
			Triangle(p[0], p[4], p[5]), //8
			Triangle(p[0], p[1], p[5]), //9
			Triangle(p[2], p[6], p[7]), //10
			Triangle(p[2], p[3], p[7]) //11
		};
		return trs;
	}

	std::vector<Quadrangle> splitQr() const {
		std::vector<Quadrangle> qrs = {
			Quadrangle(p[2], p[3], p[1], p[0]),
			Quadrangle(p[0], p[4], p[6], p[2]),
			Quadrangle(p[3], p[7], p[5], p[1]),
			Quadrangle(p[0], p[1], p[5], p[4]),
			Quadrangle(p[6], p[7], p[3], p[2]),
			Quadrangle(p[4], p[5], p[7], p[6])
		};
		for(int i = 0; i < 6; ++i)
			qrs[i].val = isExtNorm(qrs, i);
		return qrs;
	}

	std::vector<Tetrahedron> splitTh() const {
		std::vector<Tetrahedron> qrs = {
			Tetrahedron(p[0], p[4], p[5], p[6]),
			Tetrahedron(p[0], p[7], p[5], p[6]),
			Tetrahedron(p[0], p[2], p[6], p[3]),
			Tetrahedron(p[0], p[7], p[6], p[3]),
			Tetrahedron(p[0], p[3], p[1], p[7]),
			Tetrahedron(p[0], p[5], p[1], p[7])
		};
		return qrs;
	}

	double volume() const {
		double v = 0;
		for (auto &th : splitTh())
			v += th.volume();
		return v;
	}

	static int oppositeQr(const int i) {
		static int ops[] = {5, 2, 1, 4, 3, 0};
		return ops[i];
	}

	Point massCenter() const {
		Point r;
		double mass = 0;
		for (auto &th : splitTh()) {
			const double m = th.volume();
			r += th.massCenter()*m;
			mass += m;
		}
		return r / mass;
	}

	MassPoint getMassPoint() const {
		return MassPoint(massCenter(), volume()*dens);
	}

private:
	//возврещает 1 если у i-го четырехугольника нормаль внешняя, иначе -1
	static int isExtNorm(const std::vector<Quadrangle> &qrs, const int i) {
		Point n = qrs[i].normal();
		Point v = qrs[oppositeQr(i)].p3 - qrs[i].p3;
		return (n^v)>0? -1 : 1;
	}
};

#define CALC_NORMALS_ON_DEMAND

class HexahedronWid : public Hexahedron {
public:
#ifndef CALC_NORMALS_ON_DEMAND
	Point triNormals[12];
	HexahedronWid(const Hexahedron &h) : Hexahedron(h) { updNormals(); }
#else
	HexahedronWid(const Hexahedron &h) : Hexahedron(h) {}
#endif
	CUDA_HOST_DEV_FUN HexahedronWid() {}
	CUDA_HOST_DEV_FUN Triangle getTri(const int i) const { //Triangle normal() is always external WRT Hexahedron
		switch(i) {
			case 0: return Triangle(p[2], p[3], p[1]);
			case 1: return Triangle(p[2], p[1], p[0]);
			case 2: return Triangle(p[0], p[4], p[6]);
			case 3: return Triangle(p[0], p[6], p[2]);
			case 4: return Triangle(p[3], p[7], p[5]);
			case 5: return Triangle(p[3], p[5], p[1]);
			case 6: return Triangle(p[0], p[1], p[5]);
			case 7: return Triangle(p[0], p[5], p[4]);
			case 8: return Triangle(p[6], p[7], p[3]);
			case 9: return Triangle(p[6], p[3], p[2]);
			case 10: return Triangle(p[6], p[4], p[5]);
			case 11: return Triangle(p[6], p[5], p[7]);
			default: return Triangle();
		}
	}
	CUDA_HOST_DEV_FUN Point getTriNorm(const int i) const {
#ifndef CALC_NORMALS_ON_DEMAND
		return triNormals[i];
#else
		return getTri(i).normal();
#endif
	}
private:
#ifndef CALC_NORMALS_ON_DEMAND
	void updNormals() {
		for (int i = 0; i < 12; ++i)
			triNormals[i] = getTri(i).normal();
	}
#endif
};

class Ellipsoid {
public:
	const double Req, Rpl, e, e_, n, f; //n - Third flattening, f - flattening

	Ellipsoid(double const Req, double const Rpl) : Req(Req), Rpl(Rpl),
			e(sqrt(Req*Req - Rpl*Rpl)/Req), e_(sqrt(Req*Req - Rpl*Rpl)/Rpl), n((Req-Rpl)/(Req+Rpl)),
			f((Req-Rpl)/Req) {}

	Point getPoint(pair<double, double> lb, const double H = 0) const { //lambda, B
		return getPoint(Req, Rpl, lb.second, lb.first, H);
	}

	Point getPoint(double const B, double const l, const double H = 0) const {
		return getPoint(Req, Rpl, B, l, H);
	}

	static Point getNormal(double const B, double const l) {
		return Point(cos(B)*cos(l), cos(B)*sin(l), sin(B));
	}

	double Q() const { //Meridian arc unit
		const double n2 = n*n;
		return (Req/(1+n)) * (1 + n2/4 + n2*n2/64 + n2*n2*n2/256);
	}

	static Point getPoint(double const Req, double const Rpl, double const B, double const l, const double H) {
		const double L = 1/sqrt(Req*Req*cos(B)*cos(B) + Rpl*Rpl*sin(B)*sin(B));
		const double x = (Req*Req*L+H)*cos(B)*cos(l);
		const double y = (Req*Req*L+H)*cos(B)*sin(l);
		const double z = (Rpl*Rpl*L+H)*sin(B);
		return Point(x,y,z);
	}

	void map(std::function<void (const int Li, const int Bi, const int Hi)> fun,
			const limits Llim, const limits Blim, const limits Hlim = {0,1,1}) const {
		for(int Hi = 0; Hi < Hlim.n; ++Hi)
				for(int Bi = 0; Bi < Blim.n; ++Bi)
					for(int Li = 0; Li < Llim.n; ++Li)
						fun(Li, Bi, Hi);
	}

	void mapPB(std::function<void (const int Li, const int Bi, const int Hi)> fun,
			const limits Llim, const limits Blim, const limits Hlim = {0,1,1}) const {
		const unsigned long max = Hlim.n*Blim.n*Llim.n;
		unsigned long cur = 0;
		std::cout << "  0.00%" << std::flush;
		for(int Hi = 0; Hi < Hlim.n; ++Hi)
				for(int Bi = 0; Bi < Blim.n; ++Bi)
					for(int Li = 0; Li < Llim.n; ++Li) {
						fun(Li, Bi, Hi);
						++cur;
						printf("\r%3u.%02u%%", (unsigned)(cur*100UL/max), (unsigned)((cur*10000UL/max)%100));
						std::cout << std::flush;
					}
		std::cout << "\r 100.00%  " << std::endl;
	}
};

class SubCube {
	limits x;
	limits y;
	limits z;
public:
	Point p;
	SubCube() {}
	SubCube(const limits &x, const limits &y, const limits &z) : x(x), y(y), z(z) {}

	template<typename OneVarFun>
	CUDA_HOST_DEV_FUN
	void range(const double d, OneVarFun f) const {
		const int Nz = z.indu(p.z + d);
		for (int zi = z.indd(p.z - d); zi < Nz; ++zi)
			for (int yi = y.indd(p.y - d); yi < y.indu(p.y + d); ++yi)
				for (int xi = x.indd(p.x - d); xi < x.indu(p.x + d); ++xi)
					f((zi*y.n + yi)*x.n + xi);
	}
};

struct FieldPoint {
	Point p;
	Point n;
	double v;
};

void toRad(limits &l);
double toRad(const double a);
void toDeg(limits &l);
double toDeg(const double a);

std::ostream& operator<<(std::ostream& os, const limits& l);
std::ostream& operator<<(std::ostream& os, const Point& p);
std::ostream& operator<<(std::ostream& os, const Triangle& q);
std::ostream& operator<<(std::ostream& os, const Quadrangle& q);

#endif /* MOBJ_H_ */
