
#define GEOGRAPHICLIB_SHARED_LIB 0
#include <GeographicLib/TransverseMercator.hpp>

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <atomic>
#include <algorithm>
#include <initializer_list>
#include <numeric>
#include "mobj.h"
#include "calcField.h"
#include "inputLoader.h"
#include "Stopwatch.h"
#include "Quadrangles.h"
#include "MPIwrapper.h"
#include "MPIpool.h"

#include "sharedMem.h"
using gElementsShared = std::vector<HexahedronWid, ShmemAllocator<HexahedronWid>>;

using std::string;
using std::vector;
using GeographicLib::TransverseMercator;

#define MU_0 1.2566370614e-6
#define G_CONST  1	//-6.67408
#define MB (1024*1024)

#define MIN_DENS 1e-8

#define R_EQ 6378.245		//equatorial radius in km
#define R_PL 6356.863		//polar radius in km

//#define R_EQ 6367.558487
//#define R_PL 6367.558487

#define Assert(exp) do { if (!(exp)) throw std::runtime_error("Assertion failed at: " + string(__FILE__) + " # line " + string(std::to_string(__LINE__))); } while (0)

//get amount of quadrangles for nx*ny*nz discretization
int getQdAm(const int nx, const int ny, const int nz) {
	const int n_2 = 3 * nx*ny*nz;
	return n_2 + nx*ny + nx*nz + ny*nz;
}

//get amount of hexahedrons for nx*ny*nz discretization
int getHexAm(const int nx, const int ny, const int nz) {
	return nx*ny*nz;
}

//x (East direction) to Gaus-Kruger
double xToGK(const double x, const double l0) {
	const int zone = int(l0 / toRad(6.)) + 1;
	return x + (zone*1e3 + 5e2);
}

//Gaus-Kruger to x (East direction)
double xFromGK(const double x, const double l0) {
	const int zone = int(l0 / toRad(6.)) + 1;
	return x - (zone*1e3 + 5e2);
}

//estimate approximate buffer size for the Hexahedrons that can't be replaced by singular source
int triBufferSize(const limits &Nlim, const limits &Elim, const limits &Hlim, const double r) {
	auto f = [&](const limits &lim)->int {return (int)ceil(1.62*r*double(lim.n) / lim.width()); };
	const int v1 = f(Nlim)*f(Elim)*f(Hlim);
	const int v2 = Nlim.n*Elim.n*Hlim.n;
	if(v2 < 40000) return v2;
	return std::min(v1, v2);
}

struct Circle {
	Point2D center;
	double r;
	double distToCenter(const Point &p) const {
		return (center - p).eqNorm();
	}
	bool isIn(const Point &p) const {
		return distToCenter(p) < r;
	}
	Point2D toBorder(const Point &p) const {
		return center + (p - center) * r / (p - center).eqNorm();
	}
};
struct Cylinder : public Circle {
	double h;
	Cylinder(const Circle& c, const double h) : Circle(c), h(h) {}
};
struct QuadrangleRef {
	Point& p1;
	Point& p2;
	Point& p3;
	Point& p4;
	Point J;
	double k;
	bool inner = false;
	bool isCrossing(const Circle& s) const {
		return s.isIn(p1) || s.isIn(p2) || s.isIn(p3) || s.isIn(p4);
	}
	bool isIn(const Circle& s) const {
		return s.isIn(p1) && s.isIn(p2) && s.isIn(p3) && s.isIn(p4);
	}
	operator Quadrangle() const {
		return { p1, p2, p3, p4 };
	}
};

template <class VAlloc>
void wellGen(const Volume &v, const Cylinder &well, const Point Hprime, const double Koutter, const vector<double> Kinner,
		vector<HexahedronWid, VAlloc> &hsi, vector<double> &K) {

	Assert(Kinner.size() == v.z.n);

	const Point Joutter = Hprime*Koutter;
	//const Point J0inner = Hprime*Kinner;

	//make flat mesh
	vector<Point> mesh((v.x.n+1) * (v.y.n+1));
	const limits xLim = { v.x.lower - v.x.dWh() / 2.,  v.x.upper + v.x.dWh() / 2., v.x.n };
	const limits yLim = { v.y.lower - v.y.dWh() / 2.,  v.y.upper + v.y.dWh() / 2., v.y.n };
	for (int yi = 0; yi < yLim.n + 1; ++yi)
		for (int xi = 0; xi < xLim.n + 1; ++xi)
			mesh[yi*(xLim.n+1) + xi] = Point(xLim.at(xi), yLim.at(yi), 0.);
	
	//make lateral Qadrangles
	vector<QuadrangleRef> qrs;
	for (int yi = 0; yi < yLim.n; ++yi)
		for (int xi = 0; xi < xLim.n; ++xi)
			qrs.push_back({ mesh[(yi + 1)*(xLim.n+1) + xi + 1], mesh[(yi + 1)*(xLim.n+1) + xi], mesh[yi*(xLim.n+1) + xi + 1], mesh[yi*(xLim.n+1) + xi], Joutter, Koutter });
	
	//make lateral circle in mesh
	for (QuadrangleRef& q : qrs) {
		if (q.isIn(well))
			q.inner = true;
	}
	for (QuadrangleRef& q : qrs) {
		if (!q.inner && q.isCrossing(well)) {
			if (well.isIn(q.p1)) q.p1 = well.toBorder(q.p1);
			if (well.isIn(q.p2)) q.p2 = well.toBorder(q.p2);
			if (well.isIn(q.p3)) q.p3 = well.toBorder(q.p3);
			if (well.isIn(q.p4)) q.p4 = well.toBorder(q.p4);
		}
	}

	//make hexahedrons
	hsi.resize(xLim.n*yLim.n*v.z.n);
	K.resize(xLim.n*yLim.n*v.z.n);
	
//#pragma omp parallel for
	for (int zi = 0; zi < v.z.n; ++zi)
		for (int yi = 0; yi < yLim.n; ++yi)
			for (int xi = 0; xi < xLim.n; ++xi) {
				const QuadrangleRef& cur = qrs[yi*xLim.n + xi];
				const double Kval = cur.inner ? Kinner[zi] : cur.k;
				const Point Jval = cur.inner ? Hprime*Kinner[zi] : cur.J;

				const int ind = (zi*yLim.n + yi)*xLim.n + xi;
				hsi[ind] = Hexahedron{
					(Quadrangle)cur + Point{0, 0, v.z.at(zi)},
					(Quadrangle)cur + Point{0, 0, v.z.at(zi + 1)},
					Jval };
				K[ind] = Kval;
			}

	auto makeOuterHex = [&Joutter](const Point& llt, const Point& rub) { //left low top, right upper bottom
		return Hexahedron{ {
			{rub.x, rub.y, llt.z},
			{ rub.x, llt.y, llt.z },
			{ llt.x, rub.y, llt.z },
			llt,
			rub,
			{ rub.x, llt.y, rub.z },
			{ llt.x, rub.y, rub.z },
			{llt.x, llt.y, rub.z}
			}, Joutter };
	};
	
	//add 4 external quadrangles
	const double l = 1e8;
	hsi.push_back(makeOuterHex({ xLim.lower, - l, v.z.upper}, { l, yLim.lower, v.z.lower }));
	hsi.push_back(makeOuterHex({ xLim.upper, yLim.lower, v.z.upper }, { l, l, v.z.lower }));
	hsi.push_back(makeOuterHex({- l, yLim.upper, v.z.upper }, { xLim.upper, l, v.z.lower }));
	hsi.push_back(makeOuterHex({- l, - l, v.z.upper }, { xLim.lower, yLim.upper, v.z.lower }));
	K.push_back(Koutter);
	K.push_back(Koutter);
	K.push_back(Koutter);
	K.push_back(Koutter);
}

template <class VAlloc>
void cubeGen(const Volume &v, const Point J, vector<HexahedronWid, VAlloc> &hsi) {
	hsi.resize(v.x.n * v.y.n * v.z.n);

	for (int zi = 0; zi < v.z.n; ++zi)
		for (int yi = 0; yi < v.y.n; ++yi)
			for (int xi = 0; xi < v.x.n; ++xi) {
				Quadrangle cur{
					Point{v.x.at(xi+1), v.y.at(yi+1), 0}, 
					Point{v.x.at(xi+1), v.y.at(yi), 0}, 
					Point{v.x.at(xi), v.y.at(yi+1), 0}, 
					Point{v.x.at(xi), v.y.at(yi), 0}, 
				};

				const int ind = (zi*v.y.n + yi)*v.x.n + xi;
				hsi[ind] = Hexahedron{
					cur + Point{0, 0, v.z.at(zi + 1)},
					cur + Point{0, 0, v.z.at(zi)},
					J };
			}	
}

template <class VAlloc>
void ellipsoidGen(const Ellipsoid &e, const int nl, const int nB, const int nR, const Point J, vector<HexahedronWid, VAlloc> &hsi) {
	limits ll{0, M_PI_2, nl}, lB{0, M_PI_2, nB}, lReq{0, e.Req, nR}, lRpl{0, e.Rpl, nR};
	// hsi.resize(ll.n * lB.n * lReq.n);
	hsi.resize(0);
	const int riStart = 0;

	const auto makeQuadrangle = [&](const Ellipsoid &e, const int Bi, const int li) {
		return Quadrangle{
			e.getPoint(lB.at(Bi+1), ll.at(li+1)),
			e.getPoint(lB.at(Bi), ll.at(li+1)),
			e.getPoint(lB.at(Bi+1), ll.at(li)),
			e.getPoint(lB.at(Bi), ll.at(li)),
		};
	};

	for(int ri = riStart; ri < lReq.n; ++ri)
		for(int li = 0; li < ll.n; ++li)
			for(int Bi = 0; Bi < lB.n; ++Bi) {
				Quadrangle internal;
				if(ri > 0) {
					Ellipsoid ei{lReq.at(ri), lRpl.at(ri)};
					internal = makeQuadrangle(ei, Bi, li);
				}
				Ellipsoid ee{lReq.at(ri+1), lRpl.at(ri+1)};
				Quadrangle external = makeQuadrangle(ee, Bi, li);
				// const int ind = (ri*ll.n + li)*lB.n + Bi;
				// hsi[ind] = Hexahedron{ external, internal, J };
				hsi.push_back(Hexahedron{ external, internal, J });
			}
			

	const auto dbl = [&](const std::function<void(Hexahedron&)> &f) {
		vector<Hexahedron> tmp(hsi.begin(), hsi.end());
		for(auto& h: tmp) f(h);
		hsi.insert(hsi.end(), tmp.begin(), tmp.end());
	};
	dbl([](Hexahedron& h){ h.mirrorX(); });
	dbl([](Hexahedron& h){ h.mirrorY(); });
	dbl([](Hexahedron& h){ h.mirrorZ(); });
}

template <class VAlloc>
void makeCloud(const vector<HexahedronWid, VAlloc> &hsi, const string datFname) {
	Dat2D<> dat;
	for (const auto& h : hsi)
		if(h.dens.x != 0)
			for (const auto& p : h.p)
				dat.es.push_back({ {p.x, p.y}, p.z });
	dat.write(datFname);
}

template <class VAlloc>
void makeBln(const vector<HexahedronWid, VAlloc> &hsi, const string datFname, const bool isIn = false) {
	Dat2D<> dat;
	for (const auto& h : hsi)
		if (isIn? (h.dens.x == 0) : (h.dens.x != 0)) {
			/*
			dat.es.push_back({ { 5, 0 }, 0 });
			dat.es.push_back({ { h.p[0].x, h.p[0].y }, 0 });
			dat.es.push_back({ { h.p[1].x, h.p[1].y }, 0 });
			dat.es.push_back({ { h.p[3].x, h.p[3].y }, 0 });
			dat.es.push_back({ { h.p[2].x, h.p[2].y }, 0 });
			dat.es.push_back({ { h.p[0].x, h.p[0].y }, 0 });
			*/
			for (int i = 0; i < 2; ++i) {
				const auto t = h.getTri(i);
				//const Triangle t = h.getTri(i);
				dat.es.push_back({ { 4, 0 }, 0 });
				dat.es.push_back({ { t.p1.x, t.p1.y }, 0 });
				dat.es.push_back({ { t.p2.x, t.p2.y }, 0 });
				dat.es.push_back({ { t.p3.x, t.p3.y }, 0 });
				dat.es.push_back({ { t.p1.x, t.p1.y }, 0 });
			}
		}
	dat.write(datFname);
}

class VolumeMod : public Volume {
public:
	template<typename... Args>
	VolumeMod(const Volume& v) : Volume(v) {}
	void save(const string& fname) const {
		Dat3D<Point> dat;
		for (int k = 0; k < z.n; ++k)
			for (int j = 0; j < y.n; ++j)
				for (int i = 0; i < x.n; ++i) {
					const Point p0{ x.atWh(i), y.atWh(j), z.atWh(k) };
					dat.es.push_back({ { p0.x, p0.y, p0.z }, Point() });
				}
	}
};

Point intHexTr__(const Point &p0, const HexahedronWid &h);

void hexTest() {
	// Hexahedron h({
	// 	{ 20, 20, 0 },{ 20, -20, 0 },{ -20, 20, 0 },{ -20, -20, 0 },
	// 	{ 20, 20, -4 },{ 20, -20, -4 },{ -20, 20, -4 },{ -20, -20, -4 }
	// }, Point{ 14, 14, 35 }*0.2);

	Hexahedron h({
		{ 2, 1, 0 },{ 2, -1, 0 },{ -2, 1, 0 },{ -2, -1, 0 },
		{ 2, 1, -2 },{ 2, -1, -2 },{ -2, 1, -2 },{ -2, -1, -2 }
	}, Point{ 14, 14, 35 }*0.02);

	auto hw = HexahedronWid(h);
	const double H = -0.25;


	cout << "Solving..." << endl;
	Dat3D<Point> dat;
	auto l = limits{ -25 + 0.00001, 25 + 0.00001, 40 };
	for (double i = 0; i < l.n; ++i) {
		for (int j = 0; j < l.n; ++j) {
			const Point p0{ l.atWh(j), l.atWh(i), H };
			const Point res = (-intHexTr__(p0, hw)
					 + (hw.isIn(p0)? hw.dens * (4.*M_PI / 3.) : Point())
				) / (4 * M_PI);
			dat.es.push_back({ { p0.x, p0.y, p0.z }, res });
		}
	}
	
	dat.write("cubeFieldSZ_IN.dat");
	cout << "Done." << endl;
}

template<class Input = double, class Acc = Input, class Result = Acc>
class Statistics {
	std::function<Acc(const Input&, const Acc&)> add;
	std::function<Result(const Acc&, const unsigned int count)> end;
	Acc init;
	Acc acc;
	unsigned int count = 0;
public:
	Statistics(
		const Acc &acc = {}, 
		const std::function<Acc(const Input&, const Acc&)> &add = [](const Input& i, const Acc& acc){ return i + acc; },
		const std::function<Result(const Acc&, const unsigned int count)> &end = [](const Acc& acc, const unsigned int count){ return acc; }
	) : acc(acc), init(acc), add(add), end(end) {}

	void next(const Input& i) {
		acc = add(i, acc);
		++count;
	}
	Result get() const {
		return end(acc, count);
	}
	void reset() {
		acc = init;
		count = 0;
	}
};

double field_Hz(const Ellipsoid &e, const Point J, const Point p) {
	const double a = e.Req, b = e.Rpl;
	const double x = p.x, y = p.y, z = p.z;

	const Point M0 = J * (4. * M_PI * a*a * b / 3.);
	const double q = sqrt(a*a - b*b);
	const double tmp = x*x + z*z + q;
	const double t = sqrt(tmp*tmp - 4. * z * z * q * q);

	const double y1 = z, x1 = sqrt(x*x + y*y);
	const double t1 = a*a + b*b - x1*x1 - y1*y1;
	const double k = ( sqrt(t1*t1 - 4. * (a*a * b*b - a*a * y1*y1 - b*b * x1*x1)) - t1 ) / 2.;
	const double c1 = sqrt(b*b + k);

	// cout << M0 << " - " << q << " " << tmp << " " << t << " - " << y1 << " " << t1 << " " << k << " " << c1 << endl;

	const double Hz = M0.eqNorm() * (3. * (z / (t * c1*c1*c1) - 1./(q*q) * (1./c1 - 1./q * atan2(q, b))));

	return Hz;
}

Point field_sphere_H(const double R, const Point J, const Point p) {
	const Point M0 = J * (4. * M_PI) * R*R*R / 3.;
	const double mu = M0.eqNorm();
	const Point n = M0.norm();
	const double r = p.eqNorm();
	const Point H = ( p * 3. * (n^p) / (r*r*r*r*r) - n / (r*r*r) ) * mu;
	return H;
}

double field_sphere_H_in_Hz(const double R, const double Hprim_z, const double K, const Point p0) {
	const double s = p0.x*p0.x + p0.y*p0.y + p0.z*p0.z;
	const double dr = sqrt(s*s*s*s*s);
	
	return K / (K+3) * Hprim_z * R*R*R * (2*p0.z*p0.z - p0.x*p0.x - p0.y*p0.y)/dr;
}

Point magnetization_J_theor_ellipsoid(const Ellipsoid &e, const Point J0, const double K) {
	if(e.Rpl == e.Req) {	// Sphere
		return J0 / (1. + K/3.);
	}

	Assert(e.Rpl > e.Req);

	const double m = e.Rpl / e.Req;
	// cout << m << endl;
	const double lnPt = log((m + sqrt(m*m - 1))/(m - sqrt(m*m - 1)));
	// cout << lnPt << endl;
	const double L = (1. / (m*m - 1.)) * (
		(m / (2. * sqrt(m*m - 1))) * lnPt - 1 
	);
	// cout << L << endl;
	const double M = (m / (2. * (m*m - 1.))) * (
		m - (1. / (2. * sqrt(m*m - 1))) * lnPt 
	);
	// cout << M << endl;
	const Point I = J0 / (Point(1.) + Point(K).cmul({M, M, L}));
	return I;
}

class WellDemagCluster : public MPIwrapper {
public:
	const int maxGPUmemMB = 5000;
	int triBufferSize = 0;
	using Qiter = gElementsShared::const_iterator;

	SharedMemBase<gElementsShared> *sharedMem = 0;

	WellDemagCluster(const vector<int> &gpuIdMap = {}) : MPIwrapper() {
		if (gridSize < 2) throw std::runtime_error("You must run at least 2 MPI processes.");
		if (root != 0) throw std::runtime_error("Root process must have rank = 0.");
		const auto lid = localId();
		const int devId = !isRoot() && std::get<1>(lid) ? std::get<0>(lid) - 1 : std::get<0>(lid);
		const int mappedDevId = devId < gpuIdMap.size()? gpuIdMap[devId] : devId;
		cuSolver::setDevice(mappedDevId);
	}

	void runEllipsoid(int argc, char *argv[]) {
		InputParser inp(argc, argv);
		const Ellipsoid e(10, 20);

		int nl = 80, nB = 80, nR = 40;

		const double K = 2;
		const Point Hprime = { 14, 14, 35 }; //~40A/m
		// const Point Hprime = { 0, 0, 50 };
		const auto I0 = Hprime * K;

		const auto ellipsoidModelGenerator = [&](vector<HexahedronWid> &hsi, vector<double> &Kmodel, vector<Point> &I0out){
			// const auto Ipres = I0 / (1. + K/3.);	// Use known precise I
			// ellipsoidGen(e, nl, nB, nR, Ipres, hsi);
			ellipsoidGen(e, nl, nB, nR, I0, hsi);
			Kmodel.assign(hsi.size(), K);
			I0out.assign(hsi.size(), I0);
		};
		const auto cubeModelGenerator = [&](vector<HexahedronWid> &hsi, vector<double> &Kmodel, vector<Point> &I0out) {
			cubeGen(Volume{{-10, 10, nl},{-5, 5, nB}, {-5, 5, nR}}, I0, hsi);
			Kmodel.assign(hsi.size(), K);
			I0out.assign(hsi.size(), I0);
		};
		const auto createCudaSolver = [&](const vector<HexahedronWid>& hsi, const bool transpose) {
			return gFieldSolver::getCUDAsolver(&*hsi.cbegin(), &*hsi.cend(), transpose);
			// const double replDist = 3;
			// return gFieldSolver::getCUDAreplacingSolver(&*hsi.cbegin(), &*hsi.cend(), replDist, nl*nB*nR*8);
		};

		Stopwatch tmr;
		tmr.start();
		vector<HexahedronWid> hsi = demagCG<HexahedronWid>(ellipsoidModelGenerator, createCudaSolver);
		if(!isRoot()) return;
		cout << "Total time: " << tmr.stop() << "sec." << endl;

		const int layersN = hsi.size() / (nl*nB*8);
		cout << "Layers: " << layersN << endl;

		// Mean dens for model
		Statistics<HexahedronWid, Point> meanStat(
			Point{},
			[](auto& h, auto& acc){ return h.dens + acc; },
			[](auto& acc, auto count){ return acc / count; }
		);
		// Mean dens by layer
		std::vector meanStatLayers(layersN, Statistics<HexahedronWid, Point>(
			Point{},
			[](auto& h, auto& acc){ return h.dens + acc; },
			[](auto& acc, auto count){ return acc / count; }
		));
		const int inSz = nl*nB;
		for(int part = 0; part < 8; ++part)
			for(int layer = 0; layer < layersN; ++layer) 
				for(int i = 0; i < inSz; ++i) {
					const int idx = ((part * layersN) + layer) * inSz + i;
					meanStatLayers[layer].next(hsi[idx]);
				}
		for(auto& h: hsi) meanStat.next(h);
		const Point mean = meanStat.get();


		// RMS dens for model
		Statistics<Point, double> rmsStat(
			0,
			[&](auto& diff, auto& acc){
				return acc + (diff^diff); 
			},
			[](auto& acc, auto count){ return std::sqrt(acc / count); }
		);
		// RMS dens by layer
		std::vector rmsStatLayers(layersN, Statistics<Point, double>(
			0,
			[&](auto& diff, auto& acc){
				return acc + (diff^diff); 
			},
			[](auto& acc, auto count){ return std::sqrt(acc / count); }
		));
		for(int part = 0; part < 8; ++part)
			for(int layer = 0; layer < layersN; ++layer) 
				for(int i = 0; i < inSz; ++i) {
					const int idx = ((part * layersN) + layer) * inSz + i;
					rmsStatLayers[layer].next(hsi[idx].dens - meanStatLayers[layer].get());
				}
		for(auto& h: hsi) rmsStat.next(h.dens - mean);
		const double rms = rmsStat.get();

		// Testing, testing 1,2,3...

		const double fieldElipHeight = 1;

		const auto Ipres = magnetization_J_theor_ellipsoid(e, I0, K);
		cout << "Ellipsoid presice I = " << Ipres << " | rel_err = " << (Ipres-mean).eqNorm()/Ipres.eqNorm()   << " | demag_rel_err = " << (Ipres-I0).eqNorm()/Ipres.eqNorm() << endl;
		cout << "Jmean= " << mean << " | pres_err= " << (Ipres-mean)/Ipres << " | rms_err= " << rms << endl << endl;

		for(int layer = 0; layer < layersN; ++layer) {
			cout << layer << ": "  << "Jmean= " << meanStatLayers[layer].get() << " | rms_err= " << rmsStatLayers[layer].get() << endl;
		}

		const Point p0{0, 0, e.Req + fieldElipHeight};
		// cout << "Sphere presice Hsnd_z = " << field_sphere_H_in_Hz(e.Req, Hprime.z, K, p0) << endl;
		// const Point Hprec = field_sphere_H(e.Req, Ipres, p0) / (4.*M_PI);
		// cout << "Sphere presice H = " << Hprec << endl;
		// const Point noDemagHprec = field_sphere_H(e.Req, I0, p0) / (4.*M_PI);
		// cout << "Sphere presice H (no demag) = " << noDemagHprec  << " | rel_err = " << (Hprec-noDemagHprec).eqNorm()/Hprec.eqNorm() << endl;


		// {
		// 	Dat3D<Point> dd;
		// 	dd.es.resize(hsi.size());
		// 	std::transform(hsi.cbegin(), hsi.cend(), dd.es.begin(), [](const HexahedronWid &h) -> Dat3D<Point>::Element {
		// 		const auto p = h.massCenter();
		// 		return {{p.x, p.y, p.z}, h.dens};
		// 	});
		// 	dd.write("elip_J_int.dat");
		// }

		const auto solver = createCudaSolver(hsi, false);


		const auto &fOnDat = [&](Dat3D<Point> &res) {
			for (auto &i : res) 
				i.val = -solver->solve({ i.p.x, i.p.y, i.p.z }) / (4 * M_PI);
		};

		// {
		// 	Dat3D<Point> dd;
		// 	const double t = 0.001;
		// 	for (double x = -15-t; x < 15; x += 0.2)
		// 		for (double y = -15-t; y < 15; y += 0.2)
		// 			dd.es.push_back({{ x, y, e.Req/2 + t}});
		// 	fOnDat(dd);
		// 	dd.write("elip_f_in_out_no_J3.dat");
		// 	return;
		// }

		// {
		// 	Dat3D<Point> dd;
		// 	const double t = 0.007;
		// 	for (double x = -15-t; x < 15; x += 0.2)
		// 		for (double y = -15-t; y < 15; y += 0.2)
		// 			dd.es.push_back({{ x, y, e.Req/2 + t}});

		// 	fOnDat(dd);
		// 	dd.write("cube_f_in_out.dat");

		// 	for(auto &el: dd.es) el.p.z = e.Req + 1;
		// 	fOnDat(dd);
		// 	dd.write("cube_f_out.dat");
		// }

		// {
		// 	Dat3D<Point> dd;
		// 	for(int part = 0; part < 4; ++part)
		// 		for(int layer = 0; layer < layersN; ++layer) 
		// 			for(int li = 0; li < nl; ++li) {
		// 				const int Bi = 0;
		// 				const int idx = ((part * layersN) + layer) * inSz + (li * nB) + Bi;
		// 				const auto c = hsi[idx].massCenter();
		// 				dd.es.push_back({{c.x, c.y, c.z}, hsi[idx].dens});
		// 			}
		// 	dd.write("elip_J_in.dat");
		// }
		
		{
			Dat3D<Point> dd;
			for(int i = 0; i < hsi.size(); ++i) {
				const auto c = hsi[i].massCenter();
				dd.es.push_back({{c.x, c.y, c.z}, hsi[i].dens});
			}
			dd.write("ball_J_all_in.dat");
		}
	}

	void run(int argc, char *argv[]) {
/*
		if(isRoot()) {
			hexTest();
		}
*/
		InputParser inp(argc, argv);
		
		double H = 0.25;
		if (inp.exists("H")) inp["H"] >> H;

		//well Model
		const auto wellModelGenerator = [&](vector<HexahedronWid> &hsi, vector<double> &K){
			Volume v{ { -2, 10, 120 },{ -2, 10, 120 },{ -8, 0, 40 } };
			if(inp.exists("xn")) inp["xn"] >> v.x.n;
			if(inp.exists("yn")) inp["yn"] >> v.y.n;
			if(inp.exists("zn")) inp["zn"] >> v.z.n;

			Cylinder well = { { { 4.1, 4.1 }, 0.25 }, 100 }; //x_cener, y_center, r, h
			if(inp.exists("rwell")) inp["rwell"] >> well.r;
			double Kouter = 0.02;
			const double Kinner0 = 0;
			const Point Hprime = { 14, 14, 35 }; //~40A/m

			//read K inner
			vector<double> Kinner(v.z.n, Kinner0);
			if (inp.exists("kin")) {
				string KinnerFname;
				inp["kin"] >> KinnerFname;
				Dat2D<> Kin(KinnerFname);
				v.z = { -Kin.xMax(), 0, (int)Kin.size() };
				Kinner.resize(v.z.n);
				double meanK = 0;
				for (int i = 0; i < Kin.size(); ++i)
					meanK += Kinner[Kin.size() - i - 1] = Kin[i].p.y / 1e5;
				meanK /= Kin.size();
				Kouter = meanK;
				cout << "Kouter: " << Kouter << endl;
			}
		};

		//const limits fieldDim = { 0, 8.1, 60 };
		//calcDemag(wellModelGenerator, fieldDim, fieldDim, H, "well");

		// Volume cubeVolume{ { -20, 20, 160 },{ -20, 20, 160 },{ 0, -4, 16 } }; //cube "A"
		// Volume cubeVolume{ { -2, 2, 50 },{ -1, 1, 25 },{ -2, 0, 25 } }; //cube "B" 2:1 - 4x2x2

		const int ratio = 2; //x2
		Volume cubeVolume{ { -2, 2, 25*ratio },{ -2./double(ratio), 2./double(ratio), 25 },{ -4./double(ratio), 0, 25 } }; //cube "B"

		if(inp.exists("xn")) inp["xn"] >> cubeVolume.x.n;
		if(inp.exists("yn")) inp["yn"] >> cubeVolume.y.n;
		if(inp.exists("zn")) inp["zn"] >> cubeVolume.z.n;

		const auto cubeModelGenerator = [&](vector<HexahedronWid> &hsi, vector<double> &Kmodel, vector<Point> &I0) {
			const double K = 0.2;
			//const Point Hprime = { 14, 14, 35 }; //~40A/m
			const Point Hprime = { 0, 0, 50 }; //~40A/m

			cubeGen(cubeVolume, Hprime*K, hsi);
			Kmodel.assign(hsi.size(), K);
			I0.assign(hsi.size(), Hprime*K);
		};

		double DPR = std::max({ cubeVolume.x.d(), cubeVolume.y.d(), cubeVolume.z.d() }) * 4;
		if(inp.exists("DPR")) inp["DPR"] >> DPR;
		if(isRoot()) cout << "DPR: " << DPR << (DPR < 0? " (disabled)" : "") << endl;
		const limits fieldDim = { -25 + 0.00001, 25 + 0.00001, 40 };
		Stopwatch tmr;
		tmr.start();
		calcDemag<HexahedronWid>(cubeModelGenerator, cubeVolume, fieldDim, fieldDim, H, "cube", DPR);
		if(isRoot()) {
			cout << "Total time: " << tmr.stop() << "sec." << endl;
		}
	}

private:
	template<typename T>
	struct CG {
		const vector<T> &b;
		const std::function<vector<T>(vector<T>&)> Op;

		bool ready = false;
		vector<T> r;
		vector<T> z;
		vector<T> x;

		CG(const vector<T> &b, const vector<T> &x0, const std::function<vector<T>(vector<T>&)> &Op): b(b), x(x0), Op(Op) {}

		// x = ax + by
		static void ax_plus_by(const double a, vector<T>& x, const double b, const vector<T>& y) {
			std::transform(x.cbegin(), x.cend(), y.cbegin(), x.begin(), [&a, &b](const auto& x, const auto &y){ return x*a + y*b; });
		}
		static double dot(const vector<T>& a, const vector<T>& b) {
			return std::transform_reduce(a.cbegin(), a.cend(), b.cbegin(), 0., 
				[](const auto& a, const auto& b){ return a + b; }, 
				[](const auto& a, const auto& b){ return a ^ b; }
			);
		}

		double getError() const {
			return sqrt(dot(r, r) / dot(b, b));
		}

		void prepare() {
			r = Op(x);	//r0 = Ax
			Assert(x.size() == r.size());
			ax_plus_by(-1, r, 1, b); //r0 = b - Ax
			z = r;
			ready = true;
		}

		void nextIter() {
			if(!ready) throw std::runtime_error("CG: call prepare() before iter()");
			const auto Az = Op(z);
			const double r_dot = dot(r, r);
			const double alpha = r_dot / dot(Az, z);
			ax_plus_by(1, x, alpha, z); 		// x = x_prv + alpha*z_prv
			ax_plus_by(1, r, -alpha, Az);		// r = r_prv - alpha*Az_prv
			const double beta = dot(r, r) / r_dot;
			ax_plus_by(beta, z, 1, r);  		// z = r + beta*z_prv
		}

	};

	template<class ClosedShape>
	vector<ClosedShape> demagCG(
		const std::function<void(vector<ClosedShape>&, vector<double>&, vector<Point>&)> &modelGenerator,
		const std::function<std::unique_ptr<gFieldSolver>(const vector<ClosedShape>&, const bool)> createCudaSolver
	) {
		if(isRoot()) cout << "Generating model..." << endl;
		vector<ClosedShape> hsi;
		vector<double> K;
		vector<Point> J0;
		modelGenerator(hsi, K, J0);
		if(isRoot()) cout << "Model size: " << hsi.size() << endl;
		// return hsi;

		//re-calculate J (magnetization) in every ClosedShape (simple iteration)
		const auto &fJn = [&createCudaSolver, &hsi, &K, &J0, this]() -> vector<Point> {
			Bcast(hsi);
			std::unique_ptr<gFieldSolver> solver = createCudaSolver(hsi, false);
			vector<Point> fieldPoints(hsi.size());
			std::transform(hsi.cbegin(), hsi.cend(), fieldPoints.begin(), [](const ClosedShape &h) {return h.massCenter();});
			vector<Point> field(hsi.size());
			MPIpool<Point, Point> pool(*this, fieldPoints, field, 1024);
			// pool.logging = true;
			int cnt = 0;
			if (!isRoot()) {
				while (1) {
					const vector<Point> task = pool.getTask();
					if (!task.size()) break;
					if(pool.logging) cout << "Task accepted " << cnt++ << " size: " << task.size() << endl;
					vector<Point> result(task.size());
					for (int i = 0; i < task.size(); ++i)
						result[i] = -solver->solve(task[i]);
					pool.submit(result);
				}
			} else {
				cout << "result gather ok" << endl;

				for (int i = 0; i < field.size(); ++i) 
					field[i] = (field[i] / (4.*M_PI)) * K[i] + J0[i];

				Point mean = 0;
				for (int i = 0; i < field.size(); ++i) 
					mean += field[i];
				mean = mean / field.size();
				cout << "Mean Jn+1= " << mean << endl;

				return field;
			}

			return {};
		};

		// Re-calculate J (magnetization) in every ClosedShape (CG)
		// Should be valid: 'hsi' - all nodes; x - root; result - root
		const auto OpCGt = [&createCudaSolver, &hsi, &K, &J0, this](vector<Point> x = {}, const bool transpose = false) -> vector<Point> {
			Bcast(x);
			auto shapes{ hsi };
			for (int i = 0; i < shapes.size(); ++i) shapes[i].dens = x[i];
			const std::unique_ptr<gFieldSolver> solver = createCudaSolver(shapes, transpose);
			vector<Point> fieldPoints(shapes.size());
			std::transform(shapes.cbegin(), shapes.cend(), fieldPoints.begin(), [](const ClosedShape &h) {return h.massCenter();});
			vector<Point> res(shapes.size());
			MPIpool<Point, Point> pool(*this, fieldPoints, res, 1024);
			// pool.logging = true;
			int cnt = 0;
			if (!isRoot()) {
				while (1) {
					const vector<Point> task = pool.getTask();
					if (!task.size()) break;
					if(pool.logging) cout << "Task accepted " << cnt++ << " size: " << task.size() << endl;
					vector<Point> result(task.size());
					for (int i = 0; i < task.size(); ++i)
						result[i] = -solver->solve(task[i]);
					pool.submit(result);
				}
			} else {
				if(pool.logging) cout << "result gather ok" << endl;
				for (int i = 0; i < res.size(); ++i) 
					res[i] = x[i] - (res[i] / (4.*M_PI)) * K[i];
				return res;
			}
			return {};
		};

		vector<Point> x0(hsi.size());
		std::transform(hsi.cbegin(), hsi.cend(), x0.begin(), [](const ClosedShape &h) {return h.dens;});
		CG<Point> cg{J0, x0, OpCGt};

		//calculate residual ||res[] - hsi.dens[]|| and copy res[] to hsi.dens[]
		const auto residualEqAndCopy = [](const vector<Point> &res, vector<ClosedShape> &hsi) {
			double sum = 0;
			for (int i = 0; i < res.size(); ++i) {
				const Point v = res[i] - hsi[i].dens;
				sum += v^v;
				hsi[i].dens = res[i];
			}
			return sqrt(sum);
		};

		
		if(isRoot()) cout << "Demag Solving..." << endl;

		// Not Root will do calculations here
		if (!isRoot()) {
			while (1) {
				bool cont = false;
				Bcast(cont);
				if (!cont) break;
				//fJn();
				OpCGt();
			}
			cout << "Done." << endl;
			return {};
		}
		
		const double eps = 1e-4;
		const double minRelativeErrChange = 0.05; //5%
		const int maxIter = 10;
		double err = 1, prvErr = err*100;

		cout << "CG: preparing first iter" << endl;
		{
			bool cont = true;
			Bcast(cont);
			cg.prepare();
		}

		for (int it = 0; it < maxIter && err > eps /* && (fabs(err-prvErr) / err) > minRelativeErrChange && err < prvErr */ ; ++it) {
			cout << "Iter: " << it << endl;
			bool cont = true;
			Bcast(cont);

			// const vector<Point> In = fJn();
			// cout << endl;
			cg.nextIter();
			
			prvErr = err;
			// err = residualEqAndCopy(In, hsi) / [&In]() {
			// 	double sum = 0;
			// 	for (auto& i: In)
			// 		sum += i ^ i;
			// 	return sqrt(sum);
			// }();
			err = cg.getError();
			
			const auto ediff = err - prvErr;
			cout << "Err: " << err << "(" << (ediff>0?"+":"") << ediff << ")" << " at iter: " << it << endl;
		}

		bool cont = false;
		Bcast(cont);

		for (int i = 0; i < hsi.size(); ++i) hsi[i].dens = cg.x[i];

		return hsi;
	}

	template<class ClosedShape>
	void calcDemag(
		const std::function<void(vector<ClosedShape>&, vector<double>&, vector<Point>&)> &modelGenerator,
		const Volume& v,
		const limits& fieldDimX, const limits& fieldDimY, const double H,
		const string& filePrefix,
		const double DPR
	) {
		const auto createCudaSolver = [&DPR](const vector<ClosedShape>& hsi, const bool transpose) {
			return	DPR < 0
				? gFieldSolver::getCUDAsolver(&*hsi.cbegin(), &*hsi.cend(), transpose)
				: gFieldSolver::getCUDAreplacingSolver(&*hsi.cbegin(), &*hsi.cend(), DPR, hsi.size())
			;
		};

		vector<ClosedShape> hsi = demagCG<ClosedShape>(modelGenerator, createCudaSolver);
		if(!isRoot()) return;

		const auto &fOnDat = [&](Dat3D<Point> &res) {
			std::unique_ptr<gFieldSolver> solver = createCudaSolver(hsi, false);
			for (auto &i : res)
				i.val = -solver->solve({ i.p.x, i.p.y, i.p.z }) / (4 * M_PI);
		};
/*
		Dat3D<Point> dd("cubeFieldSZ_IN.dat");
		fOnDat(dd);
		dd.write("cubeFieldSZ_IN_DESC.dat");
		return;
*/
/*
		const auto &fInWell = [&v, &well, &fOnDat](const string fname) {
			Dat3D<Point> dat;
			for (int i = 1; i < v.z.n; ++i)
				dat.es.push_back({{ well.center.x, well.center.y, v.z.atWh(i) - v.z.d() / 2.}});
			fOnDat(dat);
			dat.write(fname);
		};
*/

		//calculate field on a plane of heigh 'H' above model and save to 'fname'_x.dat, 'fname'_y.dat, 'fname'_z.dat
		const auto &fOn = [&createCudaSolver, &hsi, &H, &fieldDimX, &fieldDimY](const string fname) {
			Field x{ fieldDimX, fieldDimY }, y{ x }, z{ x };
			std::unique_ptr<gFieldSolver> solver = createCudaSolver(hsi, false);
			int pp = -1;
			for (int i = 0; i < x.y.n; ++i) {
				for (int j = 0; j < x.x.n; ++j) {
					const Point p0{ x.x.atWh(j), x.y.atWh(i), H };
					const Point res = -solver->solve(p0) / (4*M_PI);
					x.data[i*x.x.n + j] = res.x;
					y.data[i*x.x.n + j] = res.y;
					z.data[i*x.x.n + j] = res.z;
				}
				const int p = (100 * (i + 1)) / x.y.n;
				if (p != pp) {
					pp = p;
					cout << p << "%" << endl;
				}
			}
			cout << endl;
			x.toGrid().Write(fname + "_x.grd");
			y.toGrid().Write(fname + "_y.grd");
			z.toGrid().Write(fname + "_z.grd");
		};


		const auto expJ = [&hsi, &v](const string& fname, const int layerIdx = 0) {	//dump J (magnetization)
			Field x{ v.x, v.y }, y{ x }, z{ x };
			for (int i = 0; i < v.y.n; ++i) {
				for (int j = 0; j < v.x.n; ++j) {
					const auto& h = hsi[layerIdx * (v.y.n*v.x.n) + i * v.x.n + j];
					const Point p0 = h.massCenter();
					x.data[i*x.x.n + j] = h.dens.x;
					y.data[i*x.x.n + j] = h.dens.y;
					z.data[i*x.x.n + j] = h.dens.z;
				}
			}
			x.toGrid().Write(fname + "_x.grd");
			y.toGrid().Write(fname + "_y.grd");
			z.toGrid().Write(fname + "_z.grd");
			//cout << "expJ() done." << endl;
		};
		const auto expJall = [&expJ, &v](const string& fname) {
			for(int zi = 0; zi < v.z.n; ++zi)
				expJ(fname + "_" + std::to_string(zi), zi);	
		};



		const auto expJdiff = [&expJall, &v](const string& filePrefix) {
			expJall("Jn/" + filePrefix + "Jn");

			const auto gridDiff3d = [](const string& base1, const string& base2, const string& baseDest) {
				const auto diff1d = [&](const string& axis) {
					Grid t(base1 + "_" + axis + ".grd");
					auto gdiff = t - Grid(base2 + "_" + axis + ".grd");
					gdiff.Write(baseDest + "_" + axis + ".grd");
					return std::make_tuple(gdiff.mean(), gdiff.sumOfCubes(), t.sumOfCubes());
				};
				const auto x = diff1d("x"), y = diff1d("y"), z = diff1d("z");
				return std::make_tuple( 
					Point{ std::get<0>(x), std::get<0>(y), std::get<0>(z) }, //mean
					std::get<1>(x) + std::get<1>(y) + std::get<1>(z), // sum of cubes for diff
					std::get<2>(x) + std::get<2>(y) + std::get<2>(z) // sum of cubes for total
				);
			};

			Point mean;
			double sumOfCubesDiff = 0;
			double sumOfCubesTotal = 0;
			for(int zi = 0; zi < v.z.n; ++zi) {
				const auto gr = gridDiff3d(
					"Jn/" + filePrefix + "Jn" + "_" + std::to_string(zi), 
					"J0/" + filePrefix + "J0" + "_" + std::to_string(zi), 
					"Jn_J0_diff/" + filePrefix + "Jn_J0_diff" + "_" + std::to_string(zi)
				);
				mean += std::get<0>(gr);
				sumOfCubesDiff += std::get<1>(gr);
				sumOfCubesTotal += std::get<2>(gr);
			}
			return std::make_tuple( mean, sumOfCubesDiff, sumOfCubesTotal );
		};

		const auto JdiffRes = expJdiff(filePrefix);
		cout << "Mean Jn-J0 = " << std::get<0>(JdiffRes) << endl;
		cout << "|Jn-J0|/|Jn| = " << std::sqrt(std::get<1>(JdiffRes)) / std::sqrt(std::get<2>(JdiffRes)) << endl;
		fOn(filePrefix + "Field");
		//fInWell("inWell.dat");


		cout << "Master Done." << endl;
	}
};



//split dencity model to Hexahedrons
template <class VAlloc>
void makeHexs(const double l0, const Ellipsoid &e, limits Nlim, limits Elim, limits Hlim,
	const vector<vector<Point>> &dens, vector<HexahedronWid, VAlloc> &hsi) {

	Nlim = { Nlim.lower - Nlim.dWh() / 2.,  Nlim.upper + Nlim.dWh() / 2., Nlim.n };
	Elim = { Elim.lower - Elim.dWh() / 2.,  Elim.upper + Elim.dWh() / 2., Elim.n };

	cout << "grid E (x): " << Elim << " " << Elim.d() << endl;
	cout << "grid N (y): " << Nlim << " " << Nlim.d() << endl;
	hsi.resize(Hlim.n*Nlim.n*Elim.n);


	const TransverseMercator proj(e.Req*1000., e.f, 1);

	auto GKtoGeo = [&](const double N, double E) -> std::pair<double, double> {
		std::pair<double, double> lb;
		E = xFromGK(E, l0);
		proj.Reverse(toDeg(l0), E*1000., N*1000., lb.second, lb.first);
		lb.first = toRad(lb.first);
		lb.second = toRad(lb.second);
		return lb;
	};

	auto Kr = [&](const int Hi, const int Ni, const int Ei) {
		Hexahedron h({
			e.getPoint(GKtoGeo(Nlim.at(Ni + 1), Elim.at(Ei + 1)), Hlim.at(Hi + 1)),
			e.getPoint(GKtoGeo(Nlim.at(Ni), Elim.at(Ei + 1)), Hlim.at(Hi + 1)),
			e.getPoint(GKtoGeo(Nlim.at(Ni + 1), Elim.at(Ei)), Hlim.at(Hi + 1)),
			e.getPoint(GKtoGeo(Nlim.at(Ni), Elim.at(Ei)), Hlim.at(Hi + 1)),
			e.getPoint(GKtoGeo(Nlim.at(Ni + 1), Elim.at(Ei + 1)), Hlim.at(Hi)),
			e.getPoint(GKtoGeo(Nlim.at(Ni), Elim.at(Ei + 1)), Hlim.at(Hi)),
			e.getPoint(GKtoGeo(Nlim.at(Ni + 1), Elim.at(Ei)), Hlim.at(Hi)),
			e.getPoint(GKtoGeo(Nlim.at(Ni), Elim.at(Ei)), Hlim.at(Hi))
		}, dens[Hi][Ni*Elim.n + Ei]);

		hsi[(Hi*Nlim.n + Ni)*Elim.n + Ei] = HexahedronWid(h);
	};

#pragma omp parallel for
	for (size_t Hi = 0; Hi < Hlim.n; ++Hi)
		for (size_t Ni = 0; Ni < Nlim.n; ++Ni)
			for (size_t Ei = 0; Ei < Elim.n; ++Ei)
				Kr(Hi, Ni, Ei);

}

/* -- TRANS SOLVER --
void transFieldNode(Ellipsoid e, double l0, const vector<Dat3D<>::Element> psGK,
	const vector<HexahedronWid>::const_iterator &hsBegin, const vector<HexahedronWid>::const_iterator &hsEnd, const vector<double>::iterator &resBegin,
	const double dotPotentialRad) {
	vector<FieldPoint> fps(psGK.size());
	const TransverseMercator proj(e.Req*1000., e.f, 1);
	std::transform(psGK.cbegin(), psGK.cend(), fps.begin(), [&](const Dat3D<>::Element &el) {
		double l, B;
		double x = xFromGK(el.p.x, l0);
		proj.Reverse(toDeg(l0), x*1000., el.p.y*1000., B, l);
		l = toRad(l);
		B = toRad(B);
		return FieldPoint{ e.getPoint(B, l, el.p.z), e.getNormal(B, l), el.val };
	});

	//solve
	gFieldInvSolver::getCUDAtransSolver(fps, dotPotentialRad)->solve(&*hsBegin, &*hsEnd, resBegin);
	std::transform(resBegin, resBegin + (hsEnd - hsBegin), resBegin, [](auto &v) {return G_CONST*v; });
}
*/

void calcFieldNode(const Ellipsoid &e, const double l0, std::unique_ptr<gFieldSolver> &solver,
	const vector<Dat3D<>::Point> &fp, vector<Point> &result) {
	Assert(fp.size() == result.size());
	const TransverseMercator proj(e.Req*1000., e.f, 1);

	auto Kr = [&](const Dat3D<>::Point &p, Point& res) {
		double l, B;
		double x = xFromGK(p.x, l0);
		proj.Reverse(toDeg(l0), x*1000., p.y*1000., B, l);
		l = toRad(l);
		B = toRad(B);
		const Point p0 = e.getPoint(B, l, p.z);
		//const Point n0 = e.getNormal(B, l);
		res += solver->solve(p0) * G_CONST;
		//res.val += (const Point&)res ^ n0;
	};
	for (size_t i = 0; i < fp.size(); ++i)
		Kr(fp[i], result[i]);
}

class ClusterSolver : public MPIwrapper {
public:
	const int maxGPUmemMB = 5000;
	int triBufferSize = 0;
	using Qiter = gElementsShared::const_iterator;

	SharedMemBase<gElementsShared> *sharedMem = 0;

	ClusterSolver() : MPIwrapper() {
		if (gridSize < 2) throw std::runtime_error("You must run at least 2 MPI processes.");
		if (root != 0) throw std::runtime_error("Root process must have rank = 0.");
		const auto lid = localId();
		const int devId = !isRoot() && std::get<1>(lid) ? std::get<0>(lid) - 1 : std::get<0>(lid);
		cuSolver::setDevice(devId);
	}

	gElementsShared& initShared(const size_t size) {
		if (isLocalRoot()) {
			sharedMem = new SharedMemMaster<gElementsShared>(size);
			Barrier();
		}
		else {
			Barrier();
			sharedMem = new SharedMemSlave<gElementsShared>();
		}
		return *sharedMem->data;
	}

	~ClusterSolver() {
		if (sharedMem) {
			if (!isLocalRoot()) {
				delete sharedMem;
				Barrier();
			}
			else {
				Barrier();
				delete sharedMem;
			}
		}
	}

	template<typename T>
	size_t getPartSize(const size_t sz) {
		const size_t partSize = size_t(maxGPUmemMB)*size_t(MB) / sizeof(T);
		if (sz%partSize == 0) return partSize;
		return sz / (sz / partSize + 1);
	}

	void calcField(Ellipsoid e, double l0, gElementsShared &qs, Dat3D<Point> &dat, double dotPotentialRad) {
		//cout << "Dot potential replace radius: " << dotPotentialRad << endl;
		if (isRoot())
			cout << "Computing nodes: " << gridSize - 1 << endl;

		size_t qSize = qs.size();
		if (myId == 1) send(qSize, 0);
		else if (myId == 0) recv(qSize, 1);

		size_t partSize = getPartSize<gElements::base>(qSize);

		const size_t parts = ((partSize + qSize - 1) / partSize);
		for (size_t i = 0; i < qSize; i += partSize) {
			const size_t part = (i / partSize);
			const Qiter qbegin = qs.begin() + i;
			const Qiter qend = qs.begin() + std::min(i + partSize, qs.size());

			cout << ":::part " << part + 1 << " out of " << parts << "::: size: " << qend - qbegin << endl;
			calcWithPool(e, l0, qbegin, qend, dat, dotPotentialRad);
		}
	}
/* -- TRANS SOLVER --
	void calcTrans(Ellipsoid e, double l0, vector<gElements::base> &qs, const Dat3D<> &dat, vector<double> &result, const double dotPotentialRad) {
		if (isRoot())
			cout << "Computing nodes: " << gridSize - 1 << endl;

		MPIpool<gElements::base, double> pool(*this, qs, result, 50000);
		int cnt = 0;
		if (!isRoot()) {
			while (1) {
				const vector<gElements::base> task = pool.getTask();
				if (!task.size()) break;
				cout << "Task accepted " << ++cnt << " size: " << task.size() << endl;
				vector<double> result(task.size());
				transFieldNode(e, l0, dat.es, task.cbegin(), task.cend(), result.begin(), dotPotentialRad);
				pool.submit(result);
			}
		}
		else cout << "result gather ok" << endl;
	}
*/
private:
	void calcWithPool(Ellipsoid e, double l0, const Qiter &qbegin, const Qiter &qend, Dat3D<Point> &dat, double dotPotentialRad) {
		vector<Dat3D<>::Point> fp = dat.getPoints();
		vector<Point> result;
		//blocking process with rank 0 until the result's been gathered
		MPIpool<Dat3D<>::Point, Point> pool(*this, fp, result, 1024);
		int cnt = 0;
		if (!isRoot()) {
			std::unique_ptr<gFieldSolver> solver = gFieldSolver::getCUDAreplacingSolver(&*qbegin, &*qend, dotPotentialRad, triBufferSize);
			while (1) {
				vector<Dat3D<>::Point> task = pool.getTask();
				if (!task.size()) break;
				cout << "Task accepted " << ++cnt << " size: " << task.size()+1 << endl;
				vector<Point> result(task.size());
				calcFieldNode(e, l0, solver, task, result);
				pool.submit(result);
			}
		}
		else {
			dat.add(result);
			cout << "result gather ok" << endl;
		}
	}
};

int main(int argc, char *argv[]) {
	bool isRoot = true;
	try {
		WellDemagCluster().runEllipsoid(argc, argv);
		//wellTest(argc, argv);
		return 0;

		ClusterSolver cs;
		isRoot = cs.isRoot();
		Ellipsoid Earth(R_EQ, R_PL);

		//gElements qs;
		cout << "GPUs: " << cuSolver::getGPUnum() << endl;
		Input inp(argc, argv);

		const size_t qAm = getHexAm(inp.Nlim.n, inp.Elim.n, inp.Hlim.n);
		cout << (qAm * sizeof(gElements::base) + MB - 1) / MB << "MB required. ";
		cout << "Elements: " << qAm << endl;
		gElementsShared &qss = cs.initShared(qAm);//for direct solver. Should not actually allocate memory untill resize().
		gElements qs; //for transpose solver
		cout << "Allocated." << endl;

		const auto& prep = [&](const auto& f) {
			cout << "Preprocessing started" << endl;
			Stopwatch tmr;
			tmr.start();
			f();
			const double time = tmr.stop();
			cout << "Preprocessing finished: " << time << "sec." << endl;
		};
		
		//prepocessing for direct solver
		if (cs.isLocalRoot() && !inp.transSolver)
			prep([&]() {makeHexs(inp.l0, Earth, inp.Nlim, inp.Elim, inp.Hlim, inp.dens, qss); });
		//prepocessing for transpose solver
		if (cs.isRoot() && inp.transSolver)
			prep([&]() {makeHexs(inp.l0, Earth, inp.Nlim, inp.Elim, inp.Hlim, inp.dens, qs); });
		cs.Barrier();
		const size_t qSize = inp.transSolver ? qs.size() : qss.size();
		cout << "Real size: " << (qSize * sizeof(gElements::base) + MB - 1) / MB << "MB" << endl;

		if (inp.transSolver) {
			throw std::runtime_error("Transposed Solver is not implemented");
/* -- TRANS SOLVER --			
			vector<double> res;
			if (cs.isRoot()) {
				cout << "TRANSPOSE SOLVER. Your *.grd files will be overwritten!" << endl;
				res.resize(qs.size());
			}
			Stopwatch tmr;
			cs.calcTrans(Earth, inp.l0, qs, inp.dat, res, inp.dotPotentialRad);
			if (cs.isRoot()) {
				cout << "time: " << tmr.stop() << endl;

				const size_t lsize = inp.Nlim.n*inp.Elim.n;
				for (size_t i = 0; i < inp.Hlim.n; ++i) {
					Grid g(inp.fnames[i]);
					g.data.assign(res.cbegin() + i*lsize, res.cbegin() + (i + 1)*lsize);
					g.Write();
				}

				cout << "Done" << endl;
			}
			return 0;
*/
		}
		if (cs.isRoot()) {
			cout << "Clearing output file. Do NOT stop me!" << endl;
			inp.dat.set(Point(0));
			inp.dat.write();
			cout << "Clearing done. You may ctrl+c now." << endl;
		}

		cout << "Computing started" << endl;
		Stopwatch tmr;
		tmr.start();
		//set approximate size of buffer
		cs.triBufferSize = triBufferSize(inp.Nlim, inp.Elim, inp.Hlim, inp.dotPotentialRad);
		//do the job
		cs.calcField(Earth, inp.l0, qss, inp.dat, inp.dotPotentialRad);
		const double time = tmr.stop();
		if (cs.isRoot()) {
			cout << "Computing finished: " << time << "sec." << endl << endl;
			inp.dat.write();
		}
	}
	catch (std::exception &ex) {
		if(isRoot) std::cerr << "Global exception: " << ex.what() << endl;
		return 1;
	}
	cout << "Done" << endl;
	return 0;
}


/*
template <class VAlloc>
void ellipsoidGenPyramid(const Ellipsoid &e, const int nl, const int nB, const Point J, vector<Pyramid, VAlloc> &hsi) {
	limits ll{0, M_PI_2, nl}, lB{0, M_PI_2, nB};
	hsi.resize(ll.n * lB.n);

	for(int li = 0; li < ll.n; ++li)
		for(int Bi = 0; Bi < lB.n; ++Bi)
			hsi[li * lB.n + Bi] = Pyramid{
				Point{}, 	//apex
				{			//base
					e.getPoint(lB.at(Bi), ll.at(li)),
					e.getPoint(lB.at(Bi), ll.at(li+1)),
					e.getPoint(lB.at(Bi+1), ll.at(li+1)),
					e.getPoint(lB.at(Bi+1), ll.at(li))
				},
				J
			};

	const auto dbl = [&](const std::function<Pyramid(Pyramid)> &f) {
		vector<Pyramid> tmp(hsi.begin(), hsi.end());
		std::transform(tmp.begin(), tmp.end(), tmp.begin(), f);
		hsi.insert(hsi.end(), tmp.begin(), tmp.end());
	};
	dbl([](Pyramid p){ 
		p.apex.x *= -1;
		p.base[0].x *= -1;
		p.base[1].x *= -1;
		p.base[2].x *= -1;
		p.base[3].x *= -1;
		return p;
	});
	dbl([](Pyramid p){ 
		p.apex.y *= -1;
		p.base[0].y *= -1;
		p.base[1].y *= -1;
		p.base[2].y *= -1;
		p.base[3].y *= -1;
		return p;
	});
	dbl([](Pyramid p){ 
		p.apex.z *= -1;
		p.base[0].z *= -1;
		p.base[1].z *= -1;
		p.base[2].z *= -1;
		p.base[3].z *= -1;
		return p;
	});
}
*/