
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
	Quadrangle operator+(const Point& p) const {
		return { p1 + p, p2 + p, p3 + p, p4 + p };
	}
};

template <class VAlloc>
void wellGen(const Volume &v, const Cylinder &well, const Point Joutter, const double Koutter, const Point Jinner, const double Kinner, 
		vector<HexahedronWid, VAlloc> &hsi, vector<double> &K) {
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
		if (q.isIn(well)) {
			q.J = Jinner;
			q.k = Kinner;
			q.inner = true;
		}
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
				const int ind = (zi*yLim.n + yi)*xLim.n + xi;
				hsi[ind] = Hexahedron{
					cur + Point{0, 0, v.z.at(zi)},
					cur + Point{0, 0, v.z.at(zi + 1)},
					cur.J};
				K[ind] = cur.k;
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

	/*
	//add 16 external quadrangles
	auto eq1 = makeOuterHex({ xLim.lower, yLim.lower - l, v.z.upper }, { xLim.upper + l, yLim.lower, v.z.lower }).splitTo4();
	auto eq2 = makeOuterHex({ xLim.upper, yLim.lower, v.z.upper }, { xLim.upper + l, yLim.upper + l, v.z.lower }).splitTo4();
	auto eq3 = makeOuterHex({ xLim.lower - l, yLim.upper, v.z.upper }, { xLim.upper, yLim.upper + l, v.z.lower }).splitTo4();
	auto eq4 = makeOuterHex({ xLim.lower - l, yLim.lower - l, v.z.upper }, { xLim.lower, yLim.upper, v.z.lower }).splitTo4();
	hsi.insert(hsi.end(), eq1.cbegin(), eq1.cend());
	hsi.insert(hsi.end(), eq2.cbegin(), eq2.cend());
	hsi.insert(hsi.end(), eq3.cbegin(), eq3.cend());
	hsi.insert(hsi.end(), eq4.cbegin(), eq4.cend());
	for(int i = 0; i < 16; ++i)
		K.push_back(Koutter);
	*/
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
				const Triangle t = h.getTriSafeZ(i);
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
				for (int i = 0; i < x.n; +i) {
					const Point p0{ x.atWh(i), y.atWh(j), z.atWh(k) };
					dat.es.push_back({ { p0.x, p0.y, p0.z }, Point() });
				}
	}
};

Point intHexTr(const Point &p0, const HexahedronWid &h);

void hexTest(int argc, char *argv[]) {
	Hexahedron h({
		{ 20, 20, 0 },{ -20, 20, 0 },{ 20, -20, 0 },{ -20, -20, 0 },
		{ 20, 20, -4 },{ -20, 20, -4 },{ 20, -20, -4 },{ -20, -20, -4 }
	}, Point{ 14, 14, 35 }*0.2);
	auto hw = HexahedronWid(h);
	double H = 0.25;


	cout << "Solving..." << endl;
	Dat2D<Point> dat;
	auto l = limits{ -20 + 0.00001, 20 + 0.00001, 40 };
	for (double i = 0; i < l.n; ++i) {
		for (int j = 0; j < l.n; ++j) {
			const Point p0{ l.atWh(j), l.atWh(i), H };
			const Point res = intHexTr(p0, hw) / (4 * M_PI);
			dat.es.push_back({ { p0.x, p0.y }, res });
		}
	}

	dat.write("cubeFieldSZ.dat"s);
	cout << "Done." << endl;
}

void wellTest(int argc, char *argv[]) {
	const string oFname = "wellField.dat";
	const VolumeMod v = Volume{ { 0, 4, 40 },{ 0, 4, 40 },{ -4, 0, 40 } };
	const Cylinder well = { { { 2, 2 }, 0.25 }, 100 }; //x_cener, y_center, r, h
	const double Kouter = 0.02;
	const double Kinner = 0;
	const Point Hprime = { 14, 14, 35 }; //~40A/m
	const double H = 0.25;

	cout << "Generating model..." << endl;
	const Point J0outer = Hprime*Kouter;
	const Point J0inner = Hprime*Kinner;
	vector<HexahedronWid> hsi;
	vector<double> K;
	wellGen(v, well, J0outer, Kouter, J0inner, Kinner, hsi, K);
	vector<Point> J0(hsi.size());
	transform(hsi.begin(), hsi.end(), J0.begin(), [](const auto& v) {return v.dens; });

	cout << "Solving..." << endl;
	Dat2D<Point> dat;
	unique_ptr<gFieldSolver> solver = gFieldSolver::getCUDAsolver(&*hsi.cbegin(), &*hsi.cend());
	for (int i = 0; i < v.y.n; ++i) {
		for (int j = 0; j < v.x.n; ++j) {
			const Point p0{ v.x.atWh(j), v.y.atWh(i), H };
			const Point res = solver->solve(p0) / (4 * M_PI);
			dat.es.push_back({ { p0.x, p0.y }, res });
		}
		cout << "\r" << (100*(i + 1))/ v.y.n << "%";
	}
	cout << endl;

	dat.write("wellField.dat"s);
	cout << "Done." << endl;
}

class WellDemagCluster : public MPI {
public:
	const int maxGPUmemMB = 5000;
	int triBufferSize = 0;
	using Qiter = gElementsShared::const_iterator;

	SharedMemBase<gElementsShared> *sharedMem = 0;

	WellDemagCluster() : MPI() {
		if (gridSize < 2) throw runtime_error("You must run at least 2 MPI processes.");
		if (root != 0) throw runtime_error("Root process must have rank = 0.");
		const auto lid = localId();
		const int devId = !isRoot() && get<1>(lid) ? get<0>(lid) - 1 : get<0>(lid);
		cuSolver::setDevice(devId);
	}

	void run(int argc, char *argv[]) {
		InputParser inp(argc, argv);

		const string oFname = "wellField.dat";
		VolumeMod v = Volume{ { -2, 10, 120 },{ -2, 10, 120 },{ -8, 0, 20 } };
		inp["xn"] >> v.x.n;
		inp["yn"] >> v.y.n;
		inp["zn"] >> v.z.n;

		const limits fieldDimX = { 0, 8, 120 }, fieldDimY = { 0, 8, 120 };
		const Cylinder well = { { { 4, 4 }, 0.25 }, 100 }; //x_cener, y_center, r, h
		const double Kouter = 0.02;
		const double Kinner = 0;
		const Point Hprime = { 14, 14, 35 }; //~40A/m
		const double H = 0.25;
		
		cout << "Generating model..." << endl;
		const Point J0outer = Hprime*Kouter;
		const Point J0inner = Hprime*Kinner;
		vector<HexahedronWid> hsi;
		vector<double> K;
		wellGen(v, well, J0outer, Kouter, J0inner, Kinner, hsi, K);
		vector<Point> J0(hsi.size());
		transform(hsi.begin(), hsi.end(), J0.begin(), [](const auto& v) {return v.dens; });

		const auto &fOn = [&hsi, &H, &fieldDimX, &fieldDimY](const string fname) {
			Field x{ fieldDimX, fieldDimY }, y{ x }, z{ x };
			unique_ptr<gFieldSolver> solver = gFieldSolver::getCUDAsolver(&*hsi.cbegin(), &*hsi.cend());
			int pp = -1;
			for (int i = 0; i < x.y.n; ++i) {
				for (int j = 0; j < x.x.n; ++j) {
					const Point p0{ x.x.atWh(j), x.y.atWh(i), H };
					const Point res = solver->solve(p0) / (4*M_PI);
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
		const auto &fJn = [&hsi, &K, &J0, this](vector<Point> &res) {
			Bcast(hsi);
			
			unique_ptr<gFieldSolver> solver = gFieldSolver::getCUDAsolver(&*hsi.cbegin(), &*hsi.cend());
			vector<Point> fieldPoints(res.size());
			std::transform(hsi.cbegin(), hsi.cbegin() + fieldPoints.size(), fieldPoints.begin(), [](const HexahedronWid &h) {return h.massCenter();});
			MPIpool<Point, Point> pool(*this, fieldPoints, res, 1024);
			int cnt = 0;
			if (!isRoot()) {
				while (1) {
					const vector<Point> task = pool.getTask();
					if (!task.size()) break;
					cout << "Task accepted " << cnt++ << " size: " << task.size() << endl;
					vector<Point> result(task.size());
					for (int i = 0; i < task.size(); ++i)
						result[i] = solver->solve(task[i]);
					pool.submit(result);
				}
			}
			else {
				cout << "result gather ok" << endl;
				for (int i = 0; i < res.size(); ++i)
					res[i] = - (res[i] + hsi[i].dens * (4.*M_PI / 3.)) * K[i] * (1./(4.*M_PI)) + J0[i];
			}
			
			/*
			for (int i = 0; i < hsi.size(); ++i) {
				const Point p0 = hsi[i].massCenter();
				res[i] = solver->solve(p0) * K[i] + J0[i];
				if (!i || !((100 * (i + 1)) % res.size()) || i == res.size() - 1)
					cout << "\r" << (100 * (i + 1)) / res.size() << "%";
			}
			*/
		};
		const auto& residualAndCopy = [](const vector<Point> &res, vector<HexahedronWid> &hsi) {
			double sum = 0;
			for (int i = 0; i < res.size(); ++i) {
				const Point v = res[i] - hsi[i].dens;
				sum += v^v;
				hsi[i].dens = res[i];
			}
			return sqrt(sum);
		};

		const auto &expJ = [&hsi, &v](const string fname) {	//dump J of upper layer
			Dat2D<Point> dat;
			for (int i = 0; i < v.y.n; ++i) {
				for (int j = 0; j < v.x.n; ++j) {
					const auto& h = hsi[i * v.x.n + j];
					const Point p0 = h.massCenter();
					dat.es.push_back({ { p0.x, p0.y }, h.dens });
				}
			}
			cout << endl;
			dat.write(fname);
		};

		
		cout << "Solving..." << endl;

		if (!isRoot()) {
			while (1) {
				bool cont = false;
				Bcast(cont);
				if (!cont) break;
				fJn(vector<Point>());
			}
			cout << "Done." << endl;
			return;
		}

		//expJ("J0.dat");
		fOn("wellField0"s);
		
		const double eps = 1e-4;
		const int maxIter = 10;
		double err = 1;
		for (int it = 0; it < maxIter && err > eps; ++it) {
			cout << "Iter: " << it << endl;
			bool cont = true;
			Bcast(cont);
			vector<Point> In(hsi.size());
			fJn(In);
			cout << endl;
			err = residualAndCopy(In, hsi) / [&In]() {
				double sum = 0;
				for (auto& i: In)
					sum += i ^ i;
				return sqrt(sum);
			}();
			cout << "Err: " << err << endl;
			//fOn("wellField.dat"s);
		}
		bool cont = false;
		Bcast(cont);

		//expJ("Jn.dat");
		fOn("wellField"s);

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
		WellDemagCluster().run(argc, argv);
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
