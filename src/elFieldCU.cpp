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
#include "Dat.h"
#include "inputParser.h"
#include "Stopwatch.h"
#include "Quadrangles.h"
#include "MPIwrapper.h"
#include "MPIpool.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;

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

std::array<Hexahedron, 4> extendCuboidBoundedModelWith4cuboids(const Volume &innerModelBounds, const double extensionLength, const Point &Joutter) {
	std::array<Hexahedron, 4> hsi;

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
	
	hsi[0] = makeOuterHex({ innerModelBounds.x.lower, -extensionLength, innerModelBounds.z.upper}, { extensionLength, innerModelBounds.y.lower, innerModelBounds.z.lower });
	hsi[1] = makeOuterHex({ innerModelBounds.x.upper, innerModelBounds.y.lower, innerModelBounds.z.upper }, { extensionLength, extensionLength, innerModelBounds.z.lower });
	hsi[2] = makeOuterHex({ -extensionLength, innerModelBounds.y.upper, innerModelBounds.z.upper }, { innerModelBounds.x.upper, extensionLength, innerModelBounds.z.lower });
	hsi[3] = makeOuterHex({ -extensionLength, -extensionLength, innerModelBounds.z.upper }, { innerModelBounds.x.lower, innerModelBounds.y.upper, innerModelBounds.z.lower });

	return hsi;
}

template <class VAlloc>
void wellGen(const Volume &v, const Cylinder &well, const Point Hprime, const double Koutter, const vector<double> Kinner,
		vector<HexahedronWid, VAlloc> &hsi, vector<double> &K) {

	Assert(Kinner.size() == v.z.n);

	const Point Joutter = Hprime*Koutter;
	//const Point J0inner = Hprime*Kinner;
	//make flat mesh
	const limits xLim = { v.x.lower - v.x.dWh() / 2.,  v.x.upper + v.x.dWh() / 2., v.x.n };
	const limits yLim = { v.y.lower - v.y.dWh() / 2.,  v.y.upper + v.y.dWh() / 2., v.y.n };
	vector<Point> mesh((v.x.n+1) * (v.y.n+1));
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

	//add 4 external quadrangles
	const double extensionLength = 1e3;
	const auto extension = extendCuboidBoundedModelWith4cuboids(Volume{ xLim, yLim, v.z }, extensionLength, Joutter);
	hsi.insert(hsi.end(), extension.begin(), extension.end());
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
	int triBufferSize = 0;

	WellDemagCluster(const vector<int> &gpuIdMap = {}) : MPIwrapper() {
		if (gridSize < 2) throw std::runtime_error("You must run at least 2 MPI processes.");
		if (root != 0) throw std::runtime_error("Root process must have rank = 0.");
		const auto lid = localId();
		const int devId = !isRoot() && std::get<1>(lid) ? std::get<0>(lid) - 1 : std::get<0>(lid);
		const int mappedDevId = devId < gpuIdMap.size()? gpuIdMap[devId] : devId;
		cuSolver::setDevice(mappedDevId);
	}

	void runExample(int argc, char *argv[]) {
		InputParser inp(argc, argv);

		double ellipEq = 10, ellipPol = 20;
		int nl = 80, nB = 80, nR = 40;
		double K = 2;
		double HprimeX = 14, HprimeY = 14, HprimeZ = 35; //~40A/m

		inp.parseIfExists("ellipEq", ellipEq);
		inp.parseIfExists("ellipPol", ellipPol);
		inp.parseIfExists("nl", nl);
		inp.parseIfExists("nB", nB);
		inp.parseIfExists("nR", nR);
		inp.parseIfExists("K", K);
		inp.parseIfExists("HprimeX", HprimeX);
		inp.parseIfExists("HprimeY", HprimeY);
		inp.parseIfExists("HprimeZ", HprimeZ);

		const Ellipsoid e(ellipEq, ellipPol);
		const Point Hprime = { HprimeX, HprimeY, HprimeZ };
		const auto I0 = Hprime * K;

		// Cuboid
		double Lx = 1, Ly = 0.4, Lz = 0.4;
		inp.parseIfExists("Lx", Lx);
		inp.parseIfExists("Ly", Ly);
		inp.parseIfExists("Lz", Lz);
		int Nper1 = 100;
		inp.parseIfExists("Nper1", Nper1);
		const int Nx = Lx*Nper1, Ny = Ly*Nper1, Nz = Lz*Nper1;
		
		const auto ellipsoidModelGenerator = [&](vector<HexahedronWid> &hsi, vector<double> &Kmodel, vector<Point> &I0out){
			// const auto Ipres = I0 / (1. + K/3.);	// Use known precise I
			// ellipsoidGen(e, nl, nB, nR, Ipres, hsi);
			ellipsoidGen(e, nl, nB, nR, I0, hsi);
			Kmodel.assign(hsi.size(), K);
			I0out.assign(hsi.size(), I0);
		};
		const auto cubeModelGenerator = [&](vector<HexahedronWid> &hsi, vector<double> &Kmodel, vector<Point> &I0out) {
			cubeGen(Volume{{-Lx/2., Lx/2., Nx}, {-Ly/2., Ly/2., Ny}, {-Lz/2., Lz/2., Nz}}, I0, hsi);
			Kmodel.assign(hsi.size(), K);
			I0out.assign(hsi.size(), I0);
		};
		const auto layersModelGenerator = [&](vector<HexahedronWid> &hsi, vector<double> &Kmodel, vector<Point> &I0out) {
			const vector<double> KtopToBottom = { 2, 1, 0.5 };
			const vector<double> layerLowerPlaneTopToBottom = { -1, -2, -3 };
			Assert(KtopToBottom.size() == layerLowerPlaneTopToBottom.size());
			const double extensionLength = 1e3;
			const int Nlayers = KtopToBottom.size();
			const int Nx = 100, Ny = 100, Nz = 5;

			for(int i = 0; i < Nlayers; ++i) {
				Volume vol{{-10, 10, Nx}, {-10, 10, Ny}, {layerLowerPlaneTopToBottom[i], i == 0? 0 : layerLowerPlaneTopToBottom[i-1], Nz}};
				const auto I0 = Hprime * KtopToBottom[i];
				vector<HexahedronWid> centerHex;
				cubeGen(vol, I0, centerHex);
				hsi.insert(hsi.begin(), centerHex.begin(), centerHex.end());
				
				// const auto extension = vector<Hexahedron>{};
				const auto extension = extendCuboidBoundedModelWith4cuboids(vol, extensionLength, I0);
				hsi.insert(hsi.end(), extension.begin(), extension.end());

				vector<double> kTmp;
				kTmp.assign(centerHex.size() + extension.size(), KtopToBottom[i]);
				Kmodel.insert(Kmodel.end(), kTmp.begin(), kTmp.end());
				vector<Point> I0Tmp;
				I0Tmp.assign(centerHex.size() + extension.size(), I0);
				I0out.insert(I0out.end(), I0Tmp.begin(), I0Tmp.end());
			}

			Assert(hsi.size() == Kmodel.size());
			Assert(hsi.size() == I0out.size());
		};
		const auto createCudaSolver = [&](const vector<HexahedronWid>& hsi, const bool transpose) {
			return gFieldSolver::getCUDAsolver(&*hsi.cbegin(), &*hsi.cend(), transpose);
			// const double replDist = 3;
			// return gFieldSolver::getCUDAreplacingSolver(&*hsi.cbegin(), &*hsi.cend(), replDist, nl*nB*nR*8);
		};

		Stopwatch tmr;
		tmr.start();
		vector<HexahedronWid> hsi = demagCG<HexahedronWid>(layersModelGenerator, createCudaSolver);
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

		// for(int layer = 0; layer < layersN; ++layer) {
		// 	cout << layer << ": "  << "Jmean= " << meanStatLayers[layer].get() << " | rms_err= " << rmsStatLayers[layer].get() << endl;
		// }

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

		{
			Dat3D<Point> dd;
			const double size = 20.1;
			const double step = 0.2;
			const double H = 0.001;
			for (double x = -size; x < size; x += step)
				for (double y = -size; y < size; y += step)
					dd.es.push_back({{ x, y, H }});
			fOnDat(dd);
			dd.write("layers_f.dat");
			return;
		}

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

		// {
		// 	Dat3D<Point> dd;
		// 	for(int i = 0; i < hsi.size(); ++i) {
		// 		const auto c = hsi[i].massCenter();
		// 		dd.es.push_back({{c.x, c.y, c.z}, hsi[i].dens - I0});
		// 	}
		// 	dd.write("cube_Jsnd_all_in.dat");
		// }
	}

	void runWell(int argc, char *argv[]) {
		InputParser inp(argc, argv);

		//well Model
		const auto wellModelGenerator = [&](vector<HexahedronWid> &hsi, vector<double> &Kmodel, vector<Point> &I0){
			Volume v{ { -2, 10, 120 },{ -2, 10, 120 },{ -8, 0, 40 } };
			inp.parseIfExists("xlower", v.x.lower);
			inp.parseIfExists("xupper", v.x.upper);
			inp.parseIfExists("xn", v.x.n);
			inp.parseIfExists("ylower", v.y.lower);
			inp.parseIfExists("yupper", v.y.upper);
			inp.parseIfExists("yn", v.y.n);
			inp.parseIfExists("zlower", v.z.lower);
			inp.parseIfExists("zupper", v.z.upper);
			inp.parseIfExists("zn", v.z.n);

			Cylinder well = { { { 4, 4 }, 0.25 }, 100 }; //x_cener, y_center, r, h
			inp.parseIfExists("cxwell", well.center.x);
			inp.parseIfExists("cywell", well.center.y);
			inp.parseIfExists("rwell", well.r);
			inp.parseIfExists("hwell", well.h);
			double Kouter = 0.2;
			double KinnerArg = 0;
			inp.parseIfExists("Kouter", Kouter);
			inp.parseIfExists("Kinner", KinnerArg);
			const vector<double> Kinner(v.z.n, KinnerArg);
			Point Hprime = { 14, 14, 35 }; //~40A/m
			inp.parseIfExists("HprimeX", Hprime.x);
			inp.parseIfExists("HprimeY", Hprime.y);
			inp.parseIfExists("HprimeZ", Hprime.z);
			
			wellGen(v, well, Hprime, Kouter, Kinner, hsi, Kmodel);
			I0.resize(hsi.size());
			std::transform(hsi.cbegin(), hsi.cend(), I0.begin(), [](const HexahedronWid &h) {return h.dens;});
		};

		
		double fieldH = 0.25;
		inp.parseIfExists("fieldH", fieldH);
		limits fieldDimX = { 2.5, 5.5, 150 };
		limits fieldDimY = { 3, 5.5, 125 };
		inp.parseIfExists("fieldXfrom", fieldDimX.lower);
		inp.parseIfExists("fieldXto", fieldDimX.upper);
		inp.parseIfExists("fieldXn", fieldDimX.n);
		inp.parseIfExists("fieldYfrom", fieldDimY.lower);
		inp.parseIfExists("fieldYto", fieldDimY.upper);
		inp.parseIfExists("fieldYn", fieldDimY.n);
		calcDemag<HexahedronWid>(wellModelGenerator, fieldDimX, fieldDimY, fieldH, "well", -1);
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

		//re-calculate J (magnetization) in every ClosedShape (simple iteration Mod)
		const auto &fJnMod = [&createCudaSolver, &hsi, &K, &J0, this]() -> vector<Point> {
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
					field[i] = ( 
							(field[i] / (4.*M_PI)) * 2. * K[i] + 
							J0[i] * 2. 
							+ hsi[i].dens * K[i]
						) / ( 2. + K[i] ) ;

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
				// fJnMod();
				OpCGt();
			}
			cout << "Done." << endl;
			return {};
		}
		
		const double eps = 1e-3;
		const double minRelativeErrChange = 0.05; //5%
		const int maxIter = 100;
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

			// const vector<Point> In = fJnMod();
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
		// const Volume& v,
		const limits& fieldDimX, const limits& fieldDimY, const double H,
		string filePrefix,
		const double DPR
	) {
		const auto createCudaSolver = [&DPR](const vector<ClosedShape>& hsi, const bool transpose) {
			return	DPR < 0
				? gFieldSolver::getCUDAsolver(&*hsi.cbegin(), &*hsi.cend(), transpose)
				: gFieldSolver::getCUDAreplacingSolver(&*hsi.cbegin(), &*hsi.cend(), DPR, hsi.size())
			;
		};

		vector<ClosedShape> hsi = demagCG<ClosedShape>(modelGenerator, createCudaSolver);
		// No demag solving
		// vector<ClosedShape> hsi;
		// vector<double> K;
		// vector<Point> I0;
		// modelGenerator(hsi, K, I0);
		// filePrefix += "_no_demag";

		if(!isRoot()) return;

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


		fOn(filePrefix + "Field");
		//fInWell("inWell.dat");


		cout << "Master Done." << endl;
	}
};

int main(int argc, char *argv[]) {
	bool isRoot = true;
	try {
		WellDemagCluster().runWell(argc, argv);
		return 0;
	}
	catch (std::exception &ex) {
		if(isRoot) std::cerr << "Global exception: " << ex.what() << endl;
		return 1;
	}
	cout << "Done" << endl;
	return 0;
}
