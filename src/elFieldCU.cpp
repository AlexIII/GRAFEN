
#define GEOGRAPHICLIB_SHARED_LIB 0
#define BOOST_ALL_NO_LIB
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
#include <numeric>
#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include "mobj.h"
#include "calcField.h"
#include "GrafenArgs.h"
#include "Stopwatch.h"
#include "Quadrangles.h"
#include "MPIwrapper.h"
#include "MPIpool.h"

#include "sharedMem.h"
using gElementsShared = std::vector<HexahedronWid, ShmemAllocator<HexahedronWid>>;

using std::string;
using std::vector;
using GeographicLib::TransverseMercator;

#define G_CONST -6.67408
#define MB (1024*1024)		//bytes in megabyte

#define MOVE_POINT_EPS (1e-8)

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

//estimate approximate buffer size for the Hexahedrons that can't be replaced by a singular source
int triBufferSize(const limits &Nlim, const limits &Elim, const limits &Hlim, const double r) {
	return std::max(
		size_t(75),
		std::min(Nlim.n, 1 + size_t(5 * r / Nlim.d())) * 
		std::min(Elim.n, 1 + size_t(5 * r / Elim.d())) * 
		std::min(Hlim.n, 1 + size_t(2 * r / Hlim.d()))
	);
}

//{B in deg, l in deg, z} -> {E, N, z} in km, l0 in Rad
Point GeoToGK(const TransverseMercator &proj, const GKoptions &opts, const double B, const double l, const double z) {
	Point p{ 0, 0, z };
	proj.Forward(toDeg(opts.l0), B, l, p.x, p.y);
	p.x /= 1000.;
	p.y /= 1000.;
	p.x += opts.xOffset;
	return p;
}

//construct density model from topography model (mesh nodes in degrees)
//topo - mesh in deg, height in km
//generates cuboids
template <class VAlloc>
void makeHexsTopoFlatCuboids(const GKoptions &opts, const Grid &topo, const double dens, vector<HexahedronWid, VAlloc> &hsi) {
	hsi.resize(topo.nCol * topo.nRow);

	const TransverseMercator proj(opts.e.Req*1000., opts.e.f, 1);
	auto GeoToGK3 = [&](const double B, const double l, const double z) -> Point { 
		return GeoToGK(proj, opts, B, l, z);
	};

	const double hdx = topo.xSize / 2., hdy = topo.ySize / 2.;

#pragma omp parallel for
	for (int yi = 0; yi < topo.nRow; ++yi)
		for (int xi = 0; xi < topo.nCol; ++xi) {
			const Hexahedron h({
				//above
				GeoToGK3(topo.yAt(yi) + hdy, topo.xAt(xi) + hdx, topo.at(xi, yi)),	//top right
				GeoToGK3(topo.yAt(yi) - hdy, topo.xAt(xi) + hdx, topo.at(xi, yi)), //bottom right
				GeoToGK3(topo.yAt(yi) + hdy, topo.xAt(xi) - hdx, topo.at(xi, yi)), //top left
				GeoToGK3(topo.yAt(yi) - hdy, topo.xAt(xi) - hdx, topo.at(xi, yi)), //bottom left
				//below
				GeoToGK3(topo.yAt(yi) + hdy, topo.xAt(xi) + hdx, 0),	//top right
				GeoToGK3(topo.yAt(yi) - hdy, topo.xAt(xi) + hdx, 0), //bottom right
				GeoToGK3(topo.yAt(yi) + hdy, topo.xAt(xi) - hdx, 0), //top left
				GeoToGK3(topo.yAt(yi) - hdy, topo.xAt(xi) - hdx, 0), //bottom left
			}, dens);

			hsi[yi*topo.nCol + xi] = HexahedronWid(h);
		}
}

//construct spherical density model from topography model (mesh nodes in degrees)
//topo - mesh in deg (geographic coord), height in km
template <class VAlloc>
void makeHexsTopoGeo(const Ellipsoid &e, const Grid &topo, const double dens, vector<HexahedronWid, VAlloc> &hsi) {
	const int xN = topo.nCol - 1, yN = topo.nRow - 1; //ln, Bn
	hsi.resize(xN * yN);

#pragma omp parallel for
	for(int yi = 0; yi < yN; ++yi)
		for (int xi = 0; xi < xN; ++xi) {
			double h0 = topo.at(xi + 1, yi + 1), h1 = topo.at(xi + 1, yi), h2 = topo.at(xi, yi + 1), h3 = topo.at(xi, yi);
			if(h0 <= 0) h0 = MOVE_POINT_EPS;
			if(h1 <= 0) h1 = MOVE_POINT_EPS;
			if(h2 <= 0) h2 = MOVE_POINT_EPS;
			if(h3 <= 0) h3 = MOVE_POINT_EPS;

			const Hexahedron h({
				//above
				e.getPoint(toRad(topo.yAt(yi + 1)), toRad(topo.xAt(xi + 1)), h0),	//top right
				e.getPoint(toRad(topo.yAt(yi)), toRad(topo.xAt(xi + 1)), h1), //bottom right
				e.getPoint(toRad(topo.yAt(yi + 1)), toRad(topo.xAt(xi)), h2), //top left
				e.getPoint(toRad(topo.yAt(yi)), toRad(topo.xAt(xi)), h3), //bottom left
				//below
				e.getPoint(toRad(topo.yAt(yi + 1)), toRad(topo.xAt(xi + 1)), 0),	//top right
				e.getPoint(toRad(topo.yAt(yi)), toRad(topo.xAt(xi + 1)), 0), //bottom right
				e.getPoint(toRad(topo.yAt(yi + 1)), toRad(topo.xAt(xi)), 0), //top left
				e.getPoint(toRad(topo.yAt(yi)), toRad(topo.xAt(xi)), 0), //bottom left
			}, dens);

			hsi[yi*xN + xi] = HexahedronWid(h);
		}
}

//construct spherical density model from topography model (mesh nodes in degrees)
//topo - mesh in km (GK coord), height in km
template <class VAlloc>
void makeHexsTopoGK(const GKoptions &GKopts, const Grid &topo, const double dens, vector<HexahedronWid, VAlloc> &hsi) {
	const int xN = topo.nCol - 1, yN = topo.nRow - 1; //ln, Bn
	hsi.resize(xN * yN);

	const TransverseMercator proj(GKopts.e.Req * 1000., GKopts.e.f, 1);
	auto GKtoGeo = [&](const double N, double E) -> std::pair<double, double> {
		std::pair<double, double> lb;
		E -= GKopts.xOffset;
		proj.Reverse(toDeg(GKopts.l0), E*1000., N*1000., lb.second, lb.first);
		lb.first = toRad(lb.first);
		lb.second = toRad(lb.second);
		return lb;
	};

	const auto& e = GKopts.e;

#pragma omp parallel for
	for(int yi = 0; yi < yN; ++yi)
		for (int xi = 0; xi < xN; ++xi) {
			double h0 = topo.at(xi + 1, yi + 1), h1 = topo.at(xi + 1, yi), h2 = topo.at(xi, yi + 1), h3 = topo.at(xi, yi);
			if(h0 <= 0) h0 = MOVE_POINT_EPS;
			if(h1 <= 0) h1 = MOVE_POINT_EPS;
			if(h2 <= 0) h2 = MOVE_POINT_EPS;
			if(h3 <= 0) h3 = MOVE_POINT_EPS;
			
			const Hexahedron hex({
				//above
				e.getPoint(GKtoGeo(topo.yAt(yi + 1), topo.xAt(xi + 1)), h0),	//top right
				e.getPoint(GKtoGeo(topo.yAt(yi), topo.xAt(xi + 1)), h1), 		//bottom right
				e.getPoint(GKtoGeo(topo.yAt(yi + 1), topo.xAt(xi)), h2), 		//top left
				e.getPoint(GKtoGeo(topo.yAt(yi), topo.xAt(xi)), h3), 			//bottom left
				//below
				e.getPoint(GKtoGeo(topo.yAt(yi + 1), topo.xAt(xi + 1)), 0),		//top right
				e.getPoint(GKtoGeo(topo.yAt(yi), topo.xAt(xi + 1)), 0), 		//bottom right
				e.getPoint(GKtoGeo(topo.yAt(yi + 1), topo.xAt(xi)), 0), 		//top left
				e.getPoint(GKtoGeo(topo.yAt(yi), topo.xAt(xi)), 0), 			//bottom left
			}, dens);
			
			hsi[yi*xN + xi] = HexahedronWid(hex);
		}
}

template <class VectorIterator>
size_t makeHexsTopoFlatCuboidsGK(const GKoptions &GKopts, limits Nlim, limits Elim,
		const vector<double> &heights, const vector<double> &dens, 
		VectorIterator hsiBegin, VectorIterator hsiEnd) {

	Nlim = { Nlim.lower - Nlim.dWh() / 2.,  Nlim.upper + Nlim.dWh() / 2., Nlim.n };
	Elim = { Elim.lower - Elim.dWh() / 2.,  Elim.upper + Elim.dWh() / 2., Elim.n };

	const size_t hsiSize = Nlim.n * Elim.n;
	if(hsiEnd - hsiBegin < hsiSize)
		throw std::runtime_error("makeHexsTopoFlatCuboidsGK(): small hsi buffer size");

	const TransverseMercator proj(GKopts.e.Req*1000., GKopts.e.f, 1);

	auto GKtoGeo = [&](const double N, double E) -> std::pair<double, double> {
		std::pair<double, double> lb;
		E -= GKopts.xOffset;
		proj.Reverse(toDeg(GKopts.l0), E * 1000., N * 1000., lb.second, lb.first);
		lb.first = toRad(lb.first);
		lb.second = toRad(lb.second);
		return lb;
	};

	const auto& e = GKopts.e;
	const double lowLevel = 0;

	auto Kr = [&](const int Ni, const int Ei) {
		double height = heights[Ni*Elim.n + Ei];
		if(height <= lowLevel) height = lowLevel + MOVE_POINT_EPS;

		Hexahedron h({
			e.getPoint(GKtoGeo(Nlim.at(Ni + 1), Elim.at(Ei + 1)), height),
			e.getPoint(GKtoGeo(Nlim.at(Ni), Elim.at(Ei + 1)), height),
			e.getPoint(GKtoGeo(Nlim.at(Ni + 1), Elim.at(Ei)), height),
			e.getPoint(GKtoGeo(Nlim.at(Ni), Elim.at(Ei)), height),

			e.getPoint(GKtoGeo(Nlim.at(Ni + 1), Elim.at(Ei + 1)), lowLevel),
			e.getPoint(GKtoGeo(Nlim.at(Ni), Elim.at(Ei + 1)), lowLevel),
			e.getPoint(GKtoGeo(Nlim.at(Ni + 1), Elim.at(Ei)), lowLevel),
			e.getPoint(GKtoGeo(Nlim.at(Ni), Elim.at(Ei)), lowLevel)
		}, dens[Ni*Elim.n + Ei]);

		hsiBegin[Ni*Elim.n + Ei] = HexahedronWid(h);
	};

	for (size_t Ni = 0; Ni < Nlim.n; ++Ni)
		for (size_t Ei = 0; Ei < Elim.n; ++Ei)
			Kr(Ni, Ei);

	return hsiSize;
}
		
//split dencity model to Hexahedrons
template <class VectorIterator>
size_t makeHexs(const GKoptions &GKopts, limits Nlim, limits Elim, limits Hlim,
	const vector<vector<double>> &dens, 
	VectorIterator hsiBegin, VectorIterator hsiEnd) {
	
	Nlim = { Nlim.lower - Nlim.dWh() / 2.,  Nlim.upper + Nlim.dWh() / 2., Nlim.n };
	Elim = { Elim.lower - Elim.dWh() / 2.,  Elim.upper + Elim.dWh() / 2., Elim.n };

	// cout << "grid E (x): " << Elim << " " << Elim.d() << endl;
	// cout << "grid N (y): " << Nlim << " " << Nlim.d() << endl;
	
	const size_t hsiSize = Hlim.n * Nlim.n * Elim.n;
	if(hsiEnd - hsiBegin < hsiSize)
		throw std::runtime_error("makeHexs(): small hsi buffer size");

	const TransverseMercator proj(GKopts.e.Req * 1000., GKopts.e.f, 1);

	auto GKtoGeo = [&](const double N, double E) -> std::pair<double, double> {
		std::pair<double, double> lb;
		E -= GKopts.xOffset;
		proj.Reverse(toDeg(GKopts.l0), E * 1000., N * 1000., lb.second, lb.first);
		lb.first = toRad(lb.first);
		lb.second = toRad(lb.second);
		return lb;
	};

	const auto& e = GKopts.e;

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

		hsiBegin[(Hi*Nlim.n + Ni)*Elim.n + Ei] = HexahedronWid(h);
	};

#pragma omp parallel for
	for (size_t Hi = 0; Hi < Hlim.n; ++Hi)
		for (size_t Ni = 0; Ni < Nlim.n; ++Ni)
			for (size_t Ei = 0; Ei < Elim.n; ++Ei)
				Kr(Hi, Ni, Ei);

	return hsiSize;
}

//calculate transpose gravity field operator on a single node
void transFieldNode(
	const GKoptions &GKopts, const vector<Dat3D::Element> psGK,
	const vector<HexahedronWid>::const_iterator &hsBegin, const vector<HexahedronWid>::const_iterator &hsEnd, const vector<double>::iterator &resBegin,
	const double dotPotentialRad
) {
	vector<FieldPoint> fps(psGK.size());
	const TransverseMercator proj(GKopts.e.Req * 1000., GKopts.e.f, 1);
	std::transform(psGK.cbegin(), psGK.cend(), fps.begin(), [&](const Dat3D::Element &el) {
		double l, B;
		double x = el.p.x - GKopts.xOffset;
		proj.Reverse(toDeg(GKopts.l0), x * 1000., el.p.y * 1000., B, l);
		l = toRad(l);
		B = toRad(B);
		return FieldPoint{ GKopts.e.getPoint(B, l, el.p.z), GKopts.e.getNormal(B, l), el.val };
	});

	//solve
	gFieldInvSolver::getCUDAtransSolver(fps, dotPotentialRad)->solve(&*hsBegin, &*hsEnd, resBegin);
	std::transform(resBegin, resBegin + (hsEnd - hsBegin), resBegin, [](auto &v) {return G_CONST * v; });
}

struct calcFieldNodeBase {
	Ellipsoid e;
	boost::optional<Point> n;
	calcFieldNodeBase(const Ellipsoid &e, boost::optional<Point> n = boost::optional<Point>{}) : e(e), n(n) {}
};
struct GeoNormOpts : calcFieldNodeBase {
	GeoNormOpts(const Ellipsoid &e, boost::optional<Point> n = boost::optional<Point>{}) : calcFieldNodeBase(e, n) {}
};
struct GKNormOpts : calcFieldNodeBase { 
	double l0;
	double xOffset;
	GKNormOpts(GKoptions &opts, boost::optional<Point> n = boost::optional<Point>{}) : calcFieldNodeBase(opts.e, n), l0(opts.l0), xOffset(opts.xOffset) {}
};
struct GeoNormOptsFlat : GKNormOpts {
	GeoNormOptsFlat(GKoptions &opts, boost::optional<Point> n = boost::optional<Point>{}) : GKNormOpts(opts, n) {}
};
using calcFieldNodeOpts = boost::variant<GeoNormOpts, GKNormOpts, GeoNormOptsFlat>;

//Calculate gravity field on a single node
//If GKNormOpts: fp - field points in GK, field in the normal derection to reference ellipsoid in field point
//If GeoNormOpts: fp - field points in geo (E, N, H) = (l, B, H) = (lon, lat, H) = (deg, deg, km), field in the normal derection to reference ellipsoid in field point
void calcFieldNode(const calcFieldNodeOpts &opts, std::unique_ptr<gFieldSolver> &solver,
	const vector<Dat3D::Point> &fp, vector<double> &result) {
	Assert(fp.size() == result.size());

	class CalcFieldVisitor : public boost::static_visitor<void> {
		std::unique_ptr<gFieldSolver> &solver;
		const vector<Dat3D::Point> &fp;
		vector<double> &result;
	public:
		CalcFieldVisitor(std::unique_ptr<gFieldSolver> &solver, const vector<Dat3D::Point> &fp, vector<double> &result) : solver(solver), fp(fp), result(result) {}
		void operator()(const GKNormOpts& opts) const {
			const TransverseMercator proj(opts.e.Req*1000., opts.e.f, 1);
			for (size_t i = 0; i < fp.size(); ++i) {
				Dat3D::Point p = fp[i];
				p.x += MOVE_POINT_EPS;
				p.y += MOVE_POINT_EPS;
				p.z += MOVE_POINT_EPS;
				double l, B;
				const double x = p.x - opts.xOffset;
				proj.Reverse(toDeg(opts.l0), x*1000., p.y*1000., B, l);
				l = toRad(l);
				B = toRad(B);
				const Point p0 = opts.e.getPoint(B, l, p.z);
				const Point n0 = opts.n? *opts.n : opts.e.getNormal(B, l);
				result[i] += G_CONST * solver->solve(p0, n0);
			}
		}
		void operator()(const GeoNormOpts& opts) const {
			for (size_t i = 0; i < fp.size(); ++i) {
				const Dat3D::Point &p = fp[i];
				Point p0 = opts.e.getPoint(toRad(p.y), toRad(p.x), p.z);
				p0.x += MOVE_POINT_EPS;
				p0.y += MOVE_POINT_EPS;
				p0.z += MOVE_POINT_EPS;
				const Point n0 = opts.n ? *opts.n : opts.e.getNormal(toRad(p.y), toRad(p.x));
				result[i] += G_CONST * solver->solve(p0, n0);
			}
		}
		void operator()(const GeoNormOptsFlat& opts) const {
			const TransverseMercator proj(opts.e.Req*1000., opts.e.f, 1);
			for (size_t i = 0; i < fp.size(); ++i) {
				const Dat3D::Point &p = fp[i];
				Point p0 = GeoToGK(proj, GKoptions{ opts.e, opts.l0, opts.xOffset  }, p.y, p.x, p.z);
				p0.x += MOVE_POINT_EPS;
				p0.y += MOVE_POINT_EPS;
				p0.z += MOVE_POINT_EPS;
				const Point n0 = opts.n ? *opts.n : Point{ 0, 0, 1 }; //in z direction
				result[i] += G_CONST * solver->solve(p0, n0);
			}
		}
	};

	boost::apply_visitor(CalcFieldVisitor(solver, fp, result), opts);
}

class ClusterSolver : public MPIwrapper {
public:
	const int maxGPUmemMB = 5000;
	int triBufferSize = 0;
	using Qiter = gElementsShared::const_iterator;

	SharedMemBase<gElementsShared> *sharedMem = 0;

	ClusterSolver(const vector<int> &gpuIdMap = {}) : MPIwrapper() {
		if (gridSize < 2) throw std::runtime_error("You must run at least 2 MPI processes.");
		if (root != 0) throw std::runtime_error("Root process must have rank = 0.");
		const auto lid = localId();
		const int devId = !isRoot() && std::get<1>(lid) ? std::get<0>(lid) - 1 : std::get<0>(lid);
		const int mappedDevId = devId < gpuIdMap.size()? gpuIdMap[devId] : devId;
		cuSolver::setDevice(mappedDevId);
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

	void calcField(const calcFieldNodeOpts &opts, gElementsShared &qs, Dat3D &dat, double dotPotentialRad) {
		if (isRoot()) cout << "Computing nodes: " << gridSize - 1 << endl;

		size_t qSize = qs.size();
		if (myId == 1) send(qSize, 0);
		else if (myId == 0) recv(qSize, 1);

		size_t partSize = getPartSize<gElements::base>(qSize);

		const size_t parts = ((partSize + qSize - 1) / partSize);
		for (size_t i = 0; i < qSize; i += partSize) {
			const size_t part = (i / partSize);
			const Qiter qbegin = qs.begin() + i;
			const Qiter qend = qs.begin() + std::min(i + partSize, qs.size());

			if(isRoot()) cout << ":::part " << part + 1 << " out of " << parts << "::: size: " << qend - qbegin << endl;
			calcWithPool(opts, qbegin, qend, dat, dotPotentialRad);
		}
	}

	void calcTrans(const GKoptions &GKopts, vector<gElements::base> &qs, const Dat3D &dat, vector<double> &result, const double dotPotentialRad) {
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
				transFieldNode(GKopts, dat.es, task.cbegin(), task.cend(), result.begin(), dotPotentialRad);
				pool.submit(result);
			}
		}
		else cout << "result gather ok" << endl;
	}

private:
	void calcWithPool(const calcFieldNodeOpts &opts, const Qiter &qbegin, const Qiter &qend, Dat3D &dat, double dotPotentialRad) {
		const vector<Dat3D::Point> fp = dat.getPoints();
		vector<double> result;
		//blocking process with rank 0 until the result's been gathered
		MPIpool<Dat3D::Point, double> pool(*this, fp, result, 1024);
		int cnt = 0;
		if (!isRoot()) {
			std::unique_ptr<gFieldSolver> solver = gFieldSolver::getCUDAsolver(&*qbegin, &*qend, dotPotentialRad, triBufferSize);
			while (1) {
				vector<Dat3D::Point> task = pool.getTask();
				if (!task.size()) break;
				cout << "Task accepted " << ++cnt << " size: " << task.size()+1 << endl;
				vector<double> result(task.size());
				calcFieldNode(opts, solver, task, result);
				pool.submit(result);
			}
		}
		else {
			dat.add(result);
			cout << "result gather ok" << endl;
		}
	}
};


int grafenMain(int argc, char *argv[]) {
	try {
		ClusterSolver cs;

		if(cs.isRoot()) {
			cout << "GPUs: " << cuSolver::getGPUnum() << endl;
			cout << "Workers: " << cs.gridSize - 1 << endl;
		}
		GrafenArgs inp(argc, argv, cs.isRoot());

		const size_t qAm = getHexAm(inp.Nlim.n, inp.Elim.n, inp.Hlim.n + (inp.withTopo? 1 : 0));
		if(cs.isRoot()) {
			cout << (qAm * sizeof(gElements::base) + MB - 1) / MB << "MB required. ";
			cout << "Elements: " << qAm << endl;
		}
		gElementsShared &qss = cs.initShared(qAm);//for direct solver. Should not actually allocate memory untill resize().
		gElements qs; //for transpose solver

		const auto& makeHexsPreprocess = [&inp, &qAm](auto& hsi) {
			cout << "Preprocessing started" << endl;
			Stopwatch tmr;
			tmr.start();

			hsi.resize(qAm);
			cout << "Allocated." << endl;
			size_t hsiOffest = 0;
			if(inp.withTopo) {
				hsiOffest += makeHexsTopoFlatCuboidsGK(inp.GKopts, inp.Nlim, inp.Elim, inp.topoHeights, inp.topoDens, hsi.begin() + hsiOffest, hsi.end());
			}
			hsiOffest += makeHexs(inp.GKopts, inp.Nlim, inp.Elim, inp.Hlim, inp.dens, hsi.begin() + hsiOffest, hsi.end()); 

			if(hsiOffest != hsi.size())
				throw std::runtime_error("Preprocessing: hsi was incorrect");

			const double time = tmr.stop();
			cout << "Preprocessing finished: " << time << "sec." << endl;
		};
		
		//prepocessing for direct solver
		if (!inp.transSolver && cs.isLocalRoot()) {
			makeHexsPreprocess(qss);
		}
		//prepocessing for transpose solver
		else if (inp.transSolver && cs.isRoot()) {
			makeHexsPreprocess(qs);
		}
		cs.Barrier();
		const size_t qSize = inp.transSolver ? qs.size() : qss.size();
		cout << "Real size: " << (qSize * sizeof(gElements::base) + MB - 1) / MB << "MB" << " (" << qSize << " elements)" << endl;

		if (inp.transSolver) {
			// ----- TRANSPOSE SOLVER
			vector<double> res;
			if (cs.isRoot()) {
				cout << "TRANSPOSE SOLVER. Your *.grd files will be overwritten!" << endl;
				res.resize(qs.size());
			}
			Stopwatch tmr;
			cs.calcTrans(inp.GKopts, qs, inp.dat, res, inp.dotPotentialRad);
			if (cs.isRoot()) {
				cout << "Computing finished: " << tmr.stop() << "sec." << endl << endl;
				const size_t layersOffset = inp.withTopo? 1 : 0;
				const size_t lsize = inp.Nlim.n*inp.Elim.n;
				if(layersOffset) {
					Grid g(inp.topoDensFname);
					g.data.assign(res.cbegin(), res.cbegin() + lsize);
					g.Write();
				}
				for (size_t i = 0; i < inp.Hlim.n; ++i) {
					Grid g(inp.fnames[i]);
					g.data.assign(res.cbegin() + (i + layersOffset)*lsize, res.cbegin() + (i + layersOffset + 1)*lsize);
					g.Write();
				}
			}
		} else {
			// ----- DIRECT SOLVER
			if (cs.isRoot()) {
				inp.dat.set(0);
				cout << "Clearing output file. Do NOT stop me!" << endl;
				if(!inp.grdFile) inp.dat.write(inp.dat2D);
				cout << "Clearing done. You may ctrl+c now." << endl;
			}
			cout << "DIRECT SOLVER" << endl;
			Stopwatch tmr;
			tmr.start();
			//set approximate size of buffer
			cs.triBufferSize = triBufferSize(inp.Nlim, inp.Elim, inp.withTopo? limits{inp.Hlim.lower, inp.Hlim.upper + inp.Hlim.d(), inp.Hlim.n + 1} : inp.Hlim, inp.dotPotentialRad);
			cout << "triBufferSize = " << cs.triBufferSize << " " << 100 * cs.triBufferSize / qSize << "%" << endl;
			//do the job
			cs.calcField(GKNormOpts{ inp.GKopts }, qss, inp.dat, inp.dotPotentialRad);
			if (cs.isRoot()) {
				cout << "Computing finished: " << tmr.stop() << "sec." << endl << endl;
				if (!inp.grdFile) inp.dat.write(inp.dat2D);
				else GDconv::toGrd(inp.dat, inp.grdCols, inp.grdRows).Write(inp.grdFname);
			}
		}

		cout << "Done" << endl;
	} catch (std::exception &ex) {
		cout << "Global exception: " << ex.what() << endl;
		return 1;
	}
	cout << "Done" << endl;
	return 0;
}

int main(int argc, char *argv[]) {
	return grafenMain(argc, argv);
}
