
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
#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include "mobj.h"
#include "calcField.h"
#include "GrafenArgs.h"
#include "TopogravArgs.h"
#include "Stopwatch.h"
#include "Quadrangles.h"
#include "MPIwrapper.h"
#include "MPIpool.h"

#include "sharedMem.h"
using gElementsShared = std::vector<HexahedronWid, ShmemAllocator<HexahedronWid>>;

using namespace std;
using GeographicLib::TransverseMercator;

#define G_CONST -6.67408
#define MB (1024*1024)		//bytes in megabyte

#define MIN_DENS 1e-8


#define Assert(exp) do { if (!(exp)) throw runtime_error("Assertion failed at: " + string(__FILE__) + " # line " + string(std::to_string(__LINE__))); } while (0)

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
	return min(v1, v2);
}

//construct density model from topography model (mesh nodes in degrees)
template <class VAlloc>
void makeHexsTopo(const Ellipsoid &e, const Grid &topo, const double dens, vector<HexahedronWid, VAlloc> &hsi) {
	const int xN = topo.nCol - 1, yN = topo.nRow - 1; //ln, Bn
	hsi.resize(xN * yN);

#pragma omp parallel for
	for(int yi = 0; yi < yN; ++yi)
		for (int xi = 0; xi < xN; ++xi) {
			Hexahedron h({
				//above
				e.getPoint(toRad(topo.yAt(yi + 1)), toRad(topo.xAt(xi + 1)), topo.at(xi + 1, yi + 1)),	//top right
				e.getPoint(toRad(topo.yAt(yi)), toRad(topo.xAt(xi + 1)), topo.at(xi + 1, yi)), //bottom right
				e.getPoint(toRad(topo.yAt(yi + 1)), toRad(topo.xAt(xi)), topo.at(xi, yi + 1)), //top left
				e.getPoint(toRad(topo.yAt(yi)), toRad(topo.xAt(xi)), topo.at(xi, yi)), //bottom left
				//below
				e.getPoint(toRad(topo.yAt(yi + 1)), toRad(topo.xAt(xi + 1)), 0),	//top right
				e.getPoint(toRad(topo.yAt(yi)), toRad(topo.xAt(xi + 1)), 0), //bottom right
				e.getPoint(toRad(topo.yAt(yi + 1)), toRad(topo.xAt(xi)), 0), //top left
				e.getPoint(toRad(topo.yAt(yi)), toRad(topo.xAt(xi)), 0), //bottom left
			}, dens);

			hsi[yi*xN + xi] = HexahedronWid(h);
		}
}

//split dencity model to Hexahedrons
template <class VAlloc>
void makeHexs(const double l0, const Ellipsoid &e, limits Nlim, limits Elim, limits Hlim,
	const vector<vector<double>> &dens, vector<HexahedronWid, VAlloc> &hsi) {

	Nlim = { Nlim.lower - Nlim.dWh() / 2.,  Nlim.upper + Nlim.dWh() / 2., Nlim.n };
	Elim = { Elim.lower - Elim.dWh() / 2.,  Elim.upper + Elim.dWh() / 2., Elim.n };

	cout << "grid E (x): " << Elim << " " << Elim.d() << endl;
	cout << "grid N (y): " << Nlim << " " << Nlim.d() << endl;
	hsi.resize(Hlim.n*Nlim.n*Elim.n);


	const TransverseMercator proj(e.Req*1000., e.f, 1);

	auto GKtoGeo = [&](const double N, double E) -> pair<double, double> {
		pair<double, double> lb;
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
	for (int Hi = 0; Hi < Hlim.n; ++Hi)
		for (int Ni = 0; Ni < Nlim.n; ++Ni)
			for (int Ei = 0; Ei < Elim.n; ++Ei)
				Kr(Hi, Ni, Ei);

}

//calculate transpose gravity field operator on a single node
void transFieldNode(Ellipsoid e, double l0, const vector<Dat3D::Element> psGK,
	const vector<HexahedronWid>::const_iterator &hsBegin, const vector<HexahedronWid>::const_iterator &hsEnd, const vector<double>::iterator &resBegin,
	const double dotPotentialRad) {
	vector<FieldPoint> fps(psGK.size());
	const TransverseMercator proj(e.Req*1000., e.f, 1);
	std::transform(psGK.cbegin(), psGK.cend(), fps.begin(), [&](const Dat3D::Element &el) {
		double l, B;
		double x = xFromGK(el.p.x, l0);
		proj.Reverse(toDeg(l0), x*1000., el.p.y*1000., B, l);
		l = toRad(l);
		B = toRad(B);
		return FieldPoint{ e.getPoint(B, l, el.p.z), e.getNormal(B, l), el.val };
	});

	//solve
	gFieldInvSolver::getCUDAtransSolver(fps, dotPotentialRad)
		->solve(&*hsBegin, &*hsEnd, resBegin);
	std::transform(resBegin, resBegin + (hsEnd - hsBegin), resBegin, [](auto &v) {return G_CONST*v; });
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
	GKNormOpts(const Ellipsoid &e, const double l0) : calcFieldNodeBase(e), l0(l0) {}
};
using calcFieldNodeOpts = boost::variant<GeoNormOpts, GKNormOpts>;

//Calculate gravity field on a single node
//If GKNormOpts: fp - field poinst in GK, field in the normal derection to reference ellipsoid in field point
//If GeoNormOpts: fp - field poinst in geo (E, N, H) = (l, B, H) = (lon, lat, H) = (deg, deg, km), field in the normal derection to reference ellipsoid in field point
void calcFieldNode(const calcFieldNodeOpts &opts, unique_ptr<gFieldSolver> &solver,
	const vector<Dat3D::Point> &fp, vector<double> &result) {
	Assert(fp.size() == result.size());

	class CalcFieldVisitor : public boost::static_visitor<void> {
		unique_ptr<gFieldSolver> &solver;
		const vector<Dat3D::Point> &fp;
		vector<double> &result;
	public:
		CalcFieldVisitor(unique_ptr<gFieldSolver> &solver, const vector<Dat3D::Point> &fp, vector<double> &result) : solver(solver), fp(fp), result(result) {}
		void operator()(const GKNormOpts& opts) const {
			const TransverseMercator proj(opts.e.Req*1000., opts.e.f, 1);
			for (size_t i = 0; i < fp.size(); ++i) {
				const Dat3D::Point &p = fp[i];
				double l, B;
				double x = xFromGK(p.x, opts.l0);
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
				const Point p0 = opts.e.getPoint(toRad(p.y), toRad(p.x), p.z);
				const Point n0 = opts.n ? *opts.n : opts.e.getNormal(toRad(p.y), toRad(p.x));
				result[i] += G_CONST * solver->solve(p0, n0);
			}
		}
	};

	boost::apply_visitor(CalcFieldVisitor(solver, fp, result), opts);
}

class ClusterSolver : public MPI {
public:
	const int maxGPUmemMB = 5000;
	int triBufferSize = 0;
	using Qiter = gElementsShared::const_iterator;

	SharedMemBase<gElementsShared> *sharedMem = 0;

	ClusterSolver() : MPI() {
		if (gridSize < 2) throw runtime_error("You must run at least 2 MPI processes.");
		if (root != 0) throw runtime_error("Root process must have rank = 0.");
		const auto lid = localId();
		const int devId = !isRoot() && get<1>(lid) ? get<0>(lid) - 1 : get<0>(lid);
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

	void calcField(const calcFieldNodeOpts &opts, gElementsShared &qs, Dat3D &dat, double dotPotentialRad) {
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
			const Qiter qend = qs.begin() + min(i + partSize, qs.size());

			cout << ":::part " << part + 1 << " out of " << parts << "::: size: " << qend - qbegin << endl;
			calcWithPool(opts, qbegin, qend, dat, dotPotentialRad);
		}
	}

	void calcTrans(Ellipsoid e, double l0, vector<gElements::base> &qs, const Dat3D &dat, vector<double> &result, const double dotPotentialRad) {
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

private:
	void calcWithPool(const calcFieldNodeOpts &opts, const Qiter &qbegin, const Qiter &qend, Dat3D &dat, double dotPotentialRad) {
		const vector<Dat3D::Point> fp = dat.getPoints();
		vector<double> result;
		//blocking process with rank 0 until the result's been gathered
		MPIpool<Dat3D::Point, double> pool(*this, fp, result, 1024);
		int cnt = 0;
		if (!isRoot()) {
			unique_ptr<gFieldSolver> solver = gFieldSolver::getCUDAsolver(&*qbegin, &*qend, dotPotentialRad, triBufferSize);
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
		Ellipsoid Earth(DEF_R_EQ, DEF_R_PL);

		//gElements qs;
		cout << "GPUs: " << cuSolver::getGPUnum() << endl;
		GrafenArgs inp(argc, argv);

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
		if (cs.isLocalRoot() && !inp.trasSolver)
			prep([&]() {makeHexs(inp.l0, Earth, inp.Nlim, inp.Elim, inp.Hlim, inp.dens, qss); });
		//prepocessing for transpose solver
		if (cs.isRoot() && inp.trasSolver)
			prep([&]() {makeHexs(inp.l0, Earth, inp.Nlim, inp.Elim, inp.Hlim, inp.dens, qs); });
		cs.Barrier();
		const size_t qSize = inp.trasSolver ? qs.size() : qss.size();
		cout << "Real size: " << (qSize * sizeof(gElements::base) + MB - 1) / MB << "MB" << endl;

		if (inp.trasSolver) {
			vector<double> res;
			if (cs.isRoot()) {
				cout << "TRANSPOSE SOLVER. Your *.grd files will be overwritten!" << endl;
				res.resize(qs.size());
			}
			Stopwatch tmr;
			cs.calcTrans(Earth, inp.l0, qs, inp.dat, res, inp.dotPotentialRad);
			if (cs.isRoot()) {
				cout << "time: " << tmr.stop() << endl;

				const int lsize = inp.Nlim.n*inp.Elim.n;
				for (int i = 0; i < inp.Hlim.n; ++i) {
					Grid g(inp.fnames[i]);
					g.data.assign(res.cbegin() + i*lsize, res.cbegin() + (i + 1)*lsize);
					g.Write();
				}

				cout << "Done" << endl;
			}
			return 0;
		}

		if (cs.isRoot()) {
			inp.dat.set(0);
			cout << "Clearing output file. Do NOT stop me!" << endl;
			if(!inp.grdFile) inp.dat.write(inp.dat2D);
			cout << "Clearing done. You may ctrl+c now." << endl;
		}

		cout << "Computing started" << endl;
		Stopwatch tmr;
		tmr.start();
		//set approximate size of buffer
		cs.triBufferSize = triBufferSize(inp.Nlim, inp.Elim, inp.Hlim, inp.dotPotentialRad);
		//do the job
		cs.calcField(GKNormOpts{ Earth, inp.l0 }, qss, inp.dat, inp.dotPotentialRad);
		const double time = tmr.stop();
		if (cs.isRoot()) {
			cout << "Computing finished: " << time << "sec." << endl << endl;
			if (!inp.grdFile) inp.dat.write(inp.dat2D);
			else GDconv::toGrd(inp.dat, inp.grdCols, inp.grdRows).Write(inp.grdFname);
		}
	} catch (exception &ex) {
		cout << "Global exception: " << ex.what() << endl;
		return 1;
	}
	cout << "Done" << endl;
	return 0;
}

/*
template <class VAlloc>
void saveAsDat(const vector<HexahedronWid, VAlloc> &hsi) {
	Dat3D dat;
	for (int i = 0; i < hsi.size(); ++i)
		for (int pi = 0; pi < 8; ++pi)
			dat.es.push_back({ { hsi[i].p[pi].x, hsi[i].p[pi].y, hsi[i].p[pi].z}, hsi[i].dens });
	dat.write("shape.dat");
}
*/

int topogravMain(int argc, char *argv[]) {

	const auto& sectionStopwatch = [](const string& sectionName, const std::function<void(void)>& section) {
		cout << sectionName << " started" << endl;
		Stopwatch tmr;
		tmr.start();
		section();
		const double time = tmr.stop();
		cout << sectionName << " finished: " << time << "sec." << endl;
	};

	try {
		ClusterSolver cs;
		cout << "GPUs: " << cuSolver::getGPUnum() << endl;
		TopogravArgs inp(argc, argv);
		const Grid topoGrid(inp.topoGridFname);
		const size_t qAm = (topoGrid.nCol - 1) * (topoGrid.nRow - 1);
		cout << "Elements: " << qAm << endl;
		cout << (qAm * sizeof(gElements::base) + MB - 1) / MB << "MB required. ";
		gElementsShared &qss = cs.initShared(qAm); //for direct solver. Should not actually allocate memory untill resize().
		cout << "Allocated." << endl;


		//preprocessing
		if(cs.isLocalRoot()) 
			sectionStopwatch("Preprocessing", [&](){ makeHexsTopo(inp.refEllipsoid, topoGrid, inp.dens, qss); });
		cs.Barrier();
		const size_t qSize = qss.size();
		cout << "Real size: " << (qSize * sizeof(gElements::base) + MB - 1) / MB << "MB" << endl;
		//if(cs.isRoot()) saveAsDat(qss);
		
		//computation
		Dat3D fieldDat = GDconv::toDat3D(topoGrid);
		fieldDat.forEach([&](auto &el) { el.p.z += 0.00001; }); //elivate computation point above the mass volume a little
		sectionStopwatch("Computing", [&]() {
			//set approximate size of buffer
			cs.triBufferSize = inp.pprr < 1e-6? 0 : int(qss.size() / 4);
			//do the job
			cs.calcField(GeoNormOpts{ inp.refEllipsoid, inp.normal }, qss, fieldDat, inp.pprr);
		});

		//save result
		if (cs.isRoot())
			GDconv::toGrd(fieldDat, topoGrid.nCol, topoGrid.nRow).Write(inp.gravGridFname);

	} catch (exception &ex) {
		cout << "Global exception: " << ex.what() << endl;
		return 1;
	}
	cout << "Done" << endl;
	return 0;
}

int main(int argc, char *argv[]) {
	//return grafenMain(argc, argv);
	return topogravMain(argc, argv);
}
