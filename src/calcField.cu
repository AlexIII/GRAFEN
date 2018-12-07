#include "calcField.h"
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/copy.h>
#include <thrust/inner_product.h>
#include <thrust/transform_reduce.h>
#include <thrust/execution_policy.h>
#include "cuVar.h"
#include <mutex>

using std::cout;
using std::endl;

#define BLOCK_SZ 32
#define MAX_DEV 16
static void CheckCudaErrorAux(const char *, unsigned, const char *, cudaError_t);
#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)

__host__ __device__ inline double tripleprod(const Point& o1, const Point& o2, const Point& o3) { return (o1 ^ (o2*o3)); }

__host__ __device__ Point intTrAn(const Point &p0, const Triangle &t) {
	Point res;

	const Point a1 = t.p1 - p0;
	const Point a2 = t.p2 - p0;
	const Point a3 = t.p3 - p0;
	const double a1m = a1.eqNorm();
	const double a2m = a2.eqNorm();
	const double a3m = a3.eqNorm();
	const Point a12 = t.p2 - t.p1;
	const Point a23 = t.p3 - t.p2;
	const Point a31 = t.p1 - t.p3;
	const double a12m = a12.eqNorm();
	const double a23m = a23.eqNorm();
	const double a31m = a31.eqNorm();
	//if (a3m + a1m - a31m < 1e-9 || a1m + a2m - a12m < 1e-9 || a2m + a3m - a23m < 1e-9) return{ HUGE_VAL, HUGE_VAL, HUGE_VAL };
	const Point N = t.normal();

	res = a31*(log((a3m + a1m + a31m) / (a3m + a1m - a31m)) / a31m);
	res += a12*(log((a1m + a2m + a12m) / (a1m + a2m - a12m)) / a12m);
	res += a23*(log((a2m + a3m + a23m) / (a2m + a3m - a23m)) / a23m);
	res = N*res;
	res += N*(2.0*atan2(tripleprod(a1, a2, a3), (a1m*a2m*a3m + a3m*(a1^a2) + a2m*(a1^a3) + a1m*(a2^a3))));

	return res;
}

__host__ __device__ Point intHexTr(const Point &p0, const HexahedronWid &h) {
	Point sum;
	for (int i = 0; i < 12; ++i) {
		const auto tri = h.getTri(i);
		sum += intTrAn(p0, tri) * (tri.normal() ^ h.dens);
	}
	return sum;
}

__host__ __device__ Point mMassPoint(const Point &r, const Point &m) {
	const double d = r.eqNorm();
	return r * 3 * (r^m) / (d*d*d*d*d) - m / (d*d*d);
}

bool cuSolver::isCUDAavailable() {
	int deviceCount;
	cudaError_t e = cudaGetDeviceCount(&deviceCount);
	return e == cudaSuccess && deviceCount > 0;
}
int cuSolver::getGPUnum() {
	int deviceCount;
	cudaError_t e = cudaGetDeviceCount(&deviceCount);
	return e == cudaSuccess ? deviceCount : -1;
}

void cuSolver::setDevice(const int id) {
	CUDA_CHECK_RETURN(cudaSetDevice(id));
}

/* -- TRANS SOLVER --
class TransSolver : public gFieldInvSolver {
public:
	constexpr static const int parts = 4;

	TransSolver(const std::vector<FieldPoint> &fps, const double dotPotentialRad) : dotPotentialRad(dotPotentialRad) {
		for (int d = 0; d < parts; ++d) {
			const size_t from = blockOffset(fps.size(), parts, d);
			const size_t to = blockOffset(fps.size(), parts, d + 1);
			fpCUDA[d] = new thrust::device_vector<FieldPoint>(&*fps.cbegin() + from, &*fps.cbegin() + to);
		}
	}
	~TransSolver() {
		for (int d = 0; d < parts; ++d)
			delete fpCUDA[d];
	}

	void solve(const HexahedronWid * const qbegin, const HexahedronWid * const qend, const std::vector<double>::iterator &resBegin) override {
		const size_t taskSz = qend - qbegin;
		std::vector<MassPoint> qtmp(taskSz);
		transform(qbegin, qend, qtmp.begin(), [](auto &h) {return h.getMassPoint(); });
		std::mutex resMt;

		//copy to device
		const thrust::device_vector<MassPoint> mps(&*qtmp.cbegin(), &*qtmp.cend());
		const thrust::device_vector<HexahedronWid> hs(qbegin, qend);

		int dev = 0;
		CUDA_CHECK_RETURN(cudaGetDevice(&dev));

#pragma omp parallel for num_threads(parts)
		for (int d = 0; d < parts; ++d) {
			CUDA_CHECK_RETURN(cudaSetDevice(dev));
			thrust::device_vector<double> rs(taskSz);

			//solve
			const FieldPoint* fps = fpCUDA[d]->data().get();
			const size_t pnSz = fpCUDA[d]->size();
			const double rad = dotPotentialRad;
			thrust::transform(mps.cbegin(), mps.cend(), hs.begin(), rs.begin(),
				[=] __device__(const MassPoint& mp, const HexahedronWid &h) -> double {
				double sum = 0;
				for (size_t i = 0; i < pnSz; ++i) {
					const FieldPoint &fp = fps[i];
					if ((mp - fp.p).eqNorm() <= rad)
						sum += -fp.v*intHexTr(fp.p, fp.n, h);
					else sum += gMassPoint(fp.p, fp.n, fp.v, mp);
				}
				return sum;
			});
			//copy back to host explicitly
			thrust::host_vector<double> rsHost(rs);

			//reduce result back
			resMt.lock();
			std::transform(resBegin, resBegin + taskSz, rsHost.cbegin(), resBegin, std::plus<double>());
			resMt.unlock();
		}
	}
private:
	thrust::device_vector<FieldPoint> *fpCUDA[parts];
	const double dotPotentialRad;
	static size_t blockOffset(const size_t size, const size_t am, const size_t i) { //returns start of i'th block out of am
		const size_t b = size / am;
		if (i < am) return b*i;
		return size;
	}
};

std::unique_ptr<gFieldInvSolver> gFieldInvSolver::getCUDAtransSolver(const std::vector<FieldPoint> &fps, const double dotPotentialRad) {
	return std::unique_ptr<gFieldInvSolver>(new TransSolver(fps, dotPotentialRad));
}
*/

class gFieldCUDAsolver : public gFieldSolver {
public:

	gFieldCUDAsolver(const HexahedronWid* const qbegin, const HexahedronWid* const qend,
		const double dotPotentialRad, const int tirBufSz) : dotPotentialRad(dotPotentialRad) {
		qsCUDA.assign(qbegin, qend);
		qsCUDAprec.resize(tirBufSz);
		std::vector<PointValue<Point>> tmp(qend - qbegin);
		transform(qbegin, qend, tmp.begin(), [](auto &h) {return h.getMassPoint(); });
		mpsCUDA.assign(&*tmp.cbegin(), &*tmp.cend());
	}

	~gFieldCUDAsolver() {
		cudaDeviceSynchronize();
	}

	Point solve(const Point &p0) override {
		Point res;

		//rough computing
		triSz = 0;					//reset precise elemets counter
		int* const cnt = &triSz;	//precise elemets counter
		HexahedronWid * const hPres = qsCUDAprec.data().get(); //precise elemets buffer
		const double rad = dotPotentialRad;
		res += thrust::inner_product(qsCUDA.begin(), qsCUDA.end(), mpsCUDA.begin(), Point(), thrust::plus<Point>(),
			[=] __device__(const HexahedronWid& h, const PointValue<Point> &mp)->Point {
			if ((mp - p0).eqNorm() > rad)
				return mMassPoint(p0-mp, mp.val);
			hPres[atomicAdd(cnt, 1)] = h;
			return Point();
		});
		
		//precise computing
		const int blockSize = triSz;	//precise elemets buffer size
		if (blockSize == 0) return res;
		const auto& triKr = [=] __device__ __host__(const HexahedronWid &h)->Point {
			return intHexTr(p0, h);
		};
		const auto& triClac = [&](const auto &execPol) {
			return thrust::transform_reduce(execPol, qsCUDAprec.begin(), qsCUDAprec.begin() + blockSize, triKr, Point(), thrust::plus<Point>());
		};
		res += blockSize > 100 ? triClac(thrust::device) : triClac(thrust::host);
		
		return res;
	}

private:
	thrust::device_vector<HexahedronWid> qsCUDA;			//ro
	thrust::device_vector<PointValue<Point>> mpsCUDA;				//ro
	thrust::device_vector<HexahedronWid> qsCUDAprec;		//rw
	cuVar<int> triSz;										//rw
	const double dotPotentialRad;
};

std::unique_ptr<gFieldSolver> gFieldSolver::getCUDAsolver(const HexahedronWid* const qbegin, const HexahedronWid* const qend,
	const double dotPotentialRad, const int tirBufSz) {
	return std::unique_ptr<gFieldSolver>(new gFieldCUDAsolver(qbegin, qend, dotPotentialRad, tirBufSz));
}

static void CheckCudaErrorAux(const char *file, unsigned line, const char *statement, cudaError_t err)
{
	if (err == cudaSuccess)
		return;
	std::cerr << statement << " returned " << cudaGetErrorString(err) << "(" << err << ") at " << file << ":" << line << std::endl;
	exit(1);
}
