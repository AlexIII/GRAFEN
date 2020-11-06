#include "calcField.h"
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/copy.h>
#include <thrust/inner_product.h>
#include <thrust/transform_reduce.h>
#include <thrust/execution_policy.h>
#include "cuVar.h"
#include <mutex>
#include <array>

using KernelComputeType = double; //float;

using std::cout;
using std::endl;

static void CheckCudaErrorAux(const char *, unsigned, const char *, cudaError_t);
#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)

template<typename T>
__host__ __device__ static inline double tripleprod(const Point3D<T>& o1, const Point3D<T>& o2, const Point3D<T>& o3) { return (o1 ^ (o2*o3)); }

template<typename T>
__host__ __device__ static Point intTrAn(const Point3D<T> &p0, const Triangle<T> &t) {
	Point3D<T> res;

	const Point3D<T> a1 = t.p1 - p0;
	const Point3D<T> a2 = t.p2 - p0;
	const Point3D<T> a3 = t.p3 - p0;
	const T a1m = a1.eqNorm();
	const T a2m = a2.eqNorm();
	const T a3m = a3.eqNorm();
	const Point3D<T> a12 = t.p2 - t.p1;
	const Point3D<T> a23 = t.p3 - t.p2;
	const Point3D<T> a31 = t.p1 - t.p3;
	const T a12m = a12.eqNorm();
	const T a23m = a23.eqNorm();
	const T a31m = a31.eqNorm();
	//if (a3m + a1m - a31m < 1e-9 || a1m + a2m - a12m < 1e-9 || a2m + a3m - a23m < 1e-9) return{ HUGE_VAL, HUGE_VAL, HUGE_VAL };
	const Point3D<T> N = t.normal();

	res = a31*(log((a3m + a1m + a31m) / (a3m + a1m - a31m)) / a31m);
	res += a12*(log((a1m + a2m + a12m) / (a1m + a2m - a12m)) / a12m);
	res += a23*(log((a2m + a3m + a23m) / (a2m + a3m - a23m)) / a23m);
	res = N*res;
	res += N*(2.0*atan2(tripleprod(a1, a2, a3), (a1m*a2m*a3m + a3m*(a1^a2) + a2m*(a1^a3) + a1m*(a2^a3))));

	//return res;
	return Point(res);
}

Point intHexTr__(const Point &p0, const HexahedronWid &h) {
	Point sum;
	for (int i = 0; i < 12; ++i) {
		const auto tri = h.getTri(i);
		sum += intTrAn<double>(p0, tri) * (tri.normal() ^ h.dens);
	}
	return sum;
}

template<typename T>
__host__ __device__ static Point intHexTr(const Point &p0, const HexahedronWid &h) {
	Point3D<T> sum;
	const Point3D<T> p0f(p0);
	const Point3D<T> densf(h.dens);
	for (int i = 0; i < 12; ++i) {
		const Triangle<T> tri(h.getTri(i));
		sum += intTrAn<T>(p0f, tri) * (tri.normal() ^ densf);
	}
	return sum;
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
	using dvHex = thrust::device_vector<HexahedronWid>;
	gFieldCUDAsolver(const HexahedronWid* const qbegin, const HexahedronWid* const qend) {
		qsCUDA.assign(qbegin, qend);
	}
	virtual ~gFieldCUDAsolver() {}
	virtual Point solve(const Point &p0) override {
		return solve(qsCUDA.cbegin(), qsCUDA.cend(), p0);
	}

	static Point solve(const dvHex::const_iterator &qbegin, const dvHex::const_iterator &qend, const Point &p0) {
		Point res;
		const auto& triKr = [=] __device__ __host__(const HexahedronWid &h)->Point {
			return intHexTr<KernelComputeType>(p0, h);
		};
		const auto& triClac = [&](const auto &execPol) {
			return thrust::transform_reduce(execPol, qbegin, qend, triKr, Point(), thrust::plus<Point>());
		};
		res += (qend - qbegin) > 100 ? triClac(thrust::device) : triClac(thrust::host);

		return res;
	}
protected:
	dvHex qsCUDA;			//ro
};

class gFieldCUDAreplacingSolver : public gFieldCUDAsolver {
public:
	gFieldCUDAreplacingSolver(const HexahedronWid* const qbegin, const HexahedronWid* const qend, const double dotPotentialRad, 
			const int tirBufSz) : gFieldCUDAsolver(qbegin, qend), dotPotentialRad(dotPotentialRad) {
		qsCUDAprec.resize(tirBufSz);
		std::vector<std::array<MagLine<KernelComputeType>, 3>> tmp(qend - qbegin);
		transform(qbegin, qend, tmp.begin(), [](auto &h) {return h.template getLines<KernelComputeType>(); });
		linesCUDA.assign(&*tmp.cbegin(), &*tmp.cend());
	}

	virtual ~gFieldCUDAreplacingSolver() {
		cudaDeviceSynchronize();
	}

	virtual Point solve(const Point &p0) override {
		Point res;

		//rough computing
		triSz = 0;					//reset precise elemets counter
		int* const cnt = &triSz;	//precise elemets counter
		HexahedronWid * const hPres = qsCUDAprec.data().get(); //precise elemets buffer
		const double rad = dotPotentialRad;
		const Point3D<KernelComputeType> p0f = p0;
		res += thrust::inner_product(qsCUDA.begin(), qsCUDA.end(), linesCUDA.begin(), Point(), thrust::plus<Point>(),
			[=] __device__(const HexahedronWid& h, const std::array<MagLine<KernelComputeType>, 3> &ml) -> Point {
			if ((ml[0].p1 - p0f).eqNorm() > rad)
				return ml[0].Hfield(p0f) + ml[1].Hfield(p0f) + ml[2].Hfield(p0f);
			hPres[atomicAdd(cnt, 1)] = h;
			return Point();
		});
		
		//precise computing
		const int blockSize = triSz;	//precise elemets buffer size
		if (blockSize == 0) return res;
		res += gFieldCUDAsolver::solve(qsCUDAprec.cbegin(), qsCUDAprec.cbegin() + blockSize, p0);

		return res;
	}

private:
	thrust::device_vector<std::array<MagLine<KernelComputeType>, 3>> linesCUDA;		//ro
	dvHex qsCUDAprec;										//rw
	cuVar<int> triSz;										//rw
	const double dotPotentialRad;
};

std::unique_ptr<gFieldSolver> gFieldSolver::getCUDAreplacingSolver(const HexahedronWid* const qbegin, const HexahedronWid* const qend,
	const double dotPotentialRad, const int tirBufSz) {
	return std::unique_ptr<gFieldSolver>(new gFieldCUDAreplacingSolver(qbegin, qend, dotPotentialRad, tirBufSz));
}
std::unique_ptr<gFieldSolver> gFieldSolver::getCUDAsolver(const HexahedronWid* const qbegin, const HexahedronWid* const qend) {
	return std::unique_ptr<gFieldSolver>(new gFieldCUDAsolver(qbegin, qend));
}

static void CheckCudaErrorAux(const char *file, unsigned line, const char *statement, cudaError_t err)
{
	if (err == cudaSuccess)
		return;
	std::cerr << statement << " returned " << cudaGetErrorString(err) << "(" << err << ") at " << file << ":" << line << std::endl;
	exit(1);
}
