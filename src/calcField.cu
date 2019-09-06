#include "calcField.h"
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/copy.h>
#include <thrust/inner_product.h>
#include <thrust/transform_reduce.h>
#include <thrust/execution_policy.h>
#include "cuVar.h"
#include <mutex>

#define BLOCK_SZ 32
#define MAX_DEV 16
static void CheckCudaErrorAux(const char *, unsigned, const char *, cudaError_t);
#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)

__host__ __device__ double intTrAn(const Point &q, const Triangle &t) {
	const Point a1 = t.p1 - q, a2 = t.p2 - q, a3 = t.p3 - q;
	const Point a12 = t.p2 - t.p1, a23 = t.p3 - t.p2, a31 = t.p1 - t.p3;
	const Point k = (a12*a23);
	const Point K = k.norm();

	const double pr = -M_PI_2*abs(a1^K);

	const Point A31 = a31.norm();
	const Point A12 = a12.norm();
	const Point A23 = a23.norm();

	double tmp1 = (A31^a3) + sqrt(a3^a3);
	double uln1 = ((A31^a1) + sqrt(a1^a1)) / tmp1;
	if (tmp1 > 1e-10 && uln1 >= 0) {
		uln1 = (a1 ^ (A31*K))*log(uln1);
	}
	else {
		uln1 = (-(A31^a1) + sqrt(a1^a1)) / (-(A31^a3) + sqrt(a3^a3));
		uln1 = -(a1 ^ (A31*K))*log(uln1);
	}

	double tmp2 = (A12^a1) + sqrt(a1^a1);
	double uln2 = ((A12^a2) + sqrt(a2^a2)) / tmp2;
	if (tmp2 > 1e-10 && uln2 >= 0) {
		uln2 = (a2 ^ (A12*K))*log(uln2);
	}
	else {
		uln2 = (-(A12^a2) + sqrt(a2^a2)) / (-(A12^a1) + sqrt(a1^a1));
		uln2 = -(a2 ^ (A12*K))*log(uln2);
	}

	double tmp3 = (A23^a2) + sqrt(a2^a2);
	double uln3 = ((A23^a3) + sqrt(a3^a3)) / tmp3;
	if (tmp3 > 1e-10 && uln3 >= 0) {
		uln3 = (a3 ^ (A23*K))*log(uln3);
	}
	else {
		uln3 = (-(A23^a3) + sqrt(a3^a3)) / (-(A23^a2) + sqrt(a2^a2));
		uln3 = -(a3 ^ (A23*K))*log(uln3);
	}

	if (a1^K) {
		uln1 -= (a1^K)*atan(((a31*a1) ^ (a12*a1)) / ((a1^k)*sqrt(a1^a1)));
		uln2 -= (a2^K)*atan(((a12*a2) ^ (a23*a2)) / ((a2^k)*sqrt(a2^a2)));
		uln3 -= (a3^K)*atan(((a23*a3) ^ (a31*a3)) / ((a3^k)*sqrt(a3^a3)));
	}

	return pr + uln1 + uln2 + uln3;
}

__host__ __device__ double gMassPoint(const Point &p0, const Point &n0, const double mass, const Point &mp) {
	const double d = (mp - p0).eqNorm();
	return mass*(n0 ^ (mp - p0)) / (d*d*d);
}

__host__ __device__ double intHexTr(const Point &p0, const Point &n0, const HexahedronWid &h) {
	double sum = 0;
	for (int i = 0; i < 12; ++i)
		sum += intTrAn(p0, h.getTri(i)) * (n0^h.getTriNorm(i));
	return sum;
}

double gFieldTri(const Point p0, const Point n0, const HexahedronWid &h) {
	return -h.dens*intHexTr(p0, n0, h);
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
		vector<MassPoint> qtmp(taskSz);
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
	int devAm;
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


class gFieldCUDAsolver : public gFieldSolver {
public:

	gFieldCUDAsolver(const HexahedronWid* const qbegin, const HexahedronWid* const qend,
		const double dotPotentialRad, const int tirBufSz) : dotPotentialRad(dotPotentialRad) {
		qsCUDA.assign(qbegin, qend);
		qsCUDAprec.resize(tirBufSz);
		vector<MassPoint> tmp(qend - qbegin);
		transform(qbegin, qend, tmp.begin(), [](auto &h) {return h.getMassPoint(); });
		mpsCUDA.assign(&*tmp.cbegin(), &*tmp.cend());
	}

	~gFieldCUDAsolver() {
		cudaDeviceSynchronize();
	}

	double solve(const Point &p0, const Point &n0) override {
		double res = 0;

		//rough computing
		triSz = 0;					//reset precise elemets counter
		int* const cnt = &triSz;	//precise elemets counter
		HexahedronWid * const hPres = qsCUDAprec.data().get(); //precise elemets buffer
		const double rad = dotPotentialRad;
		res += thrust::inner_product(qsCUDA.begin(), qsCUDA.end(), mpsCUDA.begin(), 0., thrust::plus<double>(),
			[=] __device__(const HexahedronWid& h, const MassPoint &mp)->double {
			if ((mp - p0).eqNorm() > rad)
				return gMassPoint(p0, n0, mp.mass, mp);
			hPres[atomicAdd(cnt, 1)] = h;
			return 0;
		});

		//precise computing
		const int blockSize = triSz;	//precise elemets buffer size
		if (blockSize == 0) return res;
		const auto& triKr = [=] __device__ __host__(const HexahedronWid &h)->double {
			return -h.dens*intHexTr(p0, n0, h);
		};
		const auto& triClac = [&](const auto &execPol) {
			return thrust::transform_reduce(execPol, qsCUDAprec.begin(), qsCUDAprec.begin() + blockSize, triKr, 0., thrust::plus<double>());
		};
		res += blockSize > 100 ? triClac(thrust::device) : triClac(thrust::host);

		return res;
	}

private:
	thrust::device_vector<HexahedronWid> qsCUDA;			//ro
	thrust::device_vector<MassPoint> mpsCUDA;				//ro
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