/*
 * calcFieldFun.h
 *
 *  Created on: 01 апр. 2016 г.
 *      Author: user
 */

#ifndef CALCFIELDFUN_H_
#define CALCFIELDFUN_H_

#include <iostream>
#include <stdlib.h>
#include "mobj.h"
#include <vector>
#include <algorithm>
#include <memory>

class cuSolver {
public:
	static bool isCUDAavailable();
	static int getGPUnum();
	static void setDevice(const int id);
};
/* -- TRANS SOLVER --
class gFieldInvSolver : public cuSolver {
public:
	virtual void solve(const HexahedronWid* const qbegin, const HexahedronWid* const qend, const std::vector<double>::iterator &resBegin) = 0;
	static std::unique_ptr<gFieldInvSolver> getCUDAtransSolver(const std::vector<FieldPoint> &fps, const double dotPotentialRad);
	virtual ~gFieldInvSolver() {}
};
*/

/*
	Template interface for ClosedShape
	struct ClosedShape {
		static constexpr int nTriangles; //number of triangles
		CUDA_HOST_DEV_FUN Triangle<double> getTri(const int i) const;
		PointValue<double> getMassPoint() const;
	};
*/


class gFieldSolver : public cuSolver {
public:
	virtual Point solve(const Point &p0) = 0;
	template<class ClosedShape>
	static std::unique_ptr<gFieldSolver> getCUDAsolver(const ClosedShape* const qbegin, const ClosedShape* const qend, const bool transpose = false);
	template<class ClosedShape>
	static std::unique_ptr<gFieldSolver> getCUDAreplacingSolver(const ClosedShape* const qbegin, const ClosedShape* const qend,
		const double dotPotentialRad, const int tirBufSz);
	virtual ~gFieldSolver() {}
protected:
	gFieldSolver(const bool transpose = false) : transpose(transpose) {}
	const bool transpose;
};


#endif /* CALCFIELDFUN_H_ */
