#pragma once

#include "Grid/Grid.h"
#include "Dat.h"

class GDconv {
public:
	static Dat2D toDat(const Grid &g) {
		Dat2D dat;
		for (int y = 0; y < g.nRow; ++y)
			for (int x = 0; x < g.nCol; ++x)
				dat.es.push_back({
					{ g.xLL + x*g.xSize, g.yLL + y*g.ySize },
					g.data[y*g.nCol+x]
				});
		return dat;
	}
	static Dat3D toDat3D(const Grid &g, const double val = 0) {
		Dat3D dat;
		for (int y = 0; y < g.nRow; ++y)
			for (int x = 0; x < g.nCol; ++x)
				dat.es.push_back({
					{ g.xLL + x*g.xSize, g.yLL + y*g.ySize, g.data[y*g.nCol+x] },
					val
				});
		return dat;
	}
	template <typename PointType>
	static Grid toGrd(const Dat<PointType>& dat, int nCol, int nRow) {
		Grid g;
		g.nRow = nRow;
		g.nCol = nCol;
		g.xLL = dat.xMin();
		g.yLL = dat.yMin();
		g.xSize = (dat.xMax() - g.xLL) / (nCol-1);
		g.ySize = (dat.yMax() - g.yLL) / (nRow - 1);
		for (int y = 0; y < g.nRow; ++y)
			for (int x = 0; x < g.nCol; ++x)
				g.data.push_back(dat[y*g.nCol + x].val);
		return g;
	}
};