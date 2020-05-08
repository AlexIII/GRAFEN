#pragma once

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <functional>
using namespace std;

#include <math.h>

#define __int32 int

class Grid
{
private:
	__int32 ReadInt32(ifstream* fs);
	double ReadDouble(ifstream* fs);
	void WriteInt32(ofstream* fs, __int32 value);
	void WriteDouble(ofstream* fs, double value);

	bool Init();

public:
		vector<double> data;
		int nRow;
		int nCol;
		double xLL;
		double yLL;
		double xSize;
		double ySize;
		double zMin;
		double zMax;
		double Rotation;
		double BlankValue;
	string fname;

public:
	Grid();
	Grid(const string &fileName);

	static Grid GenerateEmptyGrid(Grid& grid);

	double getAverage();
	double getMin();
	double getMax();

	const double& at(const int col, const int row) const {
		return data[row*nCol + col];
	}
	double& at(const int col, const int row) {
		return data[row*nCol + col];
	}
	double xAt(const int col) const {
		return xLL + col * xSize;
	}
	double yAt(const int row) const {
		return yLL + row * ySize;
	}
	double width() const {
		return (nCol - 1) * xSize;
	}
	double height() const {
		return (nRow - 1) * ySize;
	}
	double mean() const {
		int count = 0;
		double mean = 0;
		forEach([&BlankValue = BlankValue, &count, &mean](int, int, const double &v) {
			if (v == BlankValue) return;
			mean += v;
			++count;
		});
		return count? mean / count : 0;
	}
	void forEach(const function<void(int, int, double&)>& f) {
		for (int i = 0; i < nCol; ++i)
			for (int j = 0; j < nRow; ++j)
				f(i, j, at(i, j));
	}
	void forEach(const function<void(int, int, const double&)>& f) const {
		for (int i = 0; i < nCol; ++i)
			for (int j = 0; j < nRow; ++j)
				f(i, j, at(i, j));
	}
	void setBlanksTo(const double bv) {
		forEach([&BlankValue = BlankValue, &bv](int, int, double& val) {
			if (val == BlankValue) val = bv;
		});
	}

	Grid& upScale(const int xFactor, const int yFactor) {
		upScaleX(xFactor);
		upScaleY(yFactor);
		return *this;
	}
	void upScaleX(const int factor) {
		vector<double> newData(nRow*nCol*factor);
		for (int i = 0; i < nRow*nCol; ++i)
			for (int k = 0; k < factor; ++k)
				newData[i*factor+k] = data[i];
		data = newData;
		nCol *= factor;
		xSize /= double(factor);
	}
	void upScaleY(const int factor) {
		vector<double> newData(nRow*nCol*factor);
		for (int i = 0; i < nRow; ++i)
			for (int j = 0; j < nCol; ++j)
				for (int k = 0; k < factor; ++k)
					newData[i*nCol*factor + k*nCol + j] = data[i*nCol + j];
		data = newData;
		nRow *= factor;
		ySize /= double(factor);
	}

	bool Read(const string &fileName);
	bool Read();
	bool Write(const string &fileName);
	bool Write();
};

std::ostream& operator<<(std::ostream& os, const Grid& g);
