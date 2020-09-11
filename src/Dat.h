/*
 * Dat.h
 *
 *  Created on: 29 Jun 2016
 *      Author: osboxes
 */

#ifndef DAT_H_
#define DAT_H_

#include <vector>
#include <fstream>
#include <iomanip>
#include <functional>
#include <sstream>
#include <algorithm>
#include <type_traits>
#include "Stopwatch.h"

struct DatPoint2D {
	double x = 0;
	double y = 0;
};

std::ostream& operator<<(std::ostream& os, const DatPoint2D& p) {
	os << p.x << " " << p.y;
	return os;
}

std::istream& operator>>(std::istream& is, DatPoint2D& p) {
	is >> p.x >> p.y;
	return is;
}

struct DatPoint3D : public DatPoint2D {
	DatPoint3D() {}
	DatPoint3D(const double &x, const double &y, const double &z) : DatPoint2D({x, y}), z(z) {}
	double z = 0;
};

std::ostream& operator<<(std::ostream& os, const DatPoint3D& p) {
	os << p.x << " " << p.y << " " << p.z;
	return os;
}

std::istream& operator>>(std::istream& is, DatPoint3D& p) {
	is >> p.x >> p.y >> p.z;
	return is;
}

template<typename PointType, typename DataType>
class Dat {
public:
	using Point = PointType;
	struct Element {
		Point p{};
		DataType val{};
	};
	std::vector<Element> es;
	std::string fileName;

	template<typename Tval = double>
	Dat(std::enable_if_t<std::is_same_v<PointType, DatPoint3D>, const Dat<DatPoint2D, Tval>> dat2D, const double z) {
		set(dat2D, z);
	}

	template<typename Tval = double>
	void set(std::enable_if_t<std::is_same_v<PointType, DatPoint3D>, const Dat<DatPoint2D, Tval> &> dat2D, const double z) {
		es.resize(dat2D.es.size());
		transform(dat2D.es.begin(), dat2D.es.end(), es.begin(), [&z](const auto &e) {return Element{ Point{e.p.x, e.p.y, z}, e.val }; });
		fileName = dat2D.fileName;
	}

	std::vector<Point> getPoints() const {
		std::vector<Point> ps(es.size());
		std::transform(es.begin(), es.end(), ps.begin(), [](const Element &e) {return e.p; });
		return ps;
	}

	void set(const DataType val) {
		for(auto &e : es)
			e.val = val;
	}

	void set(const std::vector<DataType> val) {
		for (size_t i = 0; i < es.size(); ++i)
			es[i].val = val[i];
	}

	void add(const std::vector<DataType> val) {
		for (size_t i = 0; i < es.size(); ++i)
			es[i].val += val[i];
	}

	double xMax() const {
		return std::max_element(es.begin() ,es.end() , [](const Element &e1, const Element &e2)->bool{return e1.p.x < e2.p.x;})->p.x;
	}

	double xMin() const {
		return std::min_element(es.begin() ,es.end() , [](const Element &e1, const Element &e2)->bool{return e1.p.x < e2.p.x;})->p.x;
	}

	double yMax() const {
		return std::max_element(es.begin() ,es.end() , [](const Element &e1, const Element &e2)->bool{return e1.p.y < e2.p.y;})->p.y;
	}

	double yMin() const {
		return std::min_element(es.begin() ,es.end() , [](const Element &e1, const Element &e2)->bool{return e1.p.y < e2.p.y;})->p.y;
	}

	Dat(const size_t sz = 0) : es(sz) {
		fileName = "def.dat";
	}

	Dat(const std::string fname) {
		fileName = "def.dat";
		read(fname);
	}

	Dat(typename std::vector<Element>::iterator &from, typename std::vector<Element>::iterator &to) : es(from,to) {
		fileName = "def.dat";
	}

	typename std::vector<Element>::iterator begin() {
		return es.begin();
	}

	typename std::vector<Element>::iterator end() {
		return es.end();
	}

	void read(const std::string fname) {
		fileName = fname;

		std::ifstream file(fileName);
		if(!file.is_open()) throw std::runtime_error("Can't open *.dat file: "+fname);
		Element el;
		char buff[256];
		while(!file.eof()) {
			if(!file.getline(buff, sizeof(buff))) break;
			if(!std::string(buff).length()) break;
			std::stringstream(buff) >> el.p >> el.val;
			es.push_back(el);
		}
		file.close();
	}

	void write(const bool as2D = false) {
		write(fileName, as2D);
	}

	void write(const char* fname, const bool as2D = false) {
		write(std::string(fname), as2D);
	}
	void write(const std::string fname, const bool as2D = false) {
		fileName = fname;
		std::ofstream outDatFile(fileName);
		outDatFile << std::setprecision(12);
		if(!as2D) for(size_t i = 0; i < es.size(); ++i)
			outDatFile << es[i].p << " " << es[i].val << "\n";
		else for (size_t i = 0; i < es.size(); ++i)
			outDatFile << ((DatPoint2D&)es[i].p) << " " << es[i].val << "\n";
		outDatFile << std::flush;
		outDatFile.close();
	}

	Element& operator[](const int &i) {
		return es[i];
	}
	const Element& operator[](const int &i) const {
		return es[i];
	}

	size_t size(void) {
		return es.size();
	}

	void forEach(std::function<void (Element &)> fun) {
		for(size_t i = 0; i < es.size(); ++i)
			fun(es[i]);
	}
};

template <typename DataType = double>
using Dat2D = Dat<DatPoint2D, DataType>;
template <typename DataType = double>
using Dat3D = Dat<DatPoint3D, DataType>;


#endif /* DAT_H_ */
