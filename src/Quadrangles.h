#pragma once
#include "mobj.h"
#include <vector>
#include <string>

template<typename T, class _Alloc = std::allocator<T>>
class gElementsTyped : public std::vector<T, _Alloc> {
public:
	using base = T;
	void save(const std::string fname) {
		std::fstream file(fname, std::fstream::out | std::fstream::binary);
		if (file.bad()) {
			file.close();
			throw std::runtime_error("Can't open file \"" + fname + "\"");
		}
		const size_t sz = std::vector<T>::size();
		file.write((char*)&sz, sizeof(size_t));
		file.write((char*)std::vector<T>::data(), sz * sizeof(base));
		file.flush();
		file.close();
	}
	void load(const std::string fname) {
		std::fstream file(fname, std::fstream::in | std::fstream::binary);
		if (file.bad()) {
			file.close();
			throw std::runtime_error("Can't open file \"" + fname + "\"");
		}
		size_t sz;
		file.read((char*)&sz, sizeof(size_t));
		std::vector<T, _Alloc>::resize(sz);
		file.read((char*)std::vector<T>::data(), sz * sizeof(base));
		file.close();
	}
};

using gElements = gElementsTyped<HexahedronWid>;
