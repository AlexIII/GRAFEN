#pragma once

template<typename T>
class cuVar {
public:
	cuVar() {
		cudaMalloc(&p, sizeof(T));
	}
	~cuVar() {
		cudaFree(p);
	}
	cuVar(const T &v) : cuVar() {
		set(v);
	}
	T* raw() const { return (T*)p; }
	operator T () const {
		return get();
	}
	cuVar<T>& operator=(const cuVar<T> &v) {
		cudaMemcpy(p, &v, sizeof(T), cudaMemcpyDeviceToDevice);
		return *this;
	}
	cuVar<T>& operator=(const T &v) {
		set(v);
		return *this;
	}
	T get() const {
		T v;
		get(v);
		return v;
	}
private:
	void* p;
	void set(const T &v) const { cudaMemcpy(p, &v, sizeof(T), cudaMemcpyHostToDevice); }
	void get(T &v) const { cudaMemcpy(&v, p, sizeof(T), cudaMemcpyDeviceToHost); }
};
