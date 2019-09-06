
#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <chrono>

class Stopwatch {
    public:
		Stopwatch() { start(); }
        ~Stopwatch(){}

        void start(){
        	sp = std::chrono::steady_clock::now();
        }

        double stop(){
        	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        	const double us = (double)std::chrono::duration_cast<std::chrono::microseconds>(end - sp).count();
            return us / 1000000.;
        }

    private:
       std::chrono::steady_clock::time_point sp;
};

#endif
