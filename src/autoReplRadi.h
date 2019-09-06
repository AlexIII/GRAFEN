#pragma once
#include <cmath>
#include <cstdio>

class AutoReplRadi {
public:
	//Рассчитать расстояние от центра однородного параллелепипеда (с размерами Lx, Ly, Lz), на котором относительная погрешность 
	//в z-компоненте напряженности гравитационного поля
	//при замене параллелепипеда на точку равной массы (со сдвигом координат относительно центра px, py, pz) будет не больше eps.
	static double get(const double eps, const double Lx, const double Ly, const double Lz, const double px = 0, const double py = 0, const double pz = 0)
	{
		int i;
		double z, d, z2, gp, gpd, gc, gcd;

		const double sigma = 1.0;
		const double N = 100;
		const double EPS = 1e-6;
		const double dx = 0.5*Lx;
		const double dy = 0.5*Ly;
		const double dz = 0.5*Lz;
		const double z0base = dz + 0.00001;
		const double m = sigma * Lx*Ly*Lz;
		const double sigmad = sigma + sigma;
		const double sigmaq = sigmad + sigmad;
		const double px2py2 = px * px + py * py;
		const double dxdy = dx * dy;
		const double dx2dy2 = dx * dx + dy * dy;

		auto gp_f_der1 = [&](const double z)
		{
			const double cz = pz - z;
			const double cz2 = cz * cz;
			const double t = 1.0 / (px2py2 + cz2);
			const double t2 = t * sqrt(t);
			gp = -m * cz*t2;
			gpd = -m * (cz2 + cz2 - px2py2)*t*t2;
		};

		auto gc_f_der1 = [&](const double z)
		{
			const double czp = z - dz;
			const double rp = sqrt(dx2dy2 + czp * czp);
			const double ap = atan2(dxdy, czp*rp);
			const double czm = z + dz;
			const double rm = sqrt(dx2dy2 + czm * czm);
			const double am = atan2(dxdy, czm*rm);
			gc = sigmad * (dx*log((rp + dy)*(rm - dy) / ((rp - dy)*(rm + dy))) + dy * log((rp + dx)*(rm - dx) / ((rp - dx)*(rm + dx))) + 2.0*czm*am - 2.0*czp*ap);
			gcd = sigmaq * (am - ap);
		};

		z = z0base;
		d = EPS + 1.0;
		//printf("e=%.15g, z[0]=%.15g\n", eps, z);
		for (i = 0; i < N && d > EPS; ++i)
		{
			z2 = z;
			gp_f_der1(z);
			gc_f_der1(z);
			z += gc * (gp - (1.0 + eps)*gc) / (gp*gcd - gc*gpd);
			d = fabs(z - z2);
			//printf("z[%2d]=%.15g, diff_iter=%.15g\n", i + 1, z, d);
		}
		return z;
	}
private:
	AutoReplRadi();
};