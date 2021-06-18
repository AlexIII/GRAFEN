#ifndef TOPOGRAVARGS_H_
#define TOPOGRAVARGS_H_

#include <iostream>
#include <string>
#include <vector>
#include <boost/optional.hpp>
#include "mobj.h"
#include "inputParser.h"
#include <sstream>

using std::cout;
using std::endl;
using std::istringstream;

#define DEF_R_EQ 6378.245		//equatorial radius in km
#define DEF_R_PL 6356.863		//polar radius in km

class TopogravArgs {
public:
	//input options
	std::string topoGridFname; 						// -topoGrd7 *string*
	std::string gravGridFname; 						// -gravGrd7 *string*
	double dens; 									// -dens *number*
	Ellipsoid refEllipsoid{ DEF_R_EQ, DEF_R_PL }; 	// -Rpol *number* -Req *number*
	double pprr = -1;								// -DPR *number*
														// < 0 - exact calculations
														// == 0 - approximate calculations
														// > 0 - point-potential replacement radius
	boost::optional<Point> normal; 					// -nx *number* -ny *number* -nz *number* 
	bool flatMode = false;							// -flat don't use spherical model
	double l0 = 0;									// -l0 *in deg* (converted to rad)
														// (for flat mode and if grids in GK)
	std::vector<int> gpuIdMap{};					// -gpuIdMap 0,2,4
	bool gridsInGK = false;							// -gridsInGK

	TopogravArgs(int argc, char *argv[]) {
		InputParser ip(argc, argv);

		gridsInGK = ip.exists("gridsInGK");
		if(gridsInGK) {
			ip["l0"] >> l0;
			l0 = toRad(l0);
		}

		if(ip.exists("gpuIdMap")) {
			std::string s;
			ip["gpuIdMap"] >> s;
			std::string segment;
			std::stringstream ss(s);
			while(std::getline(ss, segment, ','))
				gpuIdMap.push_back(std::stoi(segment));
		}

		ip["topoGrd7"] >> topoGridFname;
		ip["gravGrd7"] >> gravGridFname;
		ip["dens"] >> dens;

		if(ip.exists("Req") && ip.exists("Rpol")) {
			double Req, Rpol;
			ip["Req"] >> Req;
			ip["Rpol"] >> Rpol;
			refEllipsoid = Ellipsoid(Req, Rpol);
		}
		else {
			cout << "Using default reference ellipsoid with Rpol=" << refEllipsoid.Rpl << " and Req=" << refEllipsoid.Req << endl;
		}

		if(ip.exists("DPR"))
			ip["DPR"] >> pprr;
		else {
			//dotPotentialRad = AutoReplRadi::get(1e-3, Elim.dWh(), Nlim.dWh(), Hlim.d());
			//cout << "Deduced dot potential replace radius: " << dotPotentialRad << endl;
		}

		if(ip.exists("nx") || ip.exists("ny") || ip.exists("nz")) {
			Point n;
			if (ip.exists("nx")) ip["nx"] >> n.x;
			if (ip.exists("ny")) ip["ny"] >> n.y;
			if (ip.exists("nz")) ip["nz"] >> n.z;
			cout << "Using specified normal: " << n << endl;
			normal = n;
		}

		if(ip.exists("flat")) {
			ip["l0"] >> l0;
			l0 = toRad(l0);
			flatMode = true;
		}
	}
};

#endif 
