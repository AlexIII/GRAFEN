/*
 * inputLoader.h
 *
 *  Created on: 19 марта 2016 г.
 *      Author: alex
 */

#ifndef INPUTLOADER_H_
#define INPUTLOADER_H_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>
#include "mobj.h"
#include "Grid/Grid.h"
#include "inputParser.h"
#include "Dat.h"
#include "autoReplRadi.h"
#include "GrdDatConv.h"

using std::cout;
using std::endl;
using std::istringstream;

class Input {
public:
	limits Elim;
	limits Nlim;
	limits Hlim;
	double l0;
	std::vector<std::vector<Point>> dens;
	Dat3D<Point> dat;
	bool noGK;
	double dotPotentialRad = 1e8;
	std::vector<std::string> fnames;

	std::string grdFname;
	int grdCols = 0;
	int grdRows = 0;

	bool noInvFileOrder = false;
	bool transSolver = false;

	//input: -dat3D (file.dat) [-Hf (val)] -Hfrom (val) -Hto (val) -Hn (val) -l0 (val) [-dens (directory)]
		//[-qw (filename)] [-qr (filename)] *** write or read preprocessing to or from file
		//[-DPR (val)] *** point potential replace radius
		// [-noInvFileOrder]
		// [-transposeSolver]
		// [-saveDatAs3D]
		//dens:
			// layer0_x.grd
			// layer0_y.grd
			// layer0_z.grd
			// layer1_x.grd
			// layer2_y.grd
			// layer2_z.grd
			// ...
	Input(int argc, char *argv[]) {
		InputParser ip(argc, argv);

		noInvFileOrder = ip.exists("noInvFileOrder");
		transSolver = ip.exists("transposeSolver");

		//Hlim
		ip["Hfrom"] >> Hlim.lower;
		ip["Hto"] >> Hlim.upper;
		ip["Hn"] >> Hlim.n;

		if(ip.exists("dens")) {
			std::string dirName;
			ip["dens"] >> dirName;
			LoadDens(dirName);
			Grid g(fnames[0]);
			Elim = {g.xLL, g.xLL + (g.nCol-1)*g.xSize, g.nCol};
			Nlim = {g.yLL, g.yLL + (g.nRow-1)*g.ySize, g.nRow};
			if(fnames.size() < Hlim.n*3)
				throw runtime_error("Not enough density files have been found");
		} else throw runtime_error("No density specified");

		ip["l0"] >> l0;
		l0 = toRad(l0);

		std::string datFname;
		if (ip.exists("dat3D")) {
			ip["dat3D"] >> datFname;
			dat.read(datFname);
			if (ip.exists("Hf")) {
				double Hf;
				ip["Hf"] >> Hf;
				for (auto &e : dat.es)
					e.p.z = Hf;
			}
		} else throw runtime_error("No *.dat file has been specified.");
		double Hf = 0;

		if (ip.exists("DPR"))
			ip["DPR"] >> dotPotentialRad;
		else {
			dotPotentialRad = AutoReplRadi::get(1e-3, Elim.dWh(), Nlim.dWh(), Hlim.d());
			cout << "Deduced dot potential replace radius: " << dotPotentialRad << endl;
		}

		cout << "Elim: " << Elim << endl;
		cout << "Nlim: " << Nlim << endl;
		cout << "Hlim: " << Hlim << endl;
		cout << "l0: " << l0 << endl;
		cout << "inDat x range: " << dat.xMin() << " " << dat.xMax() << endl;
		cout << "inDat y range: " << dat.yMin() << " " << dat.yMax() << endl;
	}

	void LoadDens(const std::string dirName) {
		getDensFiles(dirName);
		if(!fnames.size()) throw std::runtime_error("No density files have been found");
		cout << fnames.size() << " density files have been found" << endl;
		if(!noInvFileOrder) std::sort(fnames.begin(), fnames.end(), std::less<>());
		else std::sort(fnames.begin(), fnames.end(), std::greater<>());
		dens.clear();
		for (int i = 0; i < fnames.size(); i += 3) {
			auto x = Grid(fnames[i]).data;
			auto y = Grid(fnames[i+1]).data;
			auto z = Grid(fnames[i+2]).data;
			vector<Point> t;
			for (int i = 0; i < x.size(); ++i)
				//t.push_back({ x[i], 0, 0 });
				t.push_back({x[i], y[i], z[i]});
			dens.push_back(t);
		}
	}

private:
	const std::string ext = ".grd";


	void getDensFiles(const std::string dirName) {
		std::vector<std::string> tmp = getFileNamesInFolder(dirName);
		fnames.clear();
		for(auto &i : tmp) {
			if(i.length() > ext.length() && !i.compare(i.length()-ext.length(), std::string::npos, ext))
				fnames.push_back(dirName+"/"+i);
		}
	}

	static std::vector<std::string> getFileNamesInFolder(const std::string dirName) {
		std::vector<std::string> files;
		for (auto& p : std::filesystem::directory_iterator(dirName))
			files.push_back(p.path().filename().generic_string());
		return files;
	}

	static void error(const std::string str) {
		std::cerr << "Input error: " << str << std::endl;
		std::cerr << "Program terminated." << std::endl;
		exit(1);
	}
};

#endif /* INPUTLOADER_H_ */
