#ifndef GRAFENARGS_H_
#define GRAFENARGS_H_

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

class GrafenArgs {
public:
	limits Elim;
	limits Nlim;
	limits Hlim;
	double l0;
	std::vector<std::vector<double>> dens;
	Dat3D dat;
	bool noGK;
	double dotPotentialRad = 1e8;
	std::vector<std::string> fnames;

	std::string grdFname;
	int grdCols = 0;
	int grdRows = 0;

	bool noInvFileOrder = false;
	bool trasSolver = false;
	bool dat2D = false;
	bool grdFile = false;

	//input: -dat[2D/3D] (file.dat) [-Hf (val)] -Hfrom (val) -Hto (val) -Hn (val) -l0 (val) [-dens (directory)] [-toRel] [-densVal (val)] 
		//[-DPR (val)] *** point potential replace radius
		// if "-dens" is not specified: -Efrom -Eto -En -Nfrom -Nto -Nn [-densLayers (file.dat)]
		// [-noInvFileOrder]
		// [-transposeSolver]
		// [-saveDatAs3D]
	GrafenArgs(int argc, char *argv[]) {
		InputParser ip(argc, argv);

		noInvFileOrder = ip.exists("noInvFileOrder");
		trasSolver = ip.exists("transposeSolver");

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
			if(fnames.size() < Hlim.n)
				throw std::runtime_error("Not enough density files have been found");
			if(ip.exists("toRel")) relDens(dens);
		} else {
			ip["Efrom"] >> Elim.lower;
			ip["Eto"] >> Elim.upper;
			ip["En"] >> Elim.n;
			ip["Nfrom"] >> Nlim.lower;
			ip["Nto"] >> Nlim.upper;
			ip["Nn"] >> Nlim.n;
		}

		ip["l0"] >> l0;
		l0 = toRad(l0);

		if(ip.exists("densLayers")) {
			std::string fname;
			ip["densLayers"] >> fname;
			Dat2D d(fname);
			dens.clear();
			d.forEach([&](Dat2D::Element &el){ dens.push_back(std::vector<double>(Elim.n*Nlim.n, el.p.y)); });
		}

		if(ip.exists("densVal")) {
			double tmp;
			ip["densVal"] >> tmp;
			std::vector<double> l(Elim.n*Nlim.n, tmp);
			dens = std::vector<std::vector<double>>(Hlim.n, l);
			cout << "All densities are set to " << tmp << endl;
		}

		std::string datFname;
		if (ip.exists("dat")) {
			dat2D = true;
			ip["dat"] >> datFname;
		} else if (ip.exists("dat2D")) {
			dat2D = true;
			ip["dat2D"] >> datFname;
		} else if (ip.exists("dat3D")) {
			dat2D = false;
			ip["dat3D"] >> datFname;
		} else if (ip.exists("grd7")) {
			dat2D = true;
			grdFile = true;
			ip["grd7"] >> grdFname;
			datFname = grdFname;
		} else throw std::runtime_error("No *.dat or *.grd file has been specified for the Field.");

		double Hf = 0;
		if(grdFile) {
			ip["Hf"] >> Hf;
			Grid g(datFname);
			grdCols = g.nCol;
			grdRows = g.nRow;
			dat.set(GDconv::toDat(g), Hf);
			dat.fileName = datFname.replace(datFname.length() - 3, 3, "dat");
			dat2D = !ip.exists("saveDatAs3D");
		} else if (dat2D) {
			ip["Hf"] >> Hf;
			dat.set(Dat2D(datFname), Hf);
			dat2D = !ip.exists("saveDatAs3D");
		} else dat.read(datFname);

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
		for (auto &i : fnames)
			dens.push_back(Grid(i).data);
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

	void relDens(std::vector<std::vector<double>> &dens) {
		for(auto &v : dens) {
			double med = 0;
			for(auto &d : v)
				med += d;
			med /= (double)v.size();
			for(auto &d : v)
				d -= med;
		}
	}
};

#endif 
