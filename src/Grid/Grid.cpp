#include <float.h>
#include <stdexcept>
#include "Grid.h"

#ifdef win32
int isinf(double s)
{
  // By IEEE 754 rule, 2*Inf equals Inf
  return (2*s == s) && (s != 0);
}
#endif

__int32 Grid::ReadInt32(ifstream* fs)
{
	__int32 result;
	fs->read( reinterpret_cast<char*>(&result), sizeof result );
	return result;
}

double Grid::ReadDouble(ifstream* fs)
{
	double result;
	fs->read( reinterpret_cast<char*>(&result), sizeof result );
	return result;
}

void Grid::WriteInt32(ofstream* fs, __int32 value)
{
	fs->write(reinterpret_cast<char*>(&value), sizeof value);
}

void Grid::WriteDouble(ofstream* fs, double value)
{
	fs->write(reinterpret_cast<char*>(&value), sizeof value);
}

Grid::Grid(int nRow, int nCol, double xLL, double yLL, double xSize, double ySize) : 
		nRow(nRow), nCol(nCol), xLL(xLL), yLL(yLL), xSize(xSize), ySize(ySize) {
	data.resize(nRow*nCol);
	zMin = 0;
	zMax = 0;
	Rotation = 0;
	BlankValue = DBL_MAX;
	fname = "defGridName.grd";
}

Grid::Grid(const string &fileName)
{
	Read(fileName);
}

bool Grid::Read()
{
	return Read(fname);
}

bool Grid::Read(const string &fileName)
{
	fname = fileName;
	ifstream ifs(fileName, std::ifstream::binary);
	if(!ifs.good()) throw std::runtime_error("Cannot open file: " + fname);
	if (ReadInt32(&ifs) != 0x42525344) //grid header DSRB
		throw std::runtime_error(fname + " is not a valid Grid file!");
	//int Size =
		ReadInt32(&ifs);
	//int Version =
		ReadInt32(&ifs);

	while (true)
	{
		int ID = ReadInt32(&ifs);

		if (ID == 0x44495247)	// grid GRID
		{
			//int Size =
				ReadInt32(&ifs);

			nRow = ReadInt32(&ifs);
			nCol = ReadInt32(&ifs);
			xLL = ReadDouble(&ifs);
			yLL = ReadDouble(&ifs);
			xSize = ReadDouble(&ifs);
			ySize = ReadDouble(&ifs);
			zMin = ReadDouble(&ifs);
			zMax = ReadDouble(&ifs);
			Rotation = ReadDouble(&ifs);
			BlankValue = ReadDouble(&ifs);

			data.resize(nCol*nRow);
			continue;
		}
		if (ID == 0x41544144)	// data
		{
			//int Size =
				ReadInt32(&ifs);
			ifs.read((char*)data.data(), nCol * nRow * sizeof(double));
			break;
		}
		if (ID == 0x49544c46)	// fault
		{
			ifs.close();
			return false;
		}
	}

	ifs.close();
	return true;
}

bool Grid::Write() {
	return Write(fname);
}

bool Grid::Write(const string &fileName)
{
	fname = fileName;

	getMin();
	getMax();

	ofstream ofs(fileName, ios::binary | ios::out);
	
	WriteInt32(&ofs, 0x42525344); // header DSRB
	WriteInt32(&ofs, sizeof(__int32));
	WriteInt32(&ofs, 2); // Version
	
	WriteInt32(&ofs, 0x44495247); // grid GRID
	WriteInt32(&ofs, 2 * sizeof(__int32) + 8 * sizeof(double));

	WriteInt32(&ofs, nRow);
	WriteInt32(&ofs, nCol);
	WriteDouble(&ofs, xLL);
	WriteDouble(&ofs, yLL);
	WriteDouble(&ofs, xSize);
	WriteDouble(&ofs, ySize);
	WriteDouble(&ofs, zMin);
	WriteDouble(&ofs, zMax);
	WriteDouble(&ofs, Rotation);
	WriteDouble(&ofs, BlankValue);

	for (int i = 0; i < nCol * nRow; i++)
		if (isinf(data[i]))
			data[i] = BlankValue;

	WriteInt32(&ofs, 0x41544144); // data
	__int32 size = nCol * nRow * sizeof(double);
	WriteInt32(&ofs, size);
	ofs.write((char*)data.data(), size);
	
	return true;
}

double Grid::getAverage()
{
	int count = 0;
	double sum = 0;
	for (int i = 0; i < nCol * nRow; i++)
		if (data[i] != BlankValue)
		{
			count++;
			sum += data[i];
		}
	
	return sum / count;
}

double Grid::getMin()
{
	zMin = DBL_MAX;
	for (int i = 0; i < nCol * nRow; i++)
		if (data[i] < zMin && !isinf(data[i]))
			zMin = data[i];
	return zMin;
}

double Grid::getMax()
{
	zMax = DBL_MIN;
	for (int i = 0; i < nCol * nRow; i++)
		if (data[i] > zMax && !isinf(data[i]))
			zMax = data[i];
	return zMax;
}

std::ostream& operator<<(std::ostream& os, const Grid& g) {
	os << "GridName=" << g.fname << " | " <<
			"nCol=" << g.nCol << " | " <<
			"nRow=" << g.nRow << " | " <<
			"xLL=" << g.xLL << " | " <<
			"yLL=" << g.yLL << " | " <<
			"xSize=" << g.xSize << " | " <<
			"ySize=" << g.ySize;
	return os;
}
