#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

struct Point
{
	float _p[3];
};

typedef vector< Point > PointVec;

void writeMeshToPLYFile(const PointVec&    pointVec, const string&      outFilename)
{
	ofstream outFile(outFilename.c_str());

	////
	// Header
	////

	const int pointNum = (int)pointVec.size();

	outFile << "ply" << endl;
	outFile << "format ascii 1.0" << endl;
	outFile << "element vertex " << pointNum << endl;
	outFile << "property float x" << endl;
	outFile << "property float y" << endl;
	outFile << "property float z" << endl;
	outFile << "element face " << 0 << endl;
	outFile << "property list uchar int vertex_index" << endl;
	outFile << "end_header" << endl;

	////
	// Points
	////

	for (int pi = 0; pi < pointNum; ++pi)
	{
		const Point& point = pointVec[pi];

		for (int vi = 0; vi < 3; ++vi)
			outFile << point._p[vi] << " ";

		outFile << endl;
	}

	return;
}