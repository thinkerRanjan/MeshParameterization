#ifndef RECTANGLE_H_
#define RECTANGLE_H_
#include <vector>
#include <list>
using namespace std;

struct Patch
{
public:
	unsigned int id;
	unsigned int width;
	unsigned int height;

	// x and y represent the position of the top-left corner
	// of the smaller rectangle in the larger rectangle
	// to be calculated

	unsigned int x;
	unsigned int y;

};

class EnclosingRect
{
	public:
		vector<int> x, y;
		unsigned int numX, numY, height, width;
		vector< vector<bool> > occupied;

		EnclosingRect(int, int);
		~EnclosingRect();
};

// class Grid
// {
// 	public:

// 		int xi;
// 		int xf;
// 		bool occupied;
// 		Grid();
// 		~Grid();
	
// };

// class Row
// {
// public:

// 	int xi, xf, yi, yf;
// 	vector<Grid> box;

// 	Row(int, int);
// 	~Row();
	
// };

#endif