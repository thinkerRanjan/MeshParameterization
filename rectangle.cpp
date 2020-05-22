#include "stdafx.h"
#include "rectangle.h"
#include <vector>

// Row::Row(int maxheight, int maxwidth)
// {
// 	xi = 0;
// 	yi = 0;
// 	xf = maxwidth;
// 	yf = maxheight;
// 	box.resize(1);
// }

// Row::~Row(){}

// Grid::Grid() {	occupied = false;	}

// Grid::~Grid(){}

EnclosingRect::EnclosingRect(int maxWidth, int maxHeight)
{
	height = maxHeight;
	width = maxWidth;
	numX = 1;
	numY = 1;
	vector<bool> temp(1);
	temp[0] = false;
	occupied.resize(1);
	occupied[0] = temp;
	x.push_back(0);
	x.push_back(maxWidth);
	y.push_back(0);
	y.push_back(maxHeight);
}

EnclosingRect::~EnclosingRect(){}