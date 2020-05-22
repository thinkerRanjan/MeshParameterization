// RectSpace.h: interface for the DCELFace class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PIXEL_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_)
#define AFX_PIXEL_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Vector.h"
#include "DCELHalfEdge.h"
#include <vector>


class Pixel
{
public:
	double r, g, b;

	Pixel()
	{
		r = 0.0;
		g = 0.0;
		b = 0.0;
	}

	Pixel(double x, double y, double z)
	{
		r = x;
		g = y;
		b = z;
	}

	~Pixel();

};

#endif // !defined(AFX_PIXEL_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_)
