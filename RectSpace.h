// RectSpace.h: interface for the DCELFace class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RECTSPACE_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_)
#define AFX_RECTSPACE_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Vector.h"
#include"DCELFace.h"
#include "DCELHalfEdge.h"
#include <vector>
#include"Image.h"

using namespace std;

class DCELFace;
struct Colour;

class RectSpace
{
public:
	
	int size[2];
	int patchNo;
	Colour* pixels;
	float* minmax;
	Vector* barycentricCoords = new Vector[size[0] * size[1]];
	Vector* WorldCoords;

	std::vector <DCELFace*> facesArray;
	RectSpace(int w, int h)
	{
		size[0] = w;
		size[1] = h;
		pixels = new Colour[size[0]*size[1]];
		WorldCoords = new Vector[size[0] * size[1]];
	}

	~RectSpace()
	{
		//delete pixels;
		//delete barycentricCoords;
		//delete WorldCoords;
	}
};

#endif // !defined(AFX_RECTSPACE_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_)
