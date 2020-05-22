// DCELFace.h: interface for the DCELFace class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DCELFACE_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_)
#define AFX_DCELFACE_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Vector.h"
#include "DCELHalfEdge.h"
#include "RectSpace.h"

/*
 * DCELFace class. Part of an example DCEL implementation
 * Webpage: http://www.holmes3d.net/graphics/dcel/
 * Author: Ryan Holmes
 * E-mail: ryan <at> holmes3d <dot> net
 * Usage: Use freely. Please cite the website as the source if you
 * use it substantially unchanged. Please leave this documentation
 * in the code.
 */

class DCELFace  
{
public:
	DCELFace();
	~DCELFace();

	Vector normal;
	DCELHalfEdge* edge;
	float st[6];
	RectSpace* rect;
	int group, patchNo;

};

#endif // !defined(AFX_DCELFACE_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_)
