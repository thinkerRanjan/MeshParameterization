// RectSpace.h: interface for the DCELFace class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RECT_BTREE_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_)
#define AFX_RECT_BTREE_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Vector.h"
#include "DCELHalfEdge.h"
#include <vector>
#include"Pixel.h"
#include "RectSpace.h"

struct node
{
	RectSpace* rect;
	node* left;
	node* right;
	node* parent;
};

#endif // !defined(AFX_RECT_BTREE_H__A82CC426_6585_4204_BFC4_EAC734D283A5__INCLUDED_)
