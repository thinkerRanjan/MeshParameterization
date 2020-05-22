// DCELHalfEdge.cpp: implementation of the DCELHalfEdge class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "DCELHalfEdge.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DCELHalfEdge::DCELHalfEdge()
{
	isBoundary = false;
	isProcessed = false;
	curvature = 0;
}

DCELHalfEdge::~DCELHalfEdge()
{

}

bool isBoundary = false;