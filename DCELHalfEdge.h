// DCELHalfEdge.h: interface for the DCELHalfEdge class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DCELHALFEDGE_H__A8186B0F_19D5_48EF_BD30_EB266F9A8215__INCLUDED_)
#define AFX_DCELHALFEDGE_H__A8186B0F_19D5_48EF_BD30_EB266F9A8215__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/*
 * DCELHalfEdge class. Part of an example DCEL implementation
 * Webpage: http://www.holmes3d.net/graphics/dcel/
 * Author: Ryan Holmes
 * E-mail: ryan <at> holmes3d <dot> net
 * Usage: Use freely. Please cite the website as the source if you
 * use it substantially unchanged. Please leave this documentation
 * in the code.
 */
#include<vector>
#include<tuple>
#include<vector>

using namespace std;
typedef tuple<float, float> vec2;

class DCELFace;
class DCELVertex;


class DCELHalfEdge  
{
public:
	bool isBoundary;
	bool isProcessed = false;
	float curvature;
	DCELHalfEdge();
	
	~DCELHalfEdge();

	DCELHalfEdge* twin;
	DCELHalfEdge* next;
	DCELFace* face;
	DCELVertex* origin;

	std::vector<int> boundaryPixelArray;
	std::vector<vec2> ratios;
	// Helper functions for manipulating the displayBits property
	
};

#endif // !defined(AFX_DCELHALFEDGE_H__A8186B0F_19D5_48EF_BD30_EB266F9A8215__INCLUDED_)
