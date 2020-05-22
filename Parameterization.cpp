// Parameterization.cpp : Defines the entry point for the console application.
//
#pragma warning( disable : 4996 )
#include "stdafx.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <tuple>
#include <fstream>
#include <sstream>
#include <string>
#include "Image.h"
#include"ANN\ANN.h"
#include"bitmap_image.hpp"

//DCEL headers
#include "DCELFace.h"
#include "DCELHalfEdge.h"

#include "DCELVertex.h"
#include "Vector.h"
#include "HalfEdgeList.h"
#include "RectSpace.h"
#include "Pixel.h"
#include "rectangle.h"

//VCG headers
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clustering.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/import.h>

#define PIXEL_DENSITY 450
#define _CRT_SECURE_NO_WARNINGS


using namespace std;

//globals
float mesh_minmax[6];
float mesh_max;
typedef tuple<float, float> vec2;
typedef tuple<float, float, float> vec3;
typedef tuple<int, int, int> int_vec3;
vec2 vUV;
vec3 vPosition, vNormal, fNormal;
vec3 x_pos(1.0, 0.0, 0.0), y_pos(0.0, 1.0, 0.0), z_pos(0.0, 0.0, 1.0), x_neg(-1.0, 0.0, 0.0), y_neg(0.0, -1.0, 0.0), z_neg(0.0, 0.0, -1.0);

int_vec3 faceVertexIndex, faceVertexTextureCoordIndex, faceVertexNormalIndex;
vector<vec3> modelVertices, modelNormals, modelFaceNormals;
vector<vec2> modeluvCoords, modeluvCoords_mod;
vector<int_vec3> modelFaceVertexIndices, modelFaceVertexNormalIndices, modelFaceVertexTextureCoordIndices;
vector<int> x_pos_faces, x_neg_faces, y_pos_faces, y_neg_faces, z_pos_faces, z_neg_faces;
vector<DCELVertex*> modelVerticesDCEL;
vector<Vector*> modelNormalDCEL;
vector<DCELFace*> modelFaces;
HalfEdgeList* modelEdgeList = new HalfEdgeList();
vector<RectSpace*> rectList;
ofstream fN_file("log.txt");
ofstream fout_file("output.obj");
ofstream outFile("pointCloud.ply");
ofstream wc_file("pixel_wc_st.txt");//wc = wprld coordinates
ofstream nb_file("neighbours.txt");

//variables for rectangle packing
vector<int> checkPlace(EnclosingRect, Patch, int, int);
vector<int> findPlace(EnclosingRect, Patch&);
vector<int> putPatch(EnclosingRect&, Patch&);
vector<int> putRectangles(vector<Patch>&, int, int);

//vcg content
class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>   ::AsVertexType,
	vcg::Use<MyEdge>     ::AsEdgeType,
	vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::TexCoord2f, vcg::vertex::Curvaturef, vcg::vertex::CurvatureDirf, vcg::vertex::Color4b, vcg::vertex::VFAdj, vcg::vertex::BitFlags >{};
class MyFace : public vcg::Face< MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3f, vcg::face::Color4b, vcg::face::FFAdj, vcg::face::BitFlags, vcg::face::VFAdj > {};
class MyEdge : public vcg::Edge<MyUsedTypes>{};
class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace>, std::vector<MyEdge>  > {};

//ofstream outFile("curvature.txt");
float Mc = 0.0;
float Sd = 0.0;
float Min;

/*addition of two vectors
*/
vec3 addVec(vec3 x, vec3 y)
{
	return (vec3(get<0>(x)+get<0>(y), get<1>(x)+get<1>(y), get<2>(x)+get<2>(y)));
}


//absolute value of a vector
float absVec(vec3 x)
{
	return(pow(pow(get<0>(x), 2) + pow(get<1>(x), 2) + pow(get<2>(x), 2), 0.5));
}

float absVec(Vector x)
{
	return(pow(pow(x.x, 2) + pow(x.y, 2) + pow(x.z, 2), 0.5));
}

//dot product of two vectors
float dot(vec3 x, vec3 y)
{
	return (get<0>(x)*get<0>(y)+get<1>(x)*get<1>(y)+get<2>(x)*get<2>(y));
}

float cosVec(vec3 x, vec3 y)
{
	float t = get<0>(x)*get<0>(y)+get<1>(x)*get<1>(y)+get<2>(x)*get<2>(y);
	t = t / (absVec(x)*absVec(y));
	return t;
}

float cosVec(Vector x, Vector y)
{
	float t = x.x * y.x + x.y * y.y + x.z * y.z ;
	t = t / (absVec(x)*absVec(y));
	return t;
}


//calculation of face normal by averaging the three vertex 
//modelNormals forming a face in triangular mesh
vec3 calcFaceNormal(vec3 x, vec3 y, vec3 z)
{
	vec3 t = addVec(addVec(x, y), z);
	//cout << "t:" << get<0>(t)<<"  "<<get<1>(t)<<"   " << get<2>(t) << "  Modulus:" << absVec(t)<<"\n";
	t = vec3(get<0>(t) / absVec(t), get<1>(t) / absVec(t), get<2>(t) / absVec(t));
	return t;
}



void patchIdentify(DCELFace* face, int* patch)
{
	face->patchNo = *patch;

	if (((face->edge->twin != NULL)) && (face->edge->twin->face->group == face->group) && (face->edge->twin->face->patchNo == 0))
	{
		patchIdentify(face->edge->twin->face, patch);
	}
	
	if (((face->edge->next->twin != NULL)) && (face->edge->next->twin->face->group == face->group) && (face->edge->next->twin->face->patchNo == 0))
	{
		patchIdentify(face->edge->next->twin->face, patch);
	}

	if (((face->edge->next->next->twin != NULL)) && (face->edge->next->next->twin->face->group == face->group) && (face->edge->next->next->twin->face->patchNo == 0))
	{
		face->patchNo = *patch;
		patchIdentify(face->edge->next->next->twin->face, patch);
	}

	return;
}

static int fcall = 0;
void boundaryCheck()
{
	for (vector<DCELFace*>::iterator i = modelFaces.begin(); i != modelFaces.end(); i++)
	{
		if ((*i)->edge->twin != NULL)
		{
			if ((*i)->edge->face->patchNo != (*i)->edge->twin->face->patchNo)
			{
				(*i)->edge->isBoundary = true;
				fcall++;
			}
		}
		else
		{
			(*i)->edge->isBoundary = true;
		}

		if ((*i)->edge->next->twin != NULL)
		{
			if ((*i)->edge->next->face->patchNo != (*i)->edge->next->twin->face->patchNo)
			{
				(*i)->edge->next->isBoundary = true;
				fcall++;
			}
		}
		else
		{
			(*i)->edge->next->isBoundary = true;
		}

		if ((*i)->edge->next->next->twin != NULL)
		{
			if ((*i)->edge->next->next->face->patchNo != (*i)->edge->next->next->twin->face->patchNo)
			{
				(*i)->edge->next->next->isBoundary = true;
				fcall++;
			}
		}
		else
		{
			(*i)->edge->next->next->isBoundary = true;
		}
	}
	fN_file << "fcall: " << fcall << endl;
	return;
}
/*
Calculate the projection of face onto a plane with the given normal
*/
void projectionCoords(DCELFace* face, Vector normal)
{
	fN_file << "face in projCoords():"<< face << endl;
	DCELVertex* v1 = new DCELVertex();
	DCELVertex* v2 = new DCELVertex();
	DCELVertex* v3 = new DCELVertex();

	v1 = face->edge->origin;
	v2 = face->edge->next->origin;
	v3 = face->edge->next->next->origin;


	float t1 = -(v1->coords.x * normal.x + v1->coords.y * normal.y + v1->coords.z * normal.z) / (pow(normal.x, 2) + pow(normal.y, 2) + pow(normal.z, 2));
	float t2 = -(v2->coords.x * normal.x + v2->coords.y * normal.y + v2->coords.z * normal.z) / (pow(normal.x, 2) + pow(normal.y, 2) + pow(normal.z, 2));
	float t3 = -(v3->coords.x * normal.x + v3->coords.y * normal.y + v3->coords.z * normal.z) / (pow(normal.x, 2) + pow(normal.y, 2) + pow(normal.z, 2));
	
	if (face->group == 1)
	{
		face->st[0] = (v1->coords.z + normal.z*t1);	
		face->st[1] = v1->coords.y + normal.y*t1;
		face->st[2] = (v2->coords.z + normal.z*t2);
		face->st[3] = v2->coords.y + normal.y*t2;
		face->st[4] = (v3->coords.z + normal.z*t3);
		face->st[5] = v3->coords.y + normal.y*t3;
	}

	if (face->group == 2)
	{
		face->st[0] = v1->coords.z + normal.z*t1;
		face->st[1] = v1->coords.y + normal.y*t1;
		face->st[2] = v2->coords.z + normal.z*t2;
		face->st[3] = v2->coords.y + normal.y*t2;
		face->st[4] = v3->coords.z + normal.z*t3;
		face->st[5] = v3->coords.y + normal.y*t3;
	}

	if (face->group == 3)
	{
		face->st[0] = (v1->coords.z + normal.z*t1);
		face->st[1] = v1->coords.x + normal.x*t1;
		face->st[2] = (v2->coords.z + normal.z*t2);
		face->st[3] = v2->coords.x + normal.x*t2;
		face->st[4] = (v3->coords.z + normal.z*t3);
		face->st[5] = v3->coords.x + normal.x*t3;
	}

	if (face->group == 4)
	{
		face->st[0] = v1->coords.z + normal.z*t1;
		face->st[1] = v1->coords.x + normal.x*t1;
		face->st[2] = v2->coords.z + normal.z*t2;
		face->st[3] = v2->coords.x + normal.x*t2;
		face->st[4] = v3->coords.z + normal.z*t3;
		face->st[5] = v3->coords.x + normal.x*t3;
	}

	if (face->group == 5)
	{
		face->st[0] = (v1->coords.x + normal.x * t1);
		face->st[1] = v1->coords.y + normal.y * t1;
		face->st[2] = (v2->coords.x + normal.x * t2);
		face->st[3] = v2->coords.y + normal.y * t2;
		face->st[4] = (v3->coords.x + normal.x * t3);
		face->st[5] = v3->coords.y + normal.y * t3;
	}

	if (face->group == 6)
	{
		face->st[0] = v1->coords.x + normal.x*t1; 
		face->st[1] = v1->coords.y + normal.y*t1;
		face->st[2] = v2->coords.x + normal.x*t2;
		face->st[3] = v2->coords.y + normal.y*t2; 
		face->st[4] = v3->coords.x + normal.x*t3;
		face->st[5] = v3->coords.y + normal.y*t3;
	}

	
	

	fN_file << "projection:face st: " << face->st[0] << " " << face->st[1] << " " << face->st[2]<< " ";
	fN_file << face->st[3] << " " << face->st[4] << " " << face->st[5] << endl;

	return;
}



/*
find the bottom left and top right corner of a patch
*/
void BL_TR_PatchCorners(HalfEdgeList* walker1, int patch, float* minmax)
{
	float min_x, min_y, max_x, max_y;
	min_x = 0;
	min_y = 0;
	max_x = 0;
	max_y = 0;
	bool count = true;
	while (walker1->next)						//first loop for traversing the halfedgelist
	{
		if (walker1->edge->face->patchNo == patch)
		{
			if (count == true)
			{
				min_x = walker1->edge->face->st[0];
				min_y = walker1->edge->face->st[1];
				max_x = walker1->edge->face->st[0];
				max_y = walker1->edge->face->st[1];

				count = false;
			}

			if (walker1->edge->face->st[0] < min_x)
				min_x = walker1->edge->face->st[0];

			if (walker1->edge->face->st[2] < min_x)
				min_x = walker1->edge->face->st[2];

			if (walker1->edge->face->st[4] < min_x)
				min_x = walker1->edge->face->st[4];

			if (walker1->edge->face->st[1] < min_y)
				min_y = walker1->edge->face->st[1];

			if (walker1->edge->face->st[3] < min_y)
				min_y = walker1->edge->face->st[3];

			if (walker1->edge->face->st[5] < min_y)
				min_y = walker1->edge->face->st[5];

			if (walker1->edge->face->st[0] > max_x)
				max_x = walker1->edge->face->st[0];

			if (walker1->edge->face->st[2] > max_x)
				max_x = walker1->edge->face->st[2];

			if (walker1->edge->face->st[4] > max_x)
				max_x = walker1->edge->face->st[4];

			if (walker1->edge->face->st[1] > max_y)
				max_y = walker1->edge->face->st[1];

			if (walker1->edge->face->st[3] > max_y)
				max_y = walker1->edge->face->st[3];

			if (walker1->edge->face->st[5] > max_y)
				max_y = walker1->edge->face->st[5];

		}
		
		walker1 = walker1->next;
	}

	minmax[0] = min_x;
	minmax[1] = max_x;
	minmax[2] = min_y;
	minmax[3] = max_y;

	fN_file << "BL_()--minmax: " << minmax[0] << " " << minmax[1] << " " << minmax[2] << " " << minmax[3] << endl;
	return;
}

/*
find the extreme coordinates for the whole mesh
*/
void meshExtremes()
{
	bool flag = true;
	for (vector<DCELVertex*>::iterator i = modelVerticesDCEL.begin(); i != modelVerticesDCEL.end(); i++)
	{
		if (flag)
		{
			mesh_minmax[0] = (*i)->coords.x;
			mesh_minmax[1] = (*i)->coords.x;
			mesh_minmax[2] = (*i)->coords.y;
			mesh_minmax[3] = (*i)->coords.y;
			mesh_minmax[4] = (*i)->coords.z;
			mesh_minmax[5] = (*i)->coords.z;

			flag = false;
		}

		if ((*i)->coords.x < mesh_minmax[0])
			mesh_minmax[0] = (*i)->coords.x;
		if ((*i)->coords.x > mesh_minmax[1])
			mesh_minmax[1] = (*i)->coords.x;
		if ((*i)->coords.y < mesh_minmax[2])
			mesh_minmax[2] = (*i)->coords.y;
		if ((*i)->coords.y > mesh_minmax[3])
			mesh_minmax[3] = (*i)->coords.y;
		if ((*i)->coords.z < mesh_minmax[4])
			mesh_minmax[4] = (*i)->coords.z;
		if ((*i)->coords.z > mesh_minmax[5])
			mesh_minmax[5] = (*i)->coords.z;
	}

	fN_file << "Mesh Extremes:\n" << mesh_minmax[0] << "  " << mesh_minmax[1] << "  " << mesh_minmax[2] << "  " << mesh_minmax[3] << "  " << mesh_minmax[4] << "  " << mesh_minmax[5] << endl;
	return;
}

/*
	Calculate the projections of all faces onto their respective direction planes
*/
void getProject(HalfEdgeList* walker1, int totalPatch)
{
	while (walker1->next)				//first loop for traversing the halfedgelist
	{
		fN_file << "face address: "<<walker1->edge->face << endl;
		if (walker1->edge->face->group == 1)
		{
			fN_file << "group1" << endl;
			projectionCoords(walker1->edge->face, Vector(1.0, 0.0, 0.0));
		}
		if (walker1->edge->face->group == 2)
		{
			fN_file << "group2" << endl;
			projectionCoords(walker1->edge->face, Vector(-1.0, 0.0, 0.0));
		}
		if (walker1->edge->face->group == 3)
		{
			fN_file << "group3" << endl;
			projectionCoords(walker1->edge->face, Vector(0.0, 1.0, 0.0));
		}
		if (walker1->edge->face->group == 4)
		{
			fN_file << "group4" << endl;
			projectionCoords(walker1->edge->face, Vector(0.0, -1.0, 0.0));
		}
		if (walker1->edge->face->group == 5)
		{
			fN_file << "group5" << endl;
			projectionCoords(walker1->edge->face, Vector(0.0, 0.0, 1.0));
		}
		if (walker1->edge->face->group == 6)
		{
			fN_file << "group6" << endl;
			projectionCoords(walker1->edge->face, Vector(0.0, 0.0, -1.0));
		}
		
		walker1 = walker1->next;
	}

	return;
}


void setBLtoZero()
{
	for (vector<DCELVertex*>::iterator i = modelVerticesDCEL.begin(); i != modelVerticesDCEL.end(); i++)
	{
		fN_file << "Coords before: " << (*i)->coords.x << "   " << (*i)->coords.y << "   " << (*i)->coords.z << endl;
		
		(*i)->coords = (*i)->coords - Vector(mesh_minmax[0], mesh_minmax[2], mesh_minmax[4]);

		fN_file << "Coords after: " << (*i)->coords.x << "   " << (*i)->coords.y << "   " << (*i)->coords.z << endl;
	}

	mesh_minmax[1] = mesh_minmax[1] - mesh_minmax[0];
	mesh_minmax[3] = mesh_minmax[3] - mesh_minmax[2];
	mesh_minmax[5] = mesh_minmax[5] - mesh_minmax[4];

	mesh_minmax[0] = 0;
	mesh_minmax[2] = 0;
	mesh_minmax[4] = 0;

	fN_file << "Mesh Extremes:\n" << mesh_minmax[0] << "  " << mesh_minmax[1] << "  " << mesh_minmax[2] << "  " << mesh_minmax[3] << "  " << mesh_minmax[4] << "  " << mesh_minmax[5] << endl;
	
	return;
}


void stUnitize()
{
	mesh_max = mesh_minmax[1];
	if (mesh_minmax[3] > mesh_max)
		mesh_max = mesh_minmax[3];
	if (mesh_minmax[5] > mesh_max)
		mesh_max = mesh_minmax[5];

	fN_file << "Mesh Max: " << mesh_max << endl;
	
	for (vector<DCELVertex*>::iterator i = modelVerticesDCEL.begin(); i != modelVerticesDCEL.end(); i++)
	{
		(*i)->coords = (*i)->coords / mesh_max;
	}

	return;
}


//creates a rectangular space for color info
//corrects and unitizes the s and t coords of all the faces in a patch
void createRect(int patch)
{	
	for (int j = 1; j <= patch; j++)
	{
		RectSpace* rect;
		float* minmax = new float[4];//format: min x, max x, min y, max y
		//find the dimension of the rectangle required for the patch
		//find the corner coordinates
		BL_TR_PatchCorners(modelEdgeList, j, minmax);
		
		//set left bottom corner to zero, zero
		
		//create the rectangular space with patch size obtained from BL TR calculation
		//fN_file << "Rect size:" << (int)ceil((minmax[1] - minmax[0]) * 1.5 * PIXEL_DENSITY) << "  by  " << (int)ceil((minmax[3] - minmax[2]) * 1.5 * PIXEL_DENSITY) << endl;
		rect = new RectSpace((int)ceil((minmax[1] - minmax[0]) * PIXEL_DENSITY + 4), (int)ceil((minmax[3] - minmax[2]) * PIXEL_DENSITY + 4));
		
		rect->patchNo = j;
		rect->minmax = new float[4];
		rect->minmax = minmax;

		//fill the faces array in the rectangular patch with patchNo = patch
		HalfEdgeList* walker1 = modelEdgeList;
		while (walker1->next)
		{
			if (walker1->edge->face->patchNo == j)
			{
				walker1->edge->face->rect = rect;    //associate each face to rectangular space
				rect->facesArray.push_back(walker1->edge->face);   //each rectangle knows its faces
				fN_file << "face createRect st: " << "  " << walker1->edge->face->st[0] << "  " << walker1->edge->face->st[1] << "  " << walker1->edge->face->st[2] << "  ";
				fN_file << walker1->edge->face->st[3] << "  " << walker1->edge->face->st[4] << "  " << walker1->edge->face->st[5] << endl;
			}
			walker1 = walker1->next;
		}

		rectList.push_back(rect);

	}
	return;
}

/*
sort the rect list using the height i.e size[1] (bubble sort)
*/
void sortRectList()
{
	RectSpace* swap;
	for (int j = 0; j < rectList.size()-1; j++)
	{
		for (int i = 0; i < rectList.size() - j - 1; i++)
		{
			if (rectList[i]->size[1] < rectList[i + 1]->size[1])
			{
				swap = rectList[i];
				rectList[i] = rectList[i + 1];
				rectList[i + 1] = swap;
			}
		}
	}
	for (vector<RectSpace*>::iterator i = rectList.begin(); i != rectList.end(); i++)
	{
		std::cout << "Size: " << (*i)->size[0]<<" "<<(*i)->size[1]<<endl;
	}
	
}
//fill the color values in pixel arrays
void writeImage(RectSpace* rect, string file)
{
	TGAImage *img = new TGAImage(rect->size[0], rect->size[1]);

	//Loop through image and set all pixels to red
	//std::cout << rect->size[0] * rect->size[1] << endl;
	for (int i = 0; i < (rect->size[1]); i++)
	{
		for (int j = 0; j < (rect->size[0]); j++)
		{
			//std::cout << i << " "<<j<<endl;
			img->setPixel(rect->pixels[j + rect->size[0]*i], i, j);
		}
	}

	//write the image to disk
	img->WriteImage(file);

	//delete img;
	return;

}

void writeImage_bmp(RectSpace* rect, string file)
{
	const unsigned int dim = rect->size[0];
	bitmap_image image(dim, dim);

	for (unsigned int y = 0; y < dim; ++y)
	{
		for (unsigned int x = 0; x < dim; ++x)
		{
			int t = rect->size[0] * (dim - y - 1) + x;
			image.set_pixel(x, y, rect->pixels[t].r, rect->pixels[t].g, rect->pixels[t].b);
		}
	}

	image.save_image(file);

}

bool pointInTriangle(float x, float y, float* triPt)
{
	float denominator = ((triPt[3] - triPt[5])*(triPt[0] - triPt[4]) + (triPt[4] - triPt[2])*(triPt[1] - triPt[5]));
	float a = (((triPt[3] - triPt[5]) * (x - triPt[4]) + (triPt[4] - triPt[2]) * (y - triPt[5]))) / denominator;
	float b = (((triPt[5] - triPt[1]) * (x - triPt[4]) + (triPt[0] - triPt[4]) * (y - triPt[5]))) / denominator;
	float c = (1 - a - b);

	if ((a > 0 && a < 1) && (b > 0 && b < 1) && (c > 0 && c < 1))
		return true;
	else
		return false;
	//std::cout << "a:" << a << " b:" << b << " c:" << c << endl;
	/*
	float epsilon = 0.000001;
	if ((a > epsilon) && (1.0 - a > epsilon) && (b > epsilon) && (1.0 - b > epsilon) && (c > epsilon) && (1.0 - c > epsilon))
		return true;
	else if ((a > -epsilon && a < epsilon) || (b > -epsilon && b < epsilon) || (c > -epsilon && c < epsilon))
	{
		return true;
	}
	else
		return false;
		*/
}

//origin first then the endpoint coords
bool proximityToEdge(float s, float t, float s1, float t1, float s2, float t2)
{
	//std::cout << s << " " << t << " " << s1 << " " << t1 << " " << s2 << " " << t2 << endl;
	float d = abs((s2-s1)*(t1 - t) - (s1 - s)*(t2 - t1)) / sqrt(pow((s2 - s1), 2) + pow((t2 - t1), 2));
	if (d <= 1.0)
	{
		fN_file << "proximity value(t): "<< d <<" "<< s <<" "<< t <<" "<< s1 <<" "<< t1 <<" "<< s2 <<" "<< t2 <<endl;
		return true;
	}

	else
	{
		//std::cout << "false\n";
		fN_file << "proximity value(f): " << d << " " << s << " " << t << " " << s1 << " " << t1 <<" "<< s2 << " " << t2 << endl;
		return false;
	}
}

//bresenham
void edgeRaster(int x1, int y1, int x2, int y2, DCELHalfEdge* edge, RectSpace* rect)
{
	int x, y, dx, dy, dx1, dy1, px, py, xe, ye, i, t;
	dx = x2 - x1;
	dy = y2 - y1;
	dx1 = abs(dx);
	dy1 = abs(dy);
	px = 2 * dy1 - dx1;
	py = 2 * dx1 - dy1;
	if (dy1 <= dx1)
	{
		if (dx >= 0)
		{
			x = x1;
			y = y1;
			xe = x2;
		}
		else
		{
			x = x2;
			y = y2;
			xe = x1;
		}
		
		t = rect->size[0] * y + x;
		edge->boundaryPixelArray.push_back(t);
		

		for (i = 0; x < xe; i++)
		{
			x = x + 1;
			if (px < 0)
			{
				px = px + 2 * dy1;
			}
			else
			{
				if ((dx<0 && dy<0) || (dx>0 && dy>0))
				{
					y = y + 1;
				}
				else
				{
					y = y - 1;
				}
				px = px + 2 * (dy1 - dx1);
			}

			t = rect->size[0] * y + x;
			edge->boundaryPixelArray.push_back(t);
		}
	}
	else
	{
		if (dy >= 0)
		{
			x = x1;
			y = y1;
			ye = y2;
		}
		else
		{
			x = x2;
			y = y2;
			ye = y1;
		}

		t = rect->size[0] * y + x;
		edge->boundaryPixelArray.push_back(t);

		for (i = 0; y < ye;i++)
		{
			y = y + 1;
			if (py <= 0)
			{
				py = py + 2 * dx1;
			}
			else
			{
				if ((dx<0 && dy<0) || (dx>0 && dy>0))
				{
					x = x + 1;
				}
				else
				{
					x = x - 1;
				}

				py = py + 2 * (dx1 - dy1);
			}

			t = rect->size[0] * y + x;
			edge->boundaryPixelArray.push_back(t);
		}
	}

}

//
bool colourCompare(Colour c1, Colour c2)
{
	if ((c1.r == c2.r) && (c1.g == c2.g) && (c1.b == c2.b))
	{
		return true;
	}
	else
		return false;
}

Colour* getColour(float x)
{
	//colour 1
	static Colour* blue = new Colour();
	blue->r = 0;
	blue->g = 0;
	blue->b = 255;
	//colour 2
	static Colour* lightblue = new Colour();
	lightblue->r = 102;
	lightblue->g = 153;
	lightblue->b = 255;
	//colour 3
	static Colour* turquoise = new Colour();
	turquoise->r = 102;
	turquoise->g = 255;
	turquoise->b = 255;
	//colour 4
	static Colour* lightgreen = new Colour();
	lightgreen->r = 102;
	lightgreen->g = 255;
	lightgreen->b = 102;
	//colour 5
	static Colour* yellow = new Colour();
	yellow->r = 255;
	yellow->g = 255;
	yellow->b = 102;
	//colour 6
	static Colour* red = new Colour();
	red->r = 255;
	red->g = 0;
	red->b = 0;
	
	if (x < -0.866 && x > -1.0)
		return blue;
	else if (x < -0.5 && x > -0.866)
		return lightblue;
	else if (x < 0.0 && x > -0.5)
		return turquoise;
	else if (x < 0.5 && x > 0.0 )
		return lightgreen;
	else if (x < 0.866 && x > 0.5)
		return yellow;
	else if (x < 1.0 && x > 0.866)
		return red;
}

//find exact intersections of each boundary pixel
void pixel_cuts(double x0, double y0, double x1, double y1, double p, double q, RectSpace* rect, DCELHalfEdge* edge)
{
	//slope of edge line
	double m = (y1 - y0) / (x1 - x0);
	//equation of other four lines
	/*			L3
	(p,q+1) ------------- (p+1, q+1)
	        |           |
		 L4 |           |  L2
		    |           |
	(p, q)  |___________| (p+1, q)
				L1
	y = q, x = p+1, y = q+1, x = p
	*/

	bool endCheck_1 = false;
	bool endCheck_2 = false;
	if (x0 > p && x0 < p + 1 && y0 > p && y0 < p + 1)
	{
		endCheck_1 = true;
	}

	if (x1 > p && x1 < p + 1 && y1 > p && y1 < p + 1)
	{
		endCheck_2 = true;
	}

	//find the intersections of edge with these four lines
	double u1 = 0.0;
	double u3 = 0.0;
	double v2 = 0.0;
	double v4 = 0.0;
	//with L1
	if (m != 0.0)
	{
		u1 = (q - y0) / m + x0; //the pt = (u1, q)
		u3 = (q + 1.0 - y0) / m + x0; // the pt = (u3, q+1)
	}
	v2 = m * (p + 1.0 - x0) + y0;  //the pt = (p+1, v2)
	v4 = m * (p - x0) + y0; // the pt = (p, v4)

	fN_file << "p, q, p+1, q+1: " << p << " " << q << " " << p + 1 << " " << q + 1 << endl;
	fN_file << "u1, v2, u3, v4: " << u1 << " " << v2 << " " << u3 << " " << v4 << endl;
	int t1 = rect->size[0] * q + p;
	double r1, r2;
	r1 = -1;
	r2 = -1;
	tuple<float, float> ratio;
	get<0>(ratio) = r1;
	get<1>(ratio) = r2;
	
	double d = sqrt(pow((x1 - x0), 2.0) + pow((y1 - y0), 2.0));
	// cases of intersection
	if (m != 0 && m != 1)
	{
		if ((u1 >= p && u1 <= p + 1) && (v2 >= q && v2 <= q + 1))
		{
			r1 = sqrt(pow((u1 - x0), 2.0) + pow((q - y0), 2.0)) / d;
			r2 = sqrt(pow((p + 1 - x0), 2.0) + pow((v2 - y0), 2.0)) / d;

			fN_file << "r1, r2: " << r1 << " " << r2 << endl;

			if (r1 <= r2)
			{
				get<0>(ratio) = r1;
				get<1>(ratio) = r2;
			}
			else
			{
				get<0>(ratio) = r2;
				get<1>(ratio) = r1;
			}

		

		}
		else if ((u3 >= p && u3 <= p + 1) && (v2 >= q && v2 <= q + 1))
		{
			r1 = sqrt(pow((u3 - x0), 2) + pow((q + 1 - y0), 2)) / d;
			r2 = sqrt(pow((p + 1 - x0), 2) + pow((v2 - y0), 2)) / d;

			fN_file << "r1, r2: " << r1 << " " << r2 << endl;
			
			if (r1 <= r2)
			{
				get<0>(ratio) = r1;
				get<1>(ratio) = r2;
			}
			else
			{
				get<0>(ratio) = r2;
				get<1>(ratio) = r1;
				
			}

			
		}
		else if ((u3 >= p && u3 <= p + 1) && (v4 >= q && v4 <= q + 1))
		{
			r1 = sqrt(pow((u3 - x0), 2) + pow((q + 1 - y0), 2)) / d;
			r2 = sqrt(pow((p - x0), 2) + pow((v4 - y0), 2)) / d;
			fN_file << "r1, r2: " << r1 << " " << r2 << endl;

			if (r1 <= r2)
			{
				get<0>(ratio) = r1;
				get<1>(ratio) = r2;
			}
			else
			{
				get<0>(ratio) = r2;
				get<1>(ratio) = r1;
			}
			
		}
		else if ((u1 >= p && u1 <= p + 1) && (v4 >= q && v4 <= q + 1))
		{
			r1 = sqrt(pow((u1 - x0), 2) + pow((q - y0), 2)) / d;
			r2 = sqrt(pow((p - x0), 2) + pow((v4 - y0), 2)) / d;
			fN_file << "r1, r2: " << r1 << " " << r2 << endl;
			
			if (r1 <= r2)
			{
				get<0>(ratio) = r1;
				get<1>(ratio) = r2;
			}
			else
			{
				get<0>(ratio) = r2;
				get<1>(ratio) = r1;
			}

		
		}

		else if ((u1 >= p && u1 <= p + 1) && (u3 >= p && u3 <= p + 1))
		{
			r1 = sqrt(pow((u1 - x0), 2) + pow((q - y0), 2)) / d;
			r2 = sqrt(pow((u3 - x0), 2) + pow((q+1 - y0), 2)) / d;
			fN_file << "r1, r2: " << r1 << " " << r2 << endl;

			if (r1 <= r2)
			{
				get<0>(ratio) = r1;
				get<1>(ratio) = r2;
			}
			else
			{
				get<0>(ratio) = r2;
				get<1>(ratio) = r1;
			}

		}

		else if ((v2 >= q && v2 <= q + 1) && (v4 >= q && v4 <= q + 1))
		{
			r1 = sqrt(pow((p+1 - x0), 2) + pow((v2 - y0), 2)) / d;
			r2 = sqrt(pow((p - x0), 2) + pow((v4 - y0), 2)) / d;
			fN_file << "r1, r2: " << r1 << " " << r2 << endl;

			if (r1 <= r2)
			{
				get<0>(ratio) = r1;
				get<1>(ratio) = r2;
			}
			else
			{
				get<0>(ratio) = r2;
				get<1>(ratio) = r1;
			}

			
		}
	}
	else if (m == 0.0)
	{
		if (y0 >= q && y0 <= q + 1)
		{
			r1 = fabs(p - x0) / d;
			r2 = fabs(p + 1 - x0) / d;
			fN_file << "r1, r2: " << r1 << " " << r2 << endl;

			if (r1 <= r2)
			{
				get<0>(ratio) = r1;
				get<1>(ratio) = r2;
			}
			else
			{
				get<0>(ratio) = r2;
				get<1>(ratio) = r1;
			}

		
		}
	}
	else if (m == 1.0)
	{
		if (x0 >= p && x0 <= p + 1)
		{
			r1 = fabs(q - y0) / d;
			r2 = fabs(q + 1 - y0) / d;
			fN_file << "r1, r2: " << r1 << " " << r2 << endl;

			if (r1 <= r2)
			{
				get<0>(ratio) = r1;
				get<1>(ratio) = r2;
			}
			else
			{
				get<0>(ratio) = r2;
				get<1>(ratio) = r1;
			}

		}
	}
	
	if (endCheck_1)
		get<0>(ratio) = 0.0;
	if (endCheck_2)
		get<1>(ratio) = 1.0;

	if (r1 < 0.01)
		get<0>(ratio) = 0.0;
	if (r2 < 0.01)
		get<1>(ratio) = 0.0;

	if (r1 > 0.99)
		get<0>(ratio) = 1.0;
	if (r2 > 0.99)
		get<1>(ratio) = 1.0;

	edge->ratios.push_back(ratio);
	fN_file << "ratio before return: " << get<0>(ratio) << " " << get<1>(ratio) << endl;
}
//DDA
void DDA_line(float x0, float y0, float x1, float y1, DCELHalfEdge* edge, RectSpace* rect, Colour* col, vector<DCELFace*> &atlasPixelFaces)    // DDA subpixel -> thick
{
	int i, t, n, x, y, xx, yy;

	Colour* black = new Colour();
	black->r = 0;
	black->g = 0;
	black->b = 0;
	
	Colour* red = new Colour();
	red->r = 255;
	red->g = 0;
	red->b = 0;

	float x0_original, x1_original, y0_original, y1_original;
	x0_original = x0;
	y0_original = y0;
	x1_original = x1;
	y1_original = y1;

	// prepare data n-pixels,x1,y1 is line dx,dy step per pixel
	x1 = x1 - x0;
	i = ceil(fabs(x1));
	y1 = y1 - y0;
	n = ceil(fabs(y1));
	if (n < i)
		n = i;
	if (!n)
		n = 1;
	x1 /= float(n);
	y1 /= float(n);
	n++;

	x = x0;
	y = y0;
	// rasterize DDA line
	for (xx = x0, yy = y0, i = 0; i <= n; i++, x0 += x1, y0 += y1)
	{
		// direct pixel
		t = rect->size[0] * y + x;

		if (x >= 0 && x <rect->size[0] && y >= 0 && y < rect->size[0])
		{
			rect->pixels[t].curvature = edge->curvature;
			
			rect->pixels[t] = *col; // *getColour(edge->curvature);
			
			atlasPixelFaces[t] = edge->face;
			edge->boundaryPixelArray.push_back(t);
			pixel_cuts(x0_original, y0_original, x1_original, y1_original, x, y, rect, edge);
			
			//rect->pixels[t].r = (int)(255 * get<0>(edge->ratios[edge->ratios.size() - 1]));
			//rect->pixels[t].g = (int)(255 * get<0>(edge->ratios[edge->ratios.size() - 1]));
			//rect->pixels[t].b = 0;

			float o = (get<1>(edge->ratios[edge->ratios.size() - 1]) + get<1>(edge->ratios[edge->ratios.size() - 1])) / 2;
			rect->pixels[t].K1 = edge->origin->K1 * (1-o) + edge->next->origin->K1 * o;
			rect->pixels[t].K2 = edge->origin->K2 * (1-o) + edge->next->origin->K2 * o;
			rect->pixels[t].Kg = edge->origin->Kg * (1-o) + edge->next->origin->Kg * o;

			if (atlasPixelFaces[t]->group == 1 || atlasPixelFaces[t]->group == 2)
			{
				rect->pixels[t].pd1.x = 0;
				rect->pixels[t].pd1.y = edge->origin->pd1[1] * (1 - o) + edge->next->origin->pd1[1] * o;
				rect->pixels[t].pd1.z = edge->origin->pd1[2] * (1 - o) + edge->next->origin->pd1[2] * o;

				rect->pixels[t].pd2.x = 0;
				rect->pixels[t].pd2.y = edge->origin->pd2[1] * (1 - o) + edge->next->origin->pd2[1] * o;
				rect->pixels[t].pd2.z = edge->origin->pd2[2] * (1 - o) + edge->next->origin->pd2[2] * o;

				rect->pixels[t].pd1.normalize();
				rect->pixels[t].pd2.normalize();

				rect->pixels[t].cos_theta_pd1 = rect->pixels[t].pd1[1];
				rect->pixels[t].cos_theta_pd2 = rect->pixels[t].pd2[1];
			}
			else if (atlasPixelFaces[t]->group == 3 || atlasPixelFaces[t]->group == 4)
			{
				rect->pixels[t].pd1.x = edge->origin->pd1[0] * (1 - o) + edge->next->origin->pd1[0] * o;
				rect->pixels[t].pd1.y = 0;
				rect->pixels[t].pd1.z = edge->origin->pd1[2] * (1 - o) + edge->next->origin->pd1[2] * o;

				rect->pixels[t].pd2.x = edge->origin->pd2[0] * (1 - o) + edge->next->origin->pd2[0] * o;
				rect->pixels[t].pd2.y = 0;
				rect->pixels[t].pd2.z = edge->origin->pd2[2] * (1 - o) + edge->next->origin->pd2[2] * o;

				rect->pixels[t].pd1.normalize();
				rect->pixels[t].pd2.normalize();

				rect->pixels[t].cos_theta_pd1 = rect->pixels[t].pd1[2];
				rect->pixels[t].cos_theta_pd2 = rect->pixels[t].pd2[2];
			}
			else if (atlasPixelFaces[t]->group == 5 || atlasPixelFaces[t]->group == 6)
			{
				rect->pixels[t].pd1.x = edge->origin->pd1[0] * (1 - o) + edge->next->origin->pd1[0] * o;
				rect->pixels[t].pd1.y = edge->origin->pd1[1] * (1 - o) + edge->next->origin->pd1[1] * o;
				rect->pixels[t].pd1.z = 0;

				rect->pixels[t].pd2.x = edge->origin->pd2[0] * (1 - o) + edge->next->origin->pd2[0] * o;
				rect->pixels[t].pd2.y = edge->origin->pd2[1] * (1 - o) + edge->next->origin->pd2[1] * o;
				rect->pixels[t].pd2.z = 0;

				rect->pixels[t].pd1.normalize();
				rect->pixels[t].pd2.normalize();

				rect->pixels[t].cos_theta_pd1 = rect->pixels[t].pd1[0];
				rect->pixels[t].cos_theta_pd2 = rect->pixels[t].pd2[0];
			}

		}
		
		x = x0;
		y = y0;
		
		xx = x; yy = y;
	}

}

void DDAf_line_subpixel(float x0, float y0, float x1, float y1, DCELHalfEdge* edge, RectSpace* rect, Colour* col, vector<DCELFace*> &atlasPixelFaces)    // DDA subpixel -> thick
{
	int i, n, t; 
	float a, a0, a1, aa, b, d;
	t = 0;
	// end-points
	//pnt(x0, y0, col);
	t = rect->size[0] * y0 + x0;
	rect->pixels[t].curvature = edge->curvature;
	rect->pixels[t] = *col; // *getColour(edge->curvature);
	atlasPixelFaces[t] = edge->face;
	edge->boundaryPixelArray.push_back(t);
	//std::cout << "(x0, y0, t)::" << x0 << ", " << y0 << ", " << t << endl;
	
	//pnt(x1, y1, col);
	t = rect->size[0] * y1 + x1;
	rect->pixels[t].curvature = edge->curvature;
	rect->pixels[t] = *col; // *getColour(edge->curvature);
	atlasPixelFaces[t] = edge->face;
	edge->boundaryPixelArray.push_back(t);

	// x-axis pixel cross
	a0 = 1; 
	a1 = 0; 
	n = 0;
	if (x0 < x1) 
	{
		a0 = ceil(x0); 
		a1 = floor(x1); 
		d = (y1 - y0) / (x1 - x0);  //slope of edge
		a = a0; 
		b = y0 + (a0 - x0) * d; 
		n = fabs(a1 - a0); 
	}
	else if (x0 > x1) 
	{ 
		a0 = ceil(x1);
		a1 = floor(x0); 
		d = (y1 - y0) / (x1 - x0); 
		a = a0; 
		b = y1 + (a0 - x1) * d;
		n = fabs(a1 - a0);
	}
	
	if (a0 <= a1)
	{
		for (aa = a, i = 0; i <= n; i++, aa = a, a++, b += d)
		{
			//pnt(aa, b, col);
			if ((int)aa == 101)
			{
				std::cout << "(" << aa << ", " << b << ")" << endl;
			}

			if (aa >= 0 && aa < rect->size[0] && b >= 0 && b < rect->size[1])
			{
				t = rect->size[0] * b + aa;
				rect->pixels[t].curvature = edge->curvature;
				rect->pixels[t] = *col; // *getColour(edge->curvature);
				atlasPixelFaces[t] = edge->face;
				edge->boundaryPixelArray.push_back(t);
			}


			//pnt(a, b, col);
			if ((int)a==101)
			{
				std::cout << "(" << a << ", " << b << ")" << endl;
			}

			if (a >= 0 && a < rect->size[0] && b >= 0 && b < rect->size[1])
			{
				t = rect->size[0] * b + a;
				rect->pixels[t].curvature = edge->curvature;
				rect->pixels[t] = *col; // *getColour(edge->curvature);
				atlasPixelFaces[t] = edge->face;
				edge->boundaryPixelArray.push_back(t);
			}
			
		}
	}
	

	// y-axis pixel cross
	a0 = 1; 
	a1 = 0; 
	n = 0;
	if (y0 < y1)
	{ 
		a0 = ceil(y0); 
		a1 = floor(y1); 
		d = (x1 - x0) / (y1 - y0); 
		a = a0; 
		b = x0 + (a0 - y0)*d; 
		n = fabs(a1 - a0);
	}
	else if (y0 > y1)
	{
		a0 = ceil(y1); 
		a1 = floor(y0);
		d = (x1 - x0) / (y1 - y0);
		a = a0;
		b = x1 + (a0 - y1) * d;
		n = fabs(a1 - a0); 
	}

	if (a0 <= a1)
	{
		for (aa = a, i = 0; i <= n; i++, aa = a, a++, b += d)
		{
			//pnt(b, aa, col);
			if ((int)b == 101)
				std::cout << "(" << b << ", " << aa << ")" << endl;

			if (b >= 0 && b < rect->size[0] && aa >= 0 && aa < rect->size[1])
			{
				t = rect->size[0] * aa + b;
				rect->pixels[t].curvature = edge->curvature;
				rect->pixels[t] = *col; // *getColour(edge->curvature);
				atlasPixelFaces[t] = edge->face;
				edge->boundaryPixelArray.push_back(t);
			}
			

			//pnt(b, a, col);
			if ((int)b == 101)
				std::cout << "(" << b << ", " << a << ")" << endl;

			if (b >= 0 && b < rect->size[0] && a >= 0 && a < rect->size[1])
			{
				t = rect->size[0] * a + b;
				rect->pixels[t].curvature = edge->curvature;
				rect->pixels[t] = *col; // *getColour(edge->curvature);
				atlasPixelFaces[t] = edge->face;
				edge->boundaryPixelArray.push_back(t);
			}

		}
	}
	
}

bool compare(const Patch &p1, const Patch &p2){
	if (p1.height == p2.height){
		return p1.width > p2.width;
	}
	return p1.height > p2.height;
}

vector<int> sort(vector<Patch> &p)
{
	//returns the widest rectangle width

	int w = 0;
	int area = 0;
	int widthSum = 0;

	list<Patch> l;
	for (int i = 0; i<p.size(); i++)
	{
		if (p[i].width>w)
			w = p[i].width;
		area += p[i].width*p[i].height;
		widthSum += p[i].width;
		l.push_back(p[i]);
	}
	l.sort(compare);
	int i = 0;
	for (list<Patch>::iterator it = l.begin(); it != l.end(); ++it){
		p[i] = *it;
		i++;
	}

	vector<int> temp;
	temp.push_back(w);
	temp.push_back(area);
	temp.push_back(widthSum);
	return temp;

}

void modify(vector<Patch> &p, int minH)
{
	for (int i = 0; i<p.size(); i++){
		p[i].y = minH - p[i].height - p[i].y;	// to get the bottom-left offset
	}
}

int t_value(RectSpace* rect, int x, int y)
{
	return  rect->size[0] * y + x;
}

void growOnePixel(RectSpace* rect)
{
	Colour* green = new Colour();
	green->r = 0;
	green->g = 255;
	green->b = 0;

	Colour* magenta = new Colour();
	magenta->r = 255;
	magenta->g = 0;
	magenta->b = 255;

	Colour* white = new Colour();
	white->r = 255;
	white->g = 255;
	white->b = 255;

	int x, y, t, t_l, t_r, t_lt, t_rt, t_lb, t_rb, t_d, t_u;
	t = 0; t_l = 0; t_r = 0; t_lt = 0; t_rt = 0; t_lb = 0; t_rb = 0; t_u = 0; t_d = 0;
	
	bool l, r, u, d, lt, rt, lb, rb;

	for (int y = 0; y < rect->size[0]; y++)
	{
		for (int x = 0; x < rect->size[0]; x++)
		{
			l = false;
			r = false;
			u = false;
			d = false;
			lb = false;
			rb = false;
			lt = false;
			rt = false;

			t = t_value(rect, x, y);
			if (x != 0)
			{
				t_l = t_value(rect, x - 1, y);
				if (colourCompare(rect->pixels[t_l], *magenta))
					l = true;
			}
			if (x != rect->size[0] - 1)
			{
				t_r = t_value(rect, x + 1, y);
				if (colourCompare(rect->pixels[t_r], *magenta))
					r = true;
			}
			if ((x != 0) && (y != rect->size[1] - 1))
			{
				t_lt = t_value(rect, x - 1, y + 1);
				if (colourCompare(rect->pixels[t_lt], *magenta))
					lt = true;
				
			}
			if ((x != rect->size[0] - 1) && (y != rect->size[1] - 1) )
			{
				t_rt = t_value(rect, x + 1, y + 1);
				if (colourCompare(rect->pixels[t_rt], *magenta))
					rt = true;
			}
			if ((x != 0) && (y != 0))
			{
				t_lb = t_value(rect, x - 1, y - 1);
				if (colourCompare(rect->pixels[t_lb], *magenta))
					lb = true;
			}
			if ((x != rect->size[0] - 1) && (y != 0))
			{
				t_rb = t_value(rect, x + 1, y - 1);
				if (colourCompare(rect->pixels[t_rb], *magenta))
					rb = true;
			}
			if (y != rect->size[1] - 1)
			{
				t_u = t_value(rect, x, y + 1);
				if (colourCompare(rect->pixels[t_u], *magenta))
					u = true;
			}
			if (y != 0)
			{
				t_d = t_value(rect, x, y - 1);
				if (colourCompare(rect->pixels[t_d], *magenta))
					d = true;
			}
			
			if (colourCompare(rect->pixels[t], *white) && (l || r || u || d || lt || rt || lb || rb))
			{
				rect->pixels[t] = *green;
			}
		}
		
	}

	for (int y = 0; y < rect->size[0]; y++)
	{
		for (int x = 0; x < rect->size[0]; x++)
		{
			t = t_value(rect, x, y);
			if (colourCompare(rect->pixels[t], *green))
			{
				rect->pixels[t].r = 255;
				rect->pixels[t].g = 0;
				rect->pixels[t].b = 255;
			}
		}
	}
}


//this curvature is just the angle between two faces across an edge
void curvature(HalfEdgeList* walker1)
{
	while (walker1->next)
	{
		if (walker1->edge->curvature == 0)
		{
			if (walker1->edge->twin != NULL)
			{
				walker1->edge->curvature = (float)walker1->edge->face->normal.Dot(walker1->edge->twin->face->normal);
				walker1->edge->twin->curvature = walker1->edge->curvature;
			}
		}
		
		if (walker1->edge->next->curvature == 0)
		{
			if (walker1->edge->next->twin != NULL)
			{
				walker1->edge->next->curvature = (float)walker1->edge->next->face->normal.Dot(walker1->edge->next->twin->face->normal);
				walker1->edge->next->twin->curvature = walker1->edge->next->curvature;
			}
		}

		if (walker1->edge->next->next->curvature == 0)
		{
			if (walker1->edge->next->next->twin != NULL)
			{
				walker1->edge->next->next->curvature = (float) walker1->edge->next->next->face->normal.Dot(walker1->edge->next->next->twin->face->normal);
				walker1->edge->next->next->twin->curvature = walker1->edge->next->next->curvature;
			}
		}

		walker1 = walker1->next;
	}
}

float TriangleArea(DCELVertex P1, DCELVertex P2, DCELVertex P3)
{
	float a, b, c, s;
	a = sqrt(pow((P2.coords.x - P1.coords.x), 2) + pow((P2.coords.y - P1.coords.y), 2) + pow((P2.coords.z - P1.coords.z), 2));
	b = sqrt(pow((P3.coords.x - P2.coords.x), 2) + pow((P3.coords.y - P2.coords.y), 2) + pow((P3.coords.z - P2.coords.z), 2));
	c = sqrt(pow((P1.coords.x - P3.coords.x), 2) + pow((P1.coords.y - P3.coords.y), 2) + pow((P1.coords.z - P3.coords.z), 2));
	s = a + b + c;
	return (sqrt(s*(s - a)*(s - b)*(s - c)));
}


int _tmain(int argc, _TCHAR* argv[])
{
	//vcg operations
	MyMesh m;
	vcg::tri::UpdateCurvature<MyMesh>::VertexIterator vi;
	vcg::tri::UpdateCurvature<MyMesh>::FaceIterator fi;
	vcg::tri::io::ImporterOBJ<MyMesh>::Info oi;

	vcg::tri::io::ImporterOBJ<MyMesh>::Open(m, "sphere.obj", oi);
	//vcg::tri::Torus(m, 30, 10);
	//vcg::tri::io::ExporterOFF<MyMesh>::Save(m, "torus.off");
	vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
	vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);

	vcg::tri::UpdateCurvature<MyMesh>::PrincipalDirections(m);
	vcg::tri::UpdateCurvature<MyMesh>::MeanAndGaussian(m);


	//file reading operations
	ifstream Model("sphere.obj");
	string line, beginningToken;
	int f = 0;

	while (getline(Model, line))
	{
		istringstream iss(line);
		iss >> beginningToken;

		//the read line is vertex position
		if (beginningToken == "v")
		{
			DCELVertex* vertex = new DCELVertex();
			iss >> get<0>(vPosition) >> get<1>(vPosition) >> get<2>(vPosition);
			modelVertices.push_back(vPosition);
			vertex->coords.x = get<0>(vPosition);
			vertex->coords.y = get<1>(vPosition);
			vertex->coords.z = get<2>(vPosition);
			modelVerticesDCEL.push_back(vertex);
		}

		//if beginning token is vertex normal
		if (beginningToken == "vn")
		{
			Vector* normal = new Vector();
			iss >> get<0>(vNormal) >> get<1>(vNormal) >> get<2>(vNormal);
			normal->x = get<0>(vNormal);
			normal->y = get<1>(vNormal);
			normal->z = get<2>(vNormal);
			
			modelNormals.push_back(vNormal);
			modelNormalDCEL.push_back(normal);
		}

		//if beginning token is vertex texture coords
		if (beginningToken == "vt")
		{
			iss >> get<0>(vUV) >> get<1>(vUV);
			modeluvCoords.push_back(vUV);
		}

		//if beginning token is face
		if (beginningToken == "f")
		{
			string temp, token;
			for (int i = 0; i < 3; i++)
			{
				f++;
				iss >> temp;
				int count = 0;
				while (count < 3)
				{
					token = temp.substr(0, temp.find_first_of("/"));
					temp = temp.substr(temp.find_first_of("/") + 1);
					//cout << token << endl;
					if (i == 0)
					{
						if (count == 0)
						{
							get<0>(faceVertexIndex) = stoi(token);
							count++;
							continue;
						}
						if (count == 1)
						{
							get<0>(faceVertexTextureCoordIndex) = stoi(token);
							count++;
							continue;
						}
						if (count == 2)
						{
							get<0>(faceVertexNormalIndex) = stoi(token);
							//cout<<stoi(token)<<endl;
							count++;
							continue;
						}
					}

					if (i == 1)
					{
						if (count == 0)
						{
							get<1>(faceVertexIndex) = stoi(token);
							count++;
							continue;
						}
						if (count == 1)
						{
							get<1>(faceVertexTextureCoordIndex) = stoi(token);
							count++;
							continue;
						}
						if (count == 2)
						{
							get<1>(faceVertexNormalIndex) = stoi(token);
							count++;
							continue;
						}
					}

					if (i == 2)
					{
						if (count == 0)
						{
							get<2>(faceVertexIndex) = stoi(token);
							count++;
							continue;
						}
						if (count == 1)
						{
							get<2>(faceVertexTextureCoordIndex) = stoi(token);
							count++;
							continue;
						}
						if (count == 2)
						{
							get<2>(faceVertexNormalIndex) = stoi(token);
							count++;
							continue;
						}
					}
				}
			}
			modelFaceVertexIndices.push_back(faceVertexIndex);
			modelFaceVertexTextureCoordIndices.push_back(faceVertexTextureCoordIndex);
			modelFaceVertexNormalIndices.push_back(faceVertexNormalIndex);
		}
	}


	Model.close();

	//filling the vertices with curvature values
	int i;
	for (vi = m.vert.begin(),i = 0; vi != m.vert.end(); ++vi, i++) if (!(*vi).IsD())
	{
		modelVerticesDCEL[i]->K1 = (*vi).cK1();
		modelVerticesDCEL[i]->K2 = (*vi).cK2();
		modelVerticesDCEL[i]->Kg = (*vi).cKg();
		
		modelVerticesDCEL[i]->pd1.x = (*vi).cPD1()[0];
		modelVerticesDCEL[i]->pd1.y = (*vi).cPD1()[1];
		modelVerticesDCEL[i]->pd1.z = (*vi).cPD1()[2];

		modelVerticesDCEL[i]->pd2.x = (*vi).cPD2()[0];
		modelVerticesDCEL[i]->pd2.y = (*vi).cPD2()[1];
		modelVerticesDCEL[i]->pd2.z = (*vi).cPD2()[2];
	}


	//calculating and storing face modelNormals for the whole model
	for (vector<int_vec3>::iterator i = modelFaceVertexNormalIndices.begin(); i != modelFaceVertexNormalIndices.end(); i++)
	{
		fNormal = calcFaceNormal(modelNormals[get<0>(*i) - 1], modelNormals[get<1>(*i) - 1], modelNormals[get<2>(*i) - 1]); //indices start from zero in vector
		modelFaceNormals.push_back(fNormal);
	}

	//writing modelFaceNormals to a file
	for (vector<vec3>::iterator i = modelFaceNormals.begin(); i != modelFaceNormals.end(); i++)
	{
		fN_file << "fn " << get<0>(*i) << " " << get<1>(*i) << " " << get<2>(*i) << "\n";
	}

	//find all the faces within acos(1/sqrt(3)) deg of each of the primary axes
	int index = 0;
	for (vector<vec3>::iterator i = modelFaceNormals.begin(); i != modelFaceNormals.end(); i++)
	{
		if (cosVec(x_pos, *i) >= (pow(3, -0.5)))
		{
			x_pos_faces.push_back(index);
			index++;
		}
		else if (cosVec(x_neg, *i) >= (pow(3, -0.5)))
		{
			x_neg_faces.push_back(index);
			index++;
		}
		else if (cosVec(y_pos, *i) >= (pow(3, -0.5)))
		{
			y_pos_faces.push_back(index);
			index++;
		}
		else if (cosVec(y_neg, *i) >= (pow(3, -0.5)))
		{
			y_neg_faces.push_back(index);
			index++;
		}
		else if (cosVec(z_pos, *i) >= (pow(3, -0.5)))
		{
			z_pos_faces.push_back(index);
			index++;
		}
		else if (cosVec(z_neg, *i) >= (pow(3, -0.5)))
		{
			z_neg_faces.push_back(index);
			index++;
		}
		else
		{
			std::cout << "Does not fit anywhere\n";
		}

	}


	////using DCELFace class
	int count = 0;
	for (vector<int_vec3>::iterator i = modelFaceVertexIndices.begin(); i != modelFaceVertexIndices.end(); i++, count++)
	{
		//Create the first vertex and its corresponding edge and fill its entries
		DCELHalfEdge* edge12 = new DCELHalfEdge();
		DCELHalfEdge* edge23 = new DCELHalfEdge();
		DCELHalfEdge* edge31 = new DCELHalfEdge();
		DCELFace* face = new DCELFace();
		
		edge12->origin = modelVerticesDCEL[get<0>(*i) - 1];
		edge12->next = edge23;
		edge12->face = face;
		edge12->origin->normal = modelNormalDCEL[get<0>(modelFaceVertexNormalIndices[count]) - 1];

		edge23->origin = modelVerticesDCEL[get<1>(*i) - 1];
		edge23->next = edge31;
		edge23->face = face;
		edge23->origin->normal = modelNormalDCEL[get<1>(modelFaceVertexNormalIndices[count]) - 1];

		edge31->origin = modelVerticesDCEL[get<2>(*i) - 1];
		edge31->next = edge12;
		edge31->face = face;
		edge31->origin->normal = modelNormalDCEL[get<2>(modelFaceVertexNormalIndices[count]) - 1];

		//face definition
		face->edge = edge12;
		modelFaces.push_back(face);

		//add to list of half edges
		modelEdgeList->addToList(modelEdgeList, edge12);
	}

	//assigning normals from the previous calculation
	count = 0;
	for (vector<vec3>::iterator i = modelFaceNormals.begin(); i != modelFaceNormals.end(); i++, count++)
	{
		modelFaces[count]->normal.x = get<0>(*i);
		modelFaces[count]->normal.y = get<1>(*i);
		modelFaces[count]->normal.z = get<2>(*i);
	}

	//Walk through the modelEdgeList to set the twins for each of the half edges
	//Twins identification
	HalfEdgeList* walker1;
	DCELHalfEdge* walker2 = NULL;
	HalfEdgeList* walker3 = NULL;
	DCELHalfEdge* walker4 = NULL;

	walker1 = new HalfEdgeList();
	walker1 = modelEdgeList;

	while (walker1->next)						//first loop for traversing the halfedgelist
	{
		walker2 = walker1->edge;
		do								//this loop for traversing the edges of the face
		{
			walker3 = modelEdgeList;
			while (walker3->next)
			{
				walker4 = walker3->edge;
				do
				{
					if ((walker2->origin == walker4->next->origin) && (walker4->origin == walker2->next->origin))
					{
						walker2->twin = walker4;
					}

					walker4 = walker4->next;
				} while (walker4 != walker3->edge);
				walker3 = walker3->next;
			}

			walker2 = walker2->next;
		} while (walker2 != walker1->edge);

		walker1 = walker1->next;
	}

	//the data structure has been built by this point



	//Now to categorize the modelFaces into groups according to normal directions
	count = 0;
	for (vector<DCELFace*>::iterator i = modelFaces.begin(); i != modelFaces.end(); i++)
	{
		if (cosVec((*i)->normal, Vector(1.0, 0.0, 0.0)) >= (pow(3, -0.5)))
		{
			(*i)->group = 1;
			continue;
		}
		else if (cosVec((*i)->normal, Vector(-1.0, 0.0, 0.0)) >= (pow(3, -0.5)))
		{
			(*i)->group = 2;
			continue;
		}
		else if (cosVec((*i)->normal, Vector(0.0, 1.0, 0.0)) >= (pow(3, -0.5)))
		{
			(*i)->group = 3;
			continue;
		}
		else if (cosVec((*i)->normal, Vector(0.0, -1.0, 0.0)) >= (pow(3, -0.5)))
		{
			(*i)->group = 4;
			continue;
		}
		else if (cosVec((*i)->normal, Vector(0.0, 0.0, 1.0)) >= (pow(3, -0.5)))
		{
			(*i)->group = 5;
			continue;
		}
		else if (cosVec((*i)->normal, Vector(0.0, 0.0, -1.0)) >= (pow(3, -0.5)))
		{
			(*i)->group = 6;
			continue;
		}
		else
		{
			count++;
		}
	}

	std::cout << "These many points did not resolve into any of the groups:" << count << endl;

	
	//writing to file
	count = 0;
	for (vector<DCELFace*>::iterator i = modelFaces.begin(); i != modelFaces.end(); i++, count++)
	{
		if ((*i)->group == 1)
		{
			fN_file << "fx " << 1 << "\n";
		}
		if ((*i)->group == 2)
		{
			fN_file << "-fx " << 2 << "\n";
		}
		if ((*i)->group == 3)
		{
			fN_file << "fy " << 3 << "\n";
		}
		if ((*i)->group == 4)
		{
			fN_file << "-fy " << 4 << "\n";
		}
		if ((*i)->group == 5)
		{
			fN_file << "fz " << 5 << "\n";
		}
		if ((*i)->group == 6)
		{
			fN_file << "-fz " << 6 << "\n";
		}
	}

	//consume isolated faces in one of the groups
	for (vector<DCELFace*>::iterator i = modelFaces.begin(); i != modelFaces.end(); i++)
	{
		float theta1 = 0;
		float theta2 = 0;
		float theta3 = 0;
		float theta_min = 0;

		if ((*i)->edge->twin != NULL && (*i)->edge->next->twin != NULL && (*i)->edge->next->next->twin != NULL)
		{
			if (((*i)->group != (*i)->edge->twin->face->group) &&
				((*i)->group != (*i)->edge->next->twin->face->group) &&
				((*i)->group != (*i)->edge->next->next->twin->face->group))
			{
				theta1 = (*i)->normal.Dot((*i)->edge->twin->face->normal);
				theta2 = (*i)->normal.Dot((*i)->edge->next->twin->face->normal);
				theta3 = (*i)->normal.Dot((*i)->edge->next->next->twin->face->normal);

				if (theta1 <= theta2 && theta1 <= theta3)
				{
					theta_min = theta1;
					if (theta_min < 1.0)
						(*i)->group = (*i)->edge->twin->face->group;
				}
				else if (theta2 <= theta1 && theta2 <= theta3)
				{
					theta_min = theta2;
					if (theta_min < 1.0)
						(*i)->group = (*i)->edge->next->twin->face->group;
				}
				else if (theta3 <= theta1 && theta3 <= theta2)
				{
					theta_min = theta1;
					if (theta_min < 1.0)
						(*i)->group = (*i)->edge->next->next->twin->face->group;
				}
			
			}
		}
		else if ((*i)->edge->twin == NULL && (*i)->edge->next->twin != NULL && (*i)->edge->next->next->twin != NULL)
		{
			
			if (((*i)->group != (*i)->edge->next->twin->face->group) &&
				((*i)->group != (*i)->edge->next->next->twin->face->group))
			{
				theta2 = (*i)->normal.Dot((*i)->edge->next->twin->face->normal);
				theta3 = (*i)->normal.Dot((*i)->edge->next->next->twin->face->normal);
				
				if (theta2 <= theta3 && theta2 < 1.0)
					(*i)->group = (*i)->edge->next->twin->face->group;
				else if (theta2 > theta3 && theta3 < 1.0)
					(*i)->group = (*i)->edge->next->next->twin->face->group;
			}
		}
		else if ((*i)->edge->twin != NULL && (*i)->edge->next->twin == NULL && (*i)->edge->next->next->twin != NULL)
		{
			if (((*i)->group != (*i)->edge->twin->face->group) &&
				((*i)->group != (*i)->edge->next->next->twin->face->group))
			{
				theta1 = (*i)->normal.Dot((*i)->edge->twin->face->normal);
				theta3 = (*i)->normal.Dot((*i)->edge->next->next->twin->face->normal);
				
				if (theta1 <= theta3 && theta1 < 1.0)
					(*i)->group = (*i)->edge->twin->face->group;
				else if (theta1 > theta3 && theta3 < 1.0)
					(*i)->group = (*i)->edge->next->next->twin->face->group;
			}
		}
		else if ((*i)->edge->twin != NULL && (*i)->edge->next->twin != NULL && (*i)->edge->next->next->twin == NULL)
		{
			if (((*i)->group != (*i)->edge->twin->face->group) &&
				((*i)->group != (*i)->edge->next->twin->face->group))
			{
				theta1 = (*i)->normal.Dot((*i)->edge->twin->face->normal);
				theta2 = (*i)->normal.Dot((*i)->edge->next->twin->face->normal);

				if (theta1 <= theta2 && theta1 < 1.0)
					(*i)->group = (*i)->edge->twin->face->group;
				else if (theta1 > theta2 && theta2 < 1.0)
					(*i)->group = (*i)->edge->next->twin->face->group;
			}
		}
	}


	int patch = 0;
	//inside groups to find patches and separate them
	for (vector<DCELFace*>::iterator i = modelFaces.begin(); i != modelFaces.end(); i++)
	{
		if ((*i)->patchNo == 0)
		{
			patch++;
			patchIdentify((*i), &patch);
		}
	}
	
	
	fN_file << "Total No. of Patches in the model: " << patch<<endl;
	//calculate the projected faces of all the faces
	meshExtremes();
	setBLtoZero();
	stUnitize();
	getProject(modelEdgeList, patch);
	createRect(patch);
	sortRectList();
	boundaryCheck();
	//curvature(modelEdgeList);

	//determining atlas size for more atlas space utilization
	
	int area_sum = 0;
	for (int j = 0; j < rectList.size(); j++)
	{
		area_sum = area_sum + rectList[j]->size[0] * rectList[j]->size[1];
	}

	std::cout << "area_sum: " << area_sum << endl;
	int atlas_size = 0;
	for (int i = 0; ; i++)
	{
		if (i*i > 1.8 * area_sum && i%4 == 0)
		{
			atlas_size = i;
			break;
		}
	}
	std::cout << "atlas size: " << atlas_size << endl;
	
	//create atlas rectangle with black background
	//RectSpace* AtlasRect = new RectSpace(4 * PIXEL_DENSITY, 4 * PIXEL_DENSITY);
	RectSpace* AtlasRect = new RectSpace(atlas_size, atlas_size);
	vector<DCELFace*> atlasPixelFaces;
	Colour* black = new Colour();
	black->a = 1;
	black->b = 0;
	black->g = 0;
	black->r = 0;

	Colour* white = new Colour();
	white->r = 255;
	white->g = 255;
	white->b = 255;

	Colour* yellow = new Colour();
	yellow->r = 255;
	yellow->g = 255;
	yellow->b = 0;

	for (int j = 0; j < AtlasRect->size[1]; j++)
	{
		for (int k = 0; k < AtlasRect->size[0]; k++)
		{
			AtlasRect->pixels[AtlasRect->size[0] * j + k] = *white;
			AtlasRect->pixels[AtlasRect->size[0] * j + k].visited = false;
			atlasPixelFaces.push_back(NULL);
			AtlasRect->pixels[AtlasRect->size[0] * j + k].K1 = 0;
			AtlasRect->pixels[AtlasRect->size[0] * j + k].K2 = 0;
			AtlasRect->pixels[AtlasRect->size[0] * j + k].Kg = 0;
			AtlasRect->pixels[AtlasRect->size[0] * j + k].cos_theta_pd1 = 0;
			AtlasRect->pixels[AtlasRect->size[0] * j + k].cos_theta_pd2 = 0;
		}
	}

	////
	// Header for the ply file
	////

	//fill the atlas with the sorted rectList
	int rectListSize = rectList.size();
	int x_offset = 0;
	int y_offset = 0;
	int previousRowStart = 0;
	for (int i = 0; i < rectListSize; i++)
	{
		if (x_offset + rectList[i]->size[0] > AtlasRect->size[0])
		{
			x_offset = 0;
			y_offset += rectList[previousRowStart]->size[1];
			previousRowStart = i;
		}
		
		float temp[6];
		//updating the s and t's according to atlas
		for (vector<DCELFace*>::iterator p = rectList[i]->facesArray.begin(); p != rectList[i]->facesArray.end(); p++)
		{
			temp[0] = (((*p)->st[0] - rectList[i]->minmax[0]) * PIXEL_DENSITY + x_offset + 2);
			temp[1] = (((*p)->st[1] - rectList[i]->minmax[2]) * PIXEL_DENSITY + y_offset + 2);
			temp[2] = (((*p)->st[2] - rectList[i]->minmax[0]) * PIXEL_DENSITY + x_offset + 2);
			temp[3] = (((*p)->st[3] - rectList[i]->minmax[2]) * PIXEL_DENSITY + y_offset + 2);
			temp[4] = (((*p)->st[4] - rectList[i]->minmax[0]) * PIXEL_DENSITY + x_offset + 2);
			temp[5] = (((*p)->st[5] - rectList[i]->minmax[2]) * PIXEL_DENSITY + y_offset + 2);
			
			(*p)->st[0] =  temp[0] / AtlasRect->size[0];
			(*p)->st[1] =  temp[1] / AtlasRect->size[0];
			(*p)->st[2] =  temp[2] / AtlasRect->size[0];
			(*p)->st[3] =  temp[3] / AtlasRect->size[0];
			(*p)->st[4] =  temp[4] / AtlasRect->size[0];
			(*p)->st[5] =  temp[5] / AtlasRect->size[0];

			AtlasRect->facesArray.push_back(*p);
		}

		x_offset = x_offset + rectList[i]->size[0];

	}
	
	//----------------------------------------------------------------------------------------------
	int frame = 0;
	int t, width;
	float x_min, y_min, x_max, y_max, x, y, z;
	float bary_lambda[3];
	float denom;
	vector<vec2> correspondence;
	int dataPts_count = 0;

	Colour* green = new Colour();
	green->r = 0;// rand() % 255;
	green->g = 255;// rand() % 255;
	green->b = 0;// rand() % 255;
	green->a = 1;

	int s = 0;
	
	Colour* magenta = new Colour();
	magenta->r = 255;// rand() % 255;
	magenta->g = 0; //rand() % 255;
	magenta->b = 255; //rand() % 255;
	magenta->a = 1;

	Colour* red = new Colour();
	red->r = 255;// rand() % 255;
	red->g = 0; //rand() % 255;
	red->b = 0; //rand() % 255;
	red->a = 1;

	Colour* blue = new Colour();
	blue->r = 0;// rand() % 255;
	blue->g = 0; //rand() % 255;
	blue->b = 255; //rand() % 255;
	blue->a = 1;

	Colour* cyan = new Colour();
	cyan->r = 0;// rand() % 255;
	cyan->g = 255; //rand() % 255;
	cyan->b = 255; //rand() % 255;
	cyan->a = 1;

	std::cout << "Rasterizing...\n";
	
	for (vector<DCELFace*>::iterator j = modelFaces.begin(); j != modelFaces.end(); j++)
	{
		//triangle rasterization
		float x_min, y_min, x_max, y_max;

		x_min = AtlasRect->size[0] * (((*j)->st[0] < (*j)->st[2]) ? ((*j)->st[0] < (*j)->st[4] ? ((*j)->st[0]) : ((*j)->st[4])) : ((*j)->st[2] < (*j)->st[4] ? ((*j)->st[2]) : ((*j)->st[4])));
		y_min = AtlasRect->size[0] * (((*j)->st[1] < (*j)->st[3]) ? ((*j)->st[1] < (*j)->st[5] ? ((*j)->st[1]) : ((*j)->st[5])) : ((*j)->st[3] < (*j)->st[5] ? ((*j)->st[3]) : ((*j)->st[5])));

		x_max = AtlasRect->size[0] * (((*j)->st[0] > (*j)->st[2]) ? ((*j)->st[0] > (*j)->st[4] ? ((*j)->st[0]) : ((*j)->st[4])) : ((*j)->st[2] > (*j)->st[4] ? ((*j)->st[2]) : ((*j)->st[4])));
		y_max = AtlasRect->size[0] * (((*j)->st[1] > (*j)->st[3]) ? ((*j)->st[1] > (*j)->st[5] ? ((*j)->st[1]) : ((*j)->st[5])) : ((*j)->st[3] > (*j)->st[5] ? ((*j)->st[3]) : ((*j)->st[5])));

		float triangle[6] = { (*j)->st[0], (*j)->st[1], (*j)->st[2], (*j)->st[3], (*j)->st[4], (*j)->st[5] };

		for (int z = 0; z < 6; z++)
			triangle[z] *= AtlasRect->size[0];

		for (int p = x_min; p <= x_max; p++)
		{
			for (int q = y_min; q <= y_max; q++)
			{
				t = AtlasRect->size[0] * q + p;
				//std::cout << "t: " << t << endl;
				//AtlasRect->WorldCoords[t] = Vector(0, 0, 0);
				if (pointInTriangle(p + 0.5, q + 0.5, triangle))
				{
					atlasPixelFaces[t] = (*j);

					//color selection for each of the patches
					/*if ((*j)->patchNo == 1)
						AtlasRect->pixels[t] = *blue;
					else if ((*j)->patchNo == 2)
						AtlasRect->pixels[t] = *cyan;
					else if ((*j)->patchNo == 3)
						AtlasRect->pixels[t] = *green;
					else if ((*j)->patchNo == 4)
						AtlasRect->pixels[t] = *yellow;
					else if ((*j)->patchNo == 5)
						AtlasRect->pixels[t] = *red;
					else if ((*j)->patchNo == 6)*/
						AtlasRect->pixels[t] = *magenta;

					AtlasRect->pixels[t].n_l.push_back(t - 1);
					AtlasRect->pixels[t].n_r.push_back(t + 1);
					AtlasRect->pixels[t].n_u.push_back(t + AtlasRect->size[0]);
					AtlasRect->pixels[t].n_d.push_back(t - AtlasRect->size[0]);
					AtlasRect->pixels[t].r_l.push_back(1.0);
					AtlasRect->pixels[t].r_r.push_back(1.0);
					AtlasRect->pixels[t].r_u.push_back(1.0);
					AtlasRect->pixels[t].r_d.push_back(1.0);

					//
					correspondence.push_back(vec2(dataPts_count, t));
					dataPts_count++;
					//

					//calculate the barycentric coordinates for each pixel
					denom = 1 / ((triangle[3] - triangle[5]) * (triangle[0] - triangle[4]) + (triangle[4] - triangle[2]) * (triangle[1] - triangle[5]));
					bary_lambda[0] = ((triangle[3] - triangle[5]) * (p - triangle[4]) + (triangle[4] - triangle[2]) * (q - triangle[5])) * denom;
					bary_lambda[1] = ((triangle[5] - triangle[1]) * (p - triangle[4]) + (triangle[0] - triangle[4]) * (q - triangle[5])) * denom;
					bary_lambda[2] = 1 - bary_lambda[0] - bary_lambda[1];

					//std::cout << bary_lambda[0] << " " << bary_lambda[1]<<" " << bary_lambda[2] << endl;

					x = bary_lambda[0] * (*j)->edge->origin->coords.x + bary_lambda[1] * (*j)->edge->next->origin->coords.x + bary_lambda[2] * (*j)->edge->next->next->origin->coords.x;
					y = bary_lambda[0] * (*j)->edge->origin->coords.y + bary_lambda[1] * (*j)->edge->next->origin->coords.y + bary_lambda[2] * (*j)->edge->next->next->origin->coords.y;
					z = bary_lambda[0] * (*j)->edge->origin->coords.z + bary_lambda[1] * (*j)->edge->next->origin->coords.z + bary_lambda[2] * (*j)->edge->next->next->origin->coords.z;
					x *= mesh_max;
					y *= mesh_max;
					z *= mesh_max;

					//updating per pixel curvatures
					AtlasRect->pixels[t].K1 = bary_lambda[0] * (*j)->edge->origin->K1 + bary_lambda[1] * (*j)->edge->next->origin->K1 + bary_lambda[2] * (*j)->edge->next->next->origin->K1;

					AtlasRect->pixels[t].K2 = bary_lambda[0] * (*j)->edge->origin->K2 + bary_lambda[1] * (*j)->edge->next->origin->K2 + bary_lambda[2] * (*j)->edge->next->next->origin->K2;

					AtlasRect->pixels[t].Kg = bary_lambda[0] * (*j)->edge->origin->Kg + bary_lambda[1] * (*j)->edge->next->origin->Kg + bary_lambda[2] * (*j)->edge->next->next->origin->Kg;

					//std::cout << x << " " << y << " " << z << endl;
					AtlasRect->WorldCoords[t] = Vector(x, y, z);

					//interpolate the principal directions
					if (atlasPixelFaces[t]->group == 1 || atlasPixelFaces[t]->group == 2)
					{
						AtlasRect->pixels[t].pd1.x = 0;
						AtlasRect->pixels[t].pd1.y = bary_lambda[0] * (*j)->edge->origin->pd1.y + bary_lambda[1] * (*j)->edge->next->origin->pd1.y + bary_lambda[2] * (*j)->edge->next->next->origin->pd1.y;
						AtlasRect->pixels[t].pd1.z = bary_lambda[0] * (*j)->edge->origin->pd1.z + bary_lambda[1] * (*j)->edge->next->origin->pd1.z + bary_lambda[2] * (*j)->edge->next->next->origin->pd1.z;

						AtlasRect->pixels[t].pd2.x = 0;
						AtlasRect->pixels[t].pd2.y = bary_lambda[0] * (*j)->edge->origin->pd2.y + bary_lambda[1] * (*j)->edge->next->origin->pd2.y + bary_lambda[2] * (*j)->edge->next->next->origin->pd2.y;
						AtlasRect->pixels[t].pd2.z = bary_lambda[0] * (*j)->edge->origin->pd2.z + bary_lambda[1] * (*j)->edge->next->origin->pd2.z + bary_lambda[2] * (*j)->edge->next->next->origin->pd2.z;

						AtlasRect->pixels[t].pd1.normalize();
						AtlasRect->pixels[t].pd2.normalize();

						AtlasRect->pixels[t].cos_theta_pd1 = AtlasRect->pixels[t].pd1[1];
						AtlasRect->pixels[t].cos_theta_pd2 = AtlasRect->pixels[t].pd2[1];
					}
					else if (atlasPixelFaces[t]->group == 3 || atlasPixelFaces[t]->group == 4)
					{
						AtlasRect->pixels[t].pd1.x = bary_lambda[0] * (*j)->edge->origin->pd1.x + bary_lambda[1] * (*j)->edge->next->origin->pd1.x + bary_lambda[2] * (*j)->edge->next->next->origin->pd1.x;
						AtlasRect->pixels[t].pd1.y = 0;
						AtlasRect->pixels[t].pd1.z = bary_lambda[0] * (*j)->edge->origin->pd1.z + bary_lambda[1] * (*j)->edge->next->origin->pd1.z + bary_lambda[2] * (*j)->edge->next->next->origin->pd1.z;

						AtlasRect->pixels[t].pd2.x = bary_lambda[0] * (*j)->edge->origin->pd2.x + bary_lambda[1] * (*j)->edge->next->origin->pd2.x + bary_lambda[2] * (*j)->edge->next->next->origin->pd2.x;
						AtlasRect->pixels[t].pd2.y = 0;
						AtlasRect->pixels[t].pd2.z = bary_lambda[0] * (*j)->edge->origin->pd2.z + bary_lambda[1] * (*j)->edge->next->origin->pd2.z + bary_lambda[2] * (*j)->edge->next->next->origin->pd2.z;

						AtlasRect->pixels[t].pd1.normalize();
						AtlasRect->pixels[t].pd2.normalize();

						AtlasRect->pixels[t].cos_theta_pd1 = AtlasRect->pixels[t].pd1[2];
						AtlasRect->pixels[t].cos_theta_pd2 = AtlasRect->pixels[t].pd2[2];
					}
					else if (atlasPixelFaces[t]->group == 5 || atlasPixelFaces[t]->group == 6)
					{
						AtlasRect->pixels[t].pd1.x = bary_lambda[0] * (*j)->edge->origin->pd1.x + bary_lambda[1] * (*j)->edge->next->origin->pd1.x + bary_lambda[2] * (*j)->edge->next->next->origin->pd1.x;
						AtlasRect->pixels[t].pd1.y = bary_lambda[0] * (*j)->edge->origin->pd1.y + bary_lambda[1] * (*j)->edge->next->origin->pd1.y + bary_lambda[2] * (*j)->edge->next->next->origin->pd1.y;
						AtlasRect->pixels[t].pd1.z = 0;

						AtlasRect->pixels[t].pd2.x = bary_lambda[0] * (*j)->edge->origin->pd2.x + bary_lambda[1] * (*j)->edge->next->origin->pd2.x + bary_lambda[2] * (*j)->edge->next->next->origin->pd2.x;
						AtlasRect->pixels[t].pd2.y = bary_lambda[0] * (*j)->edge->origin->pd2.y + bary_lambda[1] * (*j)->edge->next->origin->pd2.y + bary_lambda[2] * (*j)->edge->next->next->origin->pd2.y;
						AtlasRect->pixels[t].pd2.z = 0;

						AtlasRect->pixels[t].pd1.normalize();
						AtlasRect->pixels[t].pd2.normalize();

						AtlasRect->pixels[t].cos_theta_pd1 = AtlasRect->pixels[t].pd1[0];
						AtlasRect->pixels[t].cos_theta_pd2 = AtlasRect->pixels[t].pd2[0];
					}
				}

			}
		}
	}

	std::cout << "Edge Rasterization...\n";
	for (vector<DCELFace*>::iterator j = modelFaces.begin(); j != modelFaces.end(); j++)
	{
		//edge raterization
		float x1, y1, x2, y2;
		//rasterizing the edges
		DCELHalfEdge* thisEdge = (*j)->edge;
		for (int i = 0; i < 3; i++)
		{
			if (thisEdge->isBoundary == true && thisEdge->isProcessed == false)
			{
				//comparing the edge twin with the three edges of the twin face to find the correct st's
				//for the edge on first face
				if (i == 0)
				{
					x1 = (((*j)->st[0])) * AtlasRect->size[0];
					y1 = (((*j)->st[1])) * AtlasRect->size[0];
					x2 = (((*j)->st[2])) * AtlasRect->size[0];
					y2 = (((*j)->st[3])) * AtlasRect->size[0];
					
					//color selection for each of the patches
					/*if (thisEdge->face->patchNo == 1)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, blue, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 2)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, cyan, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 3)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, green, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 4)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, yellow, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 5)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, red, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 6)*/
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, magenta, atlasPixelFaces);

					//DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, magenta, atlasPixelFaces);
					thisEdge->isProcessed = true;
					if (thisEdge->twin != NULL && thisEdge->twin->isProcessed == false)
					{
						if ((thisEdge->next->origin == thisEdge->twin->face->edge->origin) &&
							(thisEdge->origin == thisEdge->twin->face->edge->next->origin))
						{
							x1 = thisEdge->twin->face->st[0] * AtlasRect->size[0];
							y1 = thisEdge->twin->face->st[1] * AtlasRect->size[0];
							x2 = thisEdge->twin->face->st[2] * AtlasRect->size[0];
							y2 = thisEdge->twin->face->st[3] * AtlasRect->size[0];

							/*if (thisEdge->twin->face->patchNo == 1)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, blue, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 2)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, cyan, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 3)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, green, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 4)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, yellow, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 5)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, red, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 6)*/
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);

							//DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);
							thisEdge->twin->isProcessed = true;
						}

						else if ((thisEdge->next->origin == thisEdge->twin->face->edge->next->origin) &&
							(thisEdge->origin == thisEdge->twin->face->edge->next->next->origin))
						{
							x1 = thisEdge->twin->face->st[2] * AtlasRect->size[0];
							y1 = thisEdge->twin->face->st[3] * AtlasRect->size[0];
							x2 = thisEdge->twin->face->st[4] * AtlasRect->size[0];
							y2 = thisEdge->twin->face->st[5] * AtlasRect->size[0];

							/*if (thisEdge->twin->face->patchNo == 1)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, blue, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 2)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, cyan, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 3)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, green, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 4)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, yellow, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 5)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, red, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 6)*/
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);

							//DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);
							thisEdge->twin->isProcessed = true;
						}
						else if ((thisEdge->next->origin == thisEdge->twin->face->edge->next->next->origin) &&
							(thisEdge->origin == thisEdge->twin->face->edge->origin))
						{
							x1 = thisEdge->twin->face->st[4] * AtlasRect->size[0];
							y1 = thisEdge->twin->face->st[5] * AtlasRect->size[0];
							x2 = thisEdge->twin->face->st[0] * AtlasRect->size[0];
							y2 = thisEdge->twin->face->st[1] * AtlasRect->size[0];
							
							/*if (thisEdge->twin->face->patchNo == 1)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, blue, atlasPixelFaces);//--
							else if (thisEdge->twin->face->patchNo == 2)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, cyan, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 3)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, green, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 4)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, yellow, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 5)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, red, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 6)*/
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);

							//DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);
							thisEdge->twin->isProcessed = true;
						}
					}
				}
				else if (i == 1)
				{
					x1 = (((*j)->st[2])) * AtlasRect->size[0];
					y1 = (((*j)->st[3])) * AtlasRect->size[0];
					x2 = (((*j)->st[4])) * AtlasRect->size[0];
					y2 = (((*j)->st[5])) * AtlasRect->size[0];
					
					//color selection for each of the patches
					/*if (thisEdge->face->patchNo == 1)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, blue, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 2)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, cyan, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 3)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, green, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 4)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, yellow, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 5)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, red, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 6)*/
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, magenta, atlasPixelFaces);

					thisEdge->isProcessed = true;

					if (thisEdge->twin != NULL && thisEdge->twin->isProcessed == false)
					{
						if ((thisEdge->next->origin == thisEdge->twin->face->edge->origin) &&
							(thisEdge->origin == thisEdge->twin->face->edge->next->origin))
						{
							x1 = thisEdge->twin->face->st[0] * AtlasRect->size[0];
							y1 = thisEdge->twin->face->st[1] * AtlasRect->size[0];
							x2 = thisEdge->twin->face->st[2] * AtlasRect->size[0];
							y2 = thisEdge->twin->face->st[3] * AtlasRect->size[0];

							/*if (thisEdge->twin->face->patchNo == 1)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, blue, atlasPixelFaces);//---major
							else if (thisEdge->twin->face->patchNo == 2)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, cyan, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 3)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, green, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 4)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, yellow, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 5)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, red, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 6)*/
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);

							//DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);
							thisEdge->twin->isProcessed = true;
						}
						else if ((thisEdge->next->origin == thisEdge->twin->face->edge->next->origin) &&
							(thisEdge->origin == thisEdge->twin->face->edge->next->next->origin))
						{
							x1 = thisEdge->twin->face->st[2] * AtlasRect->size[0];
							y1 = thisEdge->twin->face->st[3] * AtlasRect->size[0];
							x2 = thisEdge->twin->face->st[4] * AtlasRect->size[0];
							y2 = thisEdge->twin->face->st[5] * AtlasRect->size[0];

							/*if (thisEdge->twin->face->patchNo == 1)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, blue, atlasPixelFaces);//--
							else if (thisEdge->twin->face->patchNo == 2)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, cyan, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 3)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, green, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 4)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, yellow, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 5)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, red, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 6)*/
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);
							//DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);
							thisEdge->twin->isProcessed = true;
						}
						else if ((thisEdge->next->origin == thisEdge->twin->face->edge->next->next->origin) &&
							(thisEdge->origin == thisEdge->twin->face->edge->origin))
						{
							x1 = thisEdge->twin->face->st[4] * AtlasRect->size[0];
							y1 = thisEdge->twin->face->st[5] * AtlasRect->size[0];
							x2 = thisEdge->twin->face->st[0] * AtlasRect->size[0];
							y2 = thisEdge->twin->face->st[1] * AtlasRect->size[0];
							
							/*if (thisEdge->twin->face->patchNo == 1)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, blue, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 2)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, cyan, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 3)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, green, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 4)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, yellow, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 5)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, red, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 6)*/
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);

							//DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);
							thisEdge->twin->isProcessed = true;
						}
					}
				}
				else if (i == 2)
				{
					x1 = (((*j)->st[4])) * AtlasRect->size[0];
					y1 = (((*j)->st[5])) * AtlasRect->size[0];
					x2 = (((*j)->st[0])) * AtlasRect->size[0];
					y2 = (((*j)->st[1])) * AtlasRect->size[0];

					//color selection for each of the patches
					/*if (thisEdge->face->patchNo == 1)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, blue, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 2)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, cyan, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 3)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, green, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 4)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, yellow, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 5)
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, red, atlasPixelFaces);
					else if (thisEdge->face->patchNo == 6)*/
						DDA_line(x1, y1, x2, y2, thisEdge, AtlasRect, magenta, atlasPixelFaces);

					thisEdge->isProcessed = true;

					if (thisEdge->twin != NULL && thisEdge->twin->isProcessed == false)
					{
						if ((thisEdge->next->origin == thisEdge->twin->face->edge->origin) &&
							(thisEdge->origin == thisEdge->twin->face->edge->next->origin))
						{
							x1 = thisEdge->twin->face->st[0] * AtlasRect->size[0];
							y1 = thisEdge->twin->face->st[1] * AtlasRect->size[0];
							x2 = thisEdge->twin->face->st[2] * AtlasRect->size[0];
							y2 = thisEdge->twin->face->st[3] * AtlasRect->size[0];
							
							/*if (thisEdge->twin->face->patchNo == 1)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, blue, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 2)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, cyan, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 3)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, green, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 4)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, yellow, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 5)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, red, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 6)*/
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);

							//DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);
							thisEdge->twin->isProcessed = true;
						}
						else if ((thisEdge->next->origin == thisEdge->twin->face->edge->next->origin) &&
							(thisEdge->origin == thisEdge->twin->face->edge->next->next->origin))
						{
							x1 = thisEdge->twin->face->st[2] * AtlasRect->size[0];
							y1 = thisEdge->twin->face->st[3] * AtlasRect->size[0];
							x2 = thisEdge->twin->face->st[4] * AtlasRect->size[0];
							y2 = thisEdge->twin->face->st[5] * AtlasRect->size[0];
							
							/*if (thisEdge->twin->face->patchNo == 1)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, blue, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 2)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, cyan, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 3)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, green, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 4)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, yellow, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 5)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, red, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 6)*/
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);

							//DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);
							thisEdge->twin->isProcessed = true;
						}

						else if ((thisEdge->next->origin == thisEdge->twin->face->edge->next->next->origin) &&
							(thisEdge->origin == thisEdge->twin->face->edge->origin))
						{
							x1 = thisEdge->twin->face->st[4] * AtlasRect->size[0];
							y1 = thisEdge->twin->face->st[5] * AtlasRect->size[0];
							x2 = thisEdge->twin->face->st[0] * AtlasRect->size[0];
							y2 = thisEdge->twin->face->st[1] * AtlasRect->size[0];
							
							/*if (thisEdge->twin->face->patchNo == 1)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, blue, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 2)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, cyan, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 3)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, green, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 4)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, yellow, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 5)
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, red, atlasPixelFaces);
							else if (thisEdge->twin->face->patchNo == 6)*/
								DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);

							//DDA_line(x2, y2, x1, y1, thisEdge->twin, AtlasRect, magenta, atlasPixelFaces);
							thisEdge->twin->isProcessed = true;
						}
					}
				}
			}
			thisEdge = thisEdge->next;
		}

		fN_file << "rect " << AtlasRect << ":" << "faces st: ";
		fN_file << (*j)->st[0] << "  " << (*j)->st[1] << "  " << (*j)->st[2] << "  " << (*j)->st[3] << "  " << (*j)->st[4] << "  " << (*j)->st[5] << endl;
	}

	
	/*
	for (vector<DCELFace*>::iterator j = modelFaces.begin(); j != modelFaces.end(); j++)
	{
		DCELHalfEdge* edgeWalker = (*j)->edge;

		for (int k = 0; k < 3; k++)
		{
			if (edgeWalker->isBoundary == true)
			{
				for (int i = 0; i < edgeWalker->boundaryPixelArray.size(); i++)
				{
					//AtlasRect->pixels[edgeWalker->boundaryPixelArray[i]].r = 0;
					fN_file << "boundary pixel no.: " << edgeWalker->boundaryPixelArray[i]<<" ratios are "
						<< get<0>(edgeWalker->ratios[i])
						<< ", " << get<1>(edgeWalker->ratios[i]) << endl;

					if (get<0>(edgeWalker->ratios[i]) == 0 && get<1>(edgeWalker->ratios[i]) == 0)
					{
						fN_file << "Two zeros\n";
					}

					if (edgeWalker->boundaryPixelArray[i] == 51707)
					{
						std::cout << "After DDA: "<< (int)AtlasRect->pixels[edgeWalker->boundaryPixelArray[i]].r << " "	<<
							(int)AtlasRect->pixels[edgeWalker->boundaryPixelArray[i]].g << " "
							<< (int)AtlasRect->pixels[edgeWalker->boundaryPixelArray[i]].b << endl;
					}
						
				
				}
			}
			edgeWalker = edgeWalker->next;
		}
	}
	*/
	//----------------------------------------------------------------------------------------------
	int s1,t1;
	t1 = 0;
	bool flag = false;
	//write the pixels wc and s,t of atlas pixels into a text file
	for (int i = 0; i < AtlasRect->size[0] * AtlasRect->size[1]; i++)
	{
		s1 = i % AtlasRect->size[0];
		
		if (s1==0 && flag)
		{
			t1++;
		}
		
		flag = true;
		//wc_file << AtlasRect->WorldCoords[i].x << "," << AtlasRect->WorldCoords[i].y << "," << AtlasRect->WorldCoords[i].z << ","<< s1 << "," << t1 <<endl;
	}

	//change the uv coords according to the new calculations
	count = 0;
	for (vector<DCELFace*>::iterator i = modelFaces.begin(); i != modelFaces.end(); i++, count++)
	{
		get<0>(modeluvCoords[get<0>(modelFaceVertexTextureCoordIndices[count]) - 1]) = (*i)->st[0];
		get<1>(modeluvCoords[get<0>(modelFaceVertexTextureCoordIndices[count]) - 1]) = (*i)->st[1];
		
		get<0>(modeluvCoords[get<1>(modelFaceVertexTextureCoordIndices[count]) - 1]) = (*i)->st[2];
		get<1>(modeluvCoords[get<1>(modelFaceVertexTextureCoordIndices[count]) - 1]) = (*i)->st[3];

		get<0>(modeluvCoords[get<2>(modelFaceVertexTextureCoordIndices[count]) - 1]) = (*i)->st[4];
		get<1>(modeluvCoords[get<2>(modelFaceVertexTextureCoordIndices[count]) - 1]) = (*i)->st[5];
	}

	//create an obj for viewing
	for (vector<DCELVertex*>::iterator i = modelVerticesDCEL.begin(); i != modelVerticesDCEL.end(); i++)
	{
		fout_file << "v " << (*i)->coords.x * mesh_max << " " << (*i)->coords.y * mesh_max << " " << (*i)->coords.z * mesh_max << endl;
	}
	
	/*
	for (int i = 0; i < correspondence.size(); i++)
	{
		t = get<1>(correspondence[i]);
		fout_file << "v " << AtlasRect->WorldCoords[t].x << " " << AtlasRect->WorldCoords[t].y << " " << AtlasRect->WorldCoords[t].z << endl;
	}
	*/

	for (vector<vec3>::iterator i = modelNormals.begin(); i != modelNormals.end(); i++)
	{
		fout_file << "vn " << get<0>(*i) << " " << get<1>(*i) << " " << get<2>(*i) << endl;
	}

	for (vector<DCELFace*>::iterator i = modelFaces.begin(); i != modelFaces.end(); i++)
	{
		fout_file << "vt " << (*i)->st[0] << " " << (*i)->st[1] << endl;
		fout_file << "vt " << (*i)->st[2] << " " << (*i)->st[3] << endl;
		fout_file << "vt " << (*i)->st[4] << " " << (*i)->st[5] << endl;
	}

	count = 0;
	for(vector<DCELFace*>::iterator i = modelFaces.begin(); i != modelFaces.end(); i++, count++)
	{
		fout_file << "f " << get<0>(modelFaceVertexIndices[count]) << "/" << 3 * count + 1 << "/" << get<0>(modelFaceVertexNormalIndices[count]) << " ";
		fout_file << get<1>(modelFaceVertexIndices[count]) << "/" << 3 * count + 2 << "/" << get<1>(modelFaceVertexNormalIndices[count]) << " ";
		fout_file << get<2>(modelFaceVertexIndices[count]) << "/" << 3 * count + 3 << "/" << get<2>(modelFaceVertexNormalIndices[count]) << endl;
	}
	


	//growOnePixel(AtlasRect);
	

	//finding neighbours of each pixel at the boundaries
	std::cout << "Finding boundary pixel neighbors...\n";
	//run through each half edge of the model
	//for (vector<DCELFace*>::iterator j = AtlasRect->facesArray.begin(); j != AtlasRect->facesArray.end(); j++)
	bool print = true;
	for (vector<DCELFace*>::iterator j = modelFaces.begin(); j != modelFaces.end(); j++)
	{
		DCELHalfEdge* edgeWalker = (*j)->edge;

		for (int k = 0; k < 3; k++)
		{
			if (edgeWalker->isBoundary == true)
			{
				//std::cout << "Inside\n";
				int thisSideSize = edgeWalker->boundaryPixelArray.size();
				int otherSideSize = 0;

				if (edgeWalker->twin != NULL)
					otherSideSize = edgeWalker->twin->boundaryPixelArray.size();

				bool flag = false;
				int quotient = 1;
				int mod = 0;

				if (thisSideSize == 0 || otherSideSize == 0)
					continue;
				//check for which side has more number of boundary pixels
				if (thisSideSize <= otherSideSize)
				{
					quotient = otherSideSize / thisSideSize;
					flag = true;
				}
				else
				{
					if (otherSideSize != 0)
					{
						quotient = thisSideSize / otherSideSize;
						mod = thisSideSize % otherSideSize;
					}

					flag = false;
				}
				//std::cout << "thisSideSize: " << thisSideSize << endl;
				for (int i = 0; i < thisSideSize; i++)
				{
					int index = edgeWalker->boundaryPixelArray[i];

					if (AtlasRect->pixels[index].visited == true)
						continue;

					AtlasRect->pixels[index].visited = true;


					//if (index == 0) continue;
					//if (AtlasRect->pixels[index].r1)
						//fN_file << "r1 of index: " << index <<" is "<<AtlasRect->pixels[index].r1 << endl;
					//checking the left neighbours
					//if (atlasPixelFaces[index - 1] == NULL || atlasPixelFaces[index - 1]->patchNo != atlasPixelFaces[index]->patchNo)//checking if the left pixel is in the same patch or not i.e. it can be in another patch or no patch
					if (colourCompare(AtlasRect->pixels[index - 1], *white))
					{
						if (edgeWalker->twin != NULL)
						{
							float d1, d2;
							int index_twin;
							bool repetition_check;
							//checking all the neighboring pixels on the other side
							for (int f = 0; f < otherSideSize; f++)
							{
								repetition_check = false;
								index_twin = edgeWalker->twin->boundaryPixelArray[f];
								float r1, r2, r3, r4;
								r1 = get<0>(edgeWalker->ratios[i]);
								r2 = get<1>(edgeWalker->ratios[i]);
								r3 = get<0>(edgeWalker->twin->ratios[f]);
								r4 = get<1>(edgeWalker->twin->ratios[f]);

								if (!((r3 < r1 && r4 <= r1) || (r3 >= r2 && r4 > r2)))
								{
									for (int p = 0; p < AtlasRect->pixels[index].n_l.size(); p++)
									{
										if (AtlasRect->pixels[index].n_l[p] == index_twin)
											repetition_check = true;
									}
									if (!repetition_check)
									{
										float l1, l2, r;
										l1 = r2 - r1;
										if (r3 > r1 && r4 < r2)
										{
											l2 = r4 - r3;
											if (l2 == 0)
												std::cout << "r4 - r3: " << r4 << " - " << r3 << endl;
										}
										else if (r3 < r1 && r4 < r2)
										{
											l2 = r4 - r1;
											if (l2 == 0)
												std::cout << "r4 - r1: " << r4 << " - " << r1 << endl;
										}
										else if (r3 > r1 && r4 > r2)
										{
											l2 = r2 - r3;
											if (l2 == 0)
												std::cout << "r2 - r3: " << r2 << " - " << r3 << endl;
										}
										else if (r3 < r1 && r4 > r2)
										{
											l2 = r2 - r1;
											if (l2 == 0)
												std::cout << "r2 - r1: " << r2 << " - " << r1 << endl;
										}

										r = l2 / l1;
										if (r != 0)
										{
											
											if (r < 0.01);
												
											else if (r > 0.99)
											{
												AtlasRect->pixels[index].n_l.push_back(index_twin);
												AtlasRect->pixels[index].r_l.push_back(1.0);
											}
											else
											{
												AtlasRect->pixels[index].n_l.push_back(index_twin);
												AtlasRect->pixels[index].r_l.push_back(r);
											}
												
										}

										if (r == 0)
										{
											std::cout << "r1: " << r1 << "r2: " << r2 << "r3: " << r3 << "r4: " << r4 << endl;
										}
									}
								}
							}
						}
					}

					else
					{
						bool repetition_check = false;
						for (int p = 0; p < AtlasRect->pixels[index].n_l.size(); p++)
						{
							if (AtlasRect->pixels[index].n_l[p] == index - 1)
								repetition_check = true;
						}
						if (!repetition_check)
						{
							AtlasRect->pixels[index].n_l.push_back(index - 1);
							AtlasRect->pixels[index].r_l.push_back(1);
						}
					}
						

					//checking the right neighbours
					//if (atlasPixelFaces[index + 1] == NULL || atlasPixelFaces[index + 1]->patchNo != atlasPixelFaces[index]->patchNo)//checking if the right pixel is in the same patch or not i.e. it can be in another patch or no patch
					if (colourCompare(AtlasRect->pixels[index + 1], *white))
					{
						//other method
						if (edgeWalker->twin != NULL)
						{
							float d1, d2;
							int index_twin;
							bool repetition_check;
							//checking the closest pixel on other side
							for (int f = 0; f < otherSideSize; f++)
							{
								repetition_check = false;
								index_twin = edgeWalker->twin->boundaryPixelArray[f];
								float r1, r2, r3, r4;
								r1 = get<0>(edgeWalker->ratios[i]);
								r2 = get<1>(edgeWalker->ratios[i]);
								r3 = get<0>(edgeWalker->twin->ratios[f]);
								r4 = get<1>(edgeWalker->twin->ratios[f]);


								if (!((r3 < r1 && r4 <= r1) || (r3 >= r2 && r4 > r2)))
								{
									for (int p = 0; p < AtlasRect->pixels[index].n_r.size(); p++)
									{
										if (AtlasRect->pixels[index].n_r[p] == index_twin)
											repetition_check = true;
									}
									if (!repetition_check)
									{
										float l1, l2, r;
										l1 = r2 - r1;
										if (r3 > r1 && r4 < r2)
										{
											l2 = r4 - r3;
										}
										else if (r3 < r1 && r4 < r2)
										{
											l2 = r4 - r1;
										}
										else if (r3 > r1 && r4 > r2)
										{
											l2 = r2 - r3;
										}
										else if (r3 < r1 && r4 > r2)
										{
											l2 = r2 - r1;
										}

										r = l2 / l1;
										
										if (r != 0)
										{
											if (r < 0.01);
											else if (r > 0.99)
											{
												AtlasRect->pixels[index].n_r.push_back(index_twin);
												AtlasRect->pixels[index].r_r.push_back(1.0);
											}
											else
											{
												AtlasRect->pixels[index].n_r.push_back(index_twin);
												AtlasRect->pixels[index].r_r.push_back(r);
											}
										}

									}
								}
							}
						}
					}
					else
					{
						bool repetition_check = false;
						for (int p = 0; p < AtlasRect->pixels[index].n_r.size(); p++)
						{
							if (AtlasRect->pixels[index].n_r[p] == index + 1)
								repetition_check = true;
						}
						if (!repetition_check)
						{
							AtlasRect->pixels[index].n_r.push_back(index + 1);
							AtlasRect->pixels[index].r_r.push_back(1);
						}
							
					}

					//checking the up neighbours
					//if (atlasPixelFaces[index + AtlasRect->size[0]] == NULL || atlasPixelFaces[index + AtlasRect->size[0]]->patchNo != atlasPixelFaces[index]->patchNo)//checking if the right pixel is in the same patch or not i.e. it can be in another patch or no patch
					if (colourCompare(AtlasRect->pixels[index + AtlasRect->size[0]], *white))
					{
						//other method
						if (edgeWalker->twin != NULL)
						{
							float d1, d2;
							int index_twin;
							bool repetition_check;
							
							//checking the closest pixel on other side
							for (int f = 0; f < otherSideSize; f++)
							{
								repetition_check = false;
								index_twin = edgeWalker->twin->boundaryPixelArray[f];
								float r1, r2, r3, r4;
								r1 = get<0>(edgeWalker->ratios[i]);
								r2 = get<1>(edgeWalker->ratios[i]);
								r3 = get<0>(edgeWalker->twin->ratios[f]);
								r4 = get<1>(edgeWalker->twin->ratios[f]);

								if (!((r3 < r1 && r4 <= r1) || (r3 >= r2 && r4 > r2)))
								{
									for (int p = 0; p < AtlasRect->pixels[index].n_u.size(); p++)
									{
										if (AtlasRect->pixels[index].n_u[p] == index_twin)
											repetition_check = true;
									}
									if (!repetition_check)
									{
										float l1, l2, r;
										l1 = r2 - r1;
										if (r3 > r1 && r4 < r2)
										{
											l2 = r4 - r3;
										}
										else if (r3 < r1 && r4 < r2)
										{
											l2 = r4 - r1;
										}
										else if (r3 > r1 && r4 > r2)
										{
											l2 = r2 - r3;
										}
										else if (r3 < r1 && r4 > r2)
										{
											l2 = r2 - r1;
										}

										r = l2 / l1;
										if (r != 0)
										{
											if (r < 0.01);
											else if (r > 0.99)
											{
												if (index == 454945)
												{
													std::cout << "step 1: \n";
													std::cout << "l2: " << l2 << "l1: " << l1 << endl;
													std::cout << "r1: " << r1 << "r2: " << r2 << "r3: " << r3 << "r4: " << r4 << endl;
													std::cout << "Index_twin: " << index_twin << endl;
													std::cout << "Edge twin: " << edgeWalker->twin<<endl;
												}
												AtlasRect->pixels[index].n_u.push_back(index_twin);
												AtlasRect->pixels[index].r_u.push_back(1.0);
											}
											else
											{
												if (index == 454945)
												{
													std::cout << "step 2: \n";
													std::cout << "l2: " << l2 << "l1: " << l1 << endl;
													std::cout << "r1: " << r1 << "r2: " << r2 << "r3: " << r3 << "r4: " << r4 << endl;
													std::cout << "Index_twin: " << index_twin << endl;
													std::cout << "Edge twin: " << edgeWalker->twin << endl;
												}
												AtlasRect->pixels[index].n_u.push_back(index_twin);
												AtlasRect->pixels[index].r_u.push_back(r);
											}
										}
										
									}
								}
							}
						}
					}
					else
					{
						bool repetition_check = false;
						for (int p = 0; p < AtlasRect->pixels[index].n_u.size(); p++)
						{
							if (AtlasRect->pixels[index].n_u[p] == index + AtlasRect->size[0])
								repetition_check = true;
						}
						if (!repetition_check)
						{
							if (index == 454945)
							{
								std::cout << "step 3: \n";
							}
							AtlasRect->pixels[index].n_u.push_back(index + AtlasRect->size[0]);
							AtlasRect->pixels[index].r_u.push_back(1);
						}
					}

					//checking the down neighbours
					//if (atlasPixelFaces[index - AtlasRect->size[0]] == NULL || atlasPixelFaces[index - AtlasRect->size[0]]->patchNo != atlasPixelFaces[index]->patchNo)//checking if the right pixel is in the same patch or not i.e. it can be in another patch or no patch
					if (colourCompare(AtlasRect->pixels[index - AtlasRect->size[0]], *white))
					{
						//other method
						if (edgeWalker->twin != NULL)
						{
							float d1, d2;
							int index_twin;
							bool repetition_check;
							//checking the closest pixel on other side
							for (int f = 0; f < otherSideSize; f++)
							{
								repetition_check = false;
								index_twin = edgeWalker->twin->boundaryPixelArray[f];
								float r1, r2, r3, r4;
								r1 = get<0>(edgeWalker->ratios[i]);
								r2 = get<1>(edgeWalker->ratios[i]);
								r3 = get<0>(edgeWalker->twin->ratios[f]);
								r4 = get<1>(edgeWalker->twin->ratios[f]);

								if (!((r3 < r1 && r4 <= r1) || (r3 >= r2 && r4 > r2)))
								{
									for (int p = 0; p < AtlasRect->pixels[index].n_d.size(); p++)
									{
										if (AtlasRect->pixels[index].n_d[p] == index_twin)
											repetition_check = true;
									}
									if (!repetition_check)
									{
										float l1, l2, r;
										l1 = r2 - r1;
										if (r3 > r1 && r4 < r2)
										{
											l2 = r4 - r3;
										}
										else if (r3 < r1 && r4 < r2)
										{
											l2 = r4 - r1;
										}
										else if (r3 > r1 && r4 > r2)
										{
											l2 = r2 - r3;
										}
										else if (r3 < r1 && r4 > r2)
										{
											l2 = r2 - r1;
										}

										r = l2 / l1;
										if (r != 0)
										{
											if (r < 0.01);
											else if (r > 0.99)
											{
												AtlasRect->pixels[index].n_d.push_back(index_twin);
												AtlasRect->pixels[index].r_d.push_back(1.0);
											}
											else
											{
												AtlasRect->pixels[index].n_d.push_back(index_twin);
												AtlasRect->pixels[index].r_d.push_back(r);
											}
										}
										
									}
								}
							}
						}
					}
					else
					{
						bool repetition_check = false;
						for (int p = 0; p < AtlasRect->pixels[index].n_d.size(); p++)
						{
							if (AtlasRect->pixels[index].n_d[p] == index - AtlasRect->size[0])
								repetition_check = true;
						}
						if (!repetition_check)
						{
							AtlasRect->pixels[index].n_d.push_back(index - AtlasRect->size[0]);
							AtlasRect->pixels[index].r_d.push_back(1);
						}
					}
				}
			}
			edgeWalker = edgeWalker->next;
		}
	}
	
	/*
	//surface area calculation
	std::cout << "Voxel Size Calculation...\n";
	float surfaceArea = 0;
	float voxelSize;
	for (vector<DCELFace*>::iterator i = modelFaces.begin(); i != modelFaces.end(); i++)
	{
		surfaceArea += TriangleArea(*(*i)->edge->origin, *(*i)->edge->next->origin, *(*i)->edge->next->next->origin);
	}
	fN_file << "Total surface area of model: " << surfaceArea << endl;
	voxelSize = sqrt(surfaceArea/correspondence.size());
	float collision_distance = 2 * 1.732 * voxelSize;
	fN_file << "Voxel Size: " << voxelSize << endl;
	fN_file << "Collision distance: " << collision_distance << endl;
	//----------------ANN operations------------------
	std::cout << "Finding voxel collisions\n";
	//variables for ANN
	int k = 1; // number of nearest neighbors
	int dim = 3; // dimension
	double eps = 0.01; // error bound
	int nPts = dataPts_count; // actual number of data points
	ANNpointArray dataPts; // data points
	ANNpoint queryPt; // query point
	ANNidxArray nnIdx; // near neighbor indices
	ANNdistArray dists; // near neighbor distances
	ANNkd_tree* kdTree; // search structure
	vector<ANNpoint> points;
	vector<int> patchesIntersect;
	int noOfNeighbors = 8;
	queryPt = annAllocPt(dim); // allocate query point
	dataPts = annAllocPts(nPts, dim); // allocate data points
	nnIdx = new ANNidx[k]; // allocate near neigh indices
	dists = new ANNdist[k]; // allocate near neighbor dists

	for (int i = 0; i < correspondence.size(); i++)
	{
		ANNcoord* coord = new ANNcoord[3];
		t = get<1>(correspondence[i]);
		//std::cout << "t:" << t << endl;
		//AtlasRect->pixels[t] = *black;
		coord[0] = AtlasRect->WorldCoords[t].x;
		coord[1] = AtlasRect->WorldCoords[t].y;
		coord[2] = AtlasRect->WorldCoords[t].z;
		//std::cout << coord[0] << " " << coord[1] << " " << coord[2] << endl;
		points.push_back(coord);
	}
	dataPts = &points[0];
	//for (int i = 0; i < correspondence.size(); i++)
	//{
		//std::cout<<*dataPts[i]<<endl;
	//}

	kdTree = new ANNkd_tree( // build search structure
		dataPts, // the data points
		nPts, // number of points
		3); // dimension of space
	//std::cout << "kdtree ban gaya\n";

	for (int i = 0; i < correspondence.size(); i++)
	{
		ANNcoord coord[3];
		t = get<1>(correspondence[i]);
		coord[0] = AtlasRect->WorldCoords[t].x;
		coord[1] = AtlasRect->WorldCoords[t].y;
		coord[2] = AtlasRect->WorldCoords[t].z;
		queryPt = coord;
		//std::cout << coord[0] << " " << coord[1] << " " << coord[2] << endl;
		//std::cout << queryPt[0] << " " << queryPt[1] << " " << queryPt[2] << endl;
		kdTree->annkSearch( // search
			queryPt, // query point
			noOfNeighbors, // number of near neighbors
			nnIdx, // nearest neighbors (returned)
			dists, // distance (returned)
			eps); // error bound
		//std::cout << dists[1] << endl;

		for (int j = 0; j < noOfNeighbors; j++)
		{
			if (dists[j] > 0.0 && dists[j] < collision_distance) // 2*sqrt(3)
			{
				//fN_file << dists[j] << endl;
				int t1 = get<1>(correspondence[nnIdx[j]]); //correspondence is b/w valid texel index and its actual index in atlas
				if (atlasPixelFaces[t]->patchNo != atlasPixelFaces[t1]->patchNo)
				{
					//fN_file << "distance: "<<dists[j] << endl;
					patchesIntersect.push_back(atlasPixelFaces[t]->patchNo);
					patchesIntersect.push_back(atlasPixelFaces[t1]->patchNo);
					//fN_file << "distances: " << dists[j] << endl;
					//AtlasRect->pixels[t] = *black;
					//AtlasRect->pixels[t1] = *black;
				}
			}
		}
	}
	*/

	std::cout << "Writing neighborhood file...\n";
	
	//writing the neighborhood file
	for (int i = 0; i < AtlasRect->size[0] * AtlasRect->size[1]; i++)
	{
		nb_file << AtlasRect->pixels[i].n_l.size();
		for (int j = 0; j < AtlasRect->pixels[i].n_l.size(); j++)
			nb_file << " " << AtlasRect->pixels[i].n_l[j] << " " << AtlasRect->pixels[i].r_l[j];

		nb_file << "\n" << AtlasRect->pixels[i].n_r.size();
		for (int j = 0; j < AtlasRect->pixels[i].n_r.size(); j++)
			nb_file << " " << AtlasRect->pixels[i].n_r[j] << " " <<AtlasRect->pixels[i].r_r[j];

		nb_file << "\n" << AtlasRect->pixels[i].n_u.size();
		for (int j = 0; j < AtlasRect->pixels[i].n_u.size(); j++)
			nb_file << " " << AtlasRect->pixels[i].n_u[j] << " " << AtlasRect->pixels[i].r_u[j];

		nb_file << "\n" << AtlasRect->pixels[i].n_d.size();
		for (int j = 0; j < AtlasRect->pixels[i].n_d.size(); j++)
			nb_file << " " << AtlasRect->pixels[i].n_d[j] << " " << AtlasRect->pixels[i].r_d[j];

		nb_file << endl;
	}

	ofstream curv("curv.txt");
	float max_K1, max_K2, max_Kg;
	max_K1 = 0;
	max_K2 = 0;
	max_Kg = 0;
	
	for (int i = 0; i < AtlasRect->size[0] * AtlasRect->size[1]; i++)
	{
		curv << AtlasRect->pixels[i].K1 << " " << AtlasRect->pixels[i].K2 << " " << AtlasRect->pixels[i].Kg << " "<<
			fabs(AtlasRect->pixels[i].cos_theta_pd1) << " " << fabs(AtlasRect->pixels[i].cos_theta_pd2) << endl;

		if (fabs(AtlasRect->pixels[i].K1) > max_K1)
			max_K1 = AtlasRect->pixels[i].K1;
		if (fabs(AtlasRect->pixels[i].K2) > max_K2)
			max_K2 = AtlasRect->pixels[i].K2;
		if (fabs(AtlasRect->pixels[i].Kg) > max_Kg)
			max_Kg = AtlasRect->pixels[i].Kg;
	}

	std::cout << "Writing Atlas...\n";
	writeImage_bmp(AtlasRect, "Atlas.bmp");

	/*
	//writing curvatures in images
	std::cout << "Writing curvature images...\n";
	RectSpace* curv1 = new RectSpace(AtlasRect->size[0], AtlasRect->size[1]);
	RectSpace* curv2 = new RectSpace(AtlasRect->size[0], AtlasRect->size[1]);
	RectSpace* curvg = new RectSpace(AtlasRect->size[0], AtlasRect->size[1]);

	for (int i = 0; i < curv1->size[0] * curv1->size[1]; i++)
	{
		if (AtlasRect->pixels[i].K1 >= 0)
		{
			curv1->pixels[i].r = fabs(AtlasRect->pixels[i].K1 * 255) / max_K1;
			curv1->pixels[i].g = 0;
			curv1->pixels[i].b = 0;
		}
		else
		{
			curv1->pixels[i].r = 0;
			curv1->pixels[i].g = fabs(AtlasRect->pixels[i].K1 * 255) / max_K1;
		}
			
		curv1->pixels[i].b = 0;
	}

	for (int i = 0; i < curv2->size[0] * curv2->size[1]; i++)
	{
		if (AtlasRect->pixels[i].K2 >= 0)
		{
			curv2->pixels[i].r = fabs(AtlasRect->pixels[i].K2 * 255) / max_K2;
			curv2->pixels[i].g = 0;
		}
		else
		{
			curv2->pixels[i].g = fabs(AtlasRect->pixels[i].K2 * 255) / max_K2;
			curv2->pixels[i].r = 0;
		}
			
		curv2->pixels[i].b = 0;
	}

	for (int i = 0; i < curvg->size[0] * curvg->size[1]; i++)
	{
		if (AtlasRect->pixels[i].Kg >= 0)
		{
			curvg->pixels[i].r = fabs(AtlasRect->pixels[i].Kg * 255) / max_Kg;
			curvg->pixels[i].g = 0;
		}
		else
		{
			curvg->pixels[i].r = 0;
			curvg->pixels[i].g = fabs(AtlasRect->pixels[i].Kg * 255) / max_Kg;
		}
			
		curvg->pixels[i].b = 0;
	}
	
	//writeImage(curv1, "curv1.tga");
	//writeImage(curv2, "curv2.tga");
	//writeImage(curvg, "curvg.tga");
	*/
	std::cout << "done :)";
	fN_file.close();
	fout_file.close();
	outFile.close();
	wc_file.close();
	nb_file.close();
	curv.close();

	getchar();

	return 0;
}