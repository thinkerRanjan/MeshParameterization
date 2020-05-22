#ifndef __IMAGE_SAVER__
#define __IMAGE_SAVER__

//includes
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include"Vector.h"

using namespace std;
class DCELFace;
//data structures
struct Colour 
{
	unsigned char r,g,b,a;
	//long unsigned int neighbor[8];
	vector<int>n_l, n_r, n_u, n_d;
	vector<float> r_l, r_r, r_u, r_d;
	bool visited;
	float curvature, K1, K2, Kg, cos_theta_pd1, cos_theta_pd2;
	Vector pd1, pd2;
};

class TGAImage {

public:

	//Constructor
	TGAImage();
	
	//Destructor
	~TGAImage();

	//Overridden Constructor
	TGAImage(int width, int height);

	//Set all pixels at once
	void setAllPixels(Colour *pixels);

	//set individual pixels
	void setPixel(Colour inputcolor, int xposition, int yposition);

	void WriteImage(string filename);

//General getters and setters

	void setWidth(int width);
	void setHeight(int height);

	int getWidth();
	int getHeight();

private:

	//store the pixels
	Colour *m_pixels;

	int m_height;
	int m_width;

	//convert 2D to 1D indexing
	int convert2dto1d(int x, int y); 

};


#endif
