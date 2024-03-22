//
// Instructor: Sai-Keung Wong
// Email:	cswingo@cs.nycu.edu.tw
//			wingo.wong@gmail.com
//
// National Yang Ming Chiao Tung University, Hsinchu, Taiwan
// Computer Science
// Date: 2024/03/03
//
#include "mySystem_Curve.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

//
// no-arg constructor
//
CURVE_FUNCTION::CURVE_FUNCTION()
{
	//xMin = -100.0;
	//xMax= 100.0;

	xMin = -100.0;
	xMax = 100.0;

	mNumSamples = 100;

	a = b = c = d = 1;
	a = 1.0;
	b = -0.5;

	mCurveType = CURVE_TYPE_EXPONENTIAL;
}

//
// set the curve type
//
void CURVE_FUNCTION::setCurveType(CURVE_TYPE type)
{
	mCurveType = type;
}

// increase parameter b
void CURVE_FUNCTION::increase(double dv)
{
	b += dv;
}

// decrease parameter b
void CURVE_FUNCTION::decrease(double dv)
{
	b -= dv;

}

// set the value to parameter c
void CURVE_FUNCTION::setC(double c)
{
	this->c = c;
}

// set the value to parameter d
void CURVE_FUNCTION::setD(double d)
{
	this->d = d;
}

// ask for input
void CURVE_FUNCTION::askForInput()
{
	//empty
}

//
// set the interval of x
// xMin <- min_X
// xMax <- max_X
//
void CURVE_FUNCTION::set_IntervalOfX(double min_X, int max_X)
{
	xMin = min_X;
	xMax = max_X;
}

//
// return the interval of x
// xMin <- this->xMin
// xMax <- this->xMax
void CURVE_FUNCTION::getIntervalOfX(double& xMin, double& xMax) const
{
	xMin = this->xMin;;
	xMax = this->xMax;

}

//
// return the number of sample points
//
int CURVE_FUNCTION::getNumOfSamples() const
{
	return mNumSamples;
}

//
// set the number of sample points
// mNumSamples <- num
//
void CURVE_FUNCTION::setNumOfSamplePoints(int num)
{

	mNumSamples = num;
}

//
// get the value of a function for x
//
double CURVE_FUNCTION::getValue(double x) const
{
	double y = 0;
	switch (mCurveType) {
	case CURVE_TYPE_EXPONENTIAL:
		y = (x / 2.5) / 10 - 1 + (c + sin(4 * d * x)) * exp(-(x / 2.5));
		break;
	case CURVE_TYPE_COSINE:
		y = c * x - d * x * cos(x);
		break;
	case CURVE_TYPE_QUADRATIC:
		y = x*x + c * x + d;
		break;
	}

	return y;
}

//
// randomly compute a value for c 
// inside the interval [u0, u1] in a uniform manner
//
void CURVE_FUNCTION::setRandom_C(double u0, double u1)
{
	c = (u1 - u0) * rand() / (RAND_MAX + 1.0) + u0;

}

//
// randomly compute a value for d 
// inside the interval [u0, u1] in a uniform manner
//
void CURVE_FUNCTION::setRandom_D(double u0, double u1)
{
	d = (u1 - u0) * rand() / (RAND_MAX + 1.0) + u0;

}

//
// get the boundary point (x, y)
// point_index = 0: return the leftmost point (x,y)
// point_index = 1: return the rightmost point (x,y)
// 
void CURVE_FUNCTION::getBoundaryPoint(int point_index, double& x, double& y) const
{
	x = xMin;
	y = -10;
	switch (point_index) {
	case 0:
		x = xMin;
		y = getValue(xMin);
		break;
	case 1:
		x = xMax;
		y = getValue(xMax);
		break;
	};
}

//
// get the x- and y-coordinates of all the extreme points 
//
// store the x- and y-coordinates to X and Y vectors, respectively
//
void CURVE_FUNCTION::getExtremePoints(vector<double>& X, vector<double>& Y) const
{
	X.clear();
	Y.clear();

	int numSamples = mNumSamples;
	vector<double> y; 
	vector<double> x;
	y.resize(3), x.resize(3);
	for (int i = 1; i < numSamples - 1; ++i) {
		for (int j = -1; j <= 1; j++) {
			x[j+1] = xMin + ((i + j) / ((double)numSamples - 1)) * (xMax - xMin);
			y[j + 1] = getValue(x[j + 1]);
		}
		if (y[1] > y[0] && y[1] > y[2])
			X.push_back(x[1]), Y.push_back(y[1]);
		if (y[1] < y[0] && y[1] < y[2])
			X.push_back(x[1]), Y.push_back(y[1]);
	}

}