//
// Instructor: Sai-Keung Wong
// Email:	cswingo@cs.nycu.edu.tw
//			wingo.wong@gmail.com
//
// National Yang Ming Chiao Tung University, Hsinchu, Taiwan
// Computer Science
// Date: 2024/03/03
//
#ifndef __MY_CURVE_FUNCTION_H_
#define __MY_CURVE_FUNCTION_H_
#include "../headers.h"
#include <string>
#include "../primitives/vector3d.h"
#include <vector>
using std::vector;

enum CURVE_TYPE {
	CURVE_TYPE_EXPONENTIAL = 0
	, CURVE_TYPE_COSINE
	, CURVE_TYPE_QUADRATIC
};
class CURVE_FUNCTION
{
protected:
	double a, b, c, d;		// parameters
	double xMin, xMax;		// interval of x, [minimum x, maximum x]
	int mNumSamples;		// number of sample points

	CURVE_TYPE mCurveType;	// the curve type

public:
	CURVE_FUNCTION( );
	void askForInput( );
	void getIntervalOfX(double& xMin, double& xMax) const;
	int getNumOfSamples( ) const;

	void setNumOfSamplePoints(int num);
	void set_IntervalOfX(double min_X, int max_X);
	//
	double getValue( double x ) const;
    void increase(double dv = 0.1);
    void decrease(double dv = 0.1);

	void setCurveType(CURVE_TYPE type);
	void setC(double c);
	void setD(double d);

    double getA() const { return a; }
    double getB() const { return b; }
    double getC() const { return c; }
    double getD() const { return d; }
	//
	void setRandom_C(double u0, double u1);
	void setRandom_D(double u0, double u1);
	void getBoundaryPoint(int point_index, double& x, double& y) const;
	void getExtremePoints(vector<double> &X, vector<double>& Y) const;
	//
};
#endif
// CODE: 2023/02/26. Do not delete this line.