//
// Instructor: Sai-Keung Wong
// Email:	cswingo@cs.nycu.edu.tw
//			wingo.wong@gmail.com
//
// National Yang Ming Chiao Tung University, Hsinchu, Taiwan
// Computer Science
// Date: 2024/03/03
//

#include "mySystem_GraphManager.h"


#define MAX_NUM_SAMPLE_POINTS 400
#define MIN_NUM_SAMPLE_POINTS 2

//
// no-arg constructor
//
GRAPH_MANAGER::GRAPH_MANAGER()
{
	mMin_X = -2;
	mMax_X = 5;


	mD_LowerBound = -10;
	mD_UpperBound = 10;

	mC_LowerBound = -2;
	mC_UpperBound = 2;

	mParam_D = -1;

	mDelta_D = 1;

	mNumOfSamplePoints = MAX_NUM_SAMPLE_POINTS / 2;

	mNumCurves = MAX_NUM_CURVES;
	mCurves.resize( MAX_NUM_CURVES );

	setCurves_IntervalOfX(mMin_X, mMax_X);
	setCurves_D(mParam_D);
	setCurves_NumOfSamplePoints(mNumOfSamplePoints);
	setCurves_Random_C( );

	mNumCurves = MIN_NUM_CURVES;

}

//
// set the curve type to all the curves
//
void GRAPH_MANAGER::setCurveType(CURVE_TYPE type)
{
	for (int i = 0; i < MAX_NUM_CURVES; ++i) {
		//
		// implement your own stuff
		//
	}
}

//
// set the interval of x, [min_X, max_X], to all the curves
//
void  GRAPH_MANAGER::setCurves_IntervalOfX(double min_X, int max_X)
{
	for (int i = 0; i < MAX_NUM_CURVES; ++i) {
		//
		// implement your own stuff
		//
	}
}

//
// set the value of parameter d of all the curves
//
void GRAPH_MANAGER::setCurves_D(double d)
{
	for (int i = 0; i < MAX_NUM_CURVES; ++i) {
		//
		// implement your own stuff
		//
	}
}

//
// increase parameter mParam_D by mDelta_D, and set it to all the curves
// 
// mParam_D must be inside [mD_LowerBound, mD_UpperBound]
//
void GRAPH_MANAGER::increase_Param_D()
{
	//
	// implement your own stuff
	//
	mParam_D = mParam_D + mDelta_D;
	
	setCurves_D(mParam_D);
}

//
// decrease parameter mParam_D by -mDelta_D, and set it to all the curves
// 
// mParam_D must be inside [mD_LowerBound, mD_UpperBound]
//
void GRAPH_MANAGER::decrease_Param_D()
{
	//
	// implement your own stuff
	//
}

//
// get the value of mParam_D 
//
double GRAPH_MANAGER::getParam_D() const
{
	//
	// implement your own stuff
	//

	return 10;
}

//
// increase the number of curves by 10 at a time
//
// mNumCurves must be inside [MAX_NUM_CURVES, MAX_NUM_CURVES]
//
void GRAPH_MANAGER::increase_CurveNumber()
{
	//
	// implement your own stuff
	//
	mNumCurves += 1;
	setCurves_Random_C();

}

//
// decrease the number of curves by 10 at a time
//
// mNumCurves must be inside [MAX_NUM_CURVES, MAX_NUM_CURVES]
//
void GRAPH_MANAGER::decrease_CurveNumber()
{
	//
	// implement your own stuff
	//
}

//
// Set the number of curves
// mNumCurves <- num
//
void GRAPH_MANAGER::setNumCurves(int num) {
	//
	// implement your own stuff
	//
}

//
// get the number of curves
// return mNumCurves
//
int GRAPH_MANAGER::getNumCurves() const
{
	//
	// implement your own stuff
	//
	return 1;
}

//
// get curve with index curve_index
//
// return mCurves[curve_index]
//
const CURVE_FUNCTION& GRAPH_MANAGER::getCurve(int curve_index) const
{
	//
	// implement your own stuff
	//
	return mCurves[0];
}

//
// increase the number of sample points by delta_num at a time
//
// mNumOfSamplePoints must be inside [MIN_NUM_SAMPLE_POINTS, MAX_NUM_SAMPLE_POINTS]
//
void GRAPH_MANAGER::increase_NumOfSamplePoints(int delta_num)
{
	//
	// implement your own stuff
	//
	mNumOfSamplePoints += delta_num;
	setCurves_NumOfSamplePoints(mNumOfSamplePoints);
}

//
// decrease the number of sample points by delta_num at a time
//
// mNumOfSamplePoints must be inside [MIN_NUM_SAMPLE_POINTS, MAX_NUM_SAMPLE_POINTS]
//
void GRAPH_MANAGER::decrease_NumOfSamplePoints(int delta_num)
{
	//
	// implement your own stuff
	//
}

//
// set the number of sample points to all the curves
//
// For each curve i
// 
// mCurves[i].setNumOfSamplePoints(num);
//
void GRAPH_MANAGER::setCurves_NumOfSamplePoints(int num)
{
	for (int i = 0; i < MAX_NUM_CURVES; ++i) {
		//
		// implement your own stuff
		//
		mCurves[i].setNumOfSamplePoints(100);
	}
}

//
// For each curve mCurves[i]
// 
// randomly generate parameter C inside 
// [mC_LowerBound, mC_UpperBound] in a uniform manner
// 
// set the parameter value of C to mCurves[i]
//
// mCurves[i].setRandom_C(mC_LowerBound, mC_UpperBound);
//
void GRAPH_MANAGER::setCurves_Random_C()
{
	for (int i = 0; i < MAX_NUM_CURVES; ++i) {
		//
		// implement your own stuff
		//
	}
}

//
// get the average value of parameter c of all the curves
//
double GRAPH_MANAGER::getC_Average() const
{
	//
	// implement your own stuff
	//
	double v = 0;
	for (int i = 0; i < mNumCurves; ++i) {
	
	}
	v = v / mNumCurves;
	return v;
}

//
// get the lowest value of parameter c of all the curves
//
double GRAPH_MANAGER::getC_Lowest() const
{
	double v = mCurves[0].getC();
	//
	// implement your own stuff
	//
	return v;
}

//
// get the largest value of parameter c of all the curves
//
double GRAPH_MANAGER::getC_Largest() const
{
	double v = mCurves[0].getC();
	//
	// implement your own stuff
	//
	return v;
}
