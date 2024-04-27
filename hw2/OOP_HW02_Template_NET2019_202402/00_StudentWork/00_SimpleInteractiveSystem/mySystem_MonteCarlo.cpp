//
// Instructor: Sai-Keung WONG
// Email:	cswingo@nycu.edu.tw
//			wingo.wong@gmail.com
//
// National Yang Ming Chiao Tung University, Taiwan
// Computer Science
// Date: 2022/03/17
//
#include "mySystem_MonteCarlo.h"
#include "myBasicTools.h"
#include <iostream>

using namespace std;

#define MIN_SAMPLES 100
#define MAX_SAMPLES 1000000
#define DFT_SAMPLES 100000
#define DELTA_SAMPLES 10000

MONTE_CARLO_SYSTEM::MONTE_CARLO_SYSTEM( )
{
    //BEGIN DO NOT CHANGE ANYTHING HERE
    mFunctionType = FUNCTION_EXPONENTIAL_SINE;
    cout.precision(8);
    //
    mMinX = -6;
    mMaxX = 6;
    mMinY = -5;
    mMaxY = 5;
    mMinNumSamples = MIN_SAMPLES;
    mMaxNumSamples = MAX_SAMPLES;
    mNumSamples = DFT_SAMPLES;
    computeSamples( );
    computeArea( );
    //END DO NOT CHANGE ANYTHING HERE
}

/*
Show the message to the user.
*/
void MONTE_CARLO_SYSTEM::askForInput( )
{
    cout << "MONTE_CARLO_SYSTEM::askForInput( )" << endl;
    cout << "Use keys to control......" << endl;
    cout << "1: FUNCTION_EXPONENTIAL_COSINE" << endl;
    cout << "2: FUNCTION_SINE_COSINE" << endl;
    cout << "3: FUNCTION_COMPLEX" << endl;
    cout << "4: FUNCTION_HYPERBOLIC_TANGENT" << endl;
    
    cout << "<: decrease the number of sample points" << endl;
    cout << ">: increase the number of sample points" << endl;
    cout << "n: set to minimum number of sample points" << endl;
    cout << "m: set to maximum number of sample points" << endl;
    cout << "r: reset" << endl;
}

/*
Handle the key events.
*/
bool MONTE_CARLO_SYSTEM::handleKeyPressedEvent( int key )
{
    cout << "MONTE_CARLO_SYSTEM::handleKeyPressedEvent" << endl;
    bool flgHandled = false;
    switch( key ) {
    case 'r':
    case 'R':
        reset();
        flgHandled = true;
        break;
    case 'n':
    case 'N':
        mNumSamples = mMinNumSamples;
        flgHandled = true;
        break;
    case 'm':
    case 'M':
        mNumSamples = mMaxNumSamples;
        
        flgHandled = true;
        break;
    case '<':
        mNumSamples -= DELTA_SAMPLES;
        flgHandled = true;

        break;
    case '>':
        mNumSamples += DELTA_SAMPLES;
        flgHandled = true;

        break;
    case '1':
        mFunctionType = FUNCTION_EXPONENTIAL_SINE;
        flgHandled = true;
        break;
    case '2':
        mFunctionType = FUNCTION_COSINE;
        flgHandled = true;
        break;
    case '3':
        mFunctionType = FUNCTION_SINE_COSINE;
        flgHandled = true;
        break;
    case '4':
        mFunctionType = FUNCTION_EXPONENTIAL;
        flgHandled = true;
        break;
    }
    if ( flgHandled ) {
        cout << "mNumSamples:" << mNumSamples << endl;

        if ( mNumSamples < mMinNumSamples ) mNumSamples = mMinNumSamples;
        if ( mNumSamples > mMaxNumSamples ) mNumSamples = mMaxNumSamples;

        computeSamples( );
        computeArea( );
    }
    return flgHandled;
}

/*
Reset the system state.
Steps:
compute the samples
compute the area
*/
void MONTE_CARLO_SYSTEM::reset( )
{
    srand(0);
    mNumSamples = DFT_SAMPLES;
    computeSamples( );
    computeArea( );
}

/*
Return the interval [mMinX, mMaxX] of x to [minX, maxX].
*/
void MONTE_CARLO_SYSTEM::getXInterval( double &minX, double &maxX ) const
{

    // fill/modify your own stuff
    minX = mMinX;
    maxX = mMaxX;
    // minX = 1;
    // maxX = 10;
}

/*
Uniformly generate all the random sample points
inside [mMinX, mMaxX] x [mMinY, mMaxY].

The number of sample points is mNumSamples.
*/
void MONTE_CARLO_SYSTEM::computeSamples( )
{
    //
    if (mNumSamples < 0) mNumSamples = 0;
    mX.clear( );
    mY.clear( );
    mX.resize( mNumSamples );
    mY.resize( mNumSamples );
    // fill/modify your own stuff
    for ( int i = 0; i < mNumSamples; ++i ) {
        // fill/modify your own stuff
        double x = getRandDouble( mMinX, mMaxX );
        mX[ i ] = x;
        mY[ i ] = 1.5;           // How to get the function output
    }
}

/*
Compute the area based on the Monte Carlo Simulation.
Steps:
1. compute the number of samples inside the area
2. compute the ratio which is equal to

number of samples inside the area 
divided 
by total number of samples

3. Based on the ratio, compute the area

Return the area value
*/
void MONTE_CARLO_SYSTEM::computeArea( ) 
{
    mArea = 0.0;
    if (mNumSamples <= 0) return;

    // fill/modify your own stuff
    int counter = 0;
    for ( int i = 0; i < mNumSamples; ++i ) {
        if (isInsideArea(mX[i], mY[i]))
            counter++;
    }

    double ratio = (counter/(double)mNumSamples);
    double area = ratio*(mMaxX - mMinX)*(mMaxY - mMinY);

    cout << "mNumSamples:" << mNumSamples << endl;
    cout << "Area:" << area << endl;
    mArea = area;

}

/*
Return the estimated area
*/
double MONTE_CARLO_SYSTEM::getArea() const {
    return mArea;
}
/*
Return the number of sample points.
*/
int MONTE_CARLO_SYSTEM::getNumSamples( ) const
{
    //Return the number of sample points
    return mNumSamples;
    // fill/modify your own stuff
    // return mNumSamples / 10;
}

/*
* Return the coordinates (x,y) of a sample whose index is sampleIndex
*/ 
void MONTE_CARLO_SYSTEM::getSample(int sampleIndex, double &x, double &y ) const
{
    // fill/modify your own stuff
    x = mX[ sampleIndex ];
    y = mY[ sampleIndex ];
    // y = 10;

}

/*
Based on the function type mFunctionType,
compute the function value for a given value x.
Return the function value.
*/
double MONTE_CARLO_SYSTEM::getValue(double x) const
{
    double y;
    // fill/modify your own stuff
    switch( mFunctionType ) {
    case FUNCTION_EXPONENTIAL_SINE:
        y = 5.0*pow(2.718281828, -x/4.0)*sin(x);
        break;
    case FUNCTION_COSINE:
        y = cos(4*x);
        break;
    case FUNCTION_SINE_COSINE:
        y = -3*cos(x*x);
        break;
    case FUNCTION_EXPONENTIAL:
        y = (pow(2.718281828, x) - pow(2.718281828, -x));
        break;
    }
    return y;
}

/*
Given a point (x,y).
Check if the point lies inside the area.
Return true if it lies inside the area.
Return false otherwise.

The area is defined as the region bounded by
the x-axis and the curve of the function f(x).

Assume that (x,y) lies inside area, then
one of the following conditions must be satisfied.

- If f(x) >= 0 and y >=0, then y <= f(x).
- If f(x) <= 0 and y <=0, then y >= f(x).

Otherwise, (x,y) does not lay inside the area.
*/
bool MONTE_CARLO_SYSTEM::isInsideArea( double x, double y ) const
{
    double realY = getValue( x );
    //
    // fill/modify your own stuff
    //
    if (realY >= 0.0 && y >= 0.0 && y <= realY ) return true;
    if (realY <= 0.0 && y <= 0.0 && y >= realY) return true;

    return false;
}