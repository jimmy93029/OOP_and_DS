//
// Instructor: Sai-Keung WONG
// Email:	cswingo@nycu.edu.tw
//			wingo.wong@gmail.com
//
// National Yang Ming Chiao Tung University, Taiwan
// Computer Science
// Date: 2022/03/17
//
#include "mySystem_SortingVisualization.h"
#include <iostream>

using namespace std;

namespace
{
#define DFT_MAX_VALUE   10
#define MIN_VALUE   1
#define MAX_VALUE   20
}

SORTING_VISUALIZATION::SORTING_VISUALIZATION( )
{
    mNumElements = 20;
    mMinValue = MIN_VALUE;
    mMaxValue = DFT_MAX_VALUE;
    reset( );
}

int SORTING_VISUALIZATION::getMinOfInterval() const
{
    return mMinValue;
}

int SORTING_VISUALIZATION::getMaxOfInterval() const
{
    return mMaxValue;

}

/*
Generate randomly the elements between mMinValue and mMaxValue.
The total number of elements is mNumElements.
*/
void SORTING_VISUALIZATION::reset( )
{
    mNumbers.resize( mNumElements );
    for ( int i = 0; i < mNumElements; ++i ) {
        int dv = rand( )%(mMaxValue - mMinValue + 1);
        int v = mMinValue + dv;
        mNumbers[ i ] = v;
    }
    mCurrentIndex = 0;
    mMaxNonHandledElementIndex = mNumElements;
}

bool SORTING_VISUALIZATION::handleKeyPressedEvent(int key)
{
    bool flgHandled = false;
    switch (key) {
    case 'r':
    case 'R':
        reset();
        flgHandled = true;
        break;
    case ' ':
        processOneSwapStep();
        flgHandled = true;
        break;
    case 'i':
    case 'I':
        askForInput();
        flgHandled = true;
        break;
    case '<':
        mMaxValue -= 1;
        if (mMaxValue < MIN_VALUE) mMaxValue = MIN_VALUE;
        reset();
        break;
    case '>':
        mMaxValue += 1;
        if (mMaxValue > MAX_VALUE) mMaxValue = MAX_VALUE;
        reset();
        break;
    }
    return flgHandled;
}

/*
Show the system title.
Show system information
Show key usage
*/
void SORTING_VISUALIZATION::askForInput( )
{
    cout << "//////////////////////////////" << endl;
    cout << "SORTING_VISUALIZATION" << endl;
    cout << "//////////////////////////////" << endl;
    cout << "Press SPACEBAR to visualize the process..." << endl;
    cout << "//////////////////////////////" << endl;
}

/*
Two elements whose
indices are (index) and (index+1).

Swap these two elements if their order is incorrect.
*/
void SORTING_VISUALIZATION::checkAndSwapTwoElements( int index )
{
    int v0 = mNumbers[ index ];
    int v1 = mNumbers[ index + 1];

    //
    // fill/modify your own stuff
    //
}

/*
Perform one swap step for two elements.
The first element is indicated by mCurrentIndex.
*/
void SORTING_VISUALIZATION::processOneSwapStep( )
{
    if ( mNumElements <= 1) return;
    if ( mMaxNonHandledElementIndex <= 1 ) return;
    //
    // fill/modify your own stuff
    //

    //checkAndSwapTwoElements( mCurrentIndex );

    //
    // fill/modify your own stuff
    //


}

/*
Return the current element that will be processed.
*/
int SORTING_VISUALIZATION::getCurrentElementIndex( ) const
{
    return mCurrentIndex;
}

/*
Return the index of the last non-handled element.
*/
int SORTING_VISUALIZATION::getMaxNonHandledElementIndex( ) const
{
    //
    // fill/modify your own stuff
    //
    //return mMaxNonHandledElementIndex;
    return 0;
}

/*
Return the number of elements.
*/
int SORTING_VISUALIZATION::getNumOfElements( ) const
{
    //
    // fill/modify your own stuff
    //
    //Return the number of elements
    //return mNumElements;
    return mNumElements/2;
}

/*
Return the element whose index is elementIndex.
*/
int SORTING_VISUALIZATION::getElement( int elementIndex ) const
{
    //
    // fill/modify your own stuff
    //
    //return mNumbers[ elementIndex ];
    return 1;
}