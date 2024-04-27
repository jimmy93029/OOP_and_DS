//
// Instructor: Sai-Keung Wong
// Email:	cswingo@nycu.edu.tw
//
// National Yang Ming Chiao Tung University, Taiwan
// Computer Science
// Date: 2024/03/20
//
/*
A simple painter system.
*/
#ifndef __MY_SIMPLE_PAINTER_H_
#define __MY_SIMPLE_PAINTER_H_
#include "../headers.h"
#include <string>
#include "../primitives/vector3d.h"
#include <vector>
#include "base_class.h"
using std::vector;

class PAINT_RECORD {
    PAINT_RECORD( ) { }
    int nx, ny;
    int x,y;
    vector<vector3> color;
};

class SIMPLE_PAINTER : public BASE_SYSTEM
{
protected:
    int mStepSize;
    vector3 mColor;     // brush color
    float mTransparency;
    vector3 *mCanvas;
    int mNX; // number of columns
    int mNY; // number of rows
    //
    int mMinBrushSize;  // minimum
    int mMaxBrushSize;  // maximum
    int mBrushSize;     // current brush size
    //
    void drawAt(int x, int y, const vector3 &color, double w );

public:
    SIMPLE_PAINTER();
    void reset( );
    bool handleKeyPressedEvent( int key );
    void clickAt(double x, double y);
    void askForInput( );
    int getBrushSize( ) const;

    void getCanvasDimension( int &nx, int &ny ) const;
    void getColorAtPixel( int x, int y, vector3 &color ) const;
    void setColor( const vector3 &color );
    void setTransparency( float v );
    //
    vector3 getBrushCurrentColor() const;
    float getTransparency() const;
};

#endif