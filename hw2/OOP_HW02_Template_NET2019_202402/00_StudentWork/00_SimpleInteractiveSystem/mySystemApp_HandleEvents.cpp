//
// Instructor: Sai-Keung Wong
// Email:	cswingo@nycu.edu.tw
//
// National Yang Ming Chiao Tung University, Taiwan
// Computer Science
// Date: 2024/03/20
//
#include <sstream>
#include "mySystemApp.h"

using namespace std;

void MY_SYSTEM_APP::askForInput( )
{
    switch( mSystemType ) {
    case SYSTEM_TYPE_MONTE_CARLO_SIMULATION:
        mMonteCarloSystem.askForInput( );
        break;
    case SYSTEM_TYPE_GALAXY_SYSTEM:
        mGalaxySystem.askForInput( );
        break;
    case SYSTEM_TYPE_SIMPLE_PAINTER:
        mSimplePainter.askForInput( );
        break;
    case SYSTEM_TYPE_SORTING_VISUALIZATION:
        mSortingVisualization.askForInput( );
        break;
    case SYSTEM_TYPE_SIMPLE_FILE_EDITOR:
        mSimpleFileEditor.askForInput( );
        break;
    }
}

void MY_SYSTEM_APP::updateMonteCarloTitle( )
{
    string title;
    string result;
    ostringstream strs;
    ostringstream strs_int;
    ostringstream strs_float;

    double value;
    int numSamples;
    title = "SYSTEM_TYPE_MONTE_CARLO_SIMULATION";
    
    //
    getStudentInfo();
    title = title + ":" + getStudentInfo() + "   ";
    //
    value = mMonteCarloSystem.getArea();
    strs << value;
    title = title + "  Area:";
    title = title + strs.str();
    numSamples = mMonteCarloSystem.getNumSamples();
    strs_int << numSamples;
    title = title + "  #Samples:";
    title = title + strs_int.str();
    //
    glutSetWindowTitle(title.data());
}

void MY_SYSTEM_APP::updateGalaxySystemTitle()
{
    string title;
    string result;
    ostringstream strs;
    ostringstream strs_int;
    ostringstream strs_float;

    double value;
    int numSamples;
    title = "SYSTEM_TYPE_GALAXY_SYSTEM";

    //
    getStudentInfo();
    title = title + ":" + getStudentInfo() + "   ";

    int numObjs = mGalaxySystem.getNumberOfActiveObjects();
    strs << numObjs;
    title = title + "  Number of All Active Objects:";
    title = title + strs.str();

    glutSetWindowTitle(title.data());
}
void MY_SYSTEM_APP::updateSimplePainterTitle()
{
    string title;
    string result;
    ostringstream strs;
    ostringstream strs_int;
    ostringstream strs_float;

    double value;
    int numSamples;
    title = "SYSTEM_TYPE_SIMPLE_PAINTER";

    //
    getStudentInfo();
    title = title + ":" + getStudentInfo() + "   ";

    int brushSize = mSimplePainter.getBrushSize();
    strs << brushSize;
    title = title + "  Brush Size:";
    title = title + strs.str();

    {
        ostringstream strs;
        vector3 color = mSimplePainter.getBrushCurrentColor();
        float transparency = mSimplePainter.getTransparency();
        strs.precision(2);
        strs << "  BrushColor:" << "(" 
            << color.r << "," << color.g << "," << color.b 
            << ")"
            << "    Transparency:" << transparency;
        title = title + strs.str();
    }

    glutSetWindowTitle(title.data());
}
void MY_SYSTEM_APP::updateSortingVisualizationTitle()
{
    string title;
    string result;
    ostringstream strs;
    ostringstream strs_int;
    ostringstream strs_float;

    double value;
    int numSamples;
    title = "SYSTEM_TYPE_SORTING_VISUALIZATION";

    //
    getStudentInfo();
    title = title + ":" + getStudentInfo() + "   ";

    int numElements = mSortingVisualization.getNumOfElements();
    strs << numElements;
    title = title + "  numElements:";
    title = title + strs.str();
    {
        int minOfInterval = mSortingVisualization.getMinOfInterval();
        int maxOfInterval = mSortingVisualization.getMaxOfInterval();
        ostringstream strs_min;
        ostringstream strs_max;
        strs_min << minOfInterval;

        strs_max << maxOfInterval;
        title = title + "  [";
        title = title + strs_min.str() + ", " + strs_max.str() + "]";

    }

    glutSetWindowTitle(title.data());
}
void MY_SYSTEM_APP::updateSimpleFileEditorTitle() {
    string title;
    string result;
    ostringstream strs;
    ostringstream strs_int;
    ostringstream strs_float;

    double value;
    int numSamples;
    title = "SYSTEM_TYPE_SIMPLE_FILE_EDITOR";

    //
    getStudentInfo();
    title = title + ":" + getStudentInfo() + "   ";

    glutSetWindowTitle(title.data());
}

void MY_SYSTEM_APP::clickAt(float x, float y )
{
    switch( mSystemType ) {
    case SYSTEM_TYPE_MONTE_CARLO_SIMULATION:
        //mMonteCarloSystem.reset( );
        break;
    case SYSTEM_TYPE_GALAXY_SYSTEM:
        //mGalaxySystem.reset( );
        break;
    case SYSTEM_TYPE_SIMPLE_PAINTER:
        if ( mColorPalette->clickAt(x, y, mMouseX, mMouseY ) ) {
            mSimplePainter.setColor( mColorPalette->getColor( ) );
            mSimplePainter.setTransparency( mColorPalette->getTransparency( ) );

        } else {
            mSimplePainter.clickAt( x, y );
        }
        break;
    case SYSTEM_TYPE_SORTING_VISUALIZATION:
        //mSortingVisualization.reset( );
        break;
    case SYSTEM_TYPE_SIMPLE_FILE_EDITOR:
        //mSimpleFileEditor.reset( );
        break;
    }
}

bool MY_SYSTEM_APP::specialFunc(int key, int x, int y)
{
    //cout << "special key:" << key << endl;

    string title;

    switch( key ) {
    case GLUT_KEY_F1:
        //title = "SYSTEM_TYPE_MONTE_CARLO_SIMULATION";
        //glutSetWindowTitle(title.data( ));
        mFlgShow_Grid = true;
        mSystemType = SYSTEM_TYPE_MONTE_CARLO_SIMULATION;
        mCamera->setZoom(89.0);
        mCamera->set( vec3(.45, .249, 0.0) );
        updateMonteCarloTitle( );
        cout << "mSystemType = SYSTEM_TYPE_MONTE_CARLO_SIMULATION" << endl;
        break;
    case GLUT_KEY_F2:
        //title = "SYSTEM_TYPE_GALAXY_SYSTEM";
        //glutSetWindowTitle(title.data( ));
        mFlgShow_Grid = false;
        mSystemType = SYSTEM_TYPE_GALAXY_SYSTEM;
        mCamera->setZoom(-280.0);
        mCamera->set( vec3(-0.0, -0.0, 0.0) );
        updateGalaxySystemTitle();
        cout << "mSystemType = SYSTEM_TYPE_GALAXY_SYSTEM" << endl;

        break;
    case GLUT_KEY_F3:
        //title = "SYSTEM_TYPE_SIMPLE_PAINTER";
        //glutSetWindowTitle(title.data( ));
        updateSimplePainterTitle();
        mFlgShow_Grid = false;
        mSystemType = SYSTEM_TYPE_SIMPLE_PAINTER;
        mCamera->setZoom(-94.0);
        mCamera->set( vec3(-188.0, -126.0, 0.0) );
        cout << "mSystemType = SYSTEM_TYPE_SIMPLE_PAINTER" << endl;
        break;
    case GLUT_KEY_F4:
        //title = "SYSTEM_TYPE_SIMPLE_FILE_EDITOR";
        //glutSetWindowTitle(title.data( ));
        updateSimpleFileEditorTitle();
        mFlgShow_Grid = false;
        mSystemType = SYSTEM_TYPE_SIMPLE_FILE_EDITOR;
        mCamera->setZoom(-37.0);
        mCamera->set( vec3(-50.05, -50.17, 0.0) );
        cout << "mSystemType = SYSTEM_TYPE_SIMPLE_FILE_EDITOR" << endl;

        break;
    case GLUT_KEY_F5:
        //title = "SYSTEM_TYPE_SORTING_VISUALIZATION";
        //glutSetWindowTitle(title.data( ));
        updateSortingVisualizationTitle();
        mFlgShow_Grid = true;
        mSystemType = SYSTEM_TYPE_SORTING_VISUALIZATION;
        mCamera->setZoom(59.0);
        mCamera->set( vec3(-13.42, -7.999, 0.0) );
        cout << "mSystemType = SYSTEM_TYPE_STUDENT_RECORD_MANAGER" << endl;

        break;
    }
    return true;
}

bool MY_SYSTEM_APP::handleKeyEvent( unsigned char key )
{
    bool flgHandled = false;
    if ( mSystemType == SYSTEM_TYPE_SIMPLE_PAINTER ) {
        if (mSimplePainter.handleKeyPressedEvent(key)) {
            updateSimplePainterTitle();
            return true;
        }
    }

    if ( mSystemType == SYSTEM_TYPE_GALAXY_SYSTEM ) {
        if (mGalaxySystem.handleKeyPressedEvent(key)) {
            updateGalaxySystemTitle();
            return true;
        }
    }

    if ( mSystemType == SYSTEM_TYPE_MONTE_CARLO_SIMULATION ) {
        if ( mMonteCarloSystem.handleKeyPressedEvent(key) ) {
            updateMonteCarloTitle( );
            return true;
        }
    }

    if ( mSystemType == SYSTEM_TYPE_SIMPLE_FILE_EDITOR ) {
        if (mSimpleFileEditor.handleKeyPressedEvent(key))
        {
            updateSimpleFileEditorTitle();
            return true;
        }
    }

    if (mSystemType == SYSTEM_TYPE_SORTING_VISUALIZATION) {
        if (mSortingVisualization.handleKeyPressedEvent(key))
        {
            updateSortingVisualizationTitle();
            return true;
        }
    }

    switch ( key ) {
    case ' ':
        if ( mSystemType == SYSTEM_TYPE_SORTING_VISUALIZATION ) {
            mSortingVisualization.processOneSwapStep( );
        }
        /*

        if ( mSystemType == SYSTEM_TYPE_MONTE_CARLO_SIMULATION ) {
            mMonteCarloSystem.computeArea( );
        }
        */

        break;
    case 'm':
    case 'M':
        showMenu( );
        break;
    case 's':
    case 'S':
        showStudentInfo_2024( );
        break;
    case 'i':
    case 'I':
        askForInput( );
        break;
    case 'g':
    case 'G':
        mFlgShow_Grid = !mFlgShow_Grid;
        break;
    case 'z':
    case 'Z':
        mCamera->setZoom(-158.0);
        break;
    case 'c':
    case 'C':
        mCamera->setZoom(1.0);
        break;
    }
    return flgHandled;
}

bool MY_SYSTEM_APP::mouseMotionFunc( int mx, int my )
{
    static vector3 color( 1.0, 0.0, 0.0 );


    float x, y, z;
    mCamera->getCoordinates( x, y, z, (float) mx, (float) my );

    mMouseX_2DSpace = x;
    mMouseY_2DSpace = z;

    if ( mMouseButton == GLUT_LEFT_BUTTON ) {
        //cout << "mouseMotionFunc GLUT_LEFT_BUTTON" << endl;


        float rf = (rand()%1000)/1000.0f*0.5f;
        color += vector3(0.1f, 0.07f, 0.11f*rf);
        //clickAt( (int) x, (int) z, color );
        clickAt( x, z );

    } else {
        //if ( mPickedObj == 0 ) {
        return FREE_CANVAS_2D::mouseMotionFunc( mx, my );
        //}
    }

    //cout << "here 5" << endl;

    //float x, y, z;
    //mCamera->getCoordinates( x, y, z, mx, my );
    //cout << "Set Position XZ:" << x << "\t" << z << endl;
    return true;
}

bool MY_SYSTEM_APP::mouseFunc( int button, int state, int mx, int my )
{
    //cout << "button:" << button << endl;
    //if ( mPickedObj == 0 || button == GLUT_MIDDLE_BUTTON ) {
    //if ( button == 3 || button == 4 ) {

    float x, y, z;
    mCamera->getCoordinates( x, y, z, (float) mx, (float) my );

    vector3 color( 1.0, 0.0, 0.0 );
    mMouseButton = button;
    mMouseX = mx;
    mMouseY = mCanvasHeight - my;

    //cout << "mMouseX:" << mMouseX << endl;

    //cout << "mMouseY:" << mMouseY << endl;

    if ( button == GLUT_LEFT_BUTTON ) {
        clickAt( x, z );
    } else {
        return FREE_CANVAS_2D::mouseFunc( button, state, mx, my );
    }
    return true;
}


bool MY_SYSTEM_APP::passiveMouseFunc( int mx, int my )
{
    float x, y, z;
    mCamera->getCoordinates( x, y, z, (float) mx, (float) my );

    mMouseX_2DSpace = x;
    mMouseY_2DSpace = z;

    return true;
}