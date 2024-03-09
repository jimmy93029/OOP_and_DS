//
// Instructor: Sai-Keung Wong
// Email:	cswingo@cs.nycu.edu.tw
//			wingo.wong@gmail.com
//
// National Yang Ming Chiao Tung University, Hsinchu, Taiwan
// Computer Science
// Date: 2024/03/03
//
#include <sstream>
#include "mySystemApp.h"

using namespace std;

void MY_SYSTEM_APP::increase()
{
    switch( mSystemType ) {
		case SYSTEM_TYPE_STUDENT_RECORD_MANAGER:
			break;
		case SYSTEM_TYPE_GRAPH:
			mGraphManager.increase_CurveNumber();
            updateGraphManagerTitle();
			break;

	}
}

void MY_SYSTEM_APP::decrease()
{
    switch( mSystemType ) {
		case SYSTEM_TYPE_MONTE_CARLO_CIRCLE:
			break;
		case SYSTEM_TYPE_SINECOSINE_FUNCTION:
			break;
		case SYSTEM_TYPE_CUBIC_FUNCTION:
			break;
		case SYSTEM_TYPE_STUDENT_RECORD_MANAGER:
			break;
		case SYSTEM_TYPE_GRAPH:
			mGraphManager.decrease_CurveNumber();
            updateGraphManagerTitle();
            break;
	}
}

void MY_SYSTEM_APP::increaseRadius()
{
    switch( mSystemType ) {
		case SYSTEM_TYPE_MONTE_CARLO_CIRCLE:
			break;
		case SYSTEM_TYPE_SINECOSINE_FUNCTION:
			break;
		case SYSTEM_TYPE_CUBIC_FUNCTION:
			break;
		case SYSTEM_TYPE_STUDENT_RECORD_MANAGER:
			break;
	}
}

void MY_SYSTEM_APP::decreaseRadius()
{
    switch( mSystemType ) {
		case SYSTEM_TYPE_MONTE_CARLO_CIRCLE:
			break;
		case SYSTEM_TYPE_SINECOSINE_FUNCTION:
			break;
		case SYSTEM_TYPE_CUBIC_FUNCTION:
			break;
		case SYSTEM_TYPE_STUDENT_RECORD_MANAGER:
			break;
	}
}

void MY_SYSTEM_APP::handle_key01()
{
        switch( mSystemType ) {
		case SYSTEM_TYPE_MONTE_CARLO_CIRCLE:
			break;
		case SYSTEM_TYPE_SINECOSINE_FUNCTION:
			break;
		case SYSTEM_TYPE_CUBIC_FUNCTION:
			break;
		case SYSTEM_TYPE_STUDENT_RECORD_MANAGER:
			break;
        case SYSTEM_TYPE_GRAPH:
            mGraphManager.setCurveType(CURVE_TYPE_EXPONENTIAL);
            break;
	}
}

void MY_SYSTEM_APP::handle_key02()
{
        switch( mSystemType ) {
		case SYSTEM_TYPE_MONTE_CARLO_CIRCLE:
			break;
		case SYSTEM_TYPE_SINECOSINE_FUNCTION:
			break;
		case SYSTEM_TYPE_CUBIC_FUNCTION:
			break;
		case SYSTEM_TYPE_STUDENT_RECORD_MANAGER:
			break;
        case SYSTEM_TYPE_GRAPH:
            mGraphManager.setCurveType(CURVE_TYPE_COSINE);
            break;
	}

}

void MY_SYSTEM_APP::handle_key03()
{
        switch( mSystemType ) {
		case SYSTEM_TYPE_MONTE_CARLO_CIRCLE:
			break;
		case SYSTEM_TYPE_SINECOSINE_FUNCTION:
			break;
		case SYSTEM_TYPE_CUBIC_FUNCTION:
			break;
		case SYSTEM_TYPE_STUDENT_RECORD_MANAGER:
			break;
        case SYSTEM_TYPE_GRAPH:
            mGraphManager.setCurveType(CURVE_TYPE_QUADRATIC);
            break;
	}
}

void MY_SYSTEM_APP::handle_key04()
{
        switch( mSystemType ) {
		case SYSTEM_TYPE_MONTE_CARLO_CIRCLE:
			break;
		case SYSTEM_TYPE_SINECOSINE_FUNCTION:
			break;
		case SYSTEM_TYPE_CUBIC_FUNCTION:
			break;
		case SYSTEM_TYPE_STUDENT_RECORD_MANAGER:
			break;
	}

}

void MY_SYSTEM_APP::askForInput( )
{
	switch( mSystemType ) {
		case SYSTEM_TYPE_MONTE_CARLO_CIRCLE:
			break;
		case SYSTEM_TYPE_SINECOSINE_FUNCTION:
			break;
		case SYSTEM_TYPE_CUBIC_FUNCTION:
			break;
		case SYSTEM_TYPE_STUDENT_RECORD_MANAGER:
			break;
	}
}


void MY_SYSTEM_APP::updateGraphManagerTitle()
{
    string title;
    string result;
    ostringstream strs;
    ostringstream strs_d;

    title = "SYSTEM_TYPE_GRAPH  ";

    {
        double numCurves = mGraphManager.getNumCurves();
        strs_d.str("");
        strs_d.clear();
        strs_d.precision(5);
        strs_d << numCurves;
        title = title + "  numCurves:";
        title = title + strs_d.str();
    }

    {
        double a = mGraphManager.getC_Average();
        strs_d.str("");
        strs_d.clear();
        strs_d.precision(5);
        strs_d << a;
        title = title + "  Average c:";
        title = title + strs_d.str();
    }

    {
        double a = mGraphManager.getC_Lowest();
        strs_d.str("");
        strs_d.clear();
        strs_d.precision(5);
        strs_d << a;
        title = title + "  Lowest c:";
        title = title + strs_d.str();
    }

    {
        double a = mGraphManager.getC_Largest();
        strs_d.str("");
        strs_d.clear();
        strs_d.precision(5);
        strs_d << a;
        title = title + "  Largest c:";
        title = title + strs_d.str();
    }

    {
        double b = mGraphManager.getParam_D();
        strs_d.str("");
        strs_d.clear();
        strs_d.precision(5);
        strs_d << b;
        title = title + "  d:";
        title = title + strs_d.str();
    }

    {
        int b = mGraphManager.getNumOfSamplePoints();
        strs_d.str("");
        strs_d.clear();
        strs_d.precision(5);
        strs_d << b;
        title = title + "  Number of Samples:";
        title = title + strs_d.str();
    }

    glutSetWindowTitle(title.data());
}

bool MY_SYSTEM_APP::specialFunc(int key, int x, int y)
{
    string title;
	//cout << "special key:" << key << endl;

	switch( key ) {
	case GLUT_KEY_F1:
        mSystemType = SYSTEM_TYPE_GRAPH;
        cout << "mSystemType = SYSTEM_TYPE_GRAPH" << endl;
        mGraphManager.setNumCurves(1);
        updateGraphManagerTitle();

		//mSystemType = SYSTEM_TYPE_MONTE_CARLO_CIRCLE;
		//updateMonteCarloTitle();
        //cout << "mSystemType = SYSTEM_TYPE_MONTE_CARLO_CIRCLE" << endl;
		break;
	case GLUT_KEY_F2:
        mSystemType = SYSTEM_TYPE_GRAPH;
        cout << "mSystemType = SYSTEM_TYPE_GRAPH" << endl;
        mGraphManager.setNumCurves();
        updateGraphManagerTitle();

		//mSystemType = SYSTEM_TYPE_SINECOSINE_FUNCTION;
		//updateSineConsineTitle();
        //cout << "mSystemType = SYSTEM_TYPE_SINECOSINE_FUNCTION" << endl;
		break;
	case GLUT_KEY_F3:
		//mSystemType = SYSTEM_TYPE_CUBIC_FUNCTION;
        //title = "SYSTEM_TYPE_CUBIC_FUNCTION";
        //glutSetWindowTitle(title.data( ));
		//cout << "mSystemType = SYSTEM_TYPE_CUBIC_FUNCTION" << endl;
		//updateCurveFunctionTitle();
		break;
    case GLUT_KEY_F4:
		//mSystemType = SYSTEM_TYPE_STUDENT_RECORD_MANAGER;
		//title = "SYSTEM_TYPE_STUDENT_RECORD_MANAGER";
        //glutSetWindowTitle(title.data( ));
        //cout << "mSystemType = SYSTEM_TYPE_STUDENT_RECORD_MANAGER" << endl;

		break;

		case GLUT_KEY_F5:
		//	mSystemType = SYSTEM_TYPE_GRAPH;
		//	cout << "mSystemType = SYSTEM_TYPE_GRAPH" << endl;
		//	updateGraphManagerTitle( );
			break;
	}
	return true;
}

bool MY_SYSTEM_APP::handleKeyEvent( unsigned char key )
{
	bool flgHandled = false;
	switch ( key ) {
		case 's':
		case 'S': //fall through
			showMyStudentInfo_2024( );
			break;
		case 'i':
		case 'I':
			//askForInput( );
			break;
		case 'g':
		case 'G':
			mFlgShow_Grid = !mFlgShow_Grid;
			break;
        case 'v':
            decreaseRadius();
            break;
        case 'b':
            increaseRadius();
            break;
        case '>':
            increase();
            break;
        case '<':
            decrease();
            break;
            case '1':
                handle_key01();
                break;
            case '2':
                handle_key02();
                break;
                case '3':
                handle_key03();
                break;
            case '4':
                handle_key04();
                break;
            case 'n':
            case 'N':
                mGraphManager.decrease_Param_D();
                updateGraphManagerTitle();
                break;
            case 'm':
            case 'M':
                mGraphManager.increase_Param_D();
                updateGraphManagerTitle();
                break;
                //
            case 'k':
            case 'K':
                mGraphManager.decrease_NumOfSamplePoints();
                updateGraphManagerTitle();
                break;
            case 'l':
            case 'L':
                mGraphManager.increase_NumOfSamplePoints();
                updateGraphManagerTitle();
                break;
	}
	return flgHandled;
}

bool MY_SYSTEM_APP::mouseMotionFunc( int mx, int my )
{
	static vector3 color( 1.0, 0.0, 0.0 );
	if ( mMouseButton == GLUT_LEFT_BUTTON ) {
		//cout << "mouseMotionFunc GLUT_LEFT_BUTTON" << endl;
		float x, y, z;
		mCamera->getCoordinates( x, y, z, (float) mx, (float) my );

		float rf = (rand()%1000)/1000.0f*0.5f;
		color += vector3(0.1f, 0.07f, 0.11f*rf);
		//clickAt( (int) x, (int) z, color );

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
	if ( button == GLUT_LEFT_BUTTON ) {
		//clickAt( (int) x, (int) z, color );
	} else {
		return FREE_CANVAS_2D::mouseFunc( button, state, mx, my );
	}
	return true;
}


bool MY_SYSTEM_APP::passiveMouseFunc( int mx, int my )
{
	float x, y, z;
	mCamera->getCoordinates( x, y, z, (float) mx, (float) my );
	return true;
}
// CODE: 2019/02/25. Do not delete this line.