//
// Instructor: Sai-Keung Wong
// Email:	cswingo@nycu.edu.tw
//
// National Yang Ming Chiao Tung University, Taiwan
// Computer Science
// Date: 2024/03/20
//
#include "mySystemApp.h"

using namespace std;



void MY_SYSTEM_APP::showStudentInfo_2024( ) const
{
	cout << "*******************************" << endl;
    cout << "Date:2024/04/27" << endl;
	cout << "Student ID 111652040:\t" << endl;
	cout << "Student Name: jimmywu\t" << endl;
	cout << "Student Email: jimmywu0229.sc11@nycu.edu.tw\t" << endl;
	cout << "*******************************" << endl;
    cout << "Spring 2024." << endl;
    cout << "Instructor: Sai-Keung Wong" << endl;
    cout << "*******************************" << endl;
}

const char* MY_SYSTEM_APP::getStudentInfo() const {
    return STUDENT_INFO;
}

void MY_SYSTEM_APP::showMenu( ) const
{
    cout << "Welcome to our mini-galaxy!" << endl;
    cout << "Menu----------------------" << endl;
    cout << "Key usage:" << endl;
    cout << "F1: Monte Carlo Simulation" << endl;
    cout << "F2: Galaxy" << endl;
    cout << "F3: Painter" << endl;
    cout << "F4: Visualization of Bubble Sort" << endl;
    cout << "F5: Simple File Reader" << endl;
    cout << endl;
    cout << "s, S: student information" << endl;

    cout << "--------------------------" << endl;

}

MY_SYSTEM_APP::MY_SYSTEM_APP( )
{
    string title = STUDENT_INFO;
    glutSetWindowTitle(title.data( ));

	//mSystemType = SYSTEM_TYPE_GALAXY_SYSTEM;
    mSystemType = SYSTEM_TYPE_MONTE_CARLO_SIMULATION;  
    mColorPalette = new COLOR_PALETTE;
    mCamera->setZoom(-67);
}

void MY_SYSTEM_APP::initApp( )
{
	mFlgShow_Grid = false;
    showMenu( );
}

void MY_SYSTEM_APP::update( )
{
    switch( mSystemType ) {
		case SYSTEM_TYPE_MONTE_CARLO_SIMULATION:
			mMonteCarloSystem.update( );
            updateMonteCarloTitle();
			break;
		case SYSTEM_TYPE_GALAXY_SYSTEM:
			mGalaxySystem.update( );
            updateGalaxySystemTitle();
			break;
		case SYSTEM_TYPE_SIMPLE_PAINTER:
			mSimplePainter.update( );
            updateSimplePainterTitle();
			break;
		case SYSTEM_TYPE_SORTING_VISUALIZATION:
			mSortingVisualization.update( );
            updateSortingVisualizationTitle();
			break;
        case SYSTEM_TYPE_SIMPLE_FILE_EDITOR:
            mSimpleFileEditor.update( );
            updateSimpleFileEditorTitle();
            break;
	}
}




