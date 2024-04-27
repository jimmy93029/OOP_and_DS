//********************************************
// Student Name			:
// Student ID			:
// Student Email Address:
//********************************************
//
// Instructor: Sai-Keung Wong
// Email:	cswingo@nycu.edu.tw
//
// National Yang Ming Chiao Tung University, Taiwan
// Computer Science
// Date: 2024/03/20
//
#ifndef __MY_SYSTEM_APP_APP_H_
#define __MY_SYSTEM_APP_APP_H_
/////////////////////////////////////////////
#include "../headers.h"
#include <string>
#include "../openGLbasic/opengl_stuff.h"
#include "../freeCanvas2D/freeCanvas2D.h"
#include "../textureManager/texture.h"
/////////////////////////////////////////////
#include "mySystem_MonteCarlo.h"
#include "mySystem_GalaxySystem.h"
#include "mySystem_SimplePainter.h"
#include "mySystem_SimpleFileEditor.h"
#include "mySystem_SortingVisualization.h"
#include "mySystem_ColorPalette.h"

#define SYSTEM_TYPE_MONTE_CARLO_SIMULATION		1
#define SYSTEM_TYPE_GALAXY_SYSTEM           2
#define SYSTEM_TYPE_SIMPLE_PAINTER      3
#define SYSTEM_TYPE_SORTING_VISUALIZATION	4
#define SYSTEM_TYPE_SIMPLE_FILE_EDITOR	    5
    
/////////////////////////////////////////////

#define STUDENT_INFO "Name:§dªl¿Î   ID:111652040"

using namespace ns_opengl;
class MY_SYSTEM_APP : public FREE_CANVAS_2D
{
protected:
	bool mFlgShow_Grid;
	int mMouseButton;
    int mMouseX, mMouseY;
    float mMouseX_2DSpace, mMouseY_2DSpace;
    //
    COLOR_PALETTE *mColorPalette;
	//
	void drawOrigin( ) const;
    void draw_System_Content( ) const;
    //

	void draw_MonteCarlo_Simulation( ) const;
	void draw_GalaxySystem( ) const;
	void draw_SimplePainter( ) const;
	void draw_SortingVisualization( ) const;
    void draw_SimpleFileEditor( ) const;
    void askForInput( );
    void reset( );
	//
	void showStudentInfo_2024( ) const;
	const char* getStudentInfo( ) const;
    void showMenu( ) const;
    //
    void updateMonteCarloTitle( );
	void updateGalaxySystemTitle();
	void updateSimplePainterTitle();
	void updateSortingVisualizationTitle();
	void updateSimpleFileEditorTitle();

public:
	MY_SYSTEM_APP( );
	virtual void update( );
	virtual void draw( ) const;
	virtual void initApp( );
	//
    void clickAt(float x, float y );
	bool specialFunc(int key, int x, int y);
	bool handleKeyEvent( unsigned char key );
	bool passiveMouseFunc( int mx, int my );
	bool mouseMotionFunc( int x, int y );
	bool mouseFunc( int button, int state, int x, int y );
	//////////////////////////////////////////////////
    
protected:
	int mSystemType;
    GALAXY_SYSTEM mGalaxySystem;
	MONTE_CARLO_SYSTEM mMonteCarloSystem;
	SIMPLE_FILE_EDITOR mSimpleFileEditor;
	SIMPLE_PAINTER mSimplePainter;
	SORTING_VISUALIZATION mSortingVisualization;
};

#endif