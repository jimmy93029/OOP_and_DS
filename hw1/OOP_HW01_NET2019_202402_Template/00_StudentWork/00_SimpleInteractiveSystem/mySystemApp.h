//********************************************
// Student Name			:
// Student ID			:
// Student Email Address:
//********************************************
//
// Instructor: Sai-Keung Wong
// Email:	cswingo@cs.nycu.edu.tw
//			wingo.wong@gmail.com
//
// National Yang Ming Chiao Tung University, Hsinchu, Taiwan
// Computer Science
// Date: 2024/03/03
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
#include "mySystem_Curve.h"
#include "mySystem_GraphManager.h"

#define SYSTEM_TYPE_MONTE_CARLO_CIRCLE		1
#define SYSTEM_TYPE_SINECOSINE_FUNCTION		2
#define SYSTEM_TYPE_CUBIC_FUNCTION			3
#define SYSTEM_TYPE_STUDENT_RECORD_MANAGER	4
#define SYSTEM_TYPE_GRAPH					5
/////////////////////////////////////////////

using namespace ns_opengl;
class MY_SYSTEM_APP : public FREE_CANVAS_2D
{
protected:
	bool mFlgShow_Grid;
	int mMouseButton;
	//
	void updateGraphManagerTitle();
    //
	void drawOrigin( ) const;
	void draw_System_Content( ) const;
	void draw_CubicFunction(const CURVE_FUNCTION& curve) const;
	void draw_BoundaryPoints(const CURVE_FUNCTION& curve) const;
	void draw_Graphs() const;
	void draw_Points(const vector<double>& X, const vector<double>& Y) const;
	void askForInput( );
    void increase();
    void decrease();
    void increaseRadius();
    void decreaseRadius();
    void handle_key01();
    void handle_key02();
    void handle_key03();
    void handle_key04();
	//
	void showMyStudentInfo_2024( ) const;
public:
	MY_SYSTEM_APP( );
	virtual void update( );
	virtual void draw( ) const;
	virtual void initApp( );
	//
	bool specialFunc(int key, int x, int y);
	bool handleKeyEvent( unsigned char key );
	bool passiveMouseFunc( int mx, int my );
	bool mouseMotionFunc( int x, int y );
	bool mouseFunc( int button, int state, int x, int y );
	//////////////////////////////////////////////////
protected:
	int mSystemType;
	GRAPH_MANAGER mGraphManager;
};

#endif
// CODE: 2019/02/25. Do not delete this line.