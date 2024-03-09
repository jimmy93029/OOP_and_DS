//
// Instructor: Sai-Keung Wong
// Email:	cswingo@cs.nycu.edu.tw
//			wingo.wong@gmail.com
//
// National Yang Ming Chiao Tung University, Hsinchu, Taiwan
// Computer Science
// Date: 2024/03/03
//
#include "mySystemApp.h"

using namespace std;

void MY_SYSTEM_APP::drawOrigin() const
{
	glLineWidth(3.0);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(-100.0, 0.0, 0.0);
	glVertex3f(100.0, 0.0, 0.0);
	glEnd();
	//
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0.0, 0.0, -100.0);
	glVertex3f(0.0, 0.0, 100.0);
	glEnd();
}

void MY_SYSTEM_APP::draw() const
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glDisable(GL_LIGHTING);
	FREE_CANVAS_2D::begin_draw();

	glLineWidth(1.0);
	glColor3f(0.5, 0.5, 0.5);



	if (mFlgShow_Grid) {
		FREE_CANVAS_2D::draw_wire_grid();
	}

	draw_System_Content();

	drawOrigin();

	FREE_CANVAS_2D::end_draw();
	//cout << "END CROWD_SIMULATION_APP::draw( )" << endl;
	glEnable(GL_LIGHTING);
}

void MY_SYSTEM_APP::draw_System_Content() const
{
	switch (mSystemType) {
	case SYSTEM_TYPE_MONTE_CARLO_CIRCLE:
		break;
	case SYSTEM_TYPE_SINECOSINE_FUNCTION:
		break;
	case SYSTEM_TYPE_CUBIC_FUNCTION:
		break;
	case SYSTEM_TYPE_STUDENT_RECORD_MANAGER:
		break;
	case SYSTEM_TYPE_GRAPH:
		draw_Graphs();
		break;


	}
}


void MY_SYSTEM_APP::draw_Points(
	const vector<double>& X
	, const vector<double>& Y) const
{
	if (X.size() != Y.size()) return;
	glPointSize(15.0);

	int numSamples = X.size();

	//cout << "numSamples:" << numSamples << endl;

	glBegin(GL_POINTS);
	for (int i = 0; i < numSamples; ++i) {
		float x, z;
		x = X[i];
		z = Y[i];
		//glColor3f(1.0, 0.5, 0.0);
		glColor3f(0.0, 0.0, 0.7);
		glVertex3f(x, 0.0, z);
	}
	glEnd();
}

void MY_SYSTEM_APP::draw_Graphs() const
{
	vector<double> X;
	vector<double> Y;
	int numCurves = mGraphManager.getNumCurves();
	for (int i = 0; i < numCurves; ++i) {
		const CURVE_FUNCTION& curve = mGraphManager.getCurve(i);
		draw_CubicFunction(curve);
		curve.getExtremePoints(X, Y);
		draw_Points(X, Y);
		//
		draw_BoundaryPoints(curve);
	}



}

void MY_SYSTEM_APP::draw_BoundaryPoints(const CURVE_FUNCTION& curve) const
{
	double x, z;

	glPointSize(15);
	glBegin(GL_POINTS);
		curve.getBoundaryPoint(0, x, z);
		glColor3f(1.0, 0.2, 0.2);
		glVertex3f(x, 0.0, z);

		curve.getBoundaryPoint(1, x, z);
		glColor3f(0.2, 1.0, 0.2);
		glVertex3f(x, 0.0, z);
	glEnd();

	

}
void MY_SYSTEM_APP::draw_CubicFunction(const CURVE_FUNCTION& curve) const
{
	glLineWidth(6.0);
	glColor3f(0.2, 0.2, 0.2);

	double x0, x1;
	curve.getIntervalOfX(x0, x1);
	glBegin(GL_LINE_STRIP);
	int numSamples = curve.getNumOfSamples();
	for (int i = 0; i < numSamples; ++i) {
		float y;
		float x;
		x = x0 + (i / ((double)numSamples - 1)) * (x1 - x0);
		y = curve.getValue(x);
		glVertex3f(x, 0.0, y);
	}
	glEnd();

	glLineWidth(3.0);
	glColor3f(0.9, 0.9, 0.0);
	glBegin(GL_LINE_STRIP);
	numSamples = curve.getNumOfSamples();
	for (int i = 0; i < numSamples; ++i) {
		float y;
		float x;
		x = x0 + (i / ((double)numSamples - 1)) * (x1 - x0);
		y = curve.getValue(x);
		glVertex3f(x, 0.0, y);
	}
	glEnd();

	glLineWidth(2.0);

}

