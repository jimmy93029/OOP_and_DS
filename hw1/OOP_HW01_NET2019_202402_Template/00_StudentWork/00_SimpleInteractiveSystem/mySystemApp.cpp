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

#define STUDENT_INFO "Name:XYZ   ID:0123456789"

void MY_SYSTEM_APP::showMyStudentInfo_2024( ) const
{
	cout << "*******************************" << endl;
    cout << "Date:2024/03/03" << endl;
	cout << "Student ID:\t" << endl;
	cout << "Student Name:\t" << endl;
	cout << "Student Email:\t" << endl;
	cout << "*******************************" << endl;
}

MY_SYSTEM_APP::MY_SYSTEM_APP( )
{
	//mSystemType = SYSTEM_TYPE_MONTE_CARLO_CIRCLE;
	mSystemType = SYSTEM_TYPE_GRAPH;
    string title = STUDENT_INFO;
        glutSetWindowTitle(title.data( ));
}

void MY_SYSTEM_APP::initApp( )
{
	mFlgShow_Grid = true;
}

void MY_SYSTEM_APP::update( )
{

}
// CODE: 2019/02/25. Do not delete this line.



