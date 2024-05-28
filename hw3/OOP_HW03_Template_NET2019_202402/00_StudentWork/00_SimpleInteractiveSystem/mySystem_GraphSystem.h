//********************************************
// Student Name			: jimmywu
// Student ID			: 111652040
// Student Email Address: jimmywu0229.sc11@nycu.edu.tw
//********************************************
//
//
// Prof. Sai-Keung Wong
// Email:	cswingo@nycu.edu.tw
//
// National Yang Ming Chiao Tung University, Taiwan
// Computer Science
// Date: 2024/05/01
//
//
// Add your own functions if necessary
//
#ifndef __MY_GRAPH_SYSTEM_H_
#define __MY_GRAPH_SYSTEM_H_
#include "../headers.h"
#include <string>
#include "../primitives/vector3d.h"
#include <vector>
#include "base_class.h"
#include "graph_basics.h"

using std::vector;

class Param {
public: 
    static int GRAPH_MAX_NUM_NODES;
    static int GRAPH_MAX_NUM_EDGES;
    static double PI;
};


class GRAPH_SYSTEM : public BASE_SYSTEM {

protected:
    void drawNodes_IndirectMethod( ) const;
    void drawEdges_IndirectMethod( ) const;
    void drawSelectedNode_Indirect( ) const;
    //
    GRAPH_NODE *GRAPH_SYSTEM::findNearestNode( double x, double z, double &cur_distance2 ) const;
    void deleteNode( int nodeID );
    void deleteSelectedNode(  );
    //
    void createDefaultGraph( );
    void createNet_Circular( int n, int num_layers );
    void createNet_Square( int n, int num_layers );
    void createNet_RadialCircular( int n );
    void createRandomGraph_DoubleCircles(int n);
    //
    // Implement your own stuff
    //
    void initMemoryPool();
    GRAPH_NODE* getFreeNode();
    GRAPH_EDGE* getFreeEdge();
    void deleteEdge(int EdgeID);

public:
    GRAPH_SYSTEM( );
    void update( );
    void stopAutoNodeDeletion();
    void reset( );
    void askForInput( );
    void clickAt(double x, double z);
    void handleKeyPressedEvent( unsigned char key );
    void handlePassiveMouseEvent( double x, double y );
    void draw( ) const;
    
    int addNode( float x, float y, float z, float r = 1.0 );
    int addEdge( int nodeID_0, int nodeID_1 );
    bool isSelectedNode( ) const;
    void getInfoOfSelectedPoint( double &r, vector3 &p ) const;
    //
    int getNumOfEdges( ) const;
    vector3 getNodePositionOfEdge( int edgeIndex, int nodeIndex ) const;

    int getNumOfNodes( ) const;
    void getNodeInfo( int nodeIndex, double &r, vector3 &p ) const;
protected:
    //
    // Implement your own stuff
    //
    //
    GRAPH_NODE *mSelectedNode;
    GRAPH_NODE *mPassiveSelectedNode;
    //
    bool mFlgAutoNodeDeletion;
    //
    GRAPH_NODE* mNodeArr_Pool;
    GRAPH_EDGE* mEdgeArr_Pool;
    //
    int* mActiveNodeArr;
    int mCurNumOfActiveNodes;
    int* mActiveEdgeArr;
    int mCurNumOfActiveEdges;
    //
    int* mFreeNodeArr;
    int mCurNumOfFreeNodes;
    int* mFreeEdgeArr;
    int mCurNumOfFreeEdges;
    //
    bool mode5;
    int Nmode5;
};

#endif