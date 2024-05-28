//********************************************
// Student Name			: jimmywu
// Student ID			: 111652040
// Student Email Address: jimmywu0229.sc11@nycu.edu.tw
//********************************************
//
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
#include <iostream>
#include "mySystem_GraphSystem.h"
#include <time.h>
#include <cmath>
#include <vector>
#include <random>
#include <cstdlib>

using namespace std;

int Param::GRAPH_MAX_NUM_NODES = 10000;
int Param::GRAPH_MAX_NUM_EDGES = 10000;
double Param::PI = acos(-1);


GRAPH_SYSTEM::GRAPH_SYSTEM( )
{
    mFlgAutoNodeDeletion = false;
    //
    // Implement your own stuff
    //
    initMemoryPool();
    //
    createDefaultGraph();
    //
}

void GRAPH_SYSTEM::initMemoryPool()
{
    mNodeArr_Pool = new GRAPH_NODE[Param::GRAPH_MAX_NUM_NODES];
    mEdgeArr_Pool = new GRAPH_EDGE[Param::GRAPH_MAX_NUM_EDGES];

    mCurNumOfActiveNodes = 0;
    mCurNumOfActiveEdges = 0;
    mActiveNodeArr = new int[Param::GRAPH_MAX_NUM_NODES];
    mActiveEdgeArr = new int[Param::GRAPH_MAX_NUM_EDGES];

    mFreeNodeArr = new int[Param::GRAPH_MAX_NUM_NODES];
    mFreeEdgeArr = new int[Param::GRAPH_MAX_NUM_EDGES];
    //
    for (int i = 0; i < Param::GRAPH_MAX_NUM_NODES; ++i) {
        mNodeArr_Pool[i].id = i; // assign a unique id
    }
    for (int i = 0; i < Param::GRAPH_MAX_NUM_EDGES; ++i) {
        mEdgeArr_Pool[i].id = i; // assign a unique id
    }
    reset();
}

void GRAPH_SYSTEM::reset( )
{
    stopAutoNodeDeletion();
    mPassiveSelectedNode = 0;
    mSelectedNode = 0;
    //
    // Implement your own stuff
    //

    mCurNumOfActiveNodes = 0;
    mCurNumOfActiveEdges = 0;

    mCurNumOfFreeNodes = Param::GRAPH_MAX_NUM_NODES;
    mCurNumOfFreeEdges = Param::GRAPH_MAX_NUM_EDGES;

    for (int i = 0; i < Param::GRAPH_MAX_NUM_NODES; ++i) {
        mFreeNodeArr[i] = i; // index is not used
    }
    for (int i = 0; i < Param::GRAPH_MAX_NUM_EDGES; ++i) {
        mFreeEdgeArr[i] = i; // index is not used
    }
}


void GRAPH_SYSTEM::createDefaultGraph( )
{
    reset( );
    //
    // Implement your own stuff
    //
    double r = 1.0;
    int id0, id1, id2;
    id0 = addNode(0.0, 0.0, 0.0, r);
    id1 = addNode(5.0, 0.0, 0.0, r);
    id2 = addNode(0.0, 0.0, 5.0, r); // not x-y plane? then try x-z plane

    addEdge(id0, id1);
    addEdge(id2, id1);
}

void GRAPH_SYSTEM::createNet_Circular( int n, int num_layers )
{
    reset( );

    float r = 5; // radius
    float d = 5; // layer distance 
    float offset_x = 15.;
    float offset_z = 15.;
    //
    // Implement your own stuff
    //
    vector<vector<int>> layers(num_layers + 1 ,vector<int>(n, 0));
    double da = 2 * Param::PI / n;
    double a0, r0, id0;

    for (int i = 0; i <= num_layers; i++) {
        for (int j = 0; j < n; j++)
        {
            a0 = j * da;
            r0 = r + d * i;
            id0 = addNode(offset_x + r0 * cos(a0), 0, offset_z + r0 * sin(a0));
            layers[i][j] = id0;

            if (i - 1 >= 0) {
                addEdge(id0, layers[i - 1][j]);
            }
            if (j - 1 >= 0 && i < num_layers) {
                addEdge(id0, layers[i][j - 1]);
            }
            if (j == n - 1 && i < num_layers) {
                addEdge(id0, layers[i][0]);
            }
        }
    }
}

void GRAPH_SYSTEM::createNet_Square(int n, int num_layers)
{
    reset();

    float d = 5; // layer distance 
    float offset_x = 5.;
    float offset_z = 5.;
    //
    // Implement your own stuff
    //
    int inner_l = (num_layers - n) / 2;
    vector<vector<bool>> valid(num_layers, vector<bool>(num_layers, true));

    for (int i = inner_l; i < inner_l + n; i++) {
        for (int j = inner_l; j < inner_l + n; j++) {
            valid[i][j] = false;
        }
    }

    int x, z, id0, id1;
    vector<vector<int>> layers(num_layers, vector<int>(num_layers, 0));

    for (int i = 0; i < num_layers; i++) {
        for (int j = 0; j < num_layers; j++) 
        {
            if (!valid[i][j])
                continue;

            x = i * d;
            z = j * d;
            id0 = addNode(offset_x + x, 0, offset_z + z);
            layers[i][j] = id0;
            
            if (i - 1 >= 0 && valid[i - 1][j]) {
                addEdge(id0, layers[i - 1][j]);
            }
            if (j - 1 >= 0 && valid[i][j - 1]) {
                addEdge(id0, layers[i][j - 1]);
            }
        }
    }
}

void GRAPH_SYSTEM::createNet_RadialCircular( int n ) {

    reset( );

    float offset_x = 15.0;
    float offset_z = 15.0;

    float r = 15; // radius
    //
    // Implement your own stuff
    //   
    double da = 2 * Param::PI / n, a;
    int idr, id1;
    idr = addNode(offset_x, 0, offset_z);

    for (int i = 0; i < n; i++) 
    {
        a = da * i;
        id1 = addNode(r * cos(a) + offset_x, 0, r * sin(a) + offset_z);
        addEdge(idr, id1);
    }
}

void GRAPH_SYSTEM::createRandomGraph_DoubleCircles(int n)
{
    reset();

    float dx = 5.0;
    float dz = 5.0;
    float r = 15; // radius
    float d = 10; // layer distance
    float offset_x = 15.;
    float offset_z = 15.;
    //
    // Implement your own stuff
    //
    vector<vector<int>> layers(2, vector<int>(n, 0));
    double da = 2 * Param::PI / n;
    double a0, r0, id0;

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < n; j++)
        {
            a0 = da * j;
            r0 = r + d * i;
            id0 = addNode(offset_x + r0 * cos(a0), 0, offset_z + r0 * sin(a0));
            layers[i][j] = id0;
        }
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    double range_a = acos(r / (d + r));
    int range = range_a / da;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        {
            if (dis(gen) <= 0.5)
                continue;
            if (i - j > range || i - j < -range)
                continue;

            addEdge(layers[0][j], layers[1][i]);
        }
    }

}


void GRAPH_SYSTEM::askForInput()
{
    //
    // Implement your own stuff
    // 
    cout << "GRAPH_SYSTEM" << endl;
    cout << "Key usage:" << endl;
    cout << "1: create a default graph" << endl;
    cout << "2: create a graph with 10x10 nodes. Connect the consecutive nodes horizontally" << endl;
    cout << "3: create a graph with 10x10 nodes. Connect the consecutive nodes vertically" << endl;
    cout << "4: create a graph with 10x10 nodes. Create 10 randomly generated edges" << endl;
    cout << "5: create a graph with 10x10 nodes. Create 10 randomly generated edges attached at a random node" << endl;
    cout << "Delete: delete a node and all the edges attached at it" << endl;
    cout << "Spacebar: unselect the selected node" << endl;
    cout << " " << endl;
    cout << "Use the mouse to select nodes and add edges" << endl;
    cout << "Click the left button to select/unselect or create an edge" << endl;
    cout << " " << endl;
    cout << "A selected node is highlighted as red." << endl;
}


// return node id
int GRAPH_SYSTEM::addNode( float x, float y, float z, float r)
{
    //
    // Implement your own stuff
    //  
    GRAPH_NODE* g;
    g = getFreeNode();
    if (g == 0) return -1; // invalid id.
    g->p = vector3(x, y, z);
    g->r = r;
    g->edgeID.clear();
    return g->id;
}

GRAPH_NODE* GRAPH_SYSTEM::getFreeNode() {
    if (mCurNumOfFreeNodes == 0) return 0;
    --mCurNumOfFreeNodes;
    int id = mFreeNodeArr[mCurNumOfFreeNodes];
    GRAPH_NODE* n = &mNodeArr_Pool[id];
    mActiveNodeArr[mCurNumOfActiveNodes] = id;
    n->dynamicID = mCurNumOfActiveNodes;
    ++mCurNumOfActiveNodes;
    return n;
}


// return edge id
int GRAPH_SYSTEM::addEdge( int nodeID_0, int nodeID_1 )
{
    //
    // Implement your own stuff
    //  
    GRAPH_EDGE* e;
    e = getFreeEdge();
    if (e == 0) return -1;
    e->nodeID[0] = nodeID_0;
    e->nodeID[1] = nodeID_1;
    mNodeArr_Pool[nodeID_0].edgeID.push_back(e->id);
    mNodeArr_Pool[nodeID_1].edgeID.push_back(e->id);
    return e->id;
}

GRAPH_EDGE* GRAPH_SYSTEM::getFreeEdge() 
{
    if (mCurNumOfFreeEdges == 0) return 0;
    --mCurNumOfFreeEdges;
    int id = mFreeEdgeArr[mCurNumOfFreeEdges];
    GRAPH_EDGE* e = &mEdgeArr_Pool[id];
    mActiveEdgeArr[mCurNumOfActiveEdges] = id;
    e->dynamicID = mCurNumOfActiveEdges;
    ++mCurNumOfActiveEdges;
    return e;
}


void GRAPH_SYSTEM::deleteNode(int nodeID) 
{
    //
    // Implement your own stuff
    // 
    
    GRAPH_NODE* node = &mNodeArr_Pool[nodeID];

    // remove edges incident to this node
    while (!node->edgeID.empty())
    {
        deleteEdge(node->edgeID[0]);
    }

    // add node back to FN
    mFreeNodeArr[mCurNumOfFreeNodes] = nodeID;
    mCurNumOfFreeNodes++;

    // remove node from AN
    int pos_in_AN = node->dynamicID;
    mCurNumOfActiveNodes--;
    mActiveNodeArr[pos_in_AN] = mActiveNodeArr[mCurNumOfActiveNodes];

    GRAPH_NODE* update = &mNodeArr_Pool[mActiveNodeArr[mCurNumOfActiveNodes]];
    update->dynamicID = node->dynamicID;

}


void GRAPH_SYSTEM::deleteEdge(int EdgeID) 
{
    // add node back to FE
    GRAPH_EDGE* edge = &mEdgeArr_Pool[EdgeID];
    mFreeEdgeArr[mCurNumOfFreeEdges] = EdgeID;
    mCurNumOfFreeEdges++;

    // remove node from AE
    int pos_in_AE = edge->dynamicID;
    mCurNumOfActiveEdges--;
    mActiveEdgeArr[pos_in_AE] = mActiveEdgeArr[mCurNumOfActiveEdges];

    GRAPH_EDGE* update = &mEdgeArr_Pool[mActiveEdgeArr[mCurNumOfActiveEdges]];
    update->dynamicID = edge->dynamicID;

    // remove edges from its incident node
    GRAPH_NODE* node;
    for (int i = 0; i < 2; i++) 
    {
        node = &mNodeArr_Pool[edge->nodeID[i]];
        vector<int>::iterator it = find(node->edgeID.begin(), node->edgeID.end(), edge->id);
        node->edgeID.erase(it);
    }

}


void GRAPH_SYSTEM::deleteSelectedNode() {
    //
    // Implement your own stuff
    // 
    if (isSelectedNode()) {
        deleteNode(mSelectedNode->id);
        mSelectedNode = 0;
    }

}


//
// get the number of nodes
//
int GRAPH_SYSTEM::getNumOfNodes() const
{
    //
    // Implement your own stuff
    // 
    return mCurNumOfActiveNodes;
}

void GRAPH_SYSTEM::getNodeInfo(int nodeIndex, double& r, vector3& p) const
{
    //
    // Implement your own stuff
    // 
    if (nodeIndex >= mCurNumOfActiveNodes) {
        return;
    }
    int id = mActiveNodeArr[nodeIndex];
    const GRAPH_NODE* g = &mNodeArr_Pool[id];
    r = g->r;
    p = g->p;

}

//
// return the number of edges
//
int GRAPH_SYSTEM::getNumOfEdges() const
{
    //
    // Implement your own stuff
    // 
    return mCurNumOfActiveEdges;
}

//
// an edge should have two nodes: index 0 and index 1
// return the position of node with nodeIndex
//
vector3 GRAPH_SYSTEM::getNodePositionOfEdge(int edgeIndex, int nodeIndex) const
{
    vector3 p;
    if (edgeIndex >= mCurNumOfActiveEdges)
        return p;
    //
    // Implement your own stuff
    // 
    int eid = mActiveEdgeArr[edgeIndex];
    const GRAPH_EDGE* e = &mEdgeArr_Pool[eid];

    int nid = e->nodeID[nodeIndex];
    const GRAPH_NODE* n = &mNodeArr_Pool[nid];

    return n->p;
}


void GRAPH_SYSTEM::handleKeyPressedEvent( unsigned char key )
{
    switch( key ) 
    {
    case 127: // delete
        mFlgAutoNodeDeletion = false;
        mode5 = false;
        deleteSelectedNode( );
        mSelectedNode = 0;
        break;
    case '1':
        mFlgAutoNodeDeletion = false;
        mode5 = false;
        createDefaultGraph( );
        mSelectedNode = 0;
        break;
    case '2':
        mFlgAutoNodeDeletion = false;
        mode5 = false;
        createNet_Circular(12, 3);
        mSelectedNode = 0;
        break;
    case '3':
        mFlgAutoNodeDeletion = false;
        mode5 = false;
        createNet_Square(3, 11); // you can modify this
        mSelectedNode = 0;
        break;
    case '4':
        mFlgAutoNodeDeletion = false;
        mode5 = false;
        createNet_RadialCircular(24);
        mSelectedNode = 0;
        break;
    case '5':
        mFlgAutoNodeDeletion = false;
        mode5 = true;
        Nmode5 = 24;
        createRandomGraph_DoubleCircles(Nmode5);
        mSelectedNode = 0;
        break;
    case '<':
        if (mode5) {
            Nmode5 = (Nmode5 == 3) ? Nmode5 : Nmode5 - 1;
            createRandomGraph_DoubleCircles(Nmode5);
        }
        break;
    case '>':
        if (mode5) {
            Nmode5 = (Nmode5 == 36) ? Nmode5 : Nmode5 + 1;
            createRandomGraph_DoubleCircles(Nmode5);
        }
        break;
    case ' ':
        mFlgAutoNodeDeletion = false;
        mode5 = false;
        mSelectedNode = 0;
        break;
    case 'd':
        mFlgAutoNodeDeletion = !mFlgAutoNodeDeletion;
        mode5 = false;
        break;
    }

}


bool GRAPH_SYSTEM::isSelectedNode() const
{
    //
    // Implement your own stuff
    // 
    if (mSelectedNode != 0)
        return true;
    return false;
}


void GRAPH_SYSTEM::getInfoOfSelectedPoint(double& r, vector3& p) const
{
    //
    // Implement your own stuff
    // 
    if (mSelectedNode != 0) {
        r = mSelectedNode->r;
        p = mSelectedNode->p;
    }

}


//
// compute mSelectedNode
//
void GRAPH_SYSTEM::clickAt(double x, double z)
{
    //
    // Implement your own stuff
    // 
    handlePassiveMouseEvent(x, z);

    if (mPassiveSelectedNode == 0) {
        return;
    }
    
    if (mSelectedNode == 0) {
        mSelectedNode = mPassiveSelectedNode;
    }
    else if (mPassiveSelectedNode == mSelectedNode) {
        mSelectedNode = 0;
    }
    else {
        bool EdgeExist = false;
        GRAPH_EDGE* e;

        for (int i = 0; i < mCurNumOfActiveEdges; i++)
        {
            e = &mEdgeArr_Pool[mActiveEdgeArr[i]];

            if (e->nodeID[0] == mPassiveSelectedNode->id &&
                e->nodeID[1] == mSelectedNode->id ) {
                EdgeExist = true;
            }
            if (e->nodeID[0] == mSelectedNode->id &&
                e->nodeID[1] == mPassiveSelectedNode->id) {
                EdgeExist = true;
            }
        }
        if (!EdgeExist) {
            addEdge(mPassiveSelectedNode->id, mSelectedNode->id);
        }

        mSelectedNode = mPassiveSelectedNode = 0;

    }

}


void GRAPH_SYSTEM::handlePassiveMouseEvent( double x, double z )
{
    double cur_d2;
    GRAPH_NODE *n = findNearestNode( x, z, cur_d2 );
    if ( n == 0 ) return;
    if ( cur_d2 > n->r*n->r ) {
        mPassiveSelectedNode = 0;
        return;
    }
    mPassiveSelectedNode = n;
}


GRAPH_NODE* GRAPH_SYSTEM::findNearestNode(double x, double z, double& cur_distance2) const
{
    GRAPH_NODE* n = 0;
    //
    // Implement your own stuff
    // 
    GRAPH_NODE* node;
    double x_, z_, mini, dist;

    for (int i = 0; i < mCurNumOfActiveNodes; i++)
    {
        node = &mNodeArr_Pool[mActiveNodeArr[i]];
        x_ = node->p[0];
        z_ = node->p[2];
        dist = pow(x - x_, 2) + pow(z - z_, 2);

        if (i == 0 || dist < mini) {
            mini = dist;
            n = node;
        }
    }
    return n;
}


void GRAPH_SYSTEM::stopAutoNodeDeletion()
{
    mFlgAutoNodeDeletion = false;
}


//
// For every frame, update( ) function is called.
// 
//
void GRAPH_SYSTEM::update( )
{
    if (!mFlgAutoNodeDeletion || mCurNumOfActiveNodes == 0) {
        return;
    }
    mSelectedNode = 0;
    mPassiveSelectedNode = 0;

    Sleep(250);
    //
    //
    // Implement your own stuff
    // 
    int random = rand() % mCurNumOfActiveNodes;
    int id = mActiveNodeArr[random];

    deleteNode(id);

}