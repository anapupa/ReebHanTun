/*
(c) 2012 Fengtao Fan
*/
#ifndef MY_psbmReebGraphElements_HEADER_H
#define MY_psbmReebGraphElements_HEADER_H

#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <queue>

#include "SimpleMesh.h"
#include <RenderVector3.h>

class psbmReebArc;

enum REEB_GRAPH_CRITICAL_TYPE {
    MINIMUM_REEB, UP_FORKING_REEB, DOWN_FORKING_REEB, MAXIMUM_REEB, UNKNOWN_REEB
};

struct DeleteObject {     // templatization and base
    // class removed here
    template<typename T>
    void operator()(const T *ptr) const {
        delete ptr;
    }
};

struct psbmArcListNode {
    int nEdgeId;    // index in mesh edges array
    class psbmReebArc *ptrReebArcId; // pointer to the list containing itself
    struct psbmArcListNode *ptrNextPos; // pointer to next arc
};

struct SimplifiedReebGraphArc {// No simplified version of Reeb graph nodes
    // These nodes are represented by nodes in original Reeb graph nodes
    int nCriticalNode0; // index of the vertex with lower value in Reeb graph nodes array
    int nCriticalNode1; // index of the vertex with higher value in Reeb graph nodes array
    int edgeLabel;
    bool InteriorIsInside;
    /*
     true --- interior of level set is inside of closed surface
     false -- interior of level set is outside of closed surface
    */
    /*
    |-----|		|-----|
    | in  |		| out |
    |-----|		|-----|
    */
};

struct auxMeshEdge {
    class psbmReebArc *ptrReebArcId; // the highest arc intersected by this
    struct psbmArcListNode *ptrPos;
};


class psbmReebNode {
public:
    psbmReebNode() {
        nVertexId = -1; // index in mesh vertices array
//		 value = 0;		// scalar value on the vertex
        ptrArcDownId = NULL;    // upward arcs
        ptrArcUpId = NULL;        // downward arcs
        ptrHArc = NULL;            // horizontal arcs , not used
        isFinalized = false;    // finalzied , not used
        isCritical = false;        // critical point flag
        //
        crType = UNKNOWN_REEB;
    }

    psbmReebNode(const int n)//, const float f)
    {
        nVertexId = n;
        //value = f;
        ptrArcDownId = new std::list<class psbmReebArc *>;
        ptrArcUpId = new std::list<class psbmReebArc *>;
        ptrHArc = new std::list<class psbmReebArc *>;
        isFinalized = false;
        isCritical = false;
        //
        crType = UNKNOWN_REEB;
    }

    psbmReebNode(const psbmReebNode &rhs) {
        nVertexId = rhs.nVertexId;
        //value = rhs.value;
        isFinalized = rhs.isFinalized;
        isCritical = rhs.isCritical;
        //
        crType = rhs.crType;
        //
        if (rhs.ptrArcDownId) {
            ptrArcDownId = new std::list<class psbmReebArc *>;
            ptrArcDownId->assign(rhs.ptrArcDownId->begin(), rhs.ptrArcDownId->end());
        } else {
            ptrArcDownId = NULL;
        }
        //
        if (rhs.ptrArcUpId) {
            ptrArcUpId = new std::list<class psbmReebArc *>;
            ptrArcUpId->assign(rhs.ptrArcUpId->begin(), rhs.ptrArcUpId->end());
        } else {
            ptrArcUpId = NULL;
        }
        //
        if (rhs.ptrHArc) {
            ptrHArc = new std::list<class psbmReebArc *>;
            ptrHArc->assign(rhs.ptrHArc->begin(), rhs.ptrHArc->end());
        } else {
            ptrHArc = NULL;
        }
        //
    }

    psbmReebNode &operator=(const psbmReebNode &rhs) {
        nVertexId = rhs.nVertexId;
        //value = rhs.value;
        isFinalized = rhs.isFinalized;
        isCritical = rhs.isCritical;
        //
        crType = rhs.crType;
        //
        if (ptrArcDownId) {// DEALLOCATE THE MEMORY OF list
            {
                std::list<class psbmReebArc *> tmpList;
                ptrArcDownId->swap(tmpList);
            }
        } else
            ptrArcDownId = new std::list<class psbmReebArc *>;

        ptrArcDownId->assign(rhs.ptrArcDownId->begin(), rhs.ptrArcDownId->end());
        //
        if (ptrArcUpId) {// DEALLOCATE THE MEMORY OF list
            {
                std::list<class psbmReebArc *> tmpList;
                ptrArcUpId->swap(tmpList);
            }
        } else
            ptrArcUpId = new std::list<class psbmReebArc *>;

        ptrArcUpId->assign(rhs.ptrArcUpId->begin(), rhs.ptrArcUpId->end());
        //
        if (ptrHArc) {// DEALLOCATE THE MEMORY OF list
            {
                std::list<class psbmReebArc *> tmpList;
                ptrHArc->swap(tmpList);
            }
        } else
            ptrHArc = new std::list<class psbmReebArc *>;

        ptrHArc->assign(rhs.ptrHArc->begin(), rhs.ptrHArc->end());
        //
        return *this;
    }

    ~psbmReebNode() {
        //
        if (ptrArcDownId) {// DEALLOCATE THE MEMORY OF list
            {
                std::list<class psbmReebArc *> tmpList;
                ptrArcDownId->swap(tmpList);
            }
            delete ptrArcDownId;
        }
        //
        if (ptrArcUpId) {// DEALLOCATE THE MEMORY OF list
            {
                std::list<class psbmReebArc *> tmpList;
                ptrArcUpId->swap(tmpList);
            }
            delete ptrArcUpId;
        }
        //
        if (ptrHArc) {// DEALLOCATE THE MEMORY OF list
            {
                std::list<class psbmReebArc *> tmpList;
                ptrHArc->swap(tmpList);
            }
            delete ptrHArc;
        }
    }

public:
    int nVertexId; // index in mesh vertices array
    //double value;   // scalar value on the vertex
    std::list<class psbmReebArc *> *ptrArcDownId;
    std::list<class psbmReebArc *> *ptrArcUpId;
    std::list<class psbmReebArc *> *ptrHArc;
    bool isFinalized;
    bool isCritical;
    // reserved for simplified reeb graph where two degree nodes are ignored
    REEB_GRAPH_CRITICAL_TYPE crType;
};

class psbmReebArc {

public:
    psbmReebArc(const int N0, const int N1, const int edgeId) {
        nNodeId0 = N0;
        nNodeId1 = N1;
        pArcList = new std::list<struct psbmArcListNode *>;
        psbmArcListNode *tmpNode = new psbmArcListNode;
        tmpNode->nEdgeId = edgeId;
        tmpNode->ptrNextPos = NULL;
        tmpNode->ptrReebArcId = this;
        pArcList->push_back(tmpNode);
        //
        clusterLabel = 0;

    }

    ~psbmReebArc() {
        if (pArcList) {
            std::for_each(pArcList->begin(), pArcList->end(), DeleteObject());
            pArcList->clear();
            delete pArcList;
        }
    }

    void ClearArcList() {
        if (pArcList) {
            std::for_each(pArcList->begin(), pArcList->end(), DeleteObject());
            pArcList->clear();
            //delete pArcList;
        }
    }

    psbmReebArc(const psbmReebArc &rhs) {
        nNodeId0 = rhs.nNodeId0;
        nNodeId1 = rhs.nNodeId1;
        pArcList = new std::list<struct psbmArcListNode *>;
        pArcList->assign(rhs.pArcList->begin(), rhs.pArcList->end());
    }

    psbmReebArc &operator=(const psbmReebArc &rhs) {
        nNodeId0 = rhs.nNodeId0;
        nNodeId1 = rhs.nNodeId1;
        if (pArcList) {
            std::for_each(pArcList->begin(), pArcList->end(), DeleteObject());
            pArcList->clear();
        }
        pArcList->assign(rhs.pArcList->begin(), rhs.pArcList->end());
        return *this;
    }

public:
    int nNodeId0; // index of the vertex with lower value in Reeb graph node array
    int nNodeId1; // index of the vertex with higher value in Reeb graph node array
    std::list<struct psbmArcListNode *> *pArcList;
    // used by simplified reeb graph
    int clusterLabel; // 0 is reserved for not assigning labels yet

};

class psbmReebPath {
public:
    psbmReebPath() {
        pathType = 0;
    }

    psbmReebPath(const psbmReebPath &rhs) {
        pathType = rhs.pathType;
        nodeIndexSet.assign(rhs.nodeIndexSet.begin(), rhs.nodeIndexSet.end());
        simpArcIndexSet.assign(rhs.simpArcIndexSet.begin(), rhs.simpArcIndexSet.end());
    }

    psbmReebPath &operator=(const psbmReebPath &rhs) {
        pathType = rhs.pathType;
        nodeIndexSet.assign(rhs.nodeIndexSet.begin(), rhs.nodeIndexSet.end());
        simpArcIndexSet.assign(rhs.simpArcIndexSet.begin(), rhs.simpArcIndexSet.end());

        return *this;
    }

    void Clear() {
        // DELLOCATE THE MEMORY OF vecPoints
        {
            std::vector<int> tmp;
            nodeIndexSet.swap(tmp);
        }
        // DELLOCATE THE MEMORY OF vecPoints
        {
            std::vector<int> tmp;
            simpArcIndexSet.swap(tmp);
        }
    }

    ~psbmReebPath() {
        // DELLOCATE THE MEMORY OF vecPoints
        {
            std::vector<int> tmp;
            nodeIndexSet.swap(tmp);
        }
        // DELLOCATE THE MEMORY OF vecPoints
        {
            std::vector<int> tmp;
            simpArcIndexSet.swap(tmp);
        }
    }

    /*
    encoded in this way
    v0--e0--v1--e1--....vn--en(--v0)
    */
    int pathType;
    //0 --- horizontal or in level set
    //1 --- vertical or orthgonal to level set
    std::vector<int> nodeIndexSet;
    std::vector<int> simpArcIndexSet;
};

class _Polygon {
public:
    _Polygon() {
    }

    ~_Polygon() {
        // DELLOCATE THE MEMORY OF vecPoints
        std::vector<Vector3> tmp;
        vecPoints.swap(tmp);
    }

    _Polygon(const _Polygon &rhs) {
        {// DELLOCATE THE MEMORY OF vecPoints
            std::vector<Vector3> tmp;
            vecPoints.swap(tmp);
        }
        vecPoints.assign(rhs.vecPoints.begin(), rhs.vecPoints.end());
    }

    _Polygon &operator=(const _Polygon &rhs) {
        if (this == &rhs)
            return *this;
        {// DELLOCATE THE MEMORY OF vecPoints
            std::vector<Vector3> tmp;
            vecPoints.swap(tmp);
        }
        vecPoints.assign(rhs.vecPoints.begin(), rhs.vecPoints.end());
        return *this;
    }

    void Clear() {
        {// DELLOCATE THE MEMORY OF vecPoints
            std::vector<Vector3> tmp;
            vecPoints.swap(tmp);
        }
    }

    double TotalLength() {
        double len = 0.0;
        Vector3 diff;
        for (unsigned int i = 0; i < vecPoints.size() - 1; i++) {
            diff = vecPoints[i + 1] - vecPoints[i];
            len += norm(diff);
        }
        return len;
    }

public:
    std::vector<Vector3> vecPoints;
};

/*************
simple graph data structure
**************************/
class _simpGraphNode {
public:
    _simpGraphNode() : selfIndex(-1), color(0) {
        value = 0.0;
        // set the value
    }

    _simpGraphNode(const int index) : selfIndex(index), color(0) {
        value = 0.0;
    }

    _simpGraphNode(const int index, const int cc) : selfIndex(index), color(cc) {
        value = 0.0;
    }

    _simpGraphNode(const int index, const int cc, const double inVal) : selfIndex(index), color(cc), value(inVal) {

    }

    _simpGraphNode(const _simpGraphNode &rhs) {
        selfIndex = rhs.selfIndex;
        color = rhs.color;
        //
        value = rhs.value;

        std::copy(rhs.adjList.begin(),
                  rhs.adjList.end(),
                  std::inserter(adjList, adjList.begin()));
    }

    _simpGraphNode &operator=(const _simpGraphNode &rhs) {
        if (this == &rhs)
            return *this;
        adjList.clear();
        selfIndex = rhs.selfIndex;
        color = rhs.color;
        //
        value = rhs.value;

        std::copy(rhs.adjList.begin(),
                  rhs.adjList.end(),
                  std::inserter(adjList, adjList.begin()));
        return *this;
    }

    ~_simpGraphNode() {
        adjList.clear();
    }

    /************************************/
    bool is_neighbor(const int u) {
        std::set<int>::iterator sIter;
        sIter = adjList.find(u);
        if (sIter == adjList.end())
            return false;
        return true;
    }

    bool is_neighbor(const int u) const {
        std::set<int>::const_iterator sIter;
        sIter = adjList.find(u);
        if (sIter == adjList.end())
            return false;
        return true;
    }

    int degree() {
        return int(adjList.size());
    }

public:
    int selfIndex;
    int color; // color is mapped to the index  in mesh vertex array
    // in this program
    std::set<int> adjList;
    //
    double value; // the height of this point
};

class _simpGraph {
public:
    _simpGraph() {
    }

    _simpGraph(const int n) {
        _simpGraphNode tmpNode;
        vecNode.reserve(n);
        for (int i = 0; i < n; i++) {
            tmpNode.selfIndex = i;
            vecNode.push_back(tmpNode);
        }
    }

    _simpGraph(const _simpGraph &rhs) {
        vecNode.assign(rhs.vecNode.begin(),
                       rhs.vecNode.end());
    }

    _simpGraph &operator=(const _simpGraph &rhs) {
        vecNode.clear();
        vecNode.assign(rhs.vecNode.begin(),
                       rhs.vecNode.end());
        return *this;
    }

    ~_simpGraph() {
        vecNode.clear();
    }

    /********************************/
    void AddEdge(const int i, const int j) {
        vecNode[i].adjList.insert(j);
        vecNode[j].adjList.insert(i);
        return;
    }

    void RemoveEdge(const int i, const int j) {// precondition: the edge <i,j> is there
        std::set<int>::iterator sIter;
        //
        sIter = vecNode[i].adjList.find(j);
        if (sIter != vecNode[i].adjList.end())
            vecNode[i].adjList.erase(sIter);
        //
        sIter = vecNode[j].adjList.find(i);
        if (sIter != vecNode[j].adjList.end())
            vecNode[j].adjList.erase(sIter);
    }

    void RemoveEdgesAdjacentToVertex(const int i) {
        std::vector<int> vecAdjVertex(vecNode[i].adjList.begin(),
                                      vecNode[i].adjList.end());
        for (int j = 0; j < int(vecAdjVertex.size()); j++) {
            RemoveEdge(i, vecAdjVertex[j]);
        }
    }

    void AddNode(const int selfIndex) {
        _simpGraphNode tmpNode(selfIndex);
        vecNode.push_back(tmpNode);
    }

    void AddNode(const int selfIndex, const int color) {
        _simpGraphNode tmpNode(selfIndex, color);
        vecNode.push_back(tmpNode);
    }

    void AddNode(const _simpGraphNode &gnode) {
        vecNode.push_back(gnode);
    }

    void InitNodes(const int n) {
        for (int selfIndex = 0; selfIndex < n; selfIndex++) {
            _simpGraphNode tmpNode(selfIndex);
            vecNode.push_back(tmpNode);
        }
    }

    void ReadFromFile(char const *pFileName);

    void WriteBackToFile(char const *pFileName);

    long int EdgeNum();

    void BreathFirstSearch(const int sourceVertex, std::vector<int> &parents) {// color : 0 --- untouched
        //			1 --- grey (touched)
        //			2 --- black (finished)
        std::vector<int> color;
        color.resize(vecNode.size());
        parents.resize(vecNode.size());
        // initializing color and parents
        for (int i = 0; i < int(color.size()); i++) {
            color[i] = 0; // untouched
            parents[i] = -1;// no parents
        }
        // set up a queue
        std::queue<int> Q;
        // push up source into queue
        // set its color as touched
        color[sourceVertex] = 1;
        Q.push(sourceVertex);
        while (!Q.empty()) {
            int curVertex = Q.front();
            Q.pop();
            for (std::set<int>::iterator sIter = vecNode[curVertex].adjList.begin();
                 sIter != vecNode[curVertex].adjList.end();
                 ++sIter) {
                if (color[*sIter] == 0) {// current vertex is not touched
                    color[*sIter] = 1;
                    parents[*sIter] = curVertex;
                    Q.push(*sIter);
                }
            }
            color[curVertex] = 2;
        }
    }

    void
    BreathFirstSearch_height_priority(const int sourceVertex, std::vector<int> &parents) {// color : 0 --- untouched
        //			1 --- grey (touched)
        //			2 --- black (finished)
        std::vector<int> color;
        color.resize(vecNode.size());
        parents.resize(vecNode.size());
        // initializing color and parents
        for (int i = 0; i < int(color.size()); i++) {
            color[i] = 0; // untouched
            parents[i] = -1;// no parents
        }
        // set up a queue
        std::queue<int> Q;
        // push up source into queue
        // set its color as touched
        color[sourceVertex] = 1;
        Q.push(sourceVertex);
        while (!Q.empty()) {
            int curVertex = Q.front();
            Q.pop();
            for (std::set<int>::iterator sIter = vecNode[curVertex].adjList.begin();
                 sIter != vecNode[curVertex].adjList.end();
                 ++sIter) {
                if (color[*sIter] == 0 &&
                    vecNode[curVertex].value >= vecNode[*sIter].value) {// current vertex is not touched
                    color[*sIter] = 1;
                    parents[*sIter] = curVertex;
                    Q.push(*sIter);
                }
            }
            color[curVertex] = 2;
        }
    }

    void Clear() {
        std::vector<_simpGraphNode> tmp;
        vecNode.swap(tmp);
    }

public:
    std::vector<_simpGraphNode> vecNode;
};

#endif /* MY_psbmReebGraphElements_HEADER_H */
