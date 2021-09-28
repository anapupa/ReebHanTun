/*
(c) 2012 Fengtao Fan
*/
#ifndef _MY_SIMPLE_GRAPH_H_
#define _MY_SIMPLE_GRAPH_H_

#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
#include <iterator>

class _simpGraphNode_vec {
public:
    _simpGraphNode_vec() : selfIndex(-1), color(0) {
        value = 0.0;
        // set the value
    }

    _simpGraphNode_vec(const int index) : selfIndex(index), color(0) {
        value = 0.0;
    }

    _simpGraphNode_vec(const int index, const int cc) : selfIndex(index), color(cc) {
        value = 0.0;
    }

    _simpGraphNode_vec(const int index, const int cc, const float inVal) : selfIndex(index), color(cc), value(inVal) {

    }

    _simpGraphNode_vec(const _simpGraphNode_vec &rhs) {
        selfIndex = rhs.selfIndex;
        color = rhs.color;
        //
        value = rhs.value;

        adjVertexIdVec.assign(rhs.adjVertexIdVec.begin(), rhs.adjVertexIdVec.end());
        adjEdgeIdVec.assign(rhs.adjEdgeIdVec.begin(), rhs.adjEdgeIdVec.end());
        //
        /*	std::copy (	rhs.adjList.begin(),
                        rhs.adjList.end(),
                        std::inserter(adjList, adjList.begin()));*/
    }

    void CopyElementsVector(std::vector<int> &dst, const std::vector<int> &inVec) {
        dst.clear();
        dst.assign(inVec.begin(), inVec.end());
    }

    _simpGraphNode_vec &operator=(const _simpGraphNode_vec &rhs) {
        if (this == &rhs)
            return *this;
        CopyElementsVector(adjVertexIdVec, rhs.adjVertexIdVec);
        //
        CopyElementsVector(adjEdgeIdVec, rhs.adjEdgeIdVec);
        //
        selfIndex = rhs.selfIndex;
        color = rhs.color;
        //
        value = rhs.value;
        //
        return *this;
    }

    ~_simpGraphNode_vec() {
        adjVertexIdVec.clear();
        adjEdgeIdVec.clear();
    }

    /************************************/
    bool is_neighbor(const int u) {
        std::set<int> sorted_ver_id(adjVertexIdVec.begin(), adjVertexIdVec.end());
        std::set<int>::const_iterator sIter;
        sIter = sorted_ver_id.find(u);
        if (sIter == sorted_ver_id.end())
            return false;
        return true;
    }

    bool is_neighbor(const int u) const {
        std::set<int> sorted_ver_id(adjVertexIdVec.begin(), adjVertexIdVec.end());
        std::set<int>::const_iterator sIter;
        sIter = sorted_ver_id.find(u);
        if (sIter == sorted_ver_id.end())
            return false;
        return true;
    }

    int degree() {
        return int(adjVertexIdVec.size());
    }

public:
    int selfIndex;
    int color; // color is mapped to the index  in mesh vertex array
    // in this program
    std::vector<int> adjVertexIdVec;
    std::vector<int> adjEdgeIdVec;
    //
    float value; // the height of this point
};

class _simpGraph_vec {
public:
    _simpGraph_vec() : edge_size(0) {
    }

    _simpGraph_vec(const int n) : edge_size(0) {
        _simpGraphNode_vec tmpNode;
        vecNode.reserve(n);
        for (int i = 0; i < n; i++) {
            tmpNode.selfIndex = i;
            vecNode.push_back(tmpNode);
        }
    }

    _simpGraph_vec(const _simpGraph_vec &rhs) {
        vecNode.assign(rhs.vecNode.begin(),
                       rhs.vecNode.end());
        edge_size = rhs.edge_size;
    }

    _simpGraph_vec &operator=(const _simpGraph_vec &rhs) {
        if (this == &rhs)
            return *this;
        vecNode.clear();
        vecNode.assign(rhs.vecNode.begin(),
                       rhs.vecNode.end());
        edge_size = rhs.edge_size;
        return *this;
    }

    ~_simpGraph_vec() {
        vecNode.clear();
    }

    /********************************/
    void AddEdge(const int i, const int j) {
        vecNode[i].adjVertexIdVec.push_back(j);
        vecNode[j].adjVertexIdVec.push_back(i);
        // use the default order for edge
        vecNode[i].adjEdgeIdVec.push_back(edge_size);
        vecNode[j].adjEdgeIdVec.push_back(edge_size);
        //
        edge_size++;
        return;
    }

    void AddEdge(const int i, const int j, const int eid) {
        vecNode[i].adjVertexIdVec.push_back(j);
        vecNode[j].adjVertexIdVec.push_back(i);
        // use the input order for edge
        vecNode[i].adjEdgeIdVec.push_back(eid);
        vecNode[j].adjEdgeIdVec.push_back(eid);
        //
        edge_size++;
        return;
    }

    void RemoveEdge(const int i, const int j) {// precondition: the edge <i,j> is there
        std::cout << "NEVER CALL THIS RemoveEdge FUNCTION" << std::endl;
        //
        edge_size--;
    }

    void RemoveEdgesAdjacentToVertex(const int i) {
        std::cout << "NEVER CALL THIS RemoveEdgesAdjacentToVertex FUNCTION" << std::endl;
    }

    void AddNode(const int selfIndex) {
        _simpGraphNode_vec tmpNode(selfIndex);
        vecNode.push_back(tmpNode);
    }

    void AddNode(const int selfIndex, const int color) {
        _simpGraphNode_vec tmpNode(selfIndex, color);
        vecNode.push_back(tmpNode);
    }

    void AddNode(const _simpGraphNode_vec &gnode) {
        vecNode.push_back(gnode);
    }

    void InitNodes(const int n) {
        for (int selfIndex = 0; selfIndex < n; selfIndex++) {
            _simpGraphNode_vec tmpNode(selfIndex, selfIndex);
            vecNode.push_back(tmpNode);
        }
    }

    void Clear() {
        std::vector<_simpGraphNode_vec> tmp;
        vecNode.swap(tmp);
    }

public:
    std::vector<_simpGraphNode_vec> vecNode;
    int edge_size;
};

#endif //_MY_SIMPLE_GRAPH_H_
