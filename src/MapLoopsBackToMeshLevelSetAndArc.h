/*
(c) 2012 Fengtao Fan
*/
#ifndef _MY_MapLoopsBackToMeshLevelSetAndArc_H
#define _MY_MapLoopsBackToMeshLevelSetAndArc_H

#include "psbmReebGraphElements.h"
#include "IntersectionOnMesh.h"

#include <vector>

class MapLoopsBackToMeshLevelSetAndArc {
public:
    // constructor and assign operator
    MapLoopsBackToMeshLevelSetAndArc() {
        Null_all();
    }

    MapLoopsBackToMeshLevelSetAndArc(const MapLoopsBackToMeshLevelSetAndArc &rhs) {
        //
        inMeshPtr = rhs.inMeshPtr;
        //
        //
        _levelCyclePointTypePtr = rhs._levelCyclePointTypePtr;
        _arcCyclePointTypePtr = rhs._arcCyclePointTypePtr;
        //
        //
        outLevelCycle_VertexOnMesh.clear();
        outLevelCycle_VertexOnMesh = rhs.outLevelCycle_VertexOnMesh;
        outLevelCycle_EdgeOnMesh.clear();
        outLevelCycle_EdgeOnMesh = rhs.outArcCycle_EdgeOnMesh;
        //
        outArcCycle_VertexOnMesh.clear();
        outArcCycle_VertexOnMesh = rhs.outArcCycle_VertexOnMesh;
        outArcCycle_EdgeOnMesh.clear();
        outArcCycle_EdgeOnMesh = rhs.outArcCycle_EdgeOnMesh;
        //
    }

    MapLoopsBackToMeshLevelSetAndArc &operator=(const MapLoopsBackToMeshLevelSetAndArc &rhs) {
        if (this == &rhs)
            return *this;
        //
        //
        inMeshPtr = rhs.inMeshPtr;
        //
        _levelCyclePointTypePtr = rhs._levelCyclePointTypePtr;
        _arcCyclePointTypePtr = rhs._arcCyclePointTypePtr;
        //
        //
        outLevelCycle_VertexOnMesh.clear();
        outLevelCycle_VertexOnMesh = rhs.outLevelCycle_VertexOnMesh;
        outLevelCycle_EdgeOnMesh.clear();
        outLevelCycle_EdgeOnMesh = rhs.outArcCycle_EdgeOnMesh;
        //
        outArcCycle_VertexOnMesh.clear();
        outArcCycle_VertexOnMesh = rhs.outArcCycle_VertexOnMesh;
        outArcCycle_EdgeOnMesh.clear();
        outArcCycle_EdgeOnMesh = rhs.outArcCycle_EdgeOnMesh;
        //
        //
        return *this;
    }

    // deconstructor and data clean functioin
    void Null_all() {
        inMeshPtr = NULL;
        //
        _levelCyclePointTypePtr = NULL;
        _arcCyclePointTypePtr = NULL;
        //
    }

    void Clear() {
        //
        outLevelCycle_VertexOnMesh.clear();
        outLevelCycle_EdgeOnMesh.clear();
        //
        outArcCycle_VertexOnMesh.clear();
        outArcCycle_EdgeOnMesh.clear();
        //
    }

    void ClearOutPut() {
        //
        outLevelCycle_VertexOnMesh.clear();
        outLevelCycle_EdgeOnMesh.clear();
        //
        outArcCycle_VertexOnMesh.clear();
        outArcCycle_EdgeOnMesh.clear();
        //
        //
    }

    ~MapLoopsBackToMeshLevelSetAndArc() {
        Null_all();
        ClearOutPut();
    }

    // real operations functions
    // Initial data
    // any pair of arc and level set needs them
    void InitMeshPtr(_SimpleMesh *inPtr) {
        inMeshPtr = inPtr;
        return;
    }

    void InitScalarPtr(double *inPtr) {
        _scalarPtr = inPtr;
    }

    //
    // particular pair of arc and levelset needs them
    void InitLevelCyclePointType(std::vector<std::pair<int, int> > *levelCyclePointTypePtr) {
        ClearDataForLevelset();
        //
        _levelCyclePointTypePtr = levelCyclePointTypePtr;
        //
    }

    void InitArcCyclePointType(std::vector<std::pair<int, int> > *arcCyclePointTypePtr) {// pathArcOnMesh is initialized
        // clean memory first
        ClearDataForArcCycle();
        //
        _arcCyclePointTypePtr = arcCyclePointTypePtr;
        //
    }

    void ClearDataForLevelset() {
        //
        outLevelCycle_VertexOnMesh.clear();
        outLevelCycle_EdgeOnMesh.clear();
        //
    }

    void ClearDataForArcCycle() {
        //
        outArcCycle_VertexOnMesh.clear();
        outArcCycle_EdgeOnMesh.clear();
        //
    }

    //
    void FindEdgeConnectingTwoVertices(const int src_x, const int dist_y, const int eid_y, int &outEdge);

    int TriangleAdjacentToOneEdgesOneVertex(const int e0, const int opp_v) {
        int triangle_id = (*inMeshPtr).vecEdge[e0].AdjTri[0];
        if ((*inMeshPtr).vecTriangle[triangle_id].v0 != opp_v &&
            (*inMeshPtr).vecTriangle[triangle_id].v1 != opp_v &&
            (*inMeshPtr).vecTriangle[triangle_id].v2 != opp_v) {//
            triangle_id = (*inMeshPtr).vecEdge[e0].AdjTri[1];
        }
        return triangle_id;
    }

    int TriangleAdjacentToTwoEdges(const int e0, const int e1) {
        int triangle_id = (*inMeshPtr).vecEdge[e0].AdjTri[0];
        if ((*inMeshPtr).vecTriangle[triangle_id].e01 != e1 &&
            (*inMeshPtr).vecTriangle[triangle_id].e02 != e1 &&
            (*inMeshPtr).vecTriangle[triangle_id].e12 != e1) {//
            triangle_id = (*inMeshPtr).vecEdge[e0].AdjTri[1];
        }
        return triangle_id;
    }

    // actual operation functions
    void ComputeLevelCycle();

    //
    void ComputeArcCycle();

    //
    void WalkAlongArcWithInitialSide_InverseArcOrder(const int init_vertex_for_side,
                                                     int &current_active_vertex, int &refEdge, int &outTriagleId,
                                                     const std::vector<std::pair<int, int> > &arcPointType,
                                                     std::vector<int> &vecArc_vertex, std::vector<int> &vecArc_edge);

    void WalkAlongArcWithInitialSide_ArcOrder(const int init_vertex_for_side,
                                              int &current_active_vertex, int &refEdge, int &outTriagleId,
                                              const std::vector<std::pair<int, int> > &arcPointType,
                                              std::vector<int> &vecArc_vertex, std::vector<int> &vecArc_edge);

    //
    bool VerifyingEmbeddedPath(std::vector<int> &vertices, std::vector<int> &edges);

public:
    //
    _SimpleMesh *inMeshPtr;
    double *_scalarPtr;
    //
    std::vector<std::pair<int, int> > *_levelCyclePointTypePtr;
    std::vector<std::pair<int, int> > *_arcCyclePointTypePtr;
    //
    // change both as the vertical loop changes
    std::vector<int> outLevelCycle_VertexOnMesh;
    std::vector<int> outLevelCycle_EdgeOnMesh;
    /*
    The closed path on mesh is encoded in this way
    v0--e0--v1--e1--v2--e2--v3--....--vn--en--v0
    */
    // alos use it as the output for closed loop cycle
    std::vector<int> outArcCycle_VertexOnMesh;
    std::vector<int> outArcCycle_EdgeOnMesh;
    //
};


#endif //MapLoopsBackToMeshLevelSetAndArc.h
