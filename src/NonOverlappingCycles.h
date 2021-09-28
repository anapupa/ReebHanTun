/*
(c) 2012 Fengtao Fan
*/
#ifndef _MY_NON_OVERLAPPING_CYCLES_H_
#define _MY_NON_OVERLAPPING_CYCLES_H_

#include <vector>
#include "SimpleMesh.h"
#include <RenderVector3.h>
#include "psbmReebGraph.h"

class NonOverlappingCycles {
    struct doubleIntLessThan {
        bool operator()(const std::pair<double, std::pair<std::pair<int, int>, int> > &lhs,
                        const std::pair<double, std::pair<std::pair<int, int>, int> > &rhs) {
            bool ret = false;
            if (lhs.first < rhs.first)
                ret = true;
            else if (lhs.first == rhs.first) {
                if (lhs.second < rhs.second)
                    ret = true;
            }
            return ret;
        }
    };

public:
    // constructor and assign operator
    NonOverlappingCycles() {
        Null_all();
    }

    NonOverlappingCycles(const NonOverlappingCycles &rhs) {
        inMeshPtr = rhs.inMeshPtr;
        //
        _pVecSimplifiedArc = rhs._pVecSimplifiedArc;
        _pathArcOnMeshPtr = rhs._pathArcOnMeshPtr;
        _offsetPathArcOnMeshPtr = rhs._offsetPathArcOnMeshPtr;
        //
        _basisLoopPtr = rhs._basisLoopPtr;
        _loopOnMeshPtr = rhs._loopOnMeshPtr;
        //
        _criticalNodeIdSetPtr = rhs._criticalNodeIdSetPtr;
        _vecReebNodePtr = rhs._vecReebNodePtr;
        //
        unitOutNormalForCritNode.clear();
        unitOutNormalForCritNode = rhs.unitOutNormalForCritNode;
        //
        outBasisLoop.vecPoints.clear();
        outBasisLoop = rhs.outBasisLoop;
        //
        outLoopOnMesh.vecPoints.clear();
        outLoopOnMesh = rhs.outLoopOnMesh;
    }

    NonOverlappingCycles &operator=(const NonOverlappingCycles &rhs) {
        if (this == &rhs)
            return *this;
        inMeshPtr = rhs.inMeshPtr;
        //
        _pVecSimplifiedArc = rhs._pVecSimplifiedArc;
        _pathArcOnMeshPtr = rhs._pathArcOnMeshPtr;
        _offsetPathArcOnMeshPtr = rhs._offsetPathArcOnMeshPtr;
        //
        _basisLoopPtr = rhs._basisLoopPtr;
        _loopOnMeshPtr = rhs._loopOnMeshPtr;
        //
        _criticalNodeIdSetPtr = rhs._criticalNodeIdSetPtr;
        _vecReebNodePtr = rhs._vecReebNodePtr;
        //
        unitOutNormalForCritNode.clear();
        unitOutNormalForCritNode = rhs.unitOutNormalForCritNode;
        //
        outBasisLoop.Clear();
        outBasisLoop = rhs.outBasisLoop;
        //
        outLoopOnMesh.Clear();
        outLoopOnMesh = rhs.outLoopOnMesh;
        return *this;
    }

    // deconstructor and data clean functioin
    void Null_all() {
        inMeshPtr = NULL;
        //
        _pVecSimplifiedArc = NULL;
        _pathArcOnMeshPtr = NULL;
        _offsetPathArcOnMeshPtr = NULL;
        //
        _basisLoopPtr = NULL;
        _loopOnMeshPtr = NULL;
        //
        _criticalNodeIdSetPtr = NULL;
        _vecReebNodePtr = NULL;
    }

    void Clear() {
        //
        _basisLoopPtr = NULL;
        _loopOnMeshPtr = NULL;
        //
    }

    ~NonOverlappingCycles() {
        Null_all();
        unitOutNormalForCritNode.clear();
        outBasisLoop.Clear();
        outLoopOnMesh.Clear();
    }

    // real operations functions
    // Initial data
    void InitMeshPtr(_SimpleMesh *inPtr) {
        inMeshPtr = inPtr;
        return;
    }

    void InitArcPolygons(std::vector<struct SimplifiedReebGraphArc> *pVecSimplifiedArc,
                         std::vector<_Polygon> *pathArcOnMesh,
                         std::vector<_Polygon> *offsetPathArcOnMesh,
                         std::vector<std::vector<std::pair<int, int> > > *pathArcOnMeshPointTypePtr,
                         std::vector<std::vector<std::pair<int, int> > > *offsetPathArcOnMeshPointTypePtr) {
        _pVecSimplifiedArc = pVecSimplifiedArc;
        _pathArcOnMeshPtr = pathArcOnMesh;
        _offsetPathArcOnMeshPtr = offsetPathArcOnMesh;
        //
        _pathArcOnMeshPointTypePtr = pathArcOnMeshPointTypePtr;
        _offsetPathArcOnMeshPointTypePtr = offsetPathArcOnMeshPointTypePtr;
    }

    void InitArcPolygons(std::vector<struct SimplifiedReebGraphArc> *pVecSimplifiedArc,
                         std::vector<_Polygon> *pathArcOnMesh,
                         std::vector<_Polygon> *offsetPathArcOnMesh
    ) {// for comparison coding
        _pVecSimplifiedArc = pVecSimplifiedArc;
        _pathArcOnMeshPtr = pathArcOnMesh;
        _offsetPathArcOnMeshPtr = offsetPathArcOnMesh;
        //
    }

    void InitPairPathes(psbmReebPath *basisLoopPtr,
                        psbmReebPath *loopOnMeshPtr) {
        _basisLoopPtr = basisLoopPtr;
        _loopOnMeshPtr = loopOnMeshPtr;
    }

    void ResetBasisLoop(psbmReebPath *basisLoopPtr) {
        outBasisLoop.Clear();
        _basisLoopPtr = basisLoopPtr;
    }

    void ResetLoopOnMesh(psbmReebPath *loopOnMeshPtr) {
        outLoopOnMesh.Clear();
        _loopOnMeshPtr = loopOnMeshPtr;
    }

    void SetEdgeIntersectArcsPtr(std::vector<std::vector<std::pair<int, int> > > *inPtr) {
        _edgeIntersectArcsPtr = inPtr;
    }

    void InitNormalInfo(std::set<int> *criticalNodeIdSetPtr,
                        std::vector<class psbmReebNode> *vecReebNodePtr) {
        _criticalNodeIdSetPtr = criticalNodeIdSetPtr;
        _vecReebNodePtr = vecReebNodePtr;
    }

    void SetBoolVectorPtr(std::vector<bool> *inLoopOnMeshPtr,
                          std::vector<bool> *inBasisLoopPtr) {
        _bVecLoopOnMeshPtr = inLoopOnMeshPtr;
        _bVecBasisLoopPtr = inBasisLoopPtr;
    }

    void SetHeightDiretion(Vector3 &in_vec) {
        heightDirection = in_vec;
    }

    //
    bool A_Pair_3D_Line_Intersection(const Vector3 *in_A_ptr,
                                     const Vector3 *in_B_ptr,
                                     const Vector3 *in_C_ptr,
                                     const Vector3 *in_D_ptr,
                                     Vector3 &resIntersection) {// preq: there exists an intersection between these two 3d lines
        /*
            L1 = P1 + a V1
            L2 = P2 + b V2
            //
            a V1 = (P2 - P1) + b V2
            //
            a (V1 X V2) = (P2 - P1) X V2
        */
        /*
        t(x_b - x_a) + x_a = s(x_d - x_c) + x_c
        t(y_b - y_a) + y_a = s(y_d - y_c) + y_c;
        common = (x_b - x_a) * (y_d - y_c) - (x_d - x_c) * (y_b - y_a)
        dividend_s =
        */
        Vector3 ab_dir = *in_B_ptr - *in_A_ptr;
        Vector3 cd_dir = *in_D_ptr - *in_C_ptr;
        //
        Vector3 CA = *in_C_ptr - *in_A_ptr;
        Vector3 ab_dir_cross_cd_dir = ab_dir ^ cd_dir;
        Vector3 pt_cross_cd_dir = CA ^ cd_dir;
        bool bIntersection = false;
        if (norm2(ab_dir_cross_cd_dir) == 0.0) {
            std::cout << "ERROR IN 3D LINES INTERSECTION COMPUTATION" << std::endl;
            exit(9);
        } else {
            double ab_ratio = -1.0;
            double max_up = pt_cross_cd_dir[0];
            double max_denom = ab_dir_cross_cd_dir[0];
            for (int i = 1; i < 3; i++) {
                if (abs(ab_dir_cross_cd_dir[i]) > abs(max_denom)) {
                    max_denom = ab_dir_cross_cd_dir[i];
                    max_up = pt_cross_cd_dir[i];
                }
            }
            ab_ratio = max_up / max_denom;
            if (ab_ratio >= 0.0 && ab_ratio <= 1.0) {
                bIntersection = true;
                resIntersection = (*in_A_ptr) + ab_ratio * ab_dir;
                //resIntersection = (*in_C_ptr) + s_cd * ((*in_D_ptr) - (*in_C_ptr));
            } else {
                std::cout << "ERROR IN 3D LINES INTERSECTION COMPUTATION in overlapping cycle" << ab_ratio << std::endl;
                exit(9);
            }
        }
        return bIntersection;
    }

    //
    void TranslatePoint_onFace(const double dir, const double distance, const Vector3 inPt,
                               const int fid, Vector3 &outPt);

    void WalkEdgeForNormal(const int eid, Vector3 &outNormal);

    void TranslatePoint_onEdge(const double dir, const double distance, const Vector3 inPt,
                               const int eid, Vector3 &outPt);

    void
    OrderPointsReferingToOneEndPoint(const Vector3 &refPt, const std::vector<Vector3> &nodes, std::vector<int> &orders);

    std::pair<int, double> smallestDisatnce(const Vector3 &refPt, const std::vector<Vector3> &nodes);

    double smallestDisatnce(const int refVer, const int eid, Vector3 &outPt);

    void NearByTwoEdges(const int pivot_v, const int eid, int &leftEdge, int &rightEdge);

    void PointBetweenTwoPoints(const double dist_to_src, const Vector3 &src, const Vector3 &dst, Vector3 &res);

    // report results
    void SetArcSharedStatus();

    //
    void ComputeNonOverlappingCycles();

    //
    void ComputeOutNormals();

    //
    //void WalkEdgeForNormal(const int eid, Vector3& outNormal);
    void WalkTrianglesAroundVertex(const int vid, Vector3 &outNormal);

    //
    void
    TranslatePoint(const double dir, const double scale, const Vector3 &curPt, const Vector3 &normal, Vector3 &outPt);

    void TranslatePoint(const double dir, const Vector3 &prePt, const Vector3 &curPt, const Vector3 &nexPt,
                        const Vector3 &normal, Vector3 &outPt);

    //
    void TranslatePoint(const int refVer, const double distance, const Vector3 opp_ver,
                        const int fid, Vector3 &outPt);

    void ComputeLoopOnMeshPolygon();

    void ComputeLoopOnMeshPolygonOnFly(std::vector<int> sharedArcs,
                                       std::vector<_Polygon> &sharedArcsNewPoints,
                                       std::vector<_Polygon> &offsetSharedArcsNewPoints);

    //
    void ComputeBasisLoopPolygon();

    void ComputeBasisLoopPolygon_1();

    void ComputeBasisLoopPolygonOnFly(std::vector<int> sharedArcs,
                                      std::vector<_Polygon> &sharedArcsNewPoints,
                                      std::vector<_Polygon> &offsetSharedArcsNewPoints);

    //
    void SharedArcs(std::vector<int> &sharedArcs);

    void ComputeNewPathesForSharedArcs(
            std::vector<int> sharedArcs,
            std::vector<_Polygon> &sharedArcsNewPoints,
            std::vector<_Polygon> &offsetSharedArcsNewPoints);

    //
    void CheckIntersectionOnTriangle(const double dir, const int critNode, const Vector3 &ptOnEdge,
                                     const int reachedEdge, const int refEdgeVertex, const int iterEdge,
                                     const int activeTriangleIdx,
                                     const std::pair<int, int> &selfPreEdgePtType, const Vector3 &selfPreEdgePt,
                                     const std::pair<int, int> &objEdgePtType, const Vector3 &objEdgePt,
                                     const std::pair<int, int> &objPreEdgePtType, const Vector3 &objPreEdgePt,
                                     std::vector<Vector3> &outPts);

    void WalkAroundVertexToMakeIntersectionAtFace(const double dir, const int critNode,
                                                  const std::pair<int, int> &selfEdgePtType, const Vector3 &selfEdgePt,
                                                  const std::pair<int, int> &selfPreEdgePtType,
                                                  const Vector3 &selfPreEdgePt,
                                                  const std::pair<int, int> &objEdgePtType, const Vector3 &objEdgePt,
                                                  const std::pair<int, int> &objPreEdgePtType,
                                                  const Vector3 &objPreEdgePt,
                                                  std::vector<Vector3> &outPts);

    void WalkAroundVertexToMakeIntersectionAtFace_edgeCase(const double dir, const int critNode,
                                                           const std::pair<int, int> &selfEdgePtType,
                                                           const Vector3 &selfEdgePt,
                                                           const std::pair<int, int> &selfPreEdgePtType,
                                                           const Vector3 &selfPreEdgePt,
                                                           const std::pair<int, int> &objEdgePtType,
                                                           const Vector3 &objEdgePt,
                                                           const std::pair<int, int> &objPreEdgePtType,
                                                           const Vector3 &objPreEdgePt,
                                                           std::vector<Vector3> &outPts);


public:
    //
    std::vector<char> arcSharedStatus; //0 not shared arc
    // 1 shared arc
    // 3 node_0 is shared
    // 4 node_1 is shared
    //
    std::vector<bool> pathArcTakenForLoopOnMesh;
    //
    std::vector<bool> *_bVecLoopOnMeshPtr;
    std::vector<bool> *_bVecBasisLoopPtr;
    //
    _SimpleMesh *inMeshPtr;
    // pointer to path infor
    psbmReebPath *_basisLoopPtr;
    psbmReebPath *_loopOnMeshPtr;
    //
    std::vector<struct SimplifiedReebGraphArc> *_pVecSimplifiedArc;
    std::vector<_Polygon> *_pathArcOnMeshPtr;
    std::vector<_Polygon> *_offsetPathArcOnMeshPtr;
    std::vector<std::vector<std::pair<int, int> > > *_pathArcOnMeshPointTypePtr;
    std::vector<std::vector<std::pair<int, int> > > *_offsetPathArcOnMeshPointTypePtr;
    //std::vector<_Polygon> *_offsetPathArcOnMeshPtr;
    // pointer for computing normals
    std::set<int> *_criticalNodeIdSetPtr;
    std::vector<class psbmReebNode> *_vecReebNodePtr;
    std::vector<std::vector<std::pair<int, int> > > *_edgeIntersectArcsPtr;
    //
    std::map<int, Vector3> unitOutNormalForCritNode;
    // output
    _Polygon outBasisLoop;
    _Polygon outLoopOnMesh;
    //
    Vector3 heightDirection;
};

#endif //_MY_NON_OVERLAPPING_CYCLES_H_
