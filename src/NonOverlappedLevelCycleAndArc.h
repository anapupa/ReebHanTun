/*
(c) 2012 Fengtao Fan
*/
#ifndef _MY_NON_OVERLAPPED_LEVEL_CYCLE_AND_ARC_H_
#define _MY_NON_OVERLAPPED_LEVEL_CYCLE_AND_ARC_H_

#include "psbmReebGraphElements.h"
#include "IntersectionOnMesh.h"

#include <vector>

class NonOverlappedLevelCycleAndArc {
public:
    // constructor and assign operator
    NonOverlappedLevelCycleAndArc() {
        Null_all();
    }

    NonOverlappedLevelCycleAndArc(const NonOverlappedLevelCycleAndArc &rhs) {
        inMeshPtr = rhs.inMeshPtr;
        //
        _pVecSimplifiedArc = rhs._pVecSimplifiedArc;
        _arcOnMeshPtr = rhs._arcOnMeshPtr;
        _levelCyclePointTypePtr = rhs._levelCyclePointTypePtr;
        _arcPointTypePtr = rhs._arcPointTypePtr;
        //
        _basisLoopPtr = rhs._basisLoopPtr;
        //
        _levelCyclePtr = rhs._levelCyclePtr;
        //
        augmentedLevelCycles.Clear();
        augmentedLevelCycles = rhs.augmentedLevelCycles;
        //
        augmentedArc.Clear();
        augmentedArc = rhs.augmentedArc;
        //
        outLevelCycle.Clear();
        outLevelCycle = rhs.outLevelCycle;
        //
        outVerticalCycle.Clear();
        outVerticalCycle = rhs.outVerticalCycle;
    }

    NonOverlappedLevelCycleAndArc &operator=(const NonOverlappedLevelCycleAndArc &rhs) {
        if (this == &rhs)
            return *this;
        inMeshPtr = rhs.inMeshPtr;
        //
        _pVecSimplifiedArc = rhs._pVecSimplifiedArc;
        _arcOnMeshPtr = rhs._arcOnMeshPtr;
        _levelCyclePointTypePtr = rhs._levelCyclePointTypePtr;
        _arcPointTypePtr = rhs._arcPointTypePtr;
        //
        _basisLoopPtr = rhs._basisLoopPtr;
        //
        _levelCyclePtr = rhs._levelCyclePtr;
        //
        augmentedLevelCycles.Clear();
        augmentedLevelCycles = rhs.augmentedLevelCycles;
        //
        augmentedArc.Clear();
        augmentedArc = rhs.augmentedArc;
        //
        outLevelCycle.Clear();
        outLevelCycle = rhs.outLevelCycle;
        //
        outVerticalCycle.Clear();
        outVerticalCycle = rhs.outVerticalCycle;

        return *this;
    }

    // deconstructor and data clean functioin
    void Null_all() {
        inMeshPtr = NULL;
        //
        _pVecSimplifiedArc = NULL;
        _arcOnMeshPtr = NULL;
        _levelCyclePointTypePtr = NULL;
        _arcPointTypePtr = NULL;
        //
        _basisLoopPtr = NULL;
        //
        _levelCyclePtr = NULL;
        //
    }

    void Clear() {
        //
        outVerticalCycle.Clear();
        _basisLoopPtr = NULL;
        //
    }

    void ClearOutPut() {
        outLevelCycle.Clear();
        outVerticalCycle.Clear();
    }

    ~NonOverlappedLevelCycleAndArc() {
        Null_all();
        augmentedLevelCycles.Clear();
        augmentedArc.Clear();
        outLevelCycle.Clear();
        outVerticalCycle.Clear();
    }

    // real operations functions
    // Initial data
    // any pair of arc and level set needs them
    void InitMeshPtr(_SimpleMesh *inPtr) {
        inMeshPtr = inPtr;
        return;
    }

    void InitBool(const bool inbool) {
        move_common_pt_on_level_cycle = inbool;
    }

    void InitHeightDirection(Vector3 &direction) {
        heightDirection = direction;
    }

    void InitArcPolygons(std::vector<struct SimplifiedReebGraphArc> *pVecSimplifiedArc,
                         std::vector<_Polygon> *arcOnMesh) {
        _pVecSimplifiedArc = pVecSimplifiedArc;
        _arcOnMeshPtr = arcOnMesh;

    }

    //
    // particular pair of arc and levelset needs them
    void InitArcLevelPointType(
            std::vector<std::pair<int, int> > *levelCyclePointTypePtr,
            std::vector<std::pair<int, int> > *arcPointTypePtr) {
        _levelCyclePointTypePtr = levelCyclePointTypePtr;
        _arcPointTypePtr = arcPointTypePtr;
    }

    void InitLevelCycleAndArc(_Polygon *inLevelCyclePtr, const int arcId) {// pathArcOnMesh is initialized
        // clean memory first
        ClearDataForParticularLevelsetAndArc();
        //
        augmentedArc = (*_arcOnMeshPtr)[arcId]; // need to be consistent with ComputeLoopOnMesh
        _levelCyclePtr = inLevelCyclePtr;
        augmentedLevelCycles = *inLevelCyclePtr;
        _arcId = arcId;
        //
    }

    void ClearDataForParticularLevelsetAndArc() {
        levelDuplicate = false;
        arcDuplicate = false;
        _arcId = -1;
        move_common_pt_on_level_cycle = false;
        pointPositionInLevelCycleVector = -1; // the position in the point vector of each level cycles
        pointPositionInArcVector = -1; // same meaning as pointPositionInLevelCycleVector
        augmentedLevelCycles.Clear(); // keep the order of input cycles as the order means orientation
        augmentedArc.Clear(); // keep the order of input arc as order means orientation
        // change both as the vertical loop changes
        outLevelCycle.Clear();
        outVerticalCycle.Clear();
    }

    // particular path needs
    void InitVerticalPath(psbmReebPath *basisLoopPtr) {
        outLevelCycle.Clear();
        outVerticalCycle.Clear();
        _basisLoopPtr = basisLoopPtr;
    }

    void ResetVerticalPath(psbmReebPath *basisLoopPtr) {
        outLevelCycle.Clear();
        outVerticalCycle.Clear();
        _basisLoopPtr = basisLoopPtr;
    }

    void ResetBool(const bool inbool) {
        //if (move_common_pt_on_level_cycle != inbool)
        //{
        outLevelCycle.Clear();
        outVerticalCycle.Clear();
        move_common_pt_on_level_cycle = inbool;
        //}
    }

    void SetArcHeightRangeIntersectingLevelsetPtr(std::vector<bool> *inPtr) {
        _arcHeightRangeIntersectingLevelsetPtr = inPtr;
    }

    //
    void SetOwnVerticalType(const bool inType) {
        ownVerticalType = inType;
    }

    // actual operation functions
    void ComputeIntersectionBetweenLevelCycleAndArc();

    //
    void ComputeOutputCycles();

    //
    void WalkTrianglesAroundVertex(const int vid, Vector3 &outNormal);

    //
    void WalkEdgeForNormal(const int eid, Vector3 &outNormal);

    //
    void
    TranslatePoint(const double dir, const double scale, const Vector3 &curPt, const Vector3 &normal, Vector3 &outPt);

    //
    void ComputeLoopOnMeshPolygon();

    void ComputeLoopOnMeshPolygon_withArcHeightInfo();

public:
    //
    _SimpleMesh *inMeshPtr;
    // pointer to path infor
    psbmReebPath *_basisLoopPtr;
    //
    std::vector<struct SimplifiedReebGraphArc> *_pVecSimplifiedArc;
    std::vector<_Polygon> *_arcOnMeshPtr;
    std::vector<std::pair<int, int> > *_levelCyclePointTypePtr;
    std::vector<std::pair<int, int> > *_arcPointTypePtr;
    // input level cycle
    _Polygon *_levelCyclePtr;
    //
    //
    // own data
    bool levelDuplicate;
    bool arcDuplicate;
    int _arcId;
    bool move_common_pt_on_level_cycle;
    Vector3 heightDirection;
    Vector3 pointNormal; // the normals of
    int pointPositionInLevelCycleVector; // the position in the point vector of each level cycles
    int pointPositionInArcVector; // same meaning as pointPositionInLevelCycleVector
    _Polygon augmentedLevelCycles; // keep the order of input cycles as the order means orientation
    _Polygon augmentedArc; // keep the order of input arc as order means orientation
    // change both as the vertical loop changes
    _Polygon outLevelCycle;
    _Polygon outVerticalCycle;
    //
    double translateScale;
    bool ownVerticalType;
    //
    std::vector<bool> *_arcHeightRangeIntersectingLevelsetPtr;
};

#endif
