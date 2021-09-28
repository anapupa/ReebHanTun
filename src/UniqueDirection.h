/*
(c) 2012 Fengtao Fan
*/
#ifndef MY_UNIQUE_DIRECTION_H
#define MY_UNIQUE_DIRECTION_H

#include <RenderVector3.h>

#include <vector>

class UniqueDirectionComputation {
public:
    UniqueDirectionComputation() {
        pointSetPtr = NULL;
        //
        _ptOnPlane[0] = 0.;
        _ptOnPlane[1] = 0.;
        _ptOnPlane[2] = 0.;
        _planeNormal[0] = 0.;
        _planeNormal[1] = 0.;
        _planeNormal[2] = 0.;
        _xAxis[0] = 0.;
        _xAxis[1] = 0.;
        _xAxis[2] = 0.;
        _yAxis[0] = 0.;
        _yAxis[1] = 0.;
        _yAxis[2] = 0.;
    }

    UniqueDirectionComputation(std::vector<Vector3> *&inPointSetPtr, std::vector<std::pair<int, int> > &inSegmentSet) {
        pointSetPtr = inPointSetPtr;
        segmentSet = inSegmentSet;

        _ptOnPlane[0] = 0.;
        _ptOnPlane[1] = 0.;
        _ptOnPlane[2] = 0.;
        _planeNormal[0] = 0.;
        _planeNormal[1] = 0.;
        _planeNormal[2] = 0.;
        _xAxis[0] = 0.;
        _xAxis[1] = 0.;
        _xAxis[2] = 0.;
        _yAxis[0] = 0.;
        _yAxis[1] = 0.;
        _yAxis[2] = 0.;
    }


// assignment operators
    UniqueDirectionComputation(const UniqueDirectionComputation &rhs) {
        //
        pointSetPtr = rhs.pointSetPtr;
        segmentSet = rhs.segmentSet;
        //
        _ptOnPlane = rhs._ptOnPlane;
        _planeNormal = rhs._planeNormal;
        _xAxis = rhs._xAxis;
        _yAxis = rhs._yAxis;
    }

    UniqueDirectionComputation &operator=(const UniqueDirectionComputation &rhs) {
        if (this == &rhs)
            return *this; // avoid self assignment
        //
        pointSetPtr = rhs.pointSetPtr;
        segmentSet = rhs.segmentSet;
        //
        _ptOnPlane = rhs._ptOnPlane;
        _planeNormal = rhs._planeNormal;
        _xAxis = rhs._xAxis;
        _yAxis = rhs._yAxis;
    }
    // assign polygon

    void AssignPointsAndSegments(std::vector<Vector3> *&inPtr,
                                 std::vector<std::pair<int, int> > &inSegmentSet) {
        pointSetPtr = inPtr;
        segmentSet = inSegmentSet;
    }

    //
    void SetIdealDirection(Vector3 &inDir) {
        idealDirection = inDir;
        return;
    }

    void AssignPointer(std::vector<Vector3> *&inPtr) {
        pointSetPtr = inPtr;
    }

    void AssignMeshNormal(std::vector<Vector3> *&inPtr) {
        meshNormalPtr = inPtr;
    }

    void CleanData() {
        {
            std::vector<std::pair<int, int> > tmp;
            segmentSet.swap(tmp);
        }
    }

/*****************
Pre-Processing the linking number computing in three steps
********************/
    void FindingProjectionPlaneWithoutDegeneratePointProjecting(int __int___);

    void FindingProjectionPlaneWithoutDegeneratePointProjecting_2(int __int___);

    void FindingProjectionPlaneWithoutDegeneratePointProjecting_MatchIdealDirection(int __int___);

    void ComputeUniqueDirection(Vector3 &resDirection);

    void ComputeUniqueDirection_MatchIdealDirection(Vector3 &resDirection);

/******************************************/

public:
    /*
    Assumption on the polygon : the first point = last point in the polygon
    */
    std::vector<Vector3> *pointSetPtr;
    std::vector<Vector3> *meshNormalPtr;
    std::vector<std::pair<int, int> > segmentSet;
    // define the project plane
    Vector3 _ptOnPlane;
    Vector3 _planeNormal;
    Vector3 _xAxis;
    Vector3 _yAxis;
    //
    Vector3 idealDirection;
};

void ProjectingToPlane(Vector3 &basePt, Vector3 &unitNormalDir, Vector3 &inPt, Vector3 &outPt);

#endif //MY_UNIQUE_DIRECTION_H
