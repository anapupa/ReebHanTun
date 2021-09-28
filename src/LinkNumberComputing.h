/*
(c) 2012 Fengtao Fan
*/
#ifndef MY_LINK_NUMBER_COMPUTING_H
#define MY_LINK_NUMBER_COMPUTING_H

#include "psbmReebGraphElements.h"
#include <RenderVector3.h>

class LinkNumberPairPolygon {
public:
    LinkNumberPairPolygon() {
        _link_number = -100000; // not computed yet
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

    LinkNumberPairPolygon(_Polygon &inLeftPlg, _Polygon &inRightPlg) {
        leftPlg = inLeftPlg;
        rightPlg = inRightPlg;
        _link_number = -100000;
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
    LinkNumberPairPolygon(const LinkNumberPairPolygon &rhs) {
        leftPlg = rhs.leftPlg;
        rightPlg = rhs.rightPlg;
        _link_number = rhs._link_number;
        _ptOnPlane = rhs._ptOnPlane;
        _planeNormal = rhs._planeNormal;
        _xAxis = rhs._xAxis;
        _yAxis = rhs._yAxis;
    }

    LinkNumberPairPolygon &operator=(const LinkNumberPairPolygon &rhs) {
        if (this == &rhs)
            return *this; // avoid self assignment
        leftPlg = rhs.leftPlg;
        rightPlg = rhs.rightPlg;
        _link_number = rhs._link_number;
        _ptOnPlane = rhs._ptOnPlane;
        _planeNormal = rhs._planeNormal;
        _xAxis = rhs._xAxis;
        _yAxis = rhs._yAxis;
    }

    // assign polygon
    void AssignPolygons(_Polygon &inLeftPlg,
                        _Polygon &inRightPlg) {
        // clean the memory
        leftPlg.Clear();
        rightPlg.Clear();
        _link_number = 0;
        //
        leftPlg = inLeftPlg;
        rightPlg = inRightPlg;
        return;
    }

/*****************
Pre-Processing the linking number computing in three steps
********************/
    void FindingProjectionPlaneWithoutDegeneratePointProjecting();

    void FindingProjectionPlaneWithoutDegeneratePointProjecting_2(int ___int___);

    void AvoidParallelOrthogonalToXYAxis();

    void AvoidOverlapping();

    //
    void FindingProjectionPlaneWithoutDegeneratePointProjecting_ANN();

    void OptimizeXYAxisMakeNoParaOrtho();

    void OptimizeNormalMakeNoOverlapping();

/******************************************/
    bool ComputeLinkNumber();

    bool ComputeLinkNumber_CGAL();

    void NewValueForLinkNumber(const double orientationFlag, const int heightFlag, const int incrementalValue);

    //
    bool ComputeLinkNumberForTwoPolygons(_Polygon &inLeftPlg,
                                         _Polygon &inRightPlg);

public:
    /*
    Assumption on the polygon : the first point = last point in the polygon
    */
    _Polygon leftPlg;
    _Polygon rightPlg;
    int _link_number;
    // define the project plane
    Vector3 _ptOnPlane;
    Vector3 _planeNormal;
    Vector3 _xAxis;
    Vector3 _yAxis;
};

void ProjectingToPlane(Vector3 &basePt, Vector3 &unitNormalDir, Vector3 &inPt, Vector3 &outPt);

#endif
