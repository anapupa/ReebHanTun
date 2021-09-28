/*
(c) 2012 Fengtao Fan
*/
#include "LinkNumberComputing.h"
#include "SegmentIntersection.h"
#include "ANNSearch.h"
#include <vector>
#include <algorithm>

#define MY_LOCAL_PI 3.14159265

/*local functions used for LinkNumberComputing*/
char SegSegInt(double a[2], double b[2],
               double c[2], double d[2],
               double p[2]);

double angle_between_projected_segment_and_X_axis_LinkNumber_simple(Vector3 &segment) {// sgement is normalized
    return acos(segment[0]);

}

double angle_between_projected_segment_and_X_axis_LinkNumber(Vector3 &segment) {
    double angle = 0.0;
    if (segment[0] == 0.0) {
        if (segment[1] == 0.0) {
            angle = 0.0;
        } else {
            if (segment[1] > 0) {
                angle = MY_LOCAL_PI * 0.5;
            } else {
                angle = MY_LOCAL_PI * 1.5;
            }
        }
    } else {//
        if (segment[0] > 0) {// first and fourth quater
            if (segment[1] == 0.0) {
                angle = 0.0;
            } else {
                if (segment[1] > 0) {
                    angle = atan(segment[1] / segment[0]);
                } else {
                    angle = 2.0 * MY_LOCAL_PI + atan(segment[1] / segment[0]);
                }
            }

        } else {// (segment[0] < 0)
            if (segment[1] == 0.0) {
                angle = MY_LOCAL_PI;
            } else {
                angle = MY_LOCAL_PI + atan(segment[1] / segment[0]);
            }
        }
    }
    return angle;
}

double LowestProjectionDistance(Vector3 &unitAxisDir,
                                _Polygon &leftPlg,
                                _Polygon &rightPlg) {
    double lowestDist = 0.0;
    double currentDist = 0.0;
    if (leftPlg.vecPoints.empty() || rightPlg.vecPoints.empty()) {
        std::cout << "NO POINTS IN ONE PLYGON" << std::endl;
    }
    // calculate the first projection distance
    lowestDist = unitAxisDir * leftPlg.vecPoints[0];
    // process the left polygon
    for (unsigned int i = 1; i < leftPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        currentDist = unitAxisDir * leftPlg.vecPoints[i];
        if (lowestDist > currentDist)
            lowestDist = currentDist;
    }
    // process the right polygon
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        currentDist = unitAxisDir * rightPlg.vecPoints[i];
        if (lowestDist > currentDist)
            lowestDist = currentDist;
    }
    return lowestDist;
}

void ProjectingToPlane_LinkNumber(Vector3 &basePt, Vector3 &unitNormalDir, Vector3 &inPt, Vector3 &outPt) {
    Vector3 segDir = inPt - basePt;
    double projDist = (segDir * unitNormalDir);
    outPt = inPt - projDist * unitNormalDir;
    return;
}

void RotateAlongAxis_LinkNumber(double u, double v, double w,
                                double x, double y, double z,
                                double theta,
                                Vector3 &outRes) {
//rotation vector along the origin
//<u, v, w> is the rotation axis, a unit vector
// theta is the rotation angle
//<x, y, z> is the input vector,
    double uvwxyz = (u * x + v * y + w * z);
    double cos_a = cos(theta);
    double sin_a = sin(theta);
    outRes[0] = u * uvwxyz * (1 - cos_a) + x * cos_a + (-w * y + v * z) * sin_a;
    outRes[1] = v * uvwxyz * (1 - cos_a) + y * cos_a + (w * x - u * z) * sin_a;
    outRes[2] = w * uvwxyz * (1 - cos_a) + z * cos_a + (-v * x + u * y) * sin_a;
    return;
}

void angle_between_projected_segment_and_X_axis_LinkNumber(double u, double v, double w,
                                                           double x, double y, double z,
                                                           Vector3 &outRes) {
    double uvw = 1.0 / sqrt(u * u + v * v + w * w);
    double uv = 1.0 / sqrt(u * u + v * v);
    outRes[0] = ((u * w * x + w * v * y) * uv - (z) / uv) * uvw;
    outRes[1] = (u * y - v * x) * uv;
    outRes[2] = (u * x + v * y + w * z) * uvw;
}

void Rotation_to_make_uvw_become_z_axis(double u, double v, double w,
                                        double x, double y, double z,
                                        Vector3 &outRes) {
    double uvw = 1.0 / sqrt(u * u + v * v + w * w);
    double uv = 1.0 / sqrt(u * u + v * v);
    outRes[0] = ((u * w * x + w * v * y) * uv - (z) / uv) * uvw;
    outRes[1] = (u * y - v * x) * uv;
    outRes[2] = (u * x + v * y + w * z) * uvw;
}
/**************************************/
/**************************************/
/**************************************/
void LinkNumberPairPolygon::FindingProjectionPlaneWithoutDegeneratePointProjecting_ANN() {//
    const int nSegmentNumber = leftPlg.vecPoints.size() - 1 + rightPlg.vecPoints.size() - 1;
    Vector3 segment;
    vector<double> theta;
    vector<double> phi;
    // process the left polygon
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = leftPlg.vecPoints[i + 1] - leftPlg.vecPoints[i];
        if (segment[2] < 0.0) {
            for (int j = 0; j < 3; j++)
                segment[j] = -segment[j];
        }
        theta.push_back(angle_between_projected_segment_and_X_axis_LinkNumber(segment));
        //
        if (norm2(segment) == 0.0)
            phi.push_back(0.0);
        else {
            unitize(segment);
            phi.push_back(MY_LOCAL_PI * 0.5 - acos(segment[2]));
        }
    }
    // process the right polygon
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = rightPlg.vecPoints[i + 1] - rightPlg.vecPoints[i];
        if (segment[2] < 0.0) {
            for (int j = 0; j < 3; j++)
                segment[j] = -segment[j];
        }
        theta.push_back(angle_between_projected_segment_and_X_axis_LinkNumber(segment));
        //
        if (norm2(segment) == 0.0)
            phi.push_back(0.0);
        else {
            unitize(segment);
            phi.push_back(MY_LOCAL_PI * 0.5 - acos(segment[2]));
        }
    }
    //

    //search the largest distance from these various existing points in the angle parametrization space
    const double thetaRange = MY_LOCAL_PI * 2.0;
    const double phiRange = MY_LOCAL_PI * 0.5;
    const int sqrtSegNumber = int(sqrt(double(nSegmentNumber))) + 1;
    const int thetaSample = 8 * sqrtSegNumber;
    const int phiSample = sqrtSegNumber;

    std::pair<double, double> maxAnglePair;
    ANNSearch::GetMaxDistancePoint(theta, phi,
                                   thetaRange, phiRange,
                                   thetaSample, phiSample,
                                   maxAnglePair.first, maxAnglePair.second);

    //calculate the direction according to the theta and phi parameter
    _planeNormal[0] = cos(maxAnglePair.second) * cos(maxAnglePair.first);
    _planeNormal[1] = cos(maxAnglePair.second) * sin(maxAnglePair.first);
    _planeNormal[2] = sin(maxAnglePair.second);
    //unitize(_planeNormal);
    // calculate the lowest projection along the normal axis for each point in current XYZ-coordinate system
    _ptOnPlane = Vector3(0., 0., 0.);//_planeNormal * (LowestProjectionDistance(_planeNormal, leftPlg, rightPlg) - 2.0);
    //
    ProjectingToPlane(_ptOnPlane, _planeNormal, leftPlg.vecPoints[0], _xAxis);
    ProjectingToPlane(_ptOnPlane, _planeNormal, leftPlg.vecPoints[1], _yAxis);
    _xAxis -= _yAxis;
    unitize(_xAxis);
    // taking the cross product
    _yAxis = _planeNormal ^ _xAxis;
    unitize(_yAxis);
    //
    return;
}

//
void LinkNumberPairPolygon::OptimizeXYAxisMakeNoParaOrtho() {
//
    // calculating the intersection angle delta
    const int nSegmentNumber = leftPlg.vecPoints.size() - 1 + rightPlg.vecPoints.size() - 1;
    std::vector<double> delta_angle(nSegmentNumber);
    // process the left polygon
    Vector3 segment;
    double xCoord = 0.0;
    double yCoord = 0.0;
    double invPI = 1.0 / MY_LOCAL_PI;
    double halfPI = 0.5 * MY_LOCAL_PI;
    //
    // process the left polygon
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = leftPlg.vecPoints[i + 1] - leftPlg.vecPoints[i];
        xCoord = segment * _xAxis;
        yCoord = segment * _yAxis;
        if (xCoord == 0.0) {
            delta_angle[i] = 0.0;
        } else {
            delta_angle[i] = MY_LOCAL_PI + atan(yCoord / xCoord);
            delta_angle[i] = delta_angle[i] - ((int) (delta_angle[i] * 2.0 * invPI)) * halfPI;
        }
    }
    int prevSize = leftPlg.vecPoints.size() - 1;
    // process the right polygon
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = rightPlg.vecPoints[i + 1] - rightPlg.vecPoints[i];
        xCoord = segment * _xAxis;
        yCoord = segment * _yAxis;
        if (xCoord == 0.0) {
            delta_angle[i + prevSize] = 0.0;
        } else {
            delta_angle[i + prevSize] = MY_LOCAL_PI + atan(yCoord / xCoord);
            delta_angle[i + prevSize] =
                    delta_angle[i + prevSize] - ((int) (delta_angle[i + prevSize] * 2.0 * invPI)) * halfPI;
        }
    }
    //

    //
    delta_angle.push_back(0.0);
    delta_angle.push_back(halfPI);
    // sorting the angle
    std::sort(delta_angle.begin(), delta_angle.end());
    //scan the array to find the max-tempty interval
    double maxIntervalBegin = 0.0;
    double maxIntervalEnd = 0.0;
    for (unsigned int i = 1; i < delta_angle.size(); i++) {
        if (delta_angle[i] - delta_angle[i - 1] > maxIntervalEnd - maxIntervalBegin) {
            maxIntervalBegin = delta_angle[i - 1];
            maxIntervalEnd = delta_angle[i];
        }
    }
    //
    double rotAngle = (maxIntervalBegin + maxIntervalEnd) * 0.5;
    //
    RotateAlongAxis_LinkNumber(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                               _yAxis[0], _yAxis[1], _yAxis[2],
                               rotAngle,
                               _yAxis);
    //
    RotateAlongAxis_LinkNumber(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                               _xAxis[0], _xAxis[1], _xAxis[2],
                               rotAngle,
                               _xAxis);
    //
    return;
}

//
void LinkNumberPairPolygon::OptimizeNormalMakeNoOverlapping() {
//
    // calculating the intersection angle delta
    const int nSegmentNumber = leftPlg.vecPoints.size() - 1 + rightPlg.vecPoints.size() - 1;
    std::vector<double> delta_angle(nSegmentNumber);
    // process the left polygon
    Vector3 segment;
    double xCoord = 0.0;
    double zCoord = 0.0;
    double invPI = 1.0 / MY_LOCAL_PI;
    double halfPI = 0.5 * MY_LOCAL_PI;
    //
    // process the left polygon
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = leftPlg.vecPoints[i + 1] - leftPlg.vecPoints[i];
        xCoord = segment * _xAxis;
        zCoord = segment * _planeNormal;
        if (zCoord == 0.0) {
            delta_angle[i] = halfPI;
        } else {
            delta_angle[i] = atan((-xCoord) / zCoord);
            if (xCoord * zCoord > 0.0) {
                delta_angle[i] = halfPI - delta_angle[i];
            }
            //delta_angle[i] = MY_LOCAL_PI + atan((-xCoord) / zCoord);
            //delta_angle[i] = delta_angle[i] - ((int)(delta_angle[i] * invPI)) * MY_LOCAL_PI;
        }
    }
    // process the right polygon
    int prevSize = leftPlg.vecPoints.size() - 1;
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = rightPlg.vecPoints[i + 1] - rightPlg.vecPoints[i];
        xCoord = segment * _xAxis;
        zCoord = segment * _planeNormal;
        if (zCoord == 0.0) {
            delta_angle[i + prevSize] = halfPI;
        } else {
            delta_angle[i] = atan((-xCoord) / zCoord);
            if (xCoord * zCoord > 0.0) {
                delta_angle[i] = halfPI - delta_angle[i];
            }
            //delta_angle[i + prevSize] = MY_LOCAL_PI + atan((-xCoord) / zCoord);
            //delta_angle[i + prevSize] = delta_angle[i] - ((int)(delta_angle[i + prevSize] * invPI)) * MY_LOCAL_PI;
        }
    }
    //

    //
    delta_angle.push_back(0.0);
    delta_angle.push_back(MY_LOCAL_PI);
    // sorting the angle
    std::sort(delta_angle.begin(), delta_angle.end());
    //scan the array to find the max-tempty interval
    double maxIntervalBegin = 0.0;
    double maxIntervalEnd = 0.0;
    for (unsigned int i = 1; i < delta_angle.size(); i++) {
        if (delta_angle[i] - delta_angle[i - 1] > maxIntervalEnd - maxIntervalBegin) {
            maxIntervalBegin = delta_angle[i - 1];
            maxIntervalEnd = delta_angle[i];
        }
    }
    //
    double rotAngle = (maxIntervalBegin + maxIntervalEnd) * 0.5;
    rotAngle = rotAngle - halfPI;
    //
    RotateAlongAxis_LinkNumber(_yAxis[0], _yAxis[1], _yAxis[2],
                               _planeNormal[0], _planeNormal[1], _planeNormal[2],
                               rotAngle,
                               _planeNormal);
//
    RotateAlongAxis_LinkNumber(_yAxis[0], _yAxis[1], _yAxis[2],
                               _xAxis[0], _xAxis[1], _xAxis[2],
                               rotAngle,
                               _xAxis);
//
    // rotate all points
    //_planeNormal[0] = 0.0; _planeNormal[1] = 0.0; _planeNormal[2] = 1.0;
    //_xAxis[0] = 1.0; _xAxis[1] = 1.0; _xAxis[2] = 0.0;
    //_yAxis[0] = -1.0; _yAxis[1] = 1.0; _yAxis[2] = 0.0;
    Vector3 tempPt;
    for (unsigned int i = 0;
         i < leftPlg.vecPoints.size() - 1; i++) {// rotate the points to make the new direction become Z-axis
        tempPt[0] = leftPlg.vecPoints[i] * _xAxis;
        tempPt[1] = leftPlg.vecPoints[i] * _yAxis;
        tempPt[2] = leftPlg.vecPoints[i] * _planeNormal;
        //
        leftPlg.vecPoints[i] = tempPt;
    }
    leftPlg.vecPoints[leftPlg.vecPoints.size() - 1] = leftPlg.vecPoints[0];
    // process the right polygon
    for (unsigned int i = 0;
         i < rightPlg.vecPoints.size() - 1; i++) {// rotate the points to make the new direction become Z-axis
        tempPt[0] = rightPlg.vecPoints[i] * _xAxis;
        tempPt[1] = rightPlg.vecPoints[i] * _yAxis;
        tempPt[2] = rightPlg.vecPoints[i] * _planeNormal;
        //
        rightPlg.vecPoints[i] = tempPt;
    }
    rightPlg.vecPoints[rightPlg.vecPoints.size() - 1] = rightPlg.vecPoints[0];
    return;
}

/*****************************/
void LinkNumberPairPolygon::FindingProjectionPlaneWithoutDegeneratePointProjecting_2(int ___int___) {//
    Vector3 segment;
    std::set<float> theta;
    std::set<float> phi;
    const int nSegmentNumber = leftPlg.vecPoints.size() - 1 + rightPlg.vecPoints.size() - 1;

    // process the left polygon
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = leftPlg.vecPoints[i + 1] - leftPlg.vecPoints[i];
        unitize(segment);
        //
        theta.insert(angle_between_projected_segment_and_X_axis_LinkNumber_simple(segment));
        //
        {
            if (segment[2] > 0)
                phi.insert(float(acos(segment[2])));
            else {
                phi.insert(float(acos(segment[2]) - MY_LOCAL_PI * 0.5));
            }
        }
    }
    // process the right polygon
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = rightPlg.vecPoints[i + 1] - rightPlg.vecPoints[i];
        unitize(segment);
        //

        {
            if (segment[2] > 0)
                phi.insert(float(acos(segment[2])));
            else {
                phi.insert(float(acos(segment[2]) - MY_LOCAL_PI * 0.5));
            }
        }
    }
    std::pair<float, float> maxAnglePair;
    float max_var_val = 0.0;
    float max_cur_val = 0.0;
    float max_pre_val = 0.0;
    float pre_val = 0.0;
    float cur_val = 0.0;
    std::set<float>::iterator sIter;
    sIter = theta.begin();
    pre_val = *sIter;
    for (sIter++; sIter != theta.end(); sIter++) {
        cur_val = *sIter;
        if (cur_val - pre_val > max_var_val) {
            max_var_val = cur_val - pre_val;
            max_cur_val = cur_val;
            max_pre_val = pre_val;
        }
        pre_val = cur_val;
    }
    //
    maxAnglePair.first = 0.5 * (max_cur_val + max_pre_val);
    //
    max_var_val = 0.0;
    sIter = phi.begin();
    pre_val = *sIter;
    for (sIter++; sIter != phi.end(); sIter++) {
        cur_val = *sIter;
        if (cur_val - pre_val > max_var_val) {
            max_var_val = cur_val - pre_val;
            max_cur_val = cur_val;
            max_pre_val = pre_val;
        }
        pre_val = cur_val;
    }
    maxAnglePair.second = 0.5 * (max_cur_val + max_pre_val);
    //search the largest distance from these various existing points in the angle parametrization space

    //std::cout << "max angle : " << maxAnglePair.first << " " << maxAnglePair.second << std::endl;
    //calculate the direction according to the theta and phi parameter
    _planeNormal[0] = cos(maxAnglePair.first) * sin(maxAnglePair.second);
    _planeNormal[1] = sin(maxAnglePair.first) * sin(maxAnglePair.second);
    _planeNormal[2] = cos(maxAnglePair.second);
    //unitize(_planeNormal);
    // calculate the lowest projection along the normal axis for each point in current XYZ-coordinate system
    _ptOnPlane = Vector3(0., 0., 0.);//_planeNormal * (LowestProjectionDistance(_planeNormal, leftPlg, rightPlg) - 2.0);
    //
    //
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        Rotation_to_make_uvw_become_z_axis(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                           leftPlg.vecPoints[i][0], leftPlg.vecPoints[i][1], leftPlg.vecPoints[i][2],
                                           leftPlg.vecPoints[i]);
    }
    leftPlg.vecPoints.back() = leftPlg.vecPoints.front();
    // process the right polygon
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        Rotation_to_make_uvw_become_z_axis(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                           rightPlg.vecPoints[i][0], rightPlg.vecPoints[i][1], rightPlg.vecPoints[i][2],
                                           rightPlg.vecPoints[i]);
    }
    rightPlg.vecPoints.back() = rightPlg.vecPoints.front();
    //
    /*Rotation_to_make_uvw_become_z_axis(_planeNormal[0], _planeNormal[1], _planeNormal[2],
        _planeNormal[0], _planeNormal[1], _planeNormal[2],
        _planeNormal);
    std::cout << _planeNormal[0] << " " << _planeNormal[1] << " " << _planeNormal[2] << " "<< std::endl;*/
    _planeNormal[0] = 0.;
    _planeNormal[1] = 0.;
    _planeNormal[2] = 1.;
    _xAxis = Vector3(1., 0., 0.);
    _yAxis = Vector3(0., 1., 0.);

    //
    return;
}

void LinkNumberPairPolygon::FindingProjectionPlaneWithoutDegeneratePointProjecting() {//
    const int nSegmentNumber = leftPlg.vecPoints.size() - 1 + rightPlg.vecPoints.size() - 1;
    Vector3 segment;
    vector<double> theta;
    vector<double> phi;
    // process the left polygon
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = leftPlg.vecPoints[i + 1] - leftPlg.vecPoints[i];
        theta.push_back(angle_between_projected_segment_and_X_axis_LinkNumber(segment));
        //
        if (norm2(segment) == 0.0)
            phi.push_back(0.0);
        else {
            unitize(segment);
            phi.push_back(MY_LOCAL_PI * 0.5 - std::min(acos(segment * Vector3(0.0, 0.0, 1.0)),
                                                       acos(segment * Vector3(0.0, 0.0, -1.0))));
        }
    }
    // process the right polygon
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = rightPlg.vecPoints[i + 1] - rightPlg.vecPoints[i];
        theta.push_back(angle_between_projected_segment_and_X_axis_LinkNumber(segment));
        //
        if (norm2(segment) == 0.0)
            phi.push_back(0.0);
        else {
            unitize(segment);
            phi.push_back(MY_LOCAL_PI * 0.5 - std::min(acos(segment * Vector3(0.0, 0.0, 1.0)),
                                                       acos(segment * Vector3(0.0, 0.0, -1.0))));
        }
    }
    //search the largest distance from these various existing points in the angle parametrization space
    const int sqrtSegNumber = int(sqrt(double(nSegmentNumber))) + 1;
    std::cout << "sqrtSegNumber " << sqrtSegNumber << std::endl;
    //
    double delta_angle = MY_LOCAL_PI / (4.0 * sqrtSegNumber);
    double max_angle_pair_dist = 0.0;
    std::pair<double, double> maxAnglePair;
    std::pair<double, double> anglePair;
    anglePair.first = 0.0;
    anglePair.second = 0.0;
    //
    maxAnglePair = anglePair;
    //
    for (int i = 0; i < 8 * sqrtSegNumber; i++) {
        anglePair.first += delta_angle;
        anglePair.second = 0.0;
        for (int j = 0; j < 2 * sqrtSegNumber; j++) {
            anglePair.second += delta_angle;
            //
            double max_dist_to_current_angle_pair = 0.0;
            double current_dist = 0.0;
            for (int k = 0; k < nSegmentNumber; k++) {
                current_dist = std::min((anglePair.first - theta[k]) * (anglePair.first - theta[k]) +
                                        (anglePair.second - phi[k]) * (anglePair.second - phi[k]),
                                        (anglePair.first - theta[k] - 2 * MY_LOCAL_PI) *
                                        (anglePair.first - theta[k] - 2 * MY_LOCAL_PI) +
                                        (anglePair.second - phi[k]) * (anglePair.second - phi[k]));

                if (current_dist > max_dist_to_current_angle_pair) {
                    max_dist_to_current_angle_pair = current_dist;
                }
            }
            //
            if (max_dist_to_current_angle_pair > max_angle_pair_dist) {
                maxAnglePair = anglePair;
                max_angle_pair_dist = max_dist_to_current_angle_pair;
            }
        }
    }//
    //calculate the direction according to the theta and phi parameter
    _planeNormal[0] = cos(maxAnglePair.first);
    _planeNormal[1] = sin(maxAnglePair.first);
    _planeNormal[2] = cos(maxAnglePair.second);
    unitize(_planeNormal);
    // calculate the lowest projection along the normal axis for each point in current XYZ-coordinate system
    _ptOnPlane = Vector3(0., 0., 0.);//_planeNormal * (LowestProjectionDistance(_planeNormal, leftPlg, rightPlg) - 2.0);
    //
    //
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        Rotation_to_make_uvw_become_z_axis(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                           leftPlg.vecPoints[i][0], leftPlg.vecPoints[i][1], leftPlg.vecPoints[i][2],
                                           leftPlg.vecPoints[i]);
    }
    leftPlg.vecPoints.back() = leftPlg.vecPoints.front();
    // process the right polygon
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        Rotation_to_make_uvw_become_z_axis(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                           rightPlg.vecPoints[i][0], rightPlg.vecPoints[i][1], rightPlg.vecPoints[i][2],
                                           rightPlg.vecPoints[i]);
    }
    rightPlg.vecPoints.back() = rightPlg.vecPoints.front();
    //
    Rotation_to_make_uvw_become_z_axis(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                       _planeNormal[0], _planeNormal[1], _planeNormal[2],
                                       _planeNormal);
    std::cout << _planeNormal[0] << " " << _planeNormal[1] << " " << _planeNormal[2] << " " << std::endl;
    _planeNormal[0] = 0.;
    _planeNormal[1] = 0.;
    _planeNormal[2] = 1.;
    _xAxis = Vector3(1., 0., 0.);
    _yAxis = Vector3(0., 1., 0.);
    return;
}

void LinkNumberPairPolygon::AvoidParallelOrthogonalToXYAxis() {
    // preq: all points are rotated such that the plane normal is Z-axis
    //		the projection is only to take X, Y coordinate

    //specify the XY axis in the projection plane first
    if (leftPlg.vecPoints.size() < 2) {
        std::cout << "NO ENOUGH POINTS" << std::endl;
        return;
    }
    _xAxis = Vector3(1., 0., 0.);
    _yAxis = Vector3(0., 1., 0.);
    //

    // computing the angle between the projected segment and its closest axis
    Vector3 orgPt, dstPt;
    std::vector<double> min_angle_to_axis;
    // process the left polygon
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {
        dstPt = leftPlg.vecPoints[i] - leftPlg.vecPoints[i + 1];
        dstPt[2] = 0.0;
        //
        unitize(dstPt);
        //
        double cos_angle = acos(dstPt[0]);
        //
        min_angle_to_axis.push_back(cos_angle);
        //
        if (cos_angle > MY_LOCAL_PI * 0.5) {
            min_angle_to_axis.push_back(cos_angle - MY_LOCAL_PI * 0.5);

        } else {
            min_angle_to_axis.push_back(cos_angle + MY_LOCAL_PI * 0.5);
        }
    }
    // process the right polygon
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        dstPt = rightPlg.vecPoints[i] - rightPlg.vecPoints[i + 1];
        dstPt[2] = 0.0;
        //
        unitize(dstPt);
        //
        double cos_angle = acos(dstPt * _xAxis);
        //
        min_angle_to_axis.push_back(cos_angle);
        //
        if (cos_angle > MY_LOCAL_PI * 0.5) {
            min_angle_to_axis.push_back(cos_angle - MY_LOCAL_PI * 0.5);

        } else {
            min_angle_to_axis.push_back(cos_angle + MY_LOCAL_PI * 0.5);
        }
    }
    //
    min_angle_to_axis.push_back(0.0);
    //min_angle_to_axis.push_back( MY_LOCAL_PI * 0.5);
    min_angle_to_axis.push_back(MY_LOCAL_PI);
    // sorting the angle
    std::sort(min_angle_to_axis.begin(), min_angle_to_axis.end());
    //scan the array to find the max-tempty interval
    double maxIntervalBegin = 0.0;
    double maxIntervalEnd = 0.0;
    for (unsigned int i = 1; i < min_angle_to_axis.size(); i++) {
        if (min_angle_to_axis[i] - min_angle_to_axis[i - 1] > maxIntervalEnd - maxIntervalBegin) {
            maxIntervalBegin = min_angle_to_axis[i - 1];
            maxIntervalEnd = min_angle_to_axis[i];
        }
    }
    //

    //
    double rotAngle = -(maxIntervalBegin + maxIntervalEnd) * 0.5;
    //
    //std::cout << "rotAngle " << rotAngle << std::endl;
/*
	//specify the XY axis in the projection plane first
	if (leftPlg.vecPoints.size() < 2)
	{
		std::cout << "NO ENOUGH POINTS" << std::endl;
		return;
	}
	ProjectingToPlane_LinkNumber(_ptOnPlane, _planeNormal, leftPlg.vecPoints[1], _xAxis);
	ProjectingToPlane_LinkNumber(_ptOnPlane, _planeNormal, leftPlg.vecPoints[0], _yAxis);
	_xAxis -= _yAxis;
	unitize(_xAxis);
	// taking the cross product
	_yAxis = _xAxis^_planeNormal;
	//

	// computing the angle between the projected segment and its closest axis
	Vector3 orgPt, dstPt;
	std::vector<double> min_angle_to_axis;
	// process the left polygon
	for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++)
	{
		ProjectingToPlane_LinkNumber(_ptOnPlane, _planeNormal,
						leftPlg.vecPoints[i], orgPt);
		ProjectingToPlane_LinkNumber(_ptOnPlane, _planeNormal,
						leftPlg.vecPoints[i+1], dstPt);
		dstPt -= orgPt;
		//
		unitize(dstPt);
		//
		double cos_angle = acos(dstPt * _xAxis);
		min_angle_to_axis.push_back(cos_angle);
		if (cos_angle > MY_LOCAL_PI * 0.5)
		{
			min_angle_to_axis.push_back(cos_angle - MY_LOCAL_PI * 0.5);

		}
		else
		{
			min_angle_to_axis.push_back(cos_angle - MY_LOCAL_PI * 0.5);
		}

		//if (abs(dstPt[0]) < 1e-10 || abs(dstPt[1]) < 1e-10)
		//{
		//	min_angle_to_axis.push_back(0.0);
		//}
		//else
		//{
		//	double initAngle = atan(dstPt[1] / dstPt[0]);
		//	initAngle += MY_LOCAL_PI;
		//	if (initAngle >= MY_LOCAL_PI)
		//		initAngle -= MY_LOCAL_PI;
		//	else
		//	{
		//		if (initAngle >= MY_LOCAL_PI * 0.5)
		//		{
		//			initAngle -= 0.5 * MY_LOCAL_PI;
		//		}
		//	}
		//	min_angle_to_axis.push_back(initAngle);
		//}
		//
		//min_angle_to_axis.push_back(std::max( std::max(dstPt * _xAxis, dstPt * _yAxis),
		//										std::max(dstPt * (-_xAxis), dstPt * (-_yAxis)) )
		//							);

	}
	// process the right polygon
	for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++)
	{// project the segment onto the XY-plane
		ProjectingToPlane_LinkNumber(_ptOnPlane, _planeNormal,
						rightPlg.vecPoints[i], orgPt);
		ProjectingToPlane_LinkNumber(_ptOnPlane, _planeNormal,
						rightPlg.vecPoints[i+1], dstPt);
		dstPt -= orgPt;
		//
		unitize(dstPt);
		//

		double cos_angle = acos(dstPt * _xAxis);
		min_angle_to_axis.push_back(cos_angle);
		if (cos_angle > MY_LOCAL_PI * 0.5)
		{
			min_angle_to_axis.push_back(cos_angle - MY_LOCAL_PI * 0.5);

		}
		else
		{
			min_angle_to_axis.push_back(cos_angle - MY_LOCAL_PI * 0.5);
		}

		//if (abs(dstPt[0]) < 1e-10 || abs(dstPt[1]) < 1e-10)
		//{
		//	min_angle_to_axis.push_back(0.0);
		//}
		//else
		//{
		//	double initAngle = atan(dstPt[1] / dstPt[0]);
		//	initAngle += MY_LOCAL_PI;
		//	if (initAngle >= MY_LOCAL_PI)
		//		initAngle -= MY_LOCAL_PI;
		//	else
		//	{
		//		if (initAngle >= MY_LOCAL_PI * 0.5)
		//		{
		//			initAngle -= 0.5 * MY_LOCAL_PI;
		//		}
		//	}
		//	min_angle_to_axis.push_back(initAngle);
		//}
		//min_angle_to_axis.push_back(std::max( std::max(dstPt * _xAxis, dstPt * _yAxis),
		//										std::max(dstPt * (-_xAxis), dstPt * (-_yAxis)) )
		//							);
	}
	//
	min_angle_to_axis.push_back(0.0);
	min_angle_to_axis.push_back(MY_LOCAL_PI);
	// sorting the angle
	std::sort(min_angle_to_axis.begin(), min_angle_to_axis.end());
	//scan the array to find the max-tempty interval
	double maxIntervalBegin = 0.0;
	double maxIntervalEnd = 0.0;
	for (unsigned int i = 1; i < min_angle_to_axis.size(); i++)
	{
		if (min_angle_to_axis[i] - min_angle_to_axis[i-1] > maxIntervalEnd - maxIntervalBegin)
		{
			maxIntervalBegin = min_angle_to_axis[i - 1];
			maxIntervalEnd   = min_angle_to_axis[i];
		}
	}
	//
	double rotAngle = (maxIntervalBegin + maxIntervalEnd) * 0.5;
	std::cout << "rotAngle " << rotAngle << std::endl;
*/
// rotate all points
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {
        RotateAlongAxis_LinkNumber(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                   leftPlg.vecPoints[i][0], leftPlg.vecPoints[i][1], leftPlg.vecPoints[i][2],
                                   rotAngle,
                                   leftPlg.vecPoints[i]);

    }
    leftPlg.vecPoints[leftPlg.vecPoints.size() - 1] = leftPlg.vecPoints[0];
    // process the right polygon
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        RotateAlongAxis_LinkNumber(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                   rightPlg.vecPoints[i][0], rightPlg.vecPoints[i][1], rightPlg.vecPoints[i][2],
                                   rotAngle,
                                   rightPlg.vecPoints[i]);
    }
    rightPlg.vecPoints[rightPlg.vecPoints.size() - 1] = rightPlg.vecPoints[0];
    //std::cout << _xAxis * _yAxis << std::endl;
/**/
    return;
}

void LinkNumberPairPolygon::AvoidOverlapping() {
//
    // calculating the intersection angle delta
    std::vector<double> delta_angle;
    // process the left polygon
    Vector3 segment;
    Vector3 intersectLineVec;
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {
        segment = leftPlg.vecPoints[i + 1] - leftPlg.vecPoints[i];
        unitize(segment);
        if (abs(segment[2]) < 1e-10) {
            delta_angle.push_back(0.0);
        } else {
            delta_angle.push_back(atan(-segment[0] / segment[2]));
        }
    }
    // process the right polygon
    for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
        segment = rightPlg.vecPoints[i + 1] - rightPlg.vecPoints[i];
        unitize(segment);
        if (abs(segment[2]) < 1e-10) {
            delta_angle.push_back(0.0);
        } else {
            delta_angle.push_back(atan(-segment[0] / segment[2]));
        }
    }
    //
    //
    delta_angle.push_back(MY_LOCAL_PI * 0.5);
    delta_angle.push_back(MY_LOCAL_PI * (-0.5));
    // sorting the angle
    std::sort(delta_angle.begin(), delta_angle.end());
    //scan the array to find the max-tempty interval
    double maxIntervalBegin = 0.0;
    double maxIntervalEnd = 0.0;
    for (unsigned int i = 1; i < delta_angle.size(); i++) {
        if (delta_angle[i] - delta_angle[i - 1] > maxIntervalEnd - maxIntervalBegin) {
            maxIntervalBegin = delta_angle[i - 1];
            maxIntervalEnd = delta_angle[i];
        }
    }
    //
    double rotAngle = (maxIntervalBegin + maxIntervalEnd) * 0.5;
    //
    if (rotAngle < 0.0) {
        rotAngle = -rotAngle;
    } else {//
        rotAngle = 0.5 * MY_LOCAL_PI - rotAngle;
        rotAngle = -rotAngle;
    }
    // rotate all points
    for (unsigned int i = 0;
         i < leftPlg.vecPoints.size() - 1; i++) {// rotate the points to make the new direction become Z-axis
        RotateAlongAxis_LinkNumber(_yAxis[0], _yAxis[1], _yAxis[2],
                                   leftPlg.vecPoints[i][0], leftPlg.vecPoints[i][1], leftPlg.vecPoints[i][2],
                                   rotAngle,
                                   leftPlg.vecPoints[i]);
    }
    leftPlg.vecPoints[leftPlg.vecPoints.size() - 1] = leftPlg.vecPoints[0];
    // process the right polygon
    for (unsigned int i = 0;
         i < rightPlg.vecPoints.size() - 1; i++) {// rotate the points to make the new direction become Z-axis
        RotateAlongAxis_LinkNumber(_yAxis[0], _yAxis[1], _yAxis[2],
                                   rightPlg.vecPoints[i][0], rightPlg.vecPoints[i][1], rightPlg.vecPoints[i][2],
                                   rotAngle,
                                   rightPlg.vecPoints[i]);
    }
    rightPlg.vecPoints[rightPlg.vecPoints.size() - 1] = rightPlg.vecPoints[0];
    //
/*
	// calculating the intersection angle delta
	std::vector<double> delta_angle;
	// process the left polygon
	Vector3 segment;
	Vector3 intersectLineVec; 
	for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++)
	{
		segment = leftPlg.vecPoints[i + 1] - leftPlg.vecPoints[i];
		unitize(segment);
		//
		if (abs(segment[2]) < 1e-10)
		{
			delta_angle.push_back(0.0);
		}
		else
		{
			delta_angle.push_back(atan(-segment[0]/segment[2]));

			//double denominator = segment[0] * _xAxis[2] - segment[2] * _xAxis[0];
			//denominator = 1.0 / denominator;
			//intersectLineVec[0] = denominator * (segment[1] * _xAxis[2] - segment[2] * _xAxis[1]); 
			//intersectLineVec[2] = denominator * (segment[1] * _xAxis[0] - segment[0] * _xAxis[1]);
			//intersectLineVec[1] = -1;
		}
		//
		//unitize(intersectLineVec);
		//
		//delta_angle.push_back(acos(_yAxis * intersectLineVec));
		//delta_angle.push_back(acos(_yAxis * (-intersectLineVec)));
	}
	// process the right polygon
	for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++)
	{// project the segment onto the XY-plane
		segment = rightPlg.vecPoints[i + 1] - rightPlg.vecPoints[i];
		unitize(segment);
		// set up the intersection directions
		if (abs(segment[2]) < 1e-10)
		{
			delta_angle.push_back(0.0);
		}
		else
		{
			delta_angle.push_back(atan(-segment[0]/segment[2]));

			//double denominator = segment[0] * _xAxis[2] - segment[2] * _xAxis[0];
			//denominator = 1.0 / denominator;
			//intersectLineVec[0] = denominator * (segment[1] * _xAxis[2] - segment[2] * _xAxis[1]); 
			//intersectLineVec[2] = denominator * (segment[1] * _xAxis[0] - segment[0] * _xAxis[1]);
			//intersectLineVec[1] = -1;
		}
		//double denominator = segment[0] * _xAxis[2] - segment[2] * _xAxis[0];
		//denominator = 1.0 / denominator;
		//intersectLineVec[0] = denominator * (segment[1] * _xAxis[2] - segment[2] * _xAxis[1]); 
		//intersectLineVec[2] = denominator * (segment[1] * _xAxis[0] - segment[0] * _xAxis[1]);
		//intersectLineVec[1] = -1;
		////
		//unitize(intersectLineVec);
		////
		//delta_angle.push_back(acos(_yAxis * intersectLineVec));
		//delta_angle.push_back(acos(_yAxis * (-intersectLineVec)));
	}
	//
	//
	//delta_angle.push_back(0.0);
	//delta_angle.push_back(MY_LOCAL_PI);
	delta_angle.push_back(MY_LOCAL_PI * (-0.5));
	delta_angle.push_back(MY_LOCAL_PI * 0.5);
	// sorting the angle
	std::sort(delta_angle.begin(), delta_angle.end());
	//scan the array to find the max-tempty interval
	double maxIntervalBegin = 0.0;
	double maxIntervalEnd = 0.0;
	for (unsigned int i = 1; i < delta_angle.size(); i++)
	{
		if (delta_angle[i] - delta_angle[i-1] > maxIntervalEnd - maxIntervalBegin)
		{
			maxIntervalBegin = delta_angle[i - 1];
			maxIntervalEnd   = delta_angle[i];
		}
	}
	//
	double rotAngle = (maxIntervalBegin + maxIntervalEnd) * 0.5;
	//
	rotAngle = 0.5 * MY_LOCAL_PI - rotAngle ;
// rotate all points
	for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++)
	{
		RotateAlongAxis_LinkNumber(_yAxis[0], _yAxis[1], _yAxis[2],
						leftPlg.vecPoints[i][0], leftPlg.vecPoints[i][1], leftPlg.vecPoints[i][2],
						rotAngle,
						leftPlg.vecPoints[i]);
		angle_between_projected_segment_and_X_axis_LinkNumber(_planeNormal[0], _planeNormal[1], _planeNormal[2],
						leftPlg.vecPoints[i][0], leftPlg.vecPoints[i][1], leftPlg.vecPoints[i][2],
						leftPlg.vecPoints[i]);

	}
	leftPlg.vecPoints[leftPlg.vecPoints.size() - 1] = leftPlg.vecPoints[0];
	// process the right polygon
	for (unsigned int i = 0; i < rightPlg.vecPoints.size() - 1; i++)
	{// project the segment onto the XY-plane
		RotateAlongAxis_LinkNumber(_yAxis[0], _yAxis[1], _yAxis[2],
						rightPlg.vecPoints[i][0], rightPlg.vecPoints[i][1], rightPlg.vecPoints[i][2],
						rotAngle,
						rightPlg.vecPoints[i]);
		angle_between_projected_segment_and_X_axis_LinkNumber (_planeNormal[0], _planeNormal[1], _planeNormal[2],
						rightPlg.vecPoints[i][0], rightPlg.vecPoints[i][1], rightPlg.vecPoints[i][2],
						rightPlg.vecPoints[i]);
	}
	rightPlg.vecPoints[rightPlg.vecPoints.size() - 1] = rightPlg.vecPoints[0];
	//std::cout << _xAxis * _yAxis << std::endl;
	angle_between_projected_segment_and_X_axis_LinkNumber(_planeNormal[0], _planeNormal[1], _planeNormal[2],
					_yAxis[0], _yAxis[1], _yAxis[2],
					_yAxis);
	angle_between_projected_segment_and_X_axis_LinkNumber(_planeNormal[0], _planeNormal[1], _planeNormal[2],
					_xAxis[0], _xAxis[1], _xAxis[2],
					_xAxis);
	angle_between_projected_segment_and_X_axis_LinkNumber(_planeNormal[0], _planeNormal[1], _planeNormal[2],
					_planeNormal[0], _planeNormal[1], _planeNormal[2],
					_planeNormal);
	// rotate the basepoint
*/
    return;
}

bool LinkNumberPairPolygon::ComputeLinkNumber_CGAL() {
    _link_number = IntersectonComputing::LinkNumberOfTow3DLines(leftPlg.vecPoints, rightPlg.vecPoints);
    return true;

}

void LinkNumberPairPolygon::NewValueForLinkNumber(const double orientationFlag,
                                                  const int heightFlag,
                                                  const int incrementalValue) {
    if (orientationFlag > 0.0) {// uv is clockwise crossing ab
        // amplify the crossing number 4 times
        // at the end it needs to divide 8 to get the link number
        _link_number += (incrementalValue * heightFlag);
    } else {// uv is counterclockwise crossing ab
        _link_number -= (incrementalValue * heightFlag);
    }
}

bool LinkNumberPairPolygon::ComputeLinkNumber() {
    //preq: the plane normal is aligned with Z-axis
    //		the projection is only to take the XY coord
    // iterate all segments in one link to find the intersection
    // and determine the crossing number
    // iterate the segments in left polygon
    bool ret = true;
    double ab_x_max = 0.0;
    double ab_y_max = 0.0;
    double ab_x_min = 0.0;
    double ab_y_min = 0.0;
    double uv_x_max = 0.0;
    double uv_y_max = 0.0;
    double uv_x_min = 0.0;
    double uv_y_min = 0.0;
    //
    double ab_a[2];
    double ab_b[2];
    double uv_u[2];
    double uv_v[2];
    double intersection_pt[2];
    _link_number = 0;
    Vector3 projPt;
    for (unsigned int i = 0; i < leftPlg.vecPoints.size() - 1; i++) {
        /*ProjectingToPlane_LinkNumber(_ptOnPlane, _planeNormal, leftPlg.vecPoints[i], projPt);
        ab_a[0] = projPt[0];
        ab_a[1] = projPt[1];*/

        ab_a[0] = leftPlg.vecPoints[i][0];
        ab_a[1] = leftPlg.vecPoints[i][1];
        //
        /*ProjectingToPlane_LinkNumber(_ptOnPlane, _planeNormal, leftPlg.vecPoints[i+1], projPt);
        ab_a[0] = projPt[0];
        ab_a[1] = projPt[1];*/

        ab_b[0] = leftPlg.vecPoints[i + 1][0];
        ab_b[1] = leftPlg.vecPoints[i + 1][1];
        //
        ab_x_max = ab_a[0] > ab_b[0] ? ab_a[0] : ab_b[0];
        ab_y_max = ab_a[1] > ab_b[1] ? ab_a[1] : ab_b[1];
        ab_x_min = ab_a[0] < ab_b[0] ? ab_a[0] : ab_b[0];
        ab_y_min = ab_a[1] < ab_b[1] ? ab_a[1] : ab_b[1];
        // check the intersection between the candidate segment and other segments in rightPlg
        for (unsigned int j = 0; j < rightPlg.vecPoints.size() - 1; j++) {
            double s = 0.0;
            double t = 0.0;
            double h_ab = 0.0;
            double h_uv = 0.0;
            char common_point_ab = 0;
            char common_point_uv = 0;
            char common_point_num = 0;
            //
            /*ProjectingToPlane_LinkNumber(_ptOnPlane, _planeNormal, rightPlg.vecPoints[j], projPt);
            uv_u[0] = projPt[0];
            uv_u[1] = projPt[1];*/
            //
            uv_u[0] = rightPlg.vecPoints[j][0];
            uv_u[1] = rightPlg.vecPoints[j][1];
            //
            /*ProjectingToPlane_LinkNumber(_ptOnPlane, _planeNormal, rightPlg.vecPoints[j+1], projPt);
            uv_u[0] = projPt[0];
            uv_u[1] = projPt[1];*/

            uv_v[0] = rightPlg.vecPoints[j + 1][0];
            uv_v[1] = rightPlg.vecPoints[j + 1][1];
            // make sure the bounding box intersect first
            uv_x_max = uv_u[0] > uv_v[0] ? uv_u[0] : uv_v[0];
            uv_y_max = uv_u[1] > uv_v[1] ? uv_u[1] : uv_v[1];
            uv_x_min = uv_u[0] < uv_v[0] ? uv_u[0] : uv_v[0];
            uv_y_min = uv_u[1] < uv_v[1] ? uv_u[1] : uv_v[1];
            if (!((uv_x_max < ab_x_min) || (ab_x_max < uv_x_min) ||
                  (uv_y_max < ab_y_min) || (ab_y_max < uv_y_min))
                    ) {
                //
                char retVal = SegSegInt(ab_a, ab_b, uv_u, uv_v, intersection_pt);
                switch (retVal) {
                    case '1':
                        // regular intersection
                        // compute value of s
                        if (ab_a[0] == ab_b[0]) {//
                            s = (intersection_pt[1] - ab_a[1]) / (ab_b[1] - ab_a[1]);
                        } else {//
                            s = (intersection_pt[0] - ab_a[0]) / (ab_b[0] - ab_a[0]);
                        }
                        // compute value of t
                        if (uv_u[0] == uv_v[0]) {//
                            t = (intersection_pt[1] - uv_u[1]) / (uv_v[1] - uv_u[1]);
                        } else {//
                            t = (intersection_pt[0] - uv_u[0]) / (uv_v[0] - uv_u[0]);
                        }
                        //
                        h_ab = leftPlg.vecPoints[i][2] + s * (leftPlg.vecPoints[i + 1][2] - leftPlg.vecPoints[i][2]);
                        h_uv = rightPlg.vecPoints[j][2] + t * (rightPlg.vecPoints[j + 1][2] - rightPlg.vecPoints[j][2]);
                        //
                        if (h_ab > h_uv)
                            NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                  (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), 1, 4);
                        else
                            NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                  (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), -1, 4);

                        break;
                    case 'v':
                        // point in segment type of intersection
                        // += 2
                        // determine number of common endpoints
                        std::cout << " <<<<<<<<<<<<<< in point in segment case >>>>>>>>>>>>>>>>>>>>" << std::endl;
                        if (ab_a[0] == intersection_pt[0] && ab_a[1] == intersection_pt[1]) {
                            common_point_num++;
                            common_point_ab = 1;
                        }
                        if (ab_b[0] == intersection_pt[0] && ab_b[1] == intersection_pt[1]) {
                            common_point_num++;
                            common_point_ab = 2;
                        }
                        if (uv_u[0] == intersection_pt[0] && uv_u[1] == intersection_pt[1]) {
                            common_point_num++;
                            common_point_uv = 1;
                        }
                        if (uv_v[0] == intersection_pt[0] && uv_v[1] == intersection_pt[1]) {
                            common_point_num++;
                            common_point_uv = 2;
                        }
                        // common_point_intersection_type == 1 or 2
                        if (common_point_num > 2 || common_point_num == 0) {
                            std::cout << "WRONG INTERSECTION NUMBERS" << std::endl;
                            ret = false;
                        }
                        if (common_point_num == 1) {// intersection point is inside of one segment
                            if (common_point_ab) {// ab is the base segment
                                h_ab = (common_point_ab - 1) * leftPlg.vecPoints[i][2]
                                       + (common_point_ab - 1) * leftPlg.vecPoints[i + 1][2];
                                // compute value of t
                                if (uv_u[0] == uv_v[0]) {//
                                    t = (intersection_pt[1] - uv_u[1]) / (uv_v[1] - uv_u[1]);
                                } else {//
                                    t = (intersection_pt[0] - uv_u[0]) / (uv_v[0] - uv_u[0]);
                                }
                                //
                                h_uv = rightPlg.vecPoints[j][2] +
                                       t * (rightPlg.vecPoints[j + 1][2] - rightPlg.vecPoints[j][2]);
                                if (h_ab > h_uv)
                                    NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                          (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), 1, 2);
                                else
                                    NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                          (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), -1, 2);
                            } else {
                                // compute value of t
                                if (ab_a[0] == ab_b[0]) {//
                                    s = (intersection_pt[1] - ab_a[1]) / (ab_b[1] - ab_a[1]);
                                } else {//
                                    s = (intersection_pt[0] - ab_a[0]) / (ab_b[0] - ab_a[0]);
                                }
                                //
                                h_ab = leftPlg.vecPoints[i][2] +
                                       s * (leftPlg.vecPoints[i + 1][2] - leftPlg.vecPoints[i][2]);

                                if (common_point_uv == 1) {// u is on ab
                                    h_uv = rightPlg.vecPoints[j][2];
                                    if (h_ab > h_uv)
                                        NewValueForLinkNumber((uv_v[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_v[1] - ab_a[1]), 1, -2);
                                    else
                                        NewValueForLinkNumber((uv_v[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_v[1] - ab_a[1]), -1, -2);
                                } else {// v is on ab
                                    h_uv = rightPlg.vecPoints[j + 1][2];
                                    if (h_ab > h_uv)
                                        NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), 1, 2);
                                    else
                                        NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), -1, 2);
                                }
                            }
                        } else {// two common intersections
                            if (common_point_ab == 1) {
                                h_ab = leftPlg.vecPoints[i][2];
                                if (common_point_uv == 1) {// u is the common point
                                    h_uv = rightPlg.vecPoints[j][2];
                                    if (h_ab > h_uv)
                                        NewValueForLinkNumber((uv_v[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_v[1] - ab_a[1]), 1, -1);
                                    else
                                        NewValueForLinkNumber((uv_v[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_v[1] - ab_a[1]), -1, -1);
                                } else {// v is the common point
                                    h_uv = rightPlg.vecPoints[j + 1][2];
                                    if (h_ab > h_uv)
                                        NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), 1, 1);
                                    else
                                        NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), -1, 1);
                                }
                            } else {//  == 2
                                h_ab = leftPlg.vecPoints[i + 1][2];
                                if (common_point_uv == 1) {// u is the common point
                                    h_uv = rightPlg.vecPoints[j][2];
                                    if (h_ab > h_uv)
                                        NewValueForLinkNumber((uv_v[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_v[1] - ab_a[1]), 1, -1);
                                    else
                                        NewValueForLinkNumber((uv_v[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_v[1] - ab_a[1]), -1, -1);
                                } else {// v is the common point
                                    h_uv = rightPlg.vecPoints[j + 1][2];
                                    if (h_ab > h_uv)
                                        NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), 1, 1);
                                    else
                                        NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                              (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), -1, 1);
                                }
                            }
                        }
                        break;
                    case 'e':
                        //collinear intersection
                        //+= 1
                        std::cout << " >>>>>>>>>>>>> in collinear case <<<<<<<<<<<<" << std::endl;
                        if (ab_a[0] == intersection_pt[0] && ab_a[1] == intersection_pt[1]) {
                            common_point_num++;
                            common_point_ab = 1;
                        }
                        if (ab_b[0] == intersection_pt[0] && ab_b[1] == intersection_pt[1]) {
                            common_point_num++;
                            common_point_ab = 2;
                        }
                        if (uv_u[0] == intersection_pt[0] && uv_u[1] == intersection_pt[1]) {
                            common_point_num++;
                            common_point_uv = 1;
                        }
                        if (uv_v[0] == intersection_pt[0] && uv_v[1] == intersection_pt[1]) {
                            common_point_num++;
                            common_point_uv = 2;
                        }
                        // common_point_intersection_type == 1 or 2
                        if (common_point_num != 2) {
                            std::cout << "WRONG INTERSECTION NUMBERS ON COLINEAR CASE" << std::endl;
                            ret = false;
                        }
                        if (common_point_ab == 1) {
                            h_ab = leftPlg.vecPoints[i][2];
                            if (common_point_uv == 1) {// u is the common point
                                h_uv = rightPlg.vecPoints[j][2];
                                if (h_ab > h_uv)
                                    NewValueForLinkNumber((uv_v[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                          (ab_b[0] - ab_a[0]) * (uv_v[1] - ab_a[1]), 1, -1);
                                else
                                    NewValueForLinkNumber((uv_v[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                          (ab_b[0] - ab_a[0]) * (uv_v[1] - ab_a[1]), -1, -1);
                            } else {// v is the common point
                                h_uv = rightPlg.vecPoints[j + 1][2];
                                if (h_ab > h_uv)
                                    NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                          (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), 1, 1);
                                else
                                    NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                          (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), -1, 1);
                            }
                        } else {//  == 2
                            h_ab = leftPlg.vecPoints[i + 1][2];
                            if (common_point_uv == 1) {// u is the common point
                                h_uv = rightPlg.vecPoints[j][2];
                                if (h_ab > h_uv)
                                    NewValueForLinkNumber((uv_v[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                          (ab_b[0] - ab_a[0]) * (uv_v[1] - ab_a[1]), 1, -1);
                                else
                                    NewValueForLinkNumber((uv_v[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                          (ab_b[0] - ab_a[0]) * (uv_v[1] - ab_a[1]), -1, -1);
                            } else {// v is the common point
                                h_uv = rightPlg.vecPoints[j + 1][2];
                                if (h_ab > h_uv)
                                    NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                          (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), 1, 1);
                                else
                                    NewValueForLinkNumber((uv_u[0] - ab_a[0]) * (ab_b[1] - ab_a[1]) -
                                                          (ab_b[0] - ab_a[0]) * (uv_u[1] - ab_a[1]), -1, 1);
                            }
                        }
                        break;
                    case '0':
                        // doesn't intersect

                        break;
                }
            }// if bounding box
        }// for second
    }// for first
    if (abs(_link_number) % 8 != 0) {
        std::cout << "WRONG LINK NUMBER COMPUTATION" << std::endl;
        ret = false;
    } else
        _link_number = _link_number / 8;
    return ret;
}

bool LinkNumberPairPolygon::ComputeLinkNumberForTwoPolygons(_Polygon &inLeftPlg,
                                                            _Polygon &inRightPlg) {
    // Initiate two polygons
    AssignPolygons(inLeftPlg, inRightPlg);

    // Preprocess the two polygons
    //FindingProjectionPlaneWithoutDegeneratePointProjecting_2(1);

    //FindingProjectionPlaneWithoutDegeneratePointProjecting();
    //AvoidParallelOrthogonalToXYAxis();
    //AvoidOverlapping();
    FindingProjectionPlaneWithoutDegeneratePointProjecting_ANN();
    OptimizeXYAxisMakeNoParaOrtho();
    OptimizeNormalMakeNoOverlapping();
    //

    // ready to compute
    //bool ret = ComputeLinkNumber();
    bool ret = ComputeLinkNumber_CGAL();
    // output the information
    //std::cout << "crossing number : " << _link_number << std::endl;
    //
    return ret;
}
