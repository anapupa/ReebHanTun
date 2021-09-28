/*
(c) 2012 Fengtao Fan
*/
#include "UniqueDirection.h"
#include <vector>
#include <set>
#include <algorithm>

#define MY_LOCAL_PI 3.14159265

/*local functions used for LinkNumberComputing*/
float angle_between_projected_segment_and_X_axis(Vector3 &segment) {
    float angle = 0.0;
    // take the x_component
    return acosf(segment[0]);

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

void ProjectingToPlane(Vector3 &basePt, Vector3 &unitNormalDir, Vector3 &inPt, Vector3 &outPt) {
    Vector3 segDir = inPt - basePt;
    double projDist = (segDir * unitNormalDir);
    outPt = inPt - projDist * unitNormalDir;
    return;
}

void RotateAlongAxis(double u, double v, double w,
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

void Rotation_to_make_uvw_become_z_axis_directionComputing(double u, double v, double w,
                                                           double x, double y, double z,
                                                           Vector3 &outRes) {
    double uvw = 1.0 / sqrt(u * u + v * v + w * w);
    double uv = 1.0 / sqrt(u * u + v * v);
    outRes[0] = ((u * w * x + w * v * y) * uv - (z) / uv) * uvw;
    outRes[1] = (u * y - v * x) * uv;
    outRes[2] = (u * x + v * y + w * z) * uvw;
}

void INV_Rotation_to_make_uvw_become_z_axis_directionComputing(double u, double v, double w,
                                                               double x, double y, double z,
                                                               Vector3 &outRes) {
    double uvw = 1.0 / sqrt(u * u + v * v + w * w);
    double uv = 1.0 / sqrt(u * u + v * v);

    outRes[0] = u * w * uv * uvw * x - v * uv * y + u * uvw * z;
    outRes[1] = v * w * uv * uvw * x + u * uv * y + v * uvw * z;
    outRes[2] = uvw / uv * (-x) + w * uvw * z;
}
/**************************************/
/**************************************/

void UniqueDirectionComputation::FindingProjectionPlaneWithoutDegeneratePointProjecting(int ___int___) {//
    const int nSegmentNumber = segmentSet.size();
    Vector3 segment;
    vector<double> theta;
    vector<double> phi;
    // process the left polygon
    for (unsigned int i = 0; i < segmentSet.size(); i++) {
        segment = (*pointSetPtr)[segmentSet[i].first] - (*pointSetPtr)[segmentSet[i].second];
        // normalize this direction
        unitize(segment);
        //
        theta.push_back(angle_between_projected_segment_and_X_axis(segment));
        //
        if (norm2(segment) == 0.0)
            phi.push_back(0.0);
        else {
            //unitize(segment);
            //phi.push_back(MY_LOCAL_PI * 0.5 - std::min(	acos(segment * Vector3(0.0, 0.0, 1.0)),
            //											acos(segment * Vector3(0.0, 0.0, -1.0))));
            phi.push_back(acos(fabs(segment[2])));
        }
    }

    //search the largest distance from these various existing points in the angle parametrization space
    const int sqrtSegNumber = int(sqrt(double(nSegmentNumber))) + 1;
    float delta_angle = MY_LOCAL_PI / (4.0 * sqrtSegNumber);
    float max_angle_pair_dist = 0.0;
    std::pair<float, float> maxAnglePair;
    std::pair<float, float> anglePair;
    anglePair.first = 0.0;
    anglePair.second = 0.0;
    //
    maxAnglePair = anglePair;
    //
    for (int i = 0; i < 4 * sqrtSegNumber; i++) {
        anglePair.first += delta_angle;
        anglePair.second = 0.0;
        for (int j = 0; j < 2 * sqrtSegNumber; j++) {
            anglePair.second += delta_angle;
            //
            double max_dist_to_current_angle_pair = 0.0;
            double current_dist = 0.0;
            for (int k = 0; k < nSegmentNumber; k++) {
                current_dist = (anglePair.first - theta[k]) * (anglePair.first - theta[k]) +
                               (anglePair.second - phi[k]) * (anglePair.second - phi[k]);

                //current_dist =	std::min(	(anglePair.first - theta[k]) * (anglePair.first - theta[k]) +
                //							(anglePair.second - phi[k]) * (anglePair.second - phi[k]),
                //							(anglePair.first - theta[k] - 2 * MY_LOCAL_PI) * (anglePair.first - theta[k]- 2 * MY_LOCAL_PI) +
                //							(anglePair.second - phi[k]) * (anglePair.second - phi[k]));

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
    std::cout << "max angle : " << maxAnglePair.first << " " << maxAnglePair.second << std::endl;
    //calculate the direction according to the theta and phi parameter
    _planeNormal[0] = cos(maxAnglePair.first);
    _planeNormal[1] = sin(maxAnglePair.first);
    _planeNormal[2] = cos(maxAnglePair.second);
    unitize(_planeNormal);
    // calculate the lowest projection along the normal axis for each point in current XYZ-coordinate system
    _ptOnPlane = Vector3(0., 0., 0.);//_planeNormal * (LowestProjectionDistance(_planeNormal, leftPlg, rightPlg) - 2.0);
    //
    //for (unsigned int i = 0; i < (*pointSetPtr).size(); i++)
    //{// project the segment onto the XY-plane
    //	Rotation_to_make_uvw_become_z_axis_directionComputing(_planeNormal[0], _planeNormal[1], _planeNormal[2],
    //		(*pointSetPtr)[i][0], (*pointSetPtr)[i][1],(*pointSetPtr)[i][2],
    //		(*pointSetPtr)[i]);
    //}
    //
    //Rotation_to_make_uvw_become_z_axis_directionComputing(_planeNormal[0], _planeNormal[1], _planeNormal[2],
    //	_planeNormal[0], _planeNormal[1], _planeNormal[2],
    //	_planeNormal);
    //std::cout << _planeNormal[0] << " " << _planeNormal[1] << " " << _planeNormal[2] << " "<< std::endl;
    //_planeNormal[0] = 0.;
    //_planeNormal[1] = 0.;
    //_planeNormal[2] = 1.;
    _xAxis = Vector3(1., 0., 0.);
    _yAxis = Vector3(0., 1., 0.);
    //
    return;
}

//
void UniqueDirectionComputation::FindingProjectionPlaneWithoutDegeneratePointProjecting_2(int ___int___) {//
    const int nSegmentNumber = segmentSet.size();
    Vector3 segment;
    std::set<float> theta;
    std::set<float> phi;
    // process the left polygon
    for (unsigned int i = 0; i < segmentSet.size(); i++) {
        segment = (*pointSetPtr)[segmentSet[i].first] - (*pointSetPtr)[segmentSet[i].second];
        // normalize this direction
        unitize(segment);
        //
        theta.insert(angle_between_projected_segment_and_X_axis(segment));
        //
        if (norm2(segment) == 0.0)
            phi.insert(0.0);
        else {
            //unitize(segment);
            //phi.push_back(MY_LOCAL_PI * 0.5 - std::min(	acos(segment * Vector3(0.0, 0.0, 1.0)),
            //											acos(segment * Vector3(0.0, 0.0, -1.0))));
            phi.insert(float(acos(fabs(segment[2]))));
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
    std::cout << max_var_val << " in pair : " << max_cur_val << " " << max_pre_val << std::endl;
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
    std::cout << max_var_val << " in pair : " << max_cur_val << " " << max_pre_val << std::endl;
    maxAnglePair.second = 0.5 * (max_cur_val + max_pre_val);
    //search the largest distance from these various existing points in the angle parametrization space
    //const int sqrtSegNumber = int (sqrt(double (nSegmentNumber))) + 1;
    //float delta_angle = MY_LOCAL_PI / (4.0 * sqrtSegNumber);
    //float max_angle_pair_dist = 0.0;
    //std::pair<float, float> maxAnglePair;
    //std::pair<float, float> anglePair;
    //anglePair.first = 0.0;
    //anglePair.second = 0.0;
    ////
    //maxAnglePair = anglePair;
    ////
    //for (int i = 0; i < 4 * sqrtSegNumber; i++)
    //{
    //	anglePair.first += delta_angle;
    //	anglePair.second = 0.0;
    //	for (int j = 0; j < 2 * sqrtSegNumber; j++)
    //	{
    //		anglePair.second += delta_angle;
    //		//
    //		double max_dist_to_current_angle_pair = 0.0;
    //		double current_dist = 0.0;
    //		for (int k = 0; k < nSegmentNumber; k++)
    //		{
    //			current_dist = (anglePair.first - theta[k]) * (anglePair.first - theta[k]) +
    //										(anglePair.second - phi[k]) * (anglePair.second - phi[k]);

    //			//current_dist =	std::min(	(anglePair.first - theta[k]) * (anglePair.first - theta[k]) +
    //			//							(anglePair.second - phi[k]) * (anglePair.second - phi[k]),
    //			//							(anglePair.first - theta[k] - 2 * MY_LOCAL_PI) * (anglePair.first - theta[k]- 2 * MY_LOCAL_PI) +
    //			//							(anglePair.second - phi[k]) * (anglePair.second - phi[k]));

    //			if (current_dist > max_dist_to_current_angle_pair )
    //			{
    //				max_dist_to_current_angle_pair = current_dist;
    //			}
    //		}
    //		//
    //		if (max_dist_to_current_angle_pair > max_angle_pair_dist)
    //		{
    //			maxAnglePair = anglePair;
    //			max_angle_pair_dist = max_dist_to_current_angle_pair;
    //		}
    //	}
    //}//
    std::cout << "max angle : " << maxAnglePair.first << " " << maxAnglePair.second << std::endl;
    //calculate the direction according to the theta and phi parameter
    _planeNormal[0] = cos(maxAnglePair.first) * sin(maxAnglePair.second);
    _planeNormal[1] = sin(maxAnglePair.first) * sin(maxAnglePair.second);
    _planeNormal[2] = cos(maxAnglePair.second);
    //unitize(_planeNormal);
    // calculate the lowest projection along the normal axis for each point in current XYZ-coordinate system
    _ptOnPlane = Vector3(0., 0., 0.);//_planeNormal * (LowestProjectionDistance(_planeNormal, leftPlg, rightPlg) - 2.0);
    std::cout << "_planeNormal " << _planeNormal[0] << " " << _planeNormal[1] << " " << _planeNormal[2] << " "
              << std::endl;
    //

    for (unsigned int i = 0; i < (*pointSetPtr).size(); i++) {// project the segment onto the XY-plane
        Rotation_to_make_uvw_become_z_axis_directionComputing(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                                              (*pointSetPtr)[i][0], (*pointSetPtr)[i][1],
                                                              (*pointSetPtr)[i][2],
                                                              (*pointSetPtr)[i]);
    }

    //Rotation_to_make_uvw_become_z_axis_directionComputing(_planeNormal[0], _planeNormal[1], _planeNormal[2],
    //	_planeNormal[0], _planeNormal[1], _planeNormal[2],
    //	_planeNormal);


    //_planeNormal[0] = 0.;
    //_planeNormal[1] = 0.;
    //_planeNormal[2] = 1.;
    _xAxis = Vector3(1., 0., 0.);
    _yAxis = Vector3(0., 1., 0.);
    //
    return;
}

float FindMidValue(std::set<float>::iterator &startIter, std::set<float> &inSet) {
    std::set<float>::iterator nxtValIter;
    //
    float retValue = 0.0;
    //
    nxtValIter = startIter;
    nxtValIter++;
    if (nxtValIter == inSet.end()) {
        nxtValIter = inSet.begin();
        for (std::set<float>::iterator sIter = inSet.begin();
             sIter != startIter; sIter++) {
            nxtValIter = sIter;
        }
    } else {
        nxtValIter = startIter;
        nxtValIter++;
    }
    return ((*startIter + *nxtValIter) * 0.5);
}

//
void UniqueDirectionComputation::FindingProjectionPlaneWithoutDegeneratePointProjecting_MatchIdealDirection(
        int ___int___) {//
    const int nSegmentNumber = segmentSet.size();
    Vector3 segment;
    std::set<float> theta;
    std::set<float> phi;
    // process the left polygon
    for (unsigned int i = 0; i < segmentSet.size(); i++) {
        segment = (*pointSetPtr)[segmentSet[i].first] - (*pointSetPtr)[segmentSet[i].second];
        // normalize this direction
        unitize(segment);
        //
        theta.insert(angle_between_projected_segment_and_X_axis(segment));
        //
        if (norm2(segment) == 0.0)
            phi.insert(0.0);
        else {
            //unitize(segment);
            //phi.push_back(MY_LOCAL_PI * 0.5 - std::min(	acos(segment * Vector3(0.0, 0.0, 1.0)),
            //											acos(segment * Vector3(0.0, 0.0, -1.0))));
            phi.insert(float(acos(fabs(segment[2]))));
        }
    }
    // assumption: idealDirction is normalized.
    float ideal_dir_theta = angle_between_projected_segment_and_X_axis(idealDirection);
    float ideal_dir_phi = acos(fabs(idealDirection[2]));
    //
    std::pair<std::set<float>::iterator, bool> thetaRetIter = theta.insert(ideal_dir_theta);
    std::pair<std::set<float>::iterator, bool> phiRetIter = phi.insert(ideal_dir_phi);
    //
    std::set<float>::iterator nxtValIter;
    //
    std::pair<float, float> maxAnglePair;
    //
    //
    maxAnglePair.first = FindMidValue(thetaRetIter.first, theta);
    //
    maxAnglePair.second = FindMidValue(phiRetIter.first, phi);
    //search the largest distance from these various existing points in the angle parametrization space

    std::cout << "max angle : " << maxAnglePair.first << " " << maxAnglePair.second << std::endl;
    //calculate the direction according to the theta and phi parameter
    _planeNormal[0] = cos(maxAnglePair.first) * sin(maxAnglePair.second);
    _planeNormal[1] = sin(maxAnglePair.first) * sin(maxAnglePair.second);
    _planeNormal[2] = cos(maxAnglePair.second);
    //unitize(_planeNormal);
    // calculate the lowest projection along the normal axis for each point in current XYZ-coordinate system
    _ptOnPlane = Vector3(0., 0., 0.);//_planeNormal * (LowestProjectionDistance(_planeNormal, leftPlg, rightPlg) - 2.0);
    //
    for (unsigned int i = 0; i < (*pointSetPtr).size(); i++) {// project the segment onto the XY-plane
        Rotation_to_make_uvw_become_z_axis_directionComputing(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                                              (*pointSetPtr)[i][0], (*pointSetPtr)[i][1],
                                                              (*pointSetPtr)[i][2],
                                                              (*pointSetPtr)[i]);
    }

    //Rotation_to_make_uvw_become_z_axis_directionComputing(_planeNormal[0], _planeNormal[1], _planeNormal[2],
    //	_planeNormal[0], _planeNormal[1], _planeNormal[2],
    //	_planeNormal);
    std::cout << "_planeNormal " << _planeNormal[0] << " " << _planeNormal[1] << " " << _planeNormal[2] << " "
              << std::endl;
    //_planeNormal[0] = 0.;
    //_planeNormal[1] = 0.;
    //_planeNormal[2] = 1.;
    _xAxis = Vector3(1., 0., 0.);
    _yAxis = Vector3(0., 1., 0.);

    //ProjectingToPlane(_ptOnPlane, _planeNormal, (*pointSetPtr)[segmentSet[0].first], _xAxis);
    //ProjectingToPlane(_ptOnPlane, _planeNormal, (*pointSetPtr)[segmentSet[0].second], _yAxis);
    //_xAxis -= _yAxis;
    //unitize(_xAxis);
    //// taking the cross product
    //_yAxis = _xAxis^_planeNormal;
    //
    return;
}

//
void UniqueDirectionComputation::ComputeUniqueDirection(Vector3 &resDirection) {
//
    // calculating the intersection angle delta
    std::vector<double> delta_angle;
    // process the left polygon
    Vector3 segment;
    Vector3 intersectLineVec;
    for (unsigned int i = 0; i < segmentSet.size(); i++) {
        segment = (*pointSetPtr)[segmentSet[i].first] - (*pointSetPtr)[segmentSet[i].second];
        unitize(segment);
        //
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

        if (fabs(segment[2]) < 1e-10) {
            delta_angle.push_back(0.0);
        } else {
            delta_angle.push_back(atan(-segment[0] / segment[2]));
        }
    }
    //
    //
    //delta_angle.push_back(0.0);
    //delta_angle.push_back(MY_LOCAL_PI);
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
    resDirection[0] = cos(rotAngle);
    resDirection[1] = 0.0;
    resDirection[2] = sin(rotAngle);
    //
    std::cout << "resDirection " << resDirection[0] << " " <<
              resDirection[1] << " " <<
              resDirection[2] << std::endl;
    //
    INV_Rotation_to_make_uvw_become_z_axis_directionComputing(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                                              resDirection[0], resDirection[1], resDirection[2],
                                                              resDirection);

    //RotateAlongAxis(_planeNormal[0], _planeNormal[1], _planeNormal[2],
    //				_yAxis[0], _yAxis[1], _yAxis[2],
    //				rotAngle,
    //				resDirection);
    //unitize(resDirection);
    return;
}

void UniqueDirectionComputation::ComputeUniqueDirection_MatchIdealDirection(Vector3 &resDirection) {
//
    // calculating the intersection angle delta
    std::set<float> delta_angle;
    // process the left polygon
    Vector3 segment;
    Vector3 intersectLineVec;
    for (unsigned int i = 0; i < segmentSet.size(); i++) {
        segment = (*pointSetPtr)[segmentSet[i].first] - (*pointSetPtr)[segmentSet[i].second];
        unitize(segment);
        //
        //double denominator = segment[0] * _xAxis[2] - segment[2] * _xAxis[0];
        //denominator = 1.0 / denominator;
        //intersectLineVec[0] = denominator * (segment[1] * _xAxis[2] - segment[2] * _xAxis[1]);
        //intersectLineVec[2] = denominator * (segment[1] * _xAxis[0] - segment[0] * _xAxis[1]);
        //intersectLineVec[1] = -1;
        ////
        //unitize(intersectLineVec);
        ////
        //delta_angle.insert(float(acos(_yAxis * intersectLineVec)));
        //delta_angle.insert(float(acos(_yAxis * (-intersectLineVec))));
        if (fabs(segment[2]) < 1e-10) {
            delta_angle.insert(0.0);
        } else {
            delta_angle.insert(atan(-segment[0] / segment[2]));
        }
    }
    //
    //
    //delta_angle.insert(0.0);
    //delta_angle.insert(MY_LOCAL_PI);
    delta_angle.insert(MY_LOCAL_PI * 0.5);
    delta_angle.insert(MY_LOCAL_PI * (-0.5));
    // sorting the angle
    //std::sort(delta_angle.begin(), delta_angle.end());
    std::pair<std::set<float>::iterator, bool> idealIter;
    // find projected angle
    {
        Vector3 projPt;
        //void ProjectingToPlane(Vector3 &basePt, Vector3 & unitNormalDir, Vector3 &inPt, Vector3 &outPt)
        Vector3 zeroPt(0., 0., 0.0);
        ProjectingToPlane(zeroPt, _planeNormal, idealDirection, projPt);
        //
        unitize(projPt);
        //
        idealIter = delta_angle.insert(float(acos(_yAxis * projPt)));
        //
    }
    //scan the array to find the max-tempty interval
    double rotAngle = FindMidValue(idealIter.first, delta_angle);
    //
    //double maxIntervalBegin = 0.0;
    //double maxIntervalEnd = 0.0;
    //for (unsigned int i = 1; i < delta_angle.size(); i++)
    //{
    //	if (delta_angle[i] - delta_angle[i-1] > maxIntervalEnd - maxIntervalBegin)
    //	{
    //		maxIntervalBegin = delta_angle[i - 1];
    //		maxIntervalEnd   = delta_angle[i];
    //	}
    //}
    ////
    //double rotAngle = (maxIntervalBegin + maxIntervalEnd) * 0.5;
    //
    resDirection[0] = cos(rotAngle);
    resDirection[1] = 0.0;
    resDirection[2] = sin(rotAngle);
    //
    std::cout << resDirection[0] << " " <<
              resDirection[1] << " " <<
              resDirection[2] << std::endl;
    //
    INV_Rotation_to_make_uvw_become_z_axis_directionComputing(_planeNormal[0], _planeNormal[1], _planeNormal[2],
                                                              resDirection[0], resDirection[1], resDirection[2],
                                                              resDirection);

    //RotateAlongAxis(_planeNormal[0], _planeNormal[1], _planeNormal[2],
    //				_yAxis[0], _yAxis[1], _yAxis[2],
    //				rotAngle,
    //				resDirection);
    //unitize(resDirection);
    return;
}
