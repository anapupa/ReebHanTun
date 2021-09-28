/*
(c) 2012 Fengtao Fan
*/
#include "NonOverlappedLevelCycleAndArc.h"
#include <iomanip>

void NonOverlappedLevelCycleAndArc::WalkTrianglesAroundVertex(const int vid, Vector3 &outNormal) {
    int eid = 0;
    outNormal = Vector3(0., 0.0, 0.0);
    //
    std::set<int> trianglesSet;
    for (unsigned int i = 0; i < inMeshPtr->vecVertex[vid].adjEdges.size(); i++) {
        eid = inMeshPtr->vecVertex[vid].adjEdges[i];
        //
        trianglesSet.insert(inMeshPtr->vecEdge[eid].AdjTri[0]);
        if (inMeshPtr->vecEdge[eid].AdjTriNum == 2)
            trianglesSet.insert(inMeshPtr->vecEdge[eid].AdjTri[1]);
    }
    //
    for (std::set<int>::iterator sIter = trianglesSet.begin();
         sIter != trianglesSet.end();
         sIter++) {
        outNormal = outNormal + (*(inMeshPtr->meshNormalPtr))[*sIter];
    }
    outNormal = outNormal / norm(outNormal);
    return;
}

void NonOverlappedLevelCycleAndArc::WalkEdgeForNormal(const int eid, Vector3 &outNormal) {
    outNormal = Vector3(0., 0.0, 0.0);
    outNormal = outNormal + (*(inMeshPtr->meshNormalPtr))[inMeshPtr->vecEdge[eid].AdjTri[0]];
    if (inMeshPtr->vecEdge[eid].AdjTriNum == 2) {
        outNormal = outNormal + (*(inMeshPtr->meshNormalPtr))[inMeshPtr->vecEdge[eid].AdjTri[1]];
        //outNormal = outNormal / norm(outNormal);

    }
    unitize(outNormal);
    return;
}

void NonOverlappedLevelCycleAndArc::TranslatePoint(const double dir, const double scale, const Vector3 &curPt,
                                                   const Vector3 &normal, Vector3 &outPt) {

    //
    double sq_dist_ab = scale * 0.5 * dir; // this is a heuristic value
    if (scale > 5.0)
        sq_dist_ab = 5.0 * dir;

    outPt = curPt + sq_dist_ab * normal;
}

void NonOverlappedLevelCycleAndArc::ComputeIntersectionBetweenLevelCycleAndArc() {
    //
    levelDuplicate = false;
    arcDuplicate = false;
    //
    IntersectionOnMeshComputation intersect_computing;
    //
    intersect_computing.InitMeshPtr(inMeshPtr);
    // find the segment where it possible intersects the level cycle
    const double levelsetHeight = heightDirection * augmentedLevelCycles.vecPoints.front();
    pointPositionInArcVector = 0;
    pointPositionInLevelCycleVector = 0; // here use it as a counter
    double curHeight = heightDirection * augmentedArc.vecPoints[0];
    double nextHeight = 0.;
    for (unsigned int i = 0; i < augmentedArc.vecPoints.size() - 1; i++) {
        nextHeight = heightDirection * augmentedArc.vecPoints[i + 1];
        if (levelsetHeight <= curHeight && levelsetHeight >= nextHeight) {
            pointPositionInArcVector = i;
            pointPositionInLevelCycleVector++;
        }
        curHeight = nextHeight;
    }
    if (pointPositionInLevelCycleVector != 1) {
        std::cout << "MORE than one intersection between level set and arc" << std::endl;
    }
    //
    Vector3 commonPt;
    bool intersection_exist = false;
    // iterate all segments in the
    for (unsigned int i = 0; i < augmentedLevelCycles.vecPoints.size() - 1; i++) {
        intersect_computing.Clear();
        intersect_computing.InitFourPoints((*_arcPointTypePtr)[pointPositionInArcVector],
                                           (*_arcPointTypePtr)[pointPositionInArcVector + 1],
                                           (*_levelCyclePointTypePtr)[i],
                                           (*_levelCyclePointTypePtr)[i + 1],
                                           &(augmentedArc.vecPoints[pointPositionInArcVector]),
                                           &(augmentedArc.vecPoints[pointPositionInArcVector + 1]),
                                           &(augmentedLevelCycles.vecPoints[i]),
                                           &(augmentedLevelCycles.vecPoints[i + 1])
        );
        //
        intersect_computing.ComputeIntersection();
        //
        if (intersect_computing.bIntersection) {// this is the intersection
            // save the normal for this intersection
            intersection_exist = true;
            switch (intersect_computing.intersect_pt_type.second) {
                case 0:
                    // intersect at the vertex
                    WalkTrianglesAroundVertex(intersect_computing.intersect_pt_type.first, pointNormal);
                    break;
                case 1:
                    // intersect at the edge
                    WalkEdgeForNormal(intersect_computing.intersect_pt_type.first, pointNormal);
                    break;
                case 2:
                    // intersect at the face
                    pointNormal = (*(inMeshPtr->meshNormalPtr))[intersect_computing.intersect_pt_type.first];
                    break;
            }
            // save the position for level set
            if (intersect_computing.CD_is_intersection) {// position is fixed
                pointPositionInLevelCycleVector = i + intersect_computing.CD_is_intersection - 1;
                levelDuplicate = true;
            } else {// between i to i+1
                pointPositionInLevelCycleVector = i;
            }
            if (intersect_computing.AB_is_intersection) {// change it ot the fixed point
                pointPositionInArcVector = pointPositionInArcVector + intersect_computing.AB_is_intersection - 1;
                arcDuplicate = true;
            }
            // save the intersection point
            commonPt = intersect_computing.resIntersection;
            break;
        }
    }
    if (!intersection_exist) {
        std::cout << "NO INTERSECTION EXISTENCE" << std::endl;
        exit(0);
    }
    // now augment the two polygon
    // augment the level cycle
    double levelScale = 10000.0;
    if (!levelDuplicate) {
        //
        int ip1 = (pointPositionInLevelCycleVector + 1) % augmentedLevelCycles.vecPoints.size();
        int ip2 = (pointPositionInLevelCycleVector + 2) % augmentedLevelCycles.vecPoints.size();
        //
        int im1 = (pointPositionInLevelCycleVector - 1 + augmentedLevelCycles.vecPoints.size()) %
                  augmentedLevelCycles.vecPoints.size();
        //
        Vector3 backVector =
                augmentedLevelCycles.vecPoints[im1] - augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector];
        Vector3 frontVector = augmentedLevelCycles.vecPoints[ip2] - augmentedLevelCycles.vecPoints[ip1];
        Vector3 segVector =
                augmentedLevelCycles.vecPoints[ip1] - augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector];
        //
        unitize(backVector);
        unitize(frontVector);
        unitize(segVector);
        //
        //Vector3 levelNormal = segVector ^ heightDirection;
        //if (levelNormal * pointNormal < 0.0)
        //	levelNormal = -levelNormal;
        ////
        //pointNormal = levelNormal;
        //
        double backAngle = backVector * segVector;
        double frontAangle = frontVector * (-segVector);
        //
        if (backAngle > 0.70710678) // sqrt(2)/2
        {
            backAngle = ((1 - backAngle * backAngle) / (backAngle * backAngle)) *
                        norm2(augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector] - commonPt);
        } else {
            backAngle = norm2(augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector] - commonPt);
        }
        if (frontAangle > 0.70710678) // sqrt(2)/2
        {
            frontAangle = ((1 - frontAangle * frontAangle) / (frontAangle * frontAangle)) *
                          norm2(augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector + 1] - commonPt);
        } else {
            frontAangle = norm2(augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector + 1] - commonPt);
        }
        //
        augmentedLevelCycles.vecPoints.push_back(Vector3(0.0, 0.0, 0.0));
        for (int i = augmentedLevelCycles.vecPoints.size() - 1; i > pointPositionInLevelCycleVector; i--) {
            augmentedLevelCycles.vecPoints[i] = augmentedLevelCycles.vecPoints[i - 1];
        }
        pointPositionInLevelCycleVector++;
        augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector] = commonPt;
        //
        levelScale = std::min(backAngle, frontAangle);
        ////
        //int i1 = pointPositionInLevelCycleVector - 1;
        //if (i1 < 0)
        //	i1 = augmentedLevelCycles.vecPoints.size() - 2;
        ////
        //levelScale = std::min(norm2(augmentedLevelCycles.vecPoints[i1] - augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector]),
        //					    norm2(augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector+1] - augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector]));
        ////
    } else {
        //
        int ip1 = (pointPositionInLevelCycleVector + 1) % augmentedLevelCycles.vecPoints.size();
        if (pointPositionInLevelCycleVector == augmentedLevelCycles.vecPoints.size() - 1) {
            ip1 = 1;
        }
        //
        int im1 = (pointPositionInLevelCycleVector - 1 + augmentedLevelCycles.vecPoints.size()) %
                  augmentedLevelCycles.vecPoints.size();
        if (pointPositionInLevelCycleVector == 0) {
            im1 = augmentedLevelCycles.vecPoints.size() - 2;
        }
        //
        Vector3 backVector =
                augmentedLevelCycles.vecPoints[im1] - augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector];
        Vector3 frontVector =
                augmentedLevelCycles.vecPoints[ip1] - augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector];
        //
        //Vector3 levelNormal = frontVector ^ heightDirection;
        //if (levelNormal * pointNormal < 0.0)
        //	levelNormal = -levelNormal;
        ////
        //pointNormal = levelNormal;
        //
        levelScale = std::min(norm2(backVector), norm2(frontVector));
        //std::cout << "levelDuplicate " << std::endl;
        //
        //std::cout << "levelDuplicate " << std::endl;
    }
    //
    double arcScale = 10000.0;
    if (!arcDuplicate) {
        augmentedArc.vecPoints.push_back(Vector3(0.0, 0.0, 0.0));
        for (int i = augmentedArc.vecPoints.size() - 1; i > pointPositionInArcVector; i--) {
            augmentedArc.vecPoints[i] = augmentedArc.vecPoints[i - 1];
        }
        pointPositionInArcVector++;
        augmentedArc.vecPoints[pointPositionInArcVector] = commonPt;
        //
        //
        int i1 = pointPositionInArcVector - 1;
        if (i1 < 0)
            i1 = augmentedArc.vecPoints.size() - 2;
        //
        arcScale = std::min(norm2(augmentedArc.vecPoints[i1] - augmentedArc.vecPoints[pointPositionInArcVector]),
                            norm2(augmentedArc.vecPoints[pointPositionInArcVector + 1] -
                                  augmentedArc.vecPoints[pointPositionInArcVector]));
        //
    } else {
        std::cout << "arcDuplicate " << std::endl;
        std::cout << "in general, this can not happen" << std::endl;
        exit(8);

    }
    translateScale = std::min(levelScale, arcScale);
    translateScale = sqrt(translateScale);
    //std::cout << " scale " << setprecision(20) << translateScale << std::endl;
    return;
}

void NonOverlappedLevelCycleAndArc::ComputeLoopOnMeshPolygon_withArcHeightInfo() {

// path is encoded in this way
// v0--e0--v1--e1--....vn--en--v(n+1)[==v0]
// each cycle are passing from high to low
// each edge is taken in this way [a...b)-[b...c)-[c...a)
    int curArcId = 0;
    int startNodeId = 0;
    int endNodeId = 0;
    Vector3 translatedPt;
    Vector3 prePt, nxtPt;
    // walk around the path in terms of simplified arcs
    for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
        curArcId = _basisLoopPtr->simpArcIndexSet[i];
        startNodeId = _basisLoopPtr->nodeIndexSet[i];
        endNodeId = _basisLoopPtr->nodeIndexSet[i + 1];
        //
        if (!(*_arcHeightRangeIntersectingLevelsetPtr)[i]) {// just use the segment connecting two critical points
            // because its height range excludes the level set value
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction
                // and ignore the first point
                outVerticalCycle.vecPoints.push_back((*_arcOnMeshPtr)[curArcId].vecPoints.back());
            } else {// ignore the last point
                outVerticalCycle.vecPoints.push_back((*_arcOnMeshPtr)[curArcId].vecPoints.front());
            }
        } else {
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction
                //
                if (!move_common_pt_on_level_cycle && curArcId == _arcId) {
                    for (int iarc = augmentedArc.vecPoints.size() - 1; iarc >
                                                                       0; iarc--) {// be aware of which arc to take like _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
                        outVerticalCycle.vecPoints.push_back(augmentedArc.vecPoints[iarc]);
                    }
                } else {
                    for (int iarc = (*_arcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc >
                                                                                     0; iarc--) {// be aware of which arc to take like _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
                        outVerticalCycle.vecPoints.push_back((*_arcOnMeshPtr)[curArcId].vecPoints[iarc]);
                    }
                }
            } else {// ignore the last point

                if (!move_common_pt_on_level_cycle && curArcId == _arcId) {
                    for (unsigned int iarc = 0; iarc < augmentedArc.vecPoints.size() - 1; iarc++) {
                        outVerticalCycle.vecPoints.push_back(augmentedArc.vecPoints[iarc]);
                    }
                } else {
                    for (unsigned int iarc = 0; iarc < (*_arcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {
                        outVerticalCycle.vecPoints.push_back((*_arcOnMeshPtr)[curArcId].vecPoints[iarc]);
                    }
                }
            }
        }
    }
    // push the first point into the loop to close it
    outVerticalCycle.vecPoints.push_back(outVerticalCycle.vecPoints[0]);
    //
    return;
}

void NonOverlappedLevelCycleAndArc::ComputeLoopOnMeshPolygon() {
// path is encoded in this way
// v0--e0--v1--e1--....vn--en--v(n+1)[==v0]
// each cycle are passing from high to low
// each edge is taken in this way [a...b)-[b...c)-[c...a)
    int curArcId = 0;
    int startNodeId = 0;
    int endNodeId = 0;
    Vector3 translatedPt;
    Vector3 prePt, nxtPt;
    // walk around the path in terms of simplified arcs
    for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
        curArcId = _basisLoopPtr->simpArcIndexSet[i];
        startNodeId = _basisLoopPtr->nodeIndexSet[i];
        endNodeId = _basisLoopPtr->nodeIndexSet[i + 1];
        //
        if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
            startNodeId) {// since the arc are traversed from high to low,
            // need to trave in the opposite direction
            //
            if (!move_common_pt_on_level_cycle && curArcId == _arcId) {
                for (int iarc = augmentedArc.vecPoints.size() - 1; iarc >
                                                                   0; iarc--) {// be aware of which arc to take like _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
                    outVerticalCycle.vecPoints.push_back(augmentedArc.vecPoints[iarc]);
                }
            } else {
                for (int iarc = (*_arcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc >
                                                                                 0; iarc--) {// be aware of which arc to take like _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
                    outVerticalCycle.vecPoints.push_back((*_arcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
            }
        } else {// ignore the last point

            if (!move_common_pt_on_level_cycle && curArcId == _arcId) {
                for (unsigned int iarc = 0; iarc < augmentedArc.vecPoints.size() - 1; iarc++) {
                    outVerticalCycle.vecPoints.push_back(augmentedArc.vecPoints[iarc]);
                }
            } else {
                for (unsigned int iarc = 0; iarc < (*_arcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {
                    outVerticalCycle.vecPoints.push_back((*_arcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
            }
        }
    }
    // push the first point into the loop to close it
    outVerticalCycle.vecPoints.push_back(outVerticalCycle.vecPoints[0]);
    //
    return;
}

void NonOverlappedLevelCycleAndArc::ComputeOutputCycles() {
    double dir = 1.0;
    int prevArcIdx = 0;
    int nextArcIdx = 0;
    double effectScale = translateScale;
    // update the translate scale by looking at the vertical loop
    for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
        if (_basisLoopPtr->simpArcIndexSet[i] == _arcId) {
            prevArcIdx = i - 1;
            nextArcIdx = i + 1;
            if (prevArcIdx < 0) {
                prevArcIdx = _basisLoopPtr->simpArcIndexSet.size() - 1;
            }
            if (nextArcIdx == _basisLoopPtr->simpArcIndexSet.size()) {
                nextArcIdx = 0;
            }
        }
    }
    //
    //prevArcIdx = _basisLoopPtr->simpArcIndexSet[prevArcIdx];
    //nextArcIdx = _basisLoopPtr->simpArcIndexSet[nextArcIdx];
    //
    int orgPtPos = pointPositionInArcVector - 1;
    if (orgPtPos != 0 && orgPtPos != (*_arcOnMeshPtr)[_arcId].vecPoints.size() - 2) {
        std::cout << "Wrong position in arc for level arc intersection computation" << std::endl;
        std::cout << "orgPtPos " << orgPtPos << " and size " << (*_arcOnMeshPtr)[_arcId].vecPoints.size() << std::endl;
        exit(9);
    } else {
        Vector3 nextPt, prevPt;
        if (orgPtPos == 0) {// intersection is between [0,1]
            nextPt = (*_arcOnMeshPtr)[_arcId].vecPoints[2];
            if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[prevArcIdx]].nCriticalNode0 ==
                _arcPointTypePtr->front().first) {
                if (!(*_arcHeightRangeIntersectingLevelsetPtr)[prevArcIdx]) {
                    prevPt = (*_arcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[prevArcIdx]].vecPoints.front();
                } else {
                    int posIndex = (*_arcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[prevArcIdx]].vecPoints.size() - 2;
                    prevPt = (*_arcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[prevArcIdx]].vecPoints[posIndex];
                }
            } else {
                if (!(*_arcHeightRangeIntersectingLevelsetPtr)[prevArcIdx]) {
                    prevPt = (*_arcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[prevArcIdx]].vecPoints.back();
                } else {
                    prevPt = (*_arcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[prevArcIdx]].vecPoints[1];
                }
            }
        } else {
            prevPt = (*_arcOnMeshPtr)[_arcId].vecPoints[orgPtPos - 1];
            if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[nextArcIdx]].nCriticalNode0 ==
                _arcPointTypePtr->front().first) {
                if (!(*_arcHeightRangeIntersectingLevelsetPtr)[nextArcIdx]) {// use only the segment
                    nextPt = (*_arcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[nextArcIdx]].vecPoints.front();
                } else {
                    int posIndex = (*_arcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[nextArcIdx]].vecPoints.size() - 2;
                    nextPt = (*_arcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[nextArcIdx]].vecPoints[posIndex];
                }
            } else {
                if (!(*_arcHeightRangeIntersectingLevelsetPtr)[nextArcIdx]) {// use only the segment
                    nextPt = (*_arcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[nextArcIdx]].vecPoints.back();
                } else {
                    nextPt = (*_arcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[nextArcIdx]].vecPoints[1];
                }
            }
        }
        //
        Vector3 backVector = prevPt - (*_arcOnMeshPtr)[_arcId].vecPoints[orgPtPos];
        Vector3 frontVector = nextPt - (*_arcOnMeshPtr)[_arcId].vecPoints[orgPtPos + 1];
        Vector3 segVector =
                (*_arcOnMeshPtr)[_arcId].vecPoints[orgPtPos + 1] - (*_arcOnMeshPtr)[_arcId].vecPoints[orgPtPos];
        //
        unitize(backVector);
        unitize(frontVector);
        unitize(segVector);
        //
        double backAngle = backVector * segVector;
        double frontAangle = frontVector * (-segVector);
        //
        Vector3 commonPt = augmentedArc.vecPoints[pointPositionInArcVector];
        if (backAngle > 0.70710678) // sqrt(2)/2
        {
            backAngle = ((1 - backAngle * backAngle) / (backAngle * backAngle)) *
                        norm2((*_arcOnMeshPtr)[_arcId].vecPoints[orgPtPos] - commonPt);
        } else {
            backAngle = norm2((*_arcOnMeshPtr)[_arcId].vecPoints[orgPtPos] - commonPt);
        }
        if (frontAangle > 0.70710678) // sqrt(2)/2
        {
            frontAangle = ((1 - frontAangle * frontAangle) / (frontAangle * frontAangle)) *
                          norm2((*_arcOnMeshPtr)[_arcId].vecPoints[orgPtPos + 1] - commonPt);
        } else {
            frontAangle = norm2((*_arcOnMeshPtr)[_arcId].vecPoints[orgPtPos + 1] - commonPt);
        }
        //
        frontAangle = std::min(frontAangle, backAngle);
        effectScale = std::min(effectScale, sqrt(frontAangle));
    }
    //
    if (move_common_pt_on_level_cycle) {// move the level cycle
        Vector3 orgPt = augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector];
        if (ownVerticalType) {// in vertical loop is wanted as vertical
            // so its dual is horizontal, go outside
            dir = 1.0;
        } else
            dir = -1.0;
        //
        //int i1 = pointPositionInLevelCycleVector - 1;
        //if (i1 < 0)
        //	i1 = augmentedLevelCycles.vecPoints.size() - 2;
        //std::cout << "effect scale " << effectScale << std::endl;
        TranslatePoint(dir, effectScale,// translateScale,
                       augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector],
                       pointNormal, augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector]);
        if (pointPositionInLevelCycleVector == 0)
            augmentedLevelCycles.vecPoints.back() = augmentedLevelCycles.vecPoints.front();
        if (pointPositionInLevelCycleVector == augmentedLevelCycles.vecPoints.size() - 1)
            augmentedLevelCycles.vecPoints.front() = augmentedLevelCycles.vecPoints.back();
        //
        outLevelCycle = augmentedLevelCycles;
        //
        augmentedLevelCycles.vecPoints[pointPositionInLevelCycleVector] = orgPt;
        if (pointPositionInLevelCycleVector == 0)
            augmentedLevelCycles.vecPoints.back() = augmentedLevelCycles.vecPoints.front();
        if (pointPositionInLevelCycleVector == augmentedLevelCycles.vecPoints.size() - 1)
            augmentedLevelCycles.vecPoints.front() = augmentedLevelCycles.vecPoints.back();
        //
        //ComputeLoopOnMeshPolygon();
        ComputeLoopOnMeshPolygon_withArcHeightInfo();
    } else {
        Vector3 orgPt = augmentedArc.vecPoints[pointPositionInArcVector];
        if (_basisLoopPtr->pathType) {// in vertical loop is wanted as vertical
            // so it go inside
            dir = -1.0;
        } else
            dir = 1.0;
        //
        //int i1 = pointPositionInArcVector - 1;
        //if (i1 < 0)
        //	i1 = augmentedArc.vecPoints.size() - 2;
        TranslatePoint(dir, effectScale, //translateScale,
                       augmentedArc.vecPoints[pointPositionInArcVector],
                       pointNormal, augmentedArc.vecPoints[pointPositionInArcVector]);
        //
        //std::cout << "effect scale " << effectScale << std::endl;
        //
        //ComputeLoopOnMeshPolygon();
        ComputeLoopOnMeshPolygon_withArcHeightInfo();
        //
        augmentedArc.vecPoints[pointPositionInArcVector] = orgPt;
        //
        outLevelCycle = *_levelCyclePtr;
    }
}
