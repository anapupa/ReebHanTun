/*
(c) 2012 Fengtao Fan
*/
#include "NonOverlappingCycles.h"

void NonOverlappingCycles::WalkTrianglesAroundVertex(const int vid, Vector3 &outNormal) {
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
    unitize(outNormal);// = outNormal / norm(outNormal);
    return;
}

void NonOverlappingCycles::ComputeOutNormals() {
    int vertexIdOnMesh = 0;
    Vector3 tmpNormal;
    std::set<int>::iterator sIter;
    for (sIter = _criticalNodeIdSetPtr->begin();
         sIter != _criticalNodeIdSetPtr->end();
         sIter++) {
        //if ((*_vecReebNodePtr)[*sIter].crType == UP_FORKING_REEB ||
        //	(*_vecReebNodePtr)[*sIter].crType == DOWN_FORKING_REEB )
        {
            vertexIdOnMesh = (*_vecReebNodePtr)[*sIter].nVertexId;
            //
            WalkTrianglesAroundVertex(vertexIdOnMesh, tmpNormal);
            //
            unitOutNormalForCritNode[*sIter] = tmpNormal;
        }
    }
    return;
}

void NonOverlappingCycles::WalkEdgeForNormal(const int eid, Vector3 &outNormal) {
    outNormal = Vector3(0., 0.0, 0.0);
    outNormal = outNormal + (*(inMeshPtr->meshNormalPtr))[inMeshPtr->vecEdge[eid].AdjTri[0]];
    if (inMeshPtr->vecEdge[eid].AdjTriNum == 2) {
        outNormal = outNormal + (*(inMeshPtr->meshNormalPtr))[inMeshPtr->vecEdge[eid].AdjTri[1]];
        //outNormal = outNormal / norm(outNormal);
        unitize(outNormal);
    }

    return;
}

void
NonOverlappingCycles::TranslatePoint(const double dir, const Vector3 &prePt, const Vector3 &curPt, const Vector3 &nexPt,
                                     const Vector3 &normal, Vector3 &outPt) {
    double sq_dist_ab = norm2(prePt - curPt);
    double sq_dist_ac = norm2(nexPt - curPt);
    sq_dist_ab = std::min(sq_dist_ab, sq_dist_ac);
    //
    sq_dist_ab = sqrt(sq_dist_ab) * 0.15 * dir; // this is a heuristic value
    //
    outPt = curPt + sq_dist_ab * normal;
}

void
NonOverlappingCycles::TranslatePoint(const double dir, const double scale, const Vector3 &curPt, const Vector3 &normal,
                                     Vector3 &outPt) {
    //
    double sq_dist_ab = scale * 0.15 * dir; // this is a heuristic value
    //
    //if (scale > 1)
    //	sq_dist_ab = 50 * dir;
    std::cout << "sq " << sq_dist_ab << std::endl;
    outPt = curPt + sq_dist_ab * normal;
}

void NonOverlappingCycles::SetArcSharedStatus() {
    if (_basisLoopPtr == _loopOnMeshPtr) {

        for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
            arcSharedStatus[_basisLoopPtr->simpArcIndexSet[i]] = 1;
        }
    } else {
        std::set<int> arcSets(_basisLoopPtr->simpArcIndexSet.begin(), _basisLoopPtr->simpArcIndexSet.end());
        std::set<int>::iterator findIter;
        for (unsigned int i = 0; i < _loopOnMeshPtr->simpArcIndexSet.size(); i++) {
            int arcIdx = _loopOnMeshPtr->simpArcIndexSet[i];
            findIter = arcSets.find(arcIdx);
            if (findIter != arcSets.end())
                arcSharedStatus[arcIdx] = 1;
        }
        //
    }
    return;

}

void NonOverlappingCycles::ComputeLoopOnMeshPolygonOnFly(std::vector<int> sharedArcs,
                                                         std::vector<_Polygon> &sharedArcsNewPoints,
                                                         std::vector<_Polygon> &offsetSharedArcsNewPoints) {
// path is encoded in this way
// v0--e0--v1--e1--....vn--en--v(n+1)[==v0]
// each cycle are passing from high to low
// each edge is taken in this way [a...b)-[b...c)-[c...a)
    int curArcId = 0;
    int startNodeId = 0;
    int endNodeId = 0;
    Vector3 translatedPt;
    Vector3 prePt, nxtPt;
    //
    int shared_arc_vec_index = 0;
    _Polygon *actualPlyPtr = NULL;
    // walk around the path in terms of simplified arcs
    for (unsigned int i = 0; i < _loopOnMeshPtr->simpArcIndexSet.size(); i++) {
        curArcId = _loopOnMeshPtr->simpArcIndexSet[i];
        startNodeId = _loopOnMeshPtr->nodeIndexSet[i];
        endNodeId = _loopOnMeshPtr->nodeIndexSet[i + 1];
        //
        actualPlyPtr = &(*_pathArcOnMeshPtr)[curArcId];
        // check it is in the shared arc or not
        shared_arc_vec_index = -1;
        for (int j = 0; j < sharedArcs.size(); j++) {
            if (sharedArcs[j] == curArcId) {
                shared_arc_vec_index = j;
                break;
            }
        }
        if (shared_arc_vec_index >= 0) {
            actualPlyPtr = &offsetSharedArcsNewPoints[shared_arc_vec_index];
        }
        //
        if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
            startNodeId) {// since the arc are traversed from high to low,
            // need to trave in the opposite direction
            //
            for (int iarc = actualPlyPtr->vecPoints.size() - 1;
                 iarc > 0; iarc--) {// be aware of which arc to take like _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
                outLoopOnMesh.vecPoints.push_back(actualPlyPtr->vecPoints[iarc]);
            }
        } else {// ignore the last point

            for (unsigned int iarc = 0; iarc < actualPlyPtr->vecPoints.size() - 1; iarc++) {
                outLoopOnMesh.vecPoints.push_back(actualPlyPtr->vecPoints[iarc]);
            }
        }
    }
    // push the first point into the loop to close it
    outLoopOnMesh.vecPoints.push_back(outLoopOnMesh.vecPoints[0]);
    //
    return;
}

void NonOverlappingCycles::ComputeLoopOnMeshPolygon() {
// path is encoded in this way
// v0--e0--v1--e1--....vn--en--v(n+1)[==v0]
// each cycle are passing from high to low
// each edge is taken in this way [a...b)-[b...c)-[c...a)
    int curArcId = 0;
    int startNodeId = 0;
    int endNodeId = 0;
    Vector3 translatedPt;
    Vector3 prePt, nxtPt;
    //
    // walk around the path in terms of simplified arcs
    for (unsigned int i = 0; i < _loopOnMeshPtr->simpArcIndexSet.size(); i++) {
        curArcId = _loopOnMeshPtr->simpArcIndexSet[i];
        startNodeId = _loopOnMeshPtr->nodeIndexSet[i];
        endNodeId = _loopOnMeshPtr->nodeIndexSet[i + 1];
        //
        if ((*_bVecLoopOnMeshPtr)[i]) {
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 == startNodeId) {// ignore the first point
                outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.back());
            } else {// ignore the last point
                outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.front());
            }
        } else {
            //
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction
                //
                for (int iarc = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc >
                                                                                           0; iarc--) {// be aware of which arc to take like _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
                    outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
            } else {// ignore the last point

                for (unsigned int iarc = 0; iarc < (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {
                    outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
            }
        }
    }
    // push the first point into the loop to close it
    outLoopOnMesh.vecPoints.push_back(outLoopOnMesh.vecPoints[0]);
    //
    return;
}

void NonOverlappingCycles::ComputeBasisLoopPolygonOnFly(std::vector<int> sharedArcs,
                                                        std::vector<_Polygon> &sharedArcsNewPoints,
                                                        std::vector<_Polygon> &offsetSharedArcsNewPoints) {
    // make a mutual decision on _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
    int curArcId = 0;
    int startNodeId = 0;
    int endNodeId = 0;
    Vector3 translatedPt;
    Vector3 prePt, nxtPt;
    //
    //
    int shared_arc_vec_index = 0;
    _Polygon *actualPlyPtr = NULL;
    //
    if (_basisLoopPtr == _loopOnMeshPtr) {// need to embed every critical node on the path
        for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
            curArcId = _basisLoopPtr->simpArcIndexSet[i];
            startNodeId = _basisLoopPtr->nodeIndexSet[i];
            endNodeId = _basisLoopPtr->nodeIndexSet[i + 1];
            //
            actualPlyPtr = &(*_pathArcOnMeshPtr)[curArcId];
            // check it is in the shared arc or not
            shared_arc_vec_index = -1;
            for (int j = 0; j < sharedArcs.size(); j++) {
                if (sharedArcs[j] == curArcId) {
                    shared_arc_vec_index = j;
                    break;
                }
            }
            if (shared_arc_vec_index >= 0) {
                actualPlyPtr = &sharedArcsNewPoints[shared_arc_vec_index];
            }
            //
            // prePt and nxtPt is used for scale information, so no need to take care
            if (i == 0) {
                if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet.back()].nCriticalNode0 == startNodeId) {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[
                            (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints.size() - 2];
                } else {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[1];
                }
            } else {
                if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[i - 1]].nCriticalNode0 == startNodeId) {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[
                            (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints.size() - 2];
                } else {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[1];
                }
            }
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction
                double dir = 1.0;
                if (_basisLoopPtr->pathType) {// vertical, so pull inside
                    dir = -1.0;
                }
                nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[(*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2];
                //

                //
                TranslatePoint(dir, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints.back(),
                               nxtPt, unitOutNormalForCritNode[startNodeId],
                               translatedPt);
                //
                outBasisLoop.vecPoints.push_back(translatedPt);
                //
                //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.back());
                //
                for (int iarc = actualPlyPtr->vecPoints.size() - 2; iarc > 0; iarc--) {// ignore the first point

                    outBasisLoop.vecPoints.push_back(actualPlyPtr->vecPoints[iarc]);
                }
                /*for ( int iarc = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 2; iarc > 0 ; iarc--)
				{
					outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
				}*/
            } else {// ignore the last point
                double dir = 1.0;
                if (_basisLoopPtr->pathType) {// vertical, so pull inside
                    dir = -1.0;
                }
                nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[1];
                //

                //
                TranslatePoint(dir, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints[0],
                               nxtPt, unitOutNormalForCritNode[startNodeId],
                               translatedPt);
                //
                outBasisLoop.vecPoints.push_back(translatedPt);
                //
                //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[0]);
                for (unsigned int iarc = 1;
                     iarc < actualPlyPtr->vecPoints.size() - 1; iarc++) {// ignore the first point

                    outBasisLoop.vecPoints.push_back(actualPlyPtr->vecPoints[iarc]);
                }/*
				for (unsigned int iarc = 1; iarc < (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 1 ; iarc++)
				{
					outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
				}*/

            }
        }
        // push the first point into the loop to close it
        //outLoopOnMesh.vecPoints.push_back(outLoopOnMesh.vecPoints[0]);
        //
        outBasisLoop.vecPoints.push_back(outBasisLoop.vecPoints[0]);
    } else {// only to embed the shared critical node

        std::set<int> shared_nodes;
        for (unsigned int i = 0; i < _basisLoopPtr->nodeIndexSet.size(); i++) {
            for (unsigned int j = 0; j < _loopOnMeshPtr->nodeIndexSet.size(); j++) {
                if (_basisLoopPtr->nodeIndexSet[i] == _loopOnMeshPtr->nodeIndexSet[j])
                    shared_nodes.insert(_loopOnMeshPtr->nodeIndexSet[j]);
            }
        }
        std::set<int>::iterator sIter;
        //
        for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
            curArcId = _basisLoopPtr->simpArcIndexSet[i];
            startNodeId = _basisLoopPtr->nodeIndexSet[i];
            endNodeId = _basisLoopPtr->nodeIndexSet[i + 1];
            //
            actualPlyPtr = &(*_pathArcOnMeshPtr)[curArcId];
            // check it is in the shared arc or not
            shared_arc_vec_index = -1;
            for (int j = 0; j < sharedArcs.size(); j++) {
                if (sharedArcs[j] == curArcId) {
                    shared_arc_vec_index = j;
                    break;
                }
            }
            if (shared_arc_vec_index >= 0) {
                actualPlyPtr = &sharedArcsNewPoints[shared_arc_vec_index];
            }
            //
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction
                sIter = shared_nodes.find(startNodeId);
                if (sIter != shared_nodes.end()) {// need to pull out or in
                    double dir = 1.0;
                    if (_basisLoopPtr->pathType) {// vertical, so pull inside
                        dir = -1.0;
                    }
                    nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[(*_pathArcOnMeshPtr)[curArcId].vecPoints.size() -
                                                                     2];
                    //
                    if (i == 0) {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet.back()].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[1];
                        }
                    } else {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[i - 1]].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[1];
                        }
                    }
                    //
                    TranslatePoint(dir, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints.back(),
                                   nxtPt, unitOutNormalForCritNode[startNodeId],
                                   translatedPt);
                    //
                    outBasisLoop.vecPoints.push_back(translatedPt);
                } else {
                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints.back());
                }
                //
                //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.back());
                //
                for (int iarc = actualPlyPtr->vecPoints.size() - 2; iarc > 0; iarc--) {// ignore the first point
                    outBasisLoop.vecPoints.push_back(actualPlyPtr->vecPoints[iarc]);

                }/*
				for ( int iarc = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 2; iarc > 0 ; iarc--)
				{
					outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
				}*/
            } else {// ignore the last point
                sIter = shared_nodes.find(startNodeId);
                if (sIter != shared_nodes.end()) {// need to pull out or in
                    double dir = 1.0;
                    if (_basisLoopPtr->pathType) {// vertical, so pull inside
                        dir = -1.0;
                    }
                    nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[1];
                    //
                    if (i == 0) {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet.back()].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[1];
                        }
                    } else {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[i - 1]].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[1];
                        }
                    }
                    //
                    TranslatePoint(dir, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints.front(),
                                   nxtPt, unitOutNormalForCritNode[startNodeId],
                                   translatedPt);
                    //
                    outBasisLoop.vecPoints.push_back(translatedPt);
                } else {
                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints.front());
                }
                //
                //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.front());
                for (unsigned int iarc = 1;
                     iarc < actualPlyPtr->vecPoints.size() - 1; iarc++) {// ignore the first point
                    outBasisLoop.vecPoints.push_back(actualPlyPtr->vecPoints[iarc]);

                }/*
				for (unsigned int iarc = 1; iarc < (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 1 ; iarc++)
				{
					outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
				}*/
            }
        }
        // push the first point into the loop to close it
        //outLoopOnMesh.vecPoints.push_back(outLoopOnMesh.vecPoints[0]);
        outBasisLoop.vecPoints.push_back(outBasisLoop.vecPoints[0]);

    }
}

std::pair<int, double> NonOverlappingCycles::smallestDisatnce(const Vector3 &refPt, const std::vector<Vector3> &nodes) {
    std::pair<int, double> dist = std::pair<int, double>(0, 0.0);
    double tempDist = 0.0;
    if (nodes.size() > 0) {
        dist.second = norm2(refPt - nodes[0]);
        for (unsigned int i = 1; i < nodes.size(); i++) {
            tempDist = norm2(refPt - nodes[i]);
            if (tempDist < dist.second) {
                dist.second = tempDist;
                dist.first = i;
            }
        }
    }
    dist.second = sqrt(dist.second);
    return dist;
}

void NonOverlappingCycles::PointBetweenTwoPoints(const double dist_to_src, const Vector3 &src, const Vector3 &dst,
                                                 Vector3 &res) {
    // src---dist_to_src--dst
    res = src + dist_to_src * (dst - src);
    return;
}

double NonOverlappingCycles::smallestDisatnce(const int refVer, const int eid,
                                              Vector3 &outPt) {// refVer---v1---v2---v....---opp_v
    double dist = 0.0;
    int opp_v = (inMeshPtr->vecEdge)[eid].v0 + (inMeshPtr->vecEdge)[eid].v1 - refVer;
    //
    if ((*_edgeIntersectArcsPtr)[eid].empty()) {
        dist = ((inMeshPtr->vecVertex)[refVer].x - (inMeshPtr->vecVertex)[opp_v].x) *
               ((inMeshPtr->vecVertex)[refVer].x - (inMeshPtr->vecVertex)[opp_v].x) +
               ((inMeshPtr->vecVertex)[refVer].y - (inMeshPtr->vecVertex)[opp_v].y) *
               ((inMeshPtr->vecVertex)[refVer].y - (inMeshPtr->vecVertex)[opp_v].y) +
               ((inMeshPtr->vecVertex)[refVer].z - (inMeshPtr->vecVertex)[opp_v].z) *
               ((inMeshPtr->vecVertex)[refVer].z - (inMeshPtr->vecVertex)[opp_v].z);
        //
        outPt[0] = (inMeshPtr->vecVertex)[opp_v].x;
        outPt[1] = (inMeshPtr->vecVertex)[opp_v].y;
        outPt[2] = (inMeshPtr->vecVertex)[opp_v].z;
    } else {
        double tempDist = 0.0;
        Vector3 refPt((inMeshPtr->vecVertex)[refVer].x,
                      (inMeshPtr->vecVertex)[refVer].y,
                      (inMeshPtr->vecVertex)[refVer].z);
        std::vector<std::pair<int, int> >::iterator vIter = (*_edgeIntersectArcsPtr)[eid].begin();
        dist = norm2((*_pathArcOnMeshPtr)[vIter->first].vecPoints[vIter->second] - refPt);
        outPt = (*_pathArcOnMeshPtr)[vIter->first].vecPoints[vIter->second];
        //
        tempDist = norm2((*_offsetPathArcOnMeshPtr)[vIter->first].vecPoints[vIter->second] - refPt);
        if (tempDist < dist) {
            dist = tempDist;
            outPt = (*_offsetPathArcOnMeshPtr)[vIter->first].vecPoints[vIter->second];
        }
        for (vIter++; vIter != (*_edgeIntersectArcsPtr)[eid].end(); vIter++) {
            tempDist = norm2((*_pathArcOnMeshPtr)[vIter->first].vecPoints[vIter->second] - refPt);
            if (tempDist < dist) {
                dist = tempDist;
                outPt = (*_pathArcOnMeshPtr)[vIter->first].vecPoints[vIter->second];
            }
            tempDist = norm2((*_offsetPathArcOnMeshPtr)[vIter->first].vecPoints[vIter->second] - refPt);
            if (tempDist < dist) {
                dist = tempDist;
                outPt = (*_offsetPathArcOnMeshPtr)[vIter->first].vecPoints[vIter->second];
            }
        }
    }
    dist = sqrt(dist);
    return dist;
}

void NonOverlappingCycles::NearByTwoEdges(const int pivot_v, const int eid, int &leftEdge, int &rightEdge) {
}

void NonOverlappingCycles::OrderPointsReferingToOneEndPoint(const Vector3 &refPt, const std::vector<Vector3> &nodes,
                                                            std::vector<int> &orders) {
    //
    if (nodes.size() == 2) {
        double dist1 = norm2(refPt - nodes[0]);
        double dist2 = norm2(refPt - nodes[1]);
        if (dist1 > dist2) {
            orders[0] = 1;
            orders[1] = 0;
        } else {
            orders[0] = 0;
            orders[1] = 1;
        }
    } else {
        if (nodes.size() == 3) {
            double dist1 = norm2(refPt - nodes[0]);
            double dist2 = norm2(refPt - nodes[1]);
            double dist3 = norm2(refPt - nodes[2]);
            double tempVal = 0.0;
            int tempInt = 0;
            for (int i = 0; i < 3; i++)
                orders[i] = i;
            //
            if (dist1 > dist2) {
                tempInt = orders[0];
                orders[0] = orders[1];
                orders[1] = tempInt;
                //
                tempVal = dist1;
                dist1 = dist2;
                dist2 = dist1;
            }
            if (dist2 > dist3) {
                tempInt = orders[1];
                orders[1] = orders[2];
                orders[2] = tempInt;
                //
                tempVal = dist2;
                dist2 = dist3;
                dist3 = dist2;
            }
            if (dist1 > dist2) {
                tempInt = orders[0];
                orders[0] = orders[1];
                orders[1] = tempInt;
                //
                tempVal = dist1;
                dist1 = dist2;
                dist2 = dist1;
            }
        }
    }
}

void NonOverlappingCycles::TranslatePoint_onFace(const double dir, const double distance, const Vector3 inPt,
                                                 const int fid, Vector3 &outPt) {
    double scale = std::min(5.0, distance);
    //
    outPt = inPt + dir * scale * (*(inMeshPtr->meshNormalPtr))[fid];
    //
    return;
}

//void NonOverlappingCycles::WalkEdgeForNormal(const int eid, Vector3& outNormal)
//{
//	outNormal = Vector3(0., 0.0, 0.0);
//	outNormal = outNormal + (*(inMeshPtr->meshNormalPtr))[inMeshPtr->vecEdge[eid].AdjTri[0]];
//	if (inMeshPtr->vecEdge[eid].AdjTriNum == 2)
//	{
//		outNormal = outNormal + (*(inMeshPtr->meshNormalPtr))[inMeshPtr->vecEdge[eid].AdjTri[1]];
//		//outNormal = outNormal / norm(outNormal);
//		unitize(outNormal);
//
//	}
//
//	return;
//}
void NonOverlappingCycles::TranslatePoint_onEdge(const double dir, const double distance, const Vector3 inPt,
                                                 const int eid, Vector3 &outPt) {
    double scale = std::min(5.0, distance);
    //
    Vector3 edgeNormal;
    WalkEdgeForNormal(eid, edgeNormal);
    outPt = inPt + dir * scale * edgeNormal;
    //
    return;
}

void NonOverlappingCycles::CheckIntersectionOnTriangle(const double dir, const int critNode, const Vector3 &ptOnEdge,
                                                       const int reachedEdge, const int refEdgeVertex,
                                                       const int iterEdge, const int activeTriangleIdx,
                                                       const std::pair<int, int> &selfEdgePtType,
                                                       const Vector3 &selfEdgePt,
                                                       const std::pair<int, int> &objEdgePtType,
                                                       const Vector3 &objEdgePt,
                                                       const std::pair<int, int> &objPreEdgePtType,
                                                       const Vector3 &objPreEdgePt,
                                                       std::vector<Vector3> &outPts) {
    //
    Vector3 ptOnFace;
    Vector3 resIntersection;
    //
    double activeDist = 0.0;
    //
    bool bObjEdgeReached = false;
    bool bObjPreEdgeReached = false;
    bool bSelfEdgeReached = false;
    /*check if other two segments following at the same triangle or not*/
    if (objEdgePtType.second) {
        if (objEdgePtType.first == reachedEdge) {// this segment is on the same triangle
            bObjEdgeReached = true;
        }
    }
    if (objPreEdgePtType.second) {
        if (objPreEdgePtType.first == reachedEdge) {// this segment is also on the same triangle
            bObjPreEdgeReached = true;
        }
    }
    if (selfEdgePtType.second) {
        if (selfEdgePtType.first == reachedEdge) {
            bSelfEdgeReached = true;
        }
    }
    if (bObjEdgeReached || bObjPreEdgeReached || bSelfEdgeReached) {
        Vector3 critPt((inMeshPtr->vecVertex)[critNode].x,
                       (inMeshPtr->vecVertex)[critNode].y,
                       (inMeshPtr->vecVertex)[critNode].z);
        //
        int third_vertex = inMeshPtr->vecEdge[reachedEdge].v0 +
                           inMeshPtr->vecEdge[reachedEdge].v1 -
                           (refEdgeVertex);
        if (bSelfEdgeReached) {
            //no matter what it is,  just pull it into inside
            activeDist = 0.5;// norm(critPt - selfEdgePt) * 0.5;
            PointBetweenTwoPoints(activeDist, critPt, selfEdgePt, ptOnFace);
            if (bObjPreEdgeReached || bSelfEdgeReached)
                TranslatePoint_onFace(dir, norm(ptOnEdge - ptOnFace), ptOnFace, activeTriangleIdx, ptOnFace);
            outPts.push_back(ptOnFace);
            //if (bObjEdgeReached || bObjPreEdgeReached)
            //{
            //	if (bObjEdgeReached && bObjPreEdgeReached)
            //	{//
            //		std::vector<Vector3> ptsOnEdge(3);
            //		std::vector<int> orders(3);
            //		ptsOnEdge[0] = selfEdgePt;
            //		ptsOnEdge[1] = objPreEdgePt;
            //		ptsOnEdge[2] = objEdgePt;
            //		//
            //		OrderPointsReferingToOneEndPoint(thridVerPt, ptsOnEdge, orders);
            //		//
            //		switch(orders[0])
            //		{
            //		case 0:
            //			break;
            //		case 2:
            //			//
            //			std::cout << "check case 2 " << std::endl;
            //			if (orders[1] == 0)
            //			{//
            //				A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objPreEdgePt, resIntersection);
            //			}
            //			else
            //			{
            //				A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objEdgePt, resIntersection);
            //			}
            //			TranslatePoint_onFace(dir, norm(resIntersection - ptOnFace), resIntersection, activeTriangleIdx, resIntersection);
            //			outPts.push_back(resIntersection);
            //			//
            //			if (orders[1] == 1)
            //			{//
            //				A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objPreEdgePt, resIntersection);
            //			}
            //			else
            //			{
            //				A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objEdgePt, resIntersection);
            //			}
            //			TranslatePoint_onFace(dir, norm(resIntersection - ptOnFace), resIntersection, activeTriangleIdx, resIntersection);
            //			outPts.push_back(resIntersection);
            //			//
            //			break;
            //		case 1:							//
            //			std::cout << "check case 1 " << std::endl;
            //			// do the intersection
            //			if (orders[1] == 0)
            //			{//
            //				A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objPreEdgePt, resIntersection);
            //			}
            //			else
            //			{
            //				A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objEdgePt, resIntersection);
            //			}
            //			//
            //			TranslatePoint_onFace(dir, norm(resIntersection - ptOnFace), resIntersection, activeTriangleIdx, resIntersection);
            //			outPts.push_back(resIntersection);
            //			//
            //				break;
            //		default:
            //			break;
            //		}
            //	}
            //	else
            //	{
            //		std::vector<Vector3> ptsOnEdge(2);
            //		std::vector<int> orders(2);
            //		ptsOnEdge[0] = selfEdgePt;
            //		ptsOnEdge[1] = objPreEdgePt;
            //		//
            //		if (bObjEdgeReached)
            //		{
            //			ptsOnEdge[1] = objEdgePt;
            //		}
            //		OrderPointsReferingToOneEndPoint(thridVerPt, ptsOnEdge, orders);
            //		//
            //		std::cout << "check case 2.1 " << std::endl;
            //		if (orders[0] == 1)
            //		{
            //			A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &ptsOnEdge[1], resIntersection);
            //			TranslatePoint_onFace(dir, norm(resIntersection - ptOnFace), resIntersection, activeTriangleIdx, resIntersection);
            //			outPts.push_back(resIntersection);
            //		}
            //	}
            //
            //}
            //outPts.push_back(ptOnFace);
        } else {// always pull only one
            if (bObjEdgeReached) {
                activeDist = 0.5;// norm(critPt - selfEdgePt) * 0.5;
                PointBetweenTwoPoints(activeDist, critPt, objEdgePt, ptOnFace);
                TranslatePoint_onFace(dir, norm(critPt - ptOnFace), ptOnFace, activeTriangleIdx, ptOnFace);
            } else {
                activeDist = 0.5;// norm(critPt - selfEdgePt) * 0.5;
                PointBetweenTwoPoints(activeDist, critPt, objPreEdgePt, ptOnFace);
                TranslatePoint_onFace(dir, norm(critPt - ptOnFace), ptOnFace, activeTriangleIdx, ptOnFace);
            }
            outPts.push_back(ptOnFace);
            //
            //if (bObjEdgeReached || bObjPreEdgeReached)
            //{// predict next edge
            //// set the values for ptOnFace
            //	//
            //	if ((!objEdgePtType.second && objEdgePtType.first == refEdgeVertex) ||
            //		(!objPreEdgePtType.second && objPreEdgePtType.first == refEdgeVertex) )
            //	{
            //		Vector3 refPt(	(inMeshPtr->vecVertex)[refEdgeVertex].x,
            //						(inMeshPtr->vecVertex)[refEdgeVertex].y,
            //						(inMeshPtr->vecVertex)[refEdgeVertex].z);
            //		ptOnFace = critPt + 0.3 * (refPt - critPt);
            //	}
            //	else
            //	{// it is not occupied edge
            //		if (!selfEdgePtType.second && selfEdgePtType.first == refEdgeVertex)
            //		{
            //			Vector3 refPt(	(inMeshPtr->vecVertex)[refEdgeVertex].x,
            //							(inMeshPtr->vecVertex)[refEdgeVertex].y,
            //							(inMeshPtr->vecVertex)[refEdgeVertex].z);
            //			ptOnFace = critPt + 0.3 * (refPt - critPt);
            //		}
            //		else
            //		{
            //			activeDist = smallestDisatnce(critNode, iterEdge, resIntersection) * 0.3;
            //			//
            //			ptOnFace = critPt + 0.3 * (resIntersection - critPt);
            //		}
            //	}
            //	//
            //	if (bObjEdgeReached && bObjPreEdgeReached)
            //	{
            //		std::vector<Vector3> ptsOnEdge(2);
            //		std::vector<int> orders(2);
            //		ptsOnEdge[0] = objPreEdgePt;
            //		ptsOnEdge[1] = objEdgePt;
            //		//
            //		OrderPointsReferingToOneEndPoint(thridVerPt, ptsOnEdge, orders);
            //		//
            //		std::cout << "check case bObjEdgeReached && bObjPreEdgeReached " << std::endl;
            //		if (orders[0] == 0)
            //		{//
            //			A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objPreEdgePt, resIntersection);
            //		}
            //		else
            //		{
            //			A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objEdgePt, resIntersection);
            //		}
            //		TranslatePoint_onFace(dir, norm(resIntersection - ptOnFace), resIntersection, activeTriangleIdx, resIntersection);
            //		outPts.push_back(resIntersection);
            //		//
            //		if (orders[0] == 1)
            //		{//
            //			A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objPreEdgePt, resIntersection);
            //		}
            //		else
            //		{
            //			A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objEdgePt, resIntersection);
            //		}
            //		TranslatePoint_onFace(dir, norm(resIntersection - ptOnFace), resIntersection, activeTriangleIdx, resIntersection);
            //		outPts.push_back(resIntersection);
            //		//
            //	}
            //	else
            //	{
            //		std::cout << "check case else cae of bObjEdgeReached && bObjPreEdgeReached " << std::endl;
            //		if (bObjEdgeReached)
            //		{
            //			A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objEdgePt, resIntersection);
            //		}
            //		else
            //		{
            //			A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objPreEdgePt, resIntersection);
            //		}
            //		TranslatePoint_onFace(dir, norm(resIntersection - ptOnFace), resIntersection, activeTriangleIdx, resIntersection);
            //		outPts.push_back(resIntersection);
            //	}
            //}
        }
    }
    return;
}

void NonOverlappingCycles::WalkAroundVertexToMakeIntersectionAtFace_edgeCase(const double dir, const int critNode,
                                                                             const std::pair<int, int> &selfEdgePtType,
                                                                             const Vector3 &selfEdgePt,
                                                                             const std::pair<int, int> &selfPreEdgePtType,
                                                                             const Vector3 &selfPreEdgePt,
                                                                             const std::pair<int, int> &objEdgePtType,
                                                                             const Vector3 &objEdgePt,
                                                                             const std::pair<int, int> &objPreEdgePtType,
                                                                             const Vector3 &objPreEdgePt,
                                                                             std::vector<Vector3> &outPts) {
    //--selfPreEdge--critNode---selfEdge
    // at least one of selfPreEdge or selfEdge is an edge on the mesh
    const int selfEdge = selfEdgePtType.first;
    const int selfPreEdge = selfPreEdgePtType.first;
    const int objEdge = objEdgePtType.first;
    const int objPreEdge = objPreEdgePtType.first;
    std::pair<int, bool> nbSharedTriangle(0, false);
    if (selfEdgePtType.second == selfPreEdgePtType.second) {// both are vertices
        if (selfEdgePtType.first == selfPreEdgePtType.first) {
            std::cout << "go and come on the same edge " << std::endl;
            exit(3);
        }
        std::pair<int, bool> tempPair = inMeshPtr->EdgeIndex(selfPreEdgePtType.first, selfEdgePtType.first);
        int tempTrinagle = 0;
        if (tempPair.second) {
            nbSharedTriangle.second = true;
            nbSharedTriangle.first = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(tempPair.first, critNode);
        }
    } else {
        std::pair<int, bool> tempPair;
        if (selfPreEdgePtType.second) {
            nbSharedTriangle.first = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(selfPreEdge, critNode);
            tempPair = inMeshPtr->EdgeIndex(selfEdge, critNode);
        } else {
            nbSharedTriangle.first = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(selfEdge, critNode);
            tempPair = inMeshPtr->EdgeIndex(selfPreEdge, critNode);
        }
        //
        if (inMeshPtr->EdgeOnTriangle(tempPair.first, nbSharedTriangle.first)) {
            nbSharedTriangle.second = true;
        }
    }
//
    Vector3 ptOnFace;
    Vector3 ptOnEdge;
    Vector3 resIntersection;
    //
    double activeDist = 0.0;
    int activeTriangleIdx = 0;
    Vector3 critPt((inMeshPtr->vecVertex)[critNode].x,
                   (inMeshPtr->vecVertex)[critNode].y,
                   (inMeshPtr->vecVertex)[critNode].z);
    //
    if (nbSharedTriangle.second) {// incoming and outgoing edges are on the same triangle
        //
        if (!selfPreEdgePtType.second) {
            Vector3 refPt((inMeshPtr->vecVertex)[selfPreEdgePtType.first].x,
                          (inMeshPtr->vecVertex)[selfPreEdgePtType.first].y,
                          (inMeshPtr->vecVertex)[selfPreEdgePtType.first].z);
            if ((!objEdgePtType.second && objEdgePtType.first == selfPreEdgePtType.first) ||
                (!objPreEdgePtType.second && objPreEdgePtType.first == selfPreEdgePtType.first)) {// need to pull
                std::pair<int, bool> edge_idx = inMeshPtr->EdgeIndex(critNode, selfPreEdgePtType.first);
                activeDist = norm(critPt - refPt) * 0.3;
                ptOnEdge = critPt + 0.3 * (refPt - critPt);
                //
                TranslatePoint_onEdge(dir, activeDist, ptOnEdge, edge_idx.first, ptOnEdge);
                //
                outPts.push_back(ptOnEdge);
            } else {
                ptOnEdge = critPt + 0.3 * (refPt - critPt);
                //
                outPts.push_back(ptOnEdge);
            }
        } else {
            activeTriangleIdx = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(selfPreEdgePtType.first, critNode);
            activeDist = norm(critPt - selfPreEdgePt) * 0.5;
            ptOnFace = critPt + 0.5 * (selfPreEdgePt - critPt);
            //
            TranslatePoint_onFace(dir, activeDist, ptOnFace, activeTriangleIdx, ptOnFace);
            //
            outPts.push_back(ptOnFace);
        }
        //
        if (!selfEdgePtType.second) {
            Vector3 refPt((inMeshPtr->vecVertex)[selfEdgePtType.first].x,
                          (inMeshPtr->vecVertex)[selfEdgePtType.first].y,
                          (inMeshPtr->vecVertex)[selfEdgePtType.first].z);
            if ((!objEdgePtType.second && objEdgePtType.first == selfEdgePtType.first) ||
                (!objPreEdgePtType.second && objPreEdgePtType.first == selfEdgePtType.first)) {// need to pull
                std::pair<int, bool> edge_idx = inMeshPtr->EdgeIndex(critNode, selfEdgePtType.first);
                activeDist = norm(critPt - refPt) * 0.3;
                ptOnEdge = critPt + 0.3 * (refPt - critPt);
                //
                TranslatePoint_onEdge(dir, activeDist, ptOnEdge, edge_idx.first, ptOnEdge);
                //
                outPts.push_back(ptOnEdge);
            } else {
                ptOnEdge = critPt + 0.3 * (refPt - critPt);
                //
                outPts.push_back(ptOnEdge);
            }
        } else {
            activeTriangleIdx = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(selfEdgePtType.first, critNode);
            activeDist = norm(critPt - selfEdgePt) * 0.5;
            ptOnFace = critPt + 0.5 * (selfEdgePt - critPt);
            //
            TranslatePoint_onFace(dir, activeDist, ptOnFace, activeTriangleIdx, ptOnFace);
            //
            outPts.push_back(ptOnFace);
        }
    } else {
        int iterEdge = 0;
        int reachedEdge = 0;
        int refEdgeVertex = 0;
        int third_vertex = 0;
        if (!selfPreEdgePtType.second) {
            refEdgeVertex = selfPreEdgePtType.first;
            std::pair<int, bool> ret = inMeshPtr->EdgeIndex(critNode, refEdgeVertex);
            iterEdge = ret.first;
            activeTriangleIdx = inMeshPtr->vecEdge[iterEdge].AdjTri[0];
            //
            third_vertex = (inMeshPtr->vecTriangle)[activeTriangleIdx].v0 +
                           (inMeshPtr->vecTriangle)[activeTriangleIdx].v1 +
                           (inMeshPtr->vecTriangle)[activeTriangleIdx].v2 -
                           ((inMeshPtr->vecEdge)[iterEdge].v0 + (inMeshPtr->vecEdge)[iterEdge].v1);
            reachedEdge = inMeshPtr->EdgeConnectingTwoVertices(third_vertex, refEdgeVertex, iterEdge,
                                                               activeTriangleIdx);
        } else {
            // preprocessing
            //
            activeTriangleIdx = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(selfPreEdgePtType.first, critNode);
            //
            int nObjPreprocessingType = 0;
            int nObjPrePreprocessingType = 0;
            /*check if other two segments following at the same triangle or not*/
            if (objEdgePtType.second) {
                if (objEdge == selfPreEdge) {// this segment is on the same triangle
                    nObjPreprocessingType = 1;
                }
            } else {
                if ((inMeshPtr->vecEdge)[selfPreEdge].v0 == objEdge ||
                    (inMeshPtr->vecEdge)[selfPreEdge].v1 == objEdge) {// this referred edge is on the active triangle
                    nObjPreprocessingType = 2;
                }
            }
            if (objPreEdgePtType.second) {
                if (objPreEdge == selfPreEdge) {// this segment is also on the same triangle
                    nObjPrePreprocessingType = 1;
                }
            } else {
                if ((inMeshPtr->vecEdge)[selfPreEdge].v0 == objPreEdge ||
                    (inMeshPtr->vecEdge)[selfPreEdge].v1 == objPreEdge) {// this referred edge is on the active triangle
                    nObjPrePreprocessingType = 2;
                }
            }
            //setting refEdgeVertex
            refEdgeVertex = (inMeshPtr->vecEdge)[selfPreEdge].v0;
            //
            if (nObjPreprocessingType + nObjPrePreprocessingType > 0) {
                activeDist = 0.5;//leftDist * 0.5;
                PointBetweenTwoPoints(activeDist, critPt, selfPreEdgePt, ptOnFace);
                //
                TranslatePoint_onFace(dir, norm(critPt - selfPreEdgePt), ptOnFace, activeTriangleIdx, ptOnFace);
                //
                outPts.push_back(ptOnFace);
                //
                inMeshPtr->EdgeConnectingTwoVertices(critNode, refEdgeVertex, selfPreEdge, activeTriangleIdx, iterEdge);
            } else {// use default direction
                activeDist = 0.5;//leftDist * 0.5;
                PointBetweenTwoPoints(activeDist, critPt, selfPreEdgePt, ptOnFace);
                outPts.push_back(ptOnFace);
                //
                inMeshPtr->EdgeConnectingTwoVertices(critNode, refEdgeVertex, selfPreEdge, activeTriangleIdx, iterEdge);
            }
            //
            third_vertex = (inMeshPtr->vecTriangle)[activeTriangleIdx].v0 +
                           (inMeshPtr->vecTriangle)[activeTriangleIdx].v1 +
                           (inMeshPtr->vecTriangle)[activeTriangleIdx].v2 -
                           ((inMeshPtr->vecEdge)[iterEdge].v0 + (inMeshPtr->vecEdge)[iterEdge].v1);
            reachedEdge = selfPreEdgePtType.first;
        }
        int ver_terminate_edge = -1;
        if (!selfEdgePtType.second) {
            std::pair<int, bool> ret = inMeshPtr->EdgeIndex(critNode, selfEdgePtType.first);
            ver_terminate_edge = ret.first;
        }
        // just walk

        while (1) {// critNode--iterEdge---refEdgeVertex
            //   \					  /
            //      \third_vertex-- reachedEdge
            // compute the a point on the iterate edge
            //
            if ((!objEdgePtType.second && objEdge == refEdgeVertex) ||
                (!objPreEdgePtType.second && objPreEdge == refEdgeVertex)) {// this is intersection of edges
                Vector3 refEdgePt((inMeshPtr->vecVertex)[refEdgeVertex].x,
                                  (inMeshPtr->vecVertex)[refEdgeVertex].y,
                                  (inMeshPtr->vecVertex)[refEdgeVertex].z);
                //
                activeDist = norm(refEdgePt - critPt) * 0.3;
                PointBetweenTwoPoints(0.3, critPt, refEdgePt, ptOnEdge);
                // translate this edge intersection
                TranslatePoint_onEdge(dir, activeDist, ptOnEdge, iterEdge, ptOnEdge);
                //
                outPts.push_back(ptOnEdge);
            } else {
                activeDist = smallestDisatnce(critNode, iterEdge, ptOnEdge);
                PointBetweenTwoPoints(0.5, critPt, ptOnEdge, ptOnEdge);
                outPts.push_back(ptOnEdge);
            }
            //
            activeTriangleIdx = (inMeshPtr->vecEdge)[iterEdge].AdjTri[0] + (inMeshPtr->vecEdge)[iterEdge].AdjTri[1] -
                                activeTriangleIdx;
            //
            third_vertex = (inMeshPtr->vecTriangle)[activeTriangleIdx].v0 +
                           (inMeshPtr->vecTriangle)[activeTriangleIdx].v1 +
                           (inMeshPtr->vecTriangle)[activeTriangleIdx].v2 -
                           ((inMeshPtr->vecEdge)[iterEdge].v0 + (inMeshPtr->vecEdge)[iterEdge].v1);
            //
            //
            inMeshPtr->EdgeConnectingTwoVertices(third_vertex, refEdgeVertex, iterEdge, activeTriangleIdx, reachedEdge);
            inMeshPtr->EdgeConnectingTwoVertices(third_vertex, critNode, iterEdge, activeTriangleIdx, iterEdge);
            //
            refEdgeVertex = third_vertex;
            //
            CheckIntersectionOnTriangle(dir, critNode, ptOnEdge, reachedEdge, refEdgeVertex, iterEdge,
                                        activeTriangleIdx,
                                        selfEdgePtType, selfEdgePt, objEdgePtType, objEdgePt, objPreEdgePtType,
                                        objPreEdgePt,
                                        outPts);
            if ((!selfEdgePtType.second && iterEdge == ver_terminate_edge) ||
                (selfEdgePtType.second && reachedEdge == selfEdgePtType.first)
                    ) {
                break;
            }
        }
        if (!selfEdgePtType.second && iterEdge == ver_terminate_edge) {
            if ((!objEdgePtType.second && objEdge == refEdgeVertex) ||
                (!objPreEdgePtType.second && objPreEdge == refEdgeVertex)) {// this is intersection of edges
                Vector3 refEdgePt((inMeshPtr->vecVertex)[refEdgeVertex].x,
                                  (inMeshPtr->vecVertex)[refEdgeVertex].y,
                                  (inMeshPtr->vecVertex)[refEdgeVertex].z);
                //
                activeDist = norm(refEdgePt - critPt) * 0.3;
                PointBetweenTwoPoints(0.3, critPt, refEdgePt, ptOnEdge);
                // translate this edge intersection
                TranslatePoint_onEdge(dir, activeDist, ptOnEdge, iterEdge, ptOnEdge);
                //
                outPts.push_back(ptOnEdge);
            } else {
                Vector3 refEdgePt((inMeshPtr->vecVertex)[refEdgeVertex].x,
                                  (inMeshPtr->vecVertex)[refEdgeVertex].y,
                                  (inMeshPtr->vecVertex)[refEdgeVertex].z);
                PointBetweenTwoPoints(0.3, critPt, refEdgePt, ptOnEdge);
                outPts.push_back(ptOnEdge);
            }

        }
    }// not shared the same triangle


    return;
}

void NonOverlappingCycles::WalkAroundVertexToMakeIntersectionAtFace(const double dir, const int critNode,
                                                                    const std::pair<int, int> &selfEdgePtType,
                                                                    const Vector3 &selfEdgePt,
                                                                    const std::pair<int, int> &selfPreEdgePtType,
                                                                    const Vector3 &selfPreEdgePt,
                                                                    const std::pair<int, int> &objEdgePtType,
                                                                    const Vector3 &objEdgePt,
                                                                    const std::pair<int, int> &objPreEdgePtType,
                                                                    const Vector3 &objPreEdgePt,
                                                                    std::vector<Vector3> &outPts) {
    //--selfPreEdge--critNode---selfEdge
    const int selfEdge = selfEdgePtType.first;
    const int selfPreEdge = selfPreEdgePtType.first;
    const int objEdge = objEdgePtType.first;
    const int objPreEdge = objPreEdgePtType.first;
    if (selfEdge == selfPreEdge) {// incoming and outgoing edges are on the same triangle
        Vector3 critPt((inMeshPtr->vecVertex)[critNode].x,
                       (inMeshPtr->vecVertex)[critNode].y,
                       (inMeshPtr->vecVertex)[critNode].z);
        int activeTriangleIdx = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(selfPreEdge, critNode);
        //
        //
        double activeDist = 0.5;//std::min(leftDist, rightDist) / 2.0;
        //
        Vector3 translatedPt;
        PointBetweenTwoPoints(activeDist, critPt, selfPreEdgePt, translatedPt);
        TranslatePoint_onFace(dir, norm(critPt - translatedPt), translatedPt, activeTriangleIdx, translatedPt);
        outPts.push_back(translatedPt);
        //
        PointBetweenTwoPoints(activeDist, critPt, selfEdgePt, translatedPt);
        TranslatePoint_onFace(dir, norm(critPt - translatedPt), translatedPt, activeTriangleIdx, translatedPt);
        outPts.push_back(translatedPt);
    } else {// need to walk around the edge

        // look at different cases
        // std::cout << "before prefpro face " << std::endl;
        Vector3 critPt((inMeshPtr->vecVertex)[critNode].x,
                       (inMeshPtr->vecVertex)[critNode].y,
                       (inMeshPtr->vecVertex)[critNode].z);
        //
        Vector3 ptOnFace;
        Vector3 ptOnEdge;
        Vector3 resIntersection;
        //
        double activeDist = 0.0;
        //
        int activeTriangleIdx = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(selfPreEdge, critNode);
        int iterEdge = 0;
        // preprocessing
        //
        int nObjPreprocessingType = 0;
        int nObjPrePreprocessingType = 0;
        /*check if other two segments following at the same triangle or not*/
        if (objEdgePtType.second) {
            if (objEdge == selfPreEdge) {// this segment is on the same triangle
                nObjPreprocessingType = 1;
            }
        } else {
            if ((inMeshPtr->vecEdge)[selfPreEdge].v0 == objEdge ||
                (inMeshPtr->vecEdge)[selfPreEdge].v1 == objEdge) {// this referred edge is on the active triangle
                nObjPreprocessingType = 2;
            }
        }
        if (objPreEdgePtType.second) {
            if (objPreEdge == selfPreEdge) {// this segment is also on the same triangle
                nObjPrePreprocessingType = 1;
            }
        } else {
            if ((inMeshPtr->vecEdge)[selfPreEdge].v0 == objPreEdge ||
                (inMeshPtr->vecEdge)[selfPreEdge].v1 == objPreEdge) {// this referred edge is on the active triangle
                nObjPrePreprocessingType = 2;
            }
        }
        int refEdgeVertex = (inMeshPtr->vecEdge)[selfPreEdge].v0;
        //Vector3 refEdgePt((	inMeshPtr->vecVertex)[refEdgeVertex].x,
        //					(inMeshPtr->vecVertex)[refEdgeVertex].y,
        //					(inMeshPtr->vecVertex)[refEdgeVertex].z);
        //
        //int leftEdgeIdx = 0;
        //int rightEdgeIdx = 0;
        ////
        //inMeshPtr->EdgeConnectingTwoVertices(critNode, (inMeshPtr->vecEdge)[selfPreEdge].v0,
        //	selfPreEdge, activeTriangleIdx, leftEdgeIdx);
        ////
        //rightEdgeIdx = (inMeshPtr->vecTriangle)[activeTriangleIdx].e01 +
        //	(inMeshPtr->vecTriangle)[activeTriangleIdx].e02 + (inMeshPtr->vecTriangle)[activeTriangleIdx].e12 -
        //	(selfPreEdge + leftEdgeIdx);
        //	//
        //if ((inMeshPtr->vecTriangle)[activeTriangleIdx].e01 != rightEdgeIdx &&
        //	(inMeshPtr->vecTriangle)[activeTriangleIdx].e02 != rightEdgeIdx &&
        //	(inMeshPtr->vecTriangle)[activeTriangleIdx].e12 != rightEdgeIdx)
        //{
        //	std::cout << "face  preprocessing; right edge is not on the face " << std::endl;
        //	exit(0);
        //}
        //if ((inMeshPtr->vecTriangle)[activeTriangleIdx].e01 != leftEdgeIdx &&
        //	(inMeshPtr->vecTriangle)[activeTriangleIdx].e02 != leftEdgeIdx &&
        //	(inMeshPtr->vecTriangle)[activeTriangleIdx].e12 != leftEdgeIdx)
        //{
        //	std::cout << "face  preprocessing; right edge is not on the face " << std::endl;
        //	exit(0);
        //}
        //if ((inMeshPtr->vecTriangle)[activeTriangleIdx].e01 != selfPreEdge &&
        //	(inMeshPtr->vecTriangle)[activeTriangleIdx].e02 != selfPreEdge &&
        //	(inMeshPtr->vecTriangle)[activeTriangleIdx].e12 != selfPreEdge)
        //{
        //	std::cout << "face  preprocessing; right edge is not on the face " << std::endl;
        //	exit(0);
        //}
        //double leftDist = 0.0;
        //Vector3 leftPt;
        //double rightDist = 0.0;
        //Vector3 rightPt;
        //
        if (nObjPreprocessingType + nObjPrePreprocessingType > 0) {
            //if (nObjPreprocessingType > 0 && nObjPrePreprocessingType > 0)
            //{// both are on the
            //	std::vector<Vector3> ptsOnEdge(3);
            //	std::vector<int> orders(3);
            //	ptsOnEdge[0] = selfPreEdgePt;
            //	ptsOnEdge[1] = objPreEdgePt;
            //	ptsOnEdge[2] = objEdgePt;
            //	//
            //	OrderPointsReferingToOneEndPoint(refEdgePt, ptsOnEdge, orders);
            //	//
            //	switch(orders[0])
            //	{
            //	case 0:
            //		refEdgeVertex =  initialRefVertex;
            //		//leftDist = norm(critPt - selfPreEdgePt);//smallestDisatnce(critNode, leftEdgeIdx, leftPt);
            //		activeDist = 0.5;
            //		PointBetweenTwoPoints(activeDist, critPt, selfPreEdgePt, ptOnFace);
            //		outPts.push_back(ptOnFace);
            //		break;
            //	case 2:
            //		refEdgeVertex =  (inMeshPtr->vecEdge)[selfPreEdge].v1;
            //		//rightDist =  norm(critPt - selfPreEdgePt);// smallestDisatnce(critNode, rightEdgeIdx, rightPt);
            //		activeDist = 0.5;//rightDist * 0.5;
            //		PointBetweenTwoPoints(activeDist, critPt, selfPreEdgePt, ptOnFace);
            //		outPts.push_back(ptOnFace);
            //		break;
            //	case 1:
            //		//
            //		//leftDist = norm(critPt - selfPreEdgePt);//smallestDisatnce(critNode, leftEdgeIdx, leftPt);
            //		activeDist = 0.5;//leftDist * 0.5;
            //		PointBetweenTwoPoints(activeDist, critPt, selfPreEdgePt, ptOnFace);
            //		outPts.push_back(ptOnFace);
            //		// do the intersection
            //		if (nObjPreprocessingType == 1 || nObjPrePreprocessingType == 1)
            //		{// prefer the surface intersection
            //			if (nObjPreprocessingType == 1)
            //			{// it is a surface intersection
            //				//if (nObjPrePreprocessingType == 2 && objPreEdge == (inMeshPtr->vecEdge)[selfPreEdge].v0)
            //				//{
            //				//	refEdgeVertex =  (inMeshPtr->vecEdge)[selfPreEdge].v1;
            //				//	rightDist = smallestDisatnce(critNode, rightEdgeIdx, rightPt);
            //				//	PointBetweenTwoPoints(0.5, critPt, rightPt, ptOnEdge);
            //				//}
            //				//else
            //				//{
            //				//	leftDist = smallestDisatnce(critNode, leftEdgeIdx, leftPt);
            //				//	PointBetweenTwoPoints(0.5, critPt, leftPt, ptOnEdge);
            //				//}
            //				//
            //				std::cout << "case I " << std::endl;
            //				if (orders[2] == 2)
            //				{
            //					refEdgeVertex =  (inMeshPtr->vecEdge)[selfPreEdge].v1;
            //					rightDist = smallestDisatnce(critNode, rightEdgeIdx, rightPt);
            //					PointBetweenTwoPoints(0.5, critPt, rightPt, ptOnEdge);
            //				}
            //				else
            //				{
            //					leftDist = smallestDisatnce(critNode, leftEdgeIdx, leftPt);
            //					PointBetweenTwoPoints(0.5, critPt, leftPt, ptOnEdge);
            //				}
            //				//
            //				A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objEdgePt, resIntersection);
            //				std::cout << "fater I " << std::endl;
            //			}
            //			else
            //			{
            //				//if (objEdge == (inMeshPtr->vecEdge)[selfPreEdge].v0)
            //				//{
            //				//	refEdgeVertex =  (inMeshPtr->vecEdge)[selfPreEdge].v1;
            //				//	rightDist = smallestDisatnce(critNode, rightEdgeIdx, rightPt);
            //				//	PointBetweenTwoPoints(0.5, critPt, rightPt, ptOnEdge);
            //				//}
            //				//else
            //				//{
            //				//	leftDist = smallestDisatnce(critNode, leftEdgeIdx, leftPt);
            //				//	PointBetweenTwoPoints(0.5, critPt, leftPt, ptOnEdge);
            //				//}
            //				//
            //				std::cout << "case II " << std::endl;
            //				if (orders[1] == 2)
            //				{
            //					refEdgeVertex =  (inMeshPtr->vecEdge)[selfPreEdge].v1;
            //					rightDist = smallestDisatnce(critNode, rightEdgeIdx, rightPt);
            //					PointBetweenTwoPoints(0.5, critPt, rightPt, ptOnEdge);
            //				}
            //				else
            //				{
            //					leftDist = smallestDisatnce(critNode, leftEdgeIdx, leftPt);
            //					PointBetweenTwoPoints(0.5, critPt, leftPt, ptOnEdge);
            //				}
            //				//
            //				A_Pair_3D_Line_Intersection(&ptOnFace, &ptOnEdge, &critPt, &objPreEdgePt, resIntersection);
            //				std::cout << "aft I " << std::endl;
            //			}
            //			TranslatePoint_onFace(dir, norm(resIntersection - ptOnFace), resIntersection, activeTriangleIdx, resIntersection);
            //			outPts.push_back(resIntersection);
            //		}
            //		// first compute the point along the segment near the critical point
            //		// compute another point on the refEdge
            //		// compute the intersection between this new segment and another edge
            //		 break;
            //	default:
            //		break;
            //	}
            //}
            //else
            //{// only has one common
            //	std::vector<Vector3> ptsOnEdge(2);
            //	std::vector<int> orders(2);
            //	ptsOnEdge[0] = selfPreEdgePt;
            //	if (nObjPreprocessingType > 0)
            //	{
            //		ptsOnEdge[1] = objEdgePt;
            //	}
            //	else
            //	{
            //		ptsOnEdge[1] = objPreEdgePt;
            //	}
            //	//
            //	//leftDist =  norm(critPt - selfPreEdgePt); //smallestDisatnce(critNode, leftEdgeIdx, leftPt);
            //	activeDist = 0.5;//leftDist *
            //	PointBetweenTwoPoints(activeDist, critPt, selfPreEdgePt, ptOnFace);
            //	outPts.push_back(ptOnFace);
            //	//
            //	OrderPointsReferingToOneEndPoint(refEdgePt, ptsOnEdge, orders);
            //	if (orders[0])
            //	{
            //		refEdgeVertex = (inMeshPtr->vecEdge)[selfPreEdge].v1;
            //	}
            //}
            //
            activeDist = 0.5;// leftDist * 0.5;
            PointBetweenTwoPoints(activeDist, critPt, selfPreEdgePt, ptOnFace);
            TranslatePoint_onFace(dir, norm(critPt - ptOnFace), ptOnFace, activeTriangleIdx, ptOnFace);
            outPts.push_back(ptOnFace);
            //
            inMeshPtr->EdgeConnectingTwoVertices(critNode, refEdgeVertex, selfPreEdge, activeTriangleIdx, iterEdge);
        } else {// use default direction
            activeDist = 0.5;// leftDist * 0.5;
            PointBetweenTwoPoints(activeDist, critPt, selfPreEdgePt, ptOnFace);
            outPts.push_back(ptOnFace);
            //
            inMeshPtr->EdgeConnectingTwoVertices(critNode, refEdgeVertex, selfPreEdge, activeTriangleIdx, iterEdge);
        }

        // just walk
        int reachedEdge = selfPreEdge;
        int third_vertex = 0;
        // std::cout << "before while " << std::endl;
        while (reachedEdge != selfEdge) {// critNode--iterEdge---refEdgeVertex
            //   \					  /
            //      \third_vertex-- reachedEdge
            // compute the a point on the iterate edge
            //
            if ((!objEdgePtType.second && objEdge == refEdgeVertex) ||
                (!objPreEdgePtType.second && objPreEdge == refEdgeVertex)) {// this is intersection of edges
                Vector3 refEdgePt((inMeshPtr->vecVertex)[refEdgeVertex].x,
                                  (inMeshPtr->vecVertex)[refEdgeVertex].y,
                                  (inMeshPtr->vecVertex)[refEdgeVertex].z);
                //
                activeDist = 0.3;
                PointBetweenTwoPoints(activeDist, critPt, refEdgePt, ptOnEdge);
                // translate this edge intersection
                TranslatePoint_onEdge(dir, activeDist, ptOnEdge, iterEdge, ptOnEdge);
                //
                outPts.push_back(ptOnEdge);
            } else {
                activeDist = smallestDisatnce(critNode, iterEdge, ptOnEdge);
                PointBetweenTwoPoints(0.5, critPt, ptOnEdge, ptOnEdge);
                outPts.push_back(ptOnEdge);
            }
            //
            activeTriangleIdx = (inMeshPtr->vecEdge)[iterEdge].AdjTri[0] + (inMeshPtr->vecEdge)[iterEdge].AdjTri[1] -
                                activeTriangleIdx;
            //
            third_vertex = (inMeshPtr->vecTriangle)[activeTriangleIdx].v0 +
                           (inMeshPtr->vecTriangle)[activeTriangleIdx].v1 +
                           (inMeshPtr->vecTriangle)[activeTriangleIdx].v2 -
                           ((inMeshPtr->vecEdge)[iterEdge].v0 + (inMeshPtr->vecEdge)[iterEdge].v1);
            //
            //
            inMeshPtr->EdgeConnectingTwoVertices(third_vertex, refEdgeVertex, iterEdge, activeTriangleIdx, reachedEdge);
            inMeshPtr->EdgeConnectingTwoVertices(third_vertex, critNode, iterEdge, activeTriangleIdx, iterEdge);
            //
            refEdgeVertex = third_vertex;
            //
            CheckIntersectionOnTriangle(dir, critNode, ptOnEdge, reachedEdge, refEdgeVertex, iterEdge,
                                        activeTriangleIdx,
                                        selfEdgePtType, selfEdgePt, objEdgePtType, objEdgePt, objPreEdgePtType,
                                        objPreEdgePt,
                                        outPts);
        }


    }
    return;
}

void NonOverlappingCycles::ComputeBasisLoopPolygon_1() {
    // make a mutual decision on _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
    int curArcId = 0;
    int startNodeId = 0;
    int endNodeId = 0;
    Vector3 translatedPt;
    Vector3 prePt, nxtPt;
    //
    std::pair<int, int> selfPreEdgePtType;
    std::pair<int, int> selfEdgePtType;
    std::pair<int, int> objPreEdgePtType;
    std::pair<int, int> objEdgePtType;
    //
    Vector3 selfPreEdgePt;
    Vector3 selfEdgePt;
    Vector3 objPreEdgePt;
    Vector3 objEdgePt;
    //
    double dir = 1.0;
    if (_basisLoopPtr->pathType) {// vertical, so pull inside
        dir = -1.0;
    }
    if (_basisLoopPtr == _loopOnMeshPtr) {// need to embed every critical node on the path
        for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
            curArcId = _basisLoopPtr->simpArcIndexSet[i];
            startNodeId = _basisLoopPtr->nodeIndexSet[i];
            endNodeId = _basisLoopPtr->nodeIndexSet[i + 1];
            //
            int prevArcId = i - 1;
            if (prevArcId < 0)
                prevArcId = _basisLoopPtr->simpArcIndexSet.size() - 1;
            //
            prevArcId = _basisLoopPtr->simpArcIndexSet[prevArcId];
            //
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction

                selfEdgePt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[(*_pathArcOnMeshPtr)[curArcId].vecPoints.size() -
                                                                      2];
                selfEdgePtType = (*_pathArcOnMeshPointTypePtr)[curArcId][
                        (*_pathArcOnMeshPointTypePtr)[curArcId].size() - 2];
                //
                objEdgePt = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[
                        (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 2];
                objEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[curArcId][
                        (*_offsetPathArcOnMeshPointTypePtr)[curArcId].size() - 2];
                //
                if ((*_pVecSimplifiedArc)[prevArcId].nCriticalNode0 == startNodeId) {
                    selfPreEdgePt = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[
                            (*_pathArcOnMeshPtr)[prevArcId].vecPoints.size() - 2];
                    selfPreEdgePtType = (*_pathArcOnMeshPointTypePtr)[prevArcId][
                            (*_pathArcOnMeshPointTypePtr)[prevArcId].size() - 2];
                    //
                    objPreEdgePt = (*_offsetPathArcOnMeshPtr)[prevArcId].vecPoints[
                            (*_offsetPathArcOnMeshPtr)[prevArcId].vecPoints.size() - 2];
                    objPreEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[prevArcId][
                            (*_offsetPathArcOnMeshPointTypePtr)[prevArcId].size() - 2];
                } else {
                    selfPreEdgePt = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[1];
                    selfPreEdgePtType = (*_pathArcOnMeshPointTypePtr)[prevArcId][1];
                    //
                    objPreEdgePt = (*_offsetPathArcOnMeshPtr)[prevArcId].vecPoints[1];
                    objPreEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[prevArcId][1];
                }
                //
                if (!selfPreEdgePtType.second || !selfEdgePtType.second) {
                    WalkAroundVertexToMakeIntersectionAtFace_edgeCase(dir, startNodeId,
                                                                      selfEdgePtType, selfEdgePt,
                                                                      selfPreEdgePtType, selfPreEdgePt,
                                                                      objEdgePtType, objEdgePt,
                                                                      objPreEdgePtType, objPreEdgePt,
                                                                      outBasisLoop.vecPoints);
                } else {
                    WalkAroundVertexToMakeIntersectionAtFace(dir, startNodeId,
                                                             selfEdgePtType, selfEdgePt,
                                                             selfPreEdgePtType, selfPreEdgePt,
                                                             objEdgePtType, objEdgePt,
                                                             objPreEdgePtType, objPreEdgePt,
                                                             outBasisLoop.vecPoints);
                }
                //
                for (int iarc = (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2;
                     iarc > 0; iarc--) {// ignore the first point
                    outBasisLoop.vecPoints.push_back(
                            (*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);//(*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
                /*for ( int iarc = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 2; iarc > 0 ; iarc--)
				{
					outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
				}*/
            } else {// ignore the last point
                selfEdgePt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[1];
                selfEdgePtType = (*_pathArcOnMeshPointTypePtr)[curArcId][1];
                //
                objEdgePt = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[1];
                objEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[curArcId][1];
                //
                if ((*_pVecSimplifiedArc)[prevArcId].nCriticalNode0 == startNodeId) {
                    selfPreEdgePt = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[
                            (*_pathArcOnMeshPtr)[prevArcId].vecPoints.size() - 2];
                    selfPreEdgePtType = (*_pathArcOnMeshPointTypePtr)[prevArcId][
                            (*_pathArcOnMeshPointTypePtr)[prevArcId].size() - 2];
                    //
                    objPreEdgePt = (*_offsetPathArcOnMeshPtr)[prevArcId].vecPoints[
                            (*_offsetPathArcOnMeshPtr)[prevArcId].vecPoints.size() - 2];
                    objPreEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[prevArcId][
                            (*_offsetPathArcOnMeshPointTypePtr)[prevArcId].size() - 2];
                } else {
                    selfPreEdgePt = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[1];
                    selfPreEdgePtType = (*_pathArcOnMeshPointTypePtr)[prevArcId][1];
                    //
                    objPreEdgePt = (*_offsetPathArcOnMeshPtr)[prevArcId].vecPoints[1];
                    objPreEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[prevArcId][1];
                }
                if (!selfPreEdgePtType.second || !selfEdgePtType.second) {
                    WalkAroundVertexToMakeIntersectionAtFace_edgeCase(dir, startNodeId,
                                                                      selfEdgePtType, selfEdgePt,
                                                                      selfPreEdgePtType, selfPreEdgePt,
                                                                      objEdgePtType, objEdgePt,
                                                                      objPreEdgePtType, objPreEdgePt,
                                                                      outBasisLoop.vecPoints);
                } else {
                    WalkAroundVertexToMakeIntersectionAtFace(dir, startNodeId,
                                                             selfEdgePtType, selfEdgePt,
                                                             selfPreEdgePtType, selfPreEdgePt,
                                                             objEdgePtType, objEdgePt,
                                                             objPreEdgePtType, objPreEdgePt,
                                                             outBasisLoop.vecPoints);
                }
                //
                //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[0]);
                for (unsigned int iarc = 1;
                     iarc < (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {// ignore the first point
                    outBasisLoop.vecPoints.push_back(
                            (*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);//(*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }/*
				for (unsigned int iarc = 1; iarc < (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 1 ; iarc++)
				{
					outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
				}*/

            }
        }
        // push the first point into the loop to close it
        //outLoopOnMesh.vecPoints.push_back(outLoopOnMesh.vecPoints[0]);
        //
        outBasisLoop.vecPoints.push_back(outBasisLoop.vecPoints[0]);
    } else {// only to embed the shared critical node

        std::set<int> shared_nodes;
        std::set<int> basis_nodes_set(_basisLoopPtr->nodeIndexSet.begin(), _basisLoopPtr->nodeIndexSet.end());
        std::set<int>::iterator findIter;
        std::map<int, int> shared_nodes_pos;
        //
        for (unsigned int j = 0; j < _loopOnMeshPtr->nodeIndexSet.size() - 1; j++) {
            findIter = basis_nodes_set.find(_loopOnMeshPtr->nodeIndexSet[j]);
            if (findIter != basis_nodes_set.end()) {
                shared_nodes.insert(_loopOnMeshPtr->nodeIndexSet[j]);
                shared_nodes_pos[_loopOnMeshPtr->nodeIndexSet[j]] = j;
            }
            //if (_basisLoopPtr->nodeIndexSet[i] == _loopOnMeshPtr->nodeIndexSet[j])
            //	shared_nodes.insert(_loopOnMeshPtr->nodeIndexSet[j]);
        }
        //for (unsigned int i = 0; i < _basisLoopPtr->nodeIndexSet.size(); i++)
        //{
        //	for (unsigned int j = 0; j < _loopOnMeshPtr->nodeIndexSet.size(); j++)
        //	{
        //		if (_basisLoopPtr->nodeIndexSet[i] == _loopOnMeshPtr->nodeIndexSet[j])
        //			shared_nodes.insert(_loopOnMeshPtr->nodeIndexSet[j]);
        //	}
        //}
        std::set<int>::iterator sIter;
        //
        for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
            curArcId = _basisLoopPtr->simpArcIndexSet[i];
            startNodeId = _basisLoopPtr->nodeIndexSet[i];
            endNodeId = _basisLoopPtr->nodeIndexSet[i + 1];
            //
            //
            int prevArcId = i - 1;
            if (prevArcId < 0)
                prevArcId = _basisLoopPtr->simpArcIndexSet.size() - 1;
            //
            prevArcId = _basisLoopPtr->simpArcIndexSet[prevArcId];
            //
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction
                sIter = shared_nodes.find(startNodeId);
                if (sIter != shared_nodes.end()) {// need to pull out or in
                    //
                    int orgArcId = shared_nodes_pos[startNodeId];
                    int orgPrevArcId = orgArcId - 1;
                    if (orgPrevArcId < 0) {
                        orgPrevArcId = _loopOnMeshPtr->simpArcIndexSet.size() - 1;
                    }
                    orgArcId = _loopOnMeshPtr->simpArcIndexSet[orgArcId];
                    orgPrevArcId = _loopOnMeshPtr->simpArcIndexSet[orgPrevArcId];
                    //
                    selfEdgePt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[
                            (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2];
                    selfEdgePtType = (*_pathArcOnMeshPointTypePtr)[curArcId][
                            (*_pathArcOnMeshPointTypePtr)[curArcId].size() - 2];
                    //
                    if ((*_pVecSimplifiedArc)[orgArcId].nCriticalNode0 == startNodeId) {
                        objEdgePt = (*_offsetPathArcOnMeshPtr)[orgArcId].vecPoints[
                                (*_offsetPathArcOnMeshPtr)[orgArcId].vecPoints.size() - 2];
                        objEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[orgArcId][
                                (*_offsetPathArcOnMeshPointTypePtr)[orgArcId].size() - 2];
                    } else {
                        objEdgePt = (*_offsetPathArcOnMeshPtr)[orgArcId].vecPoints[1];
                        objEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[orgArcId][1];
                    }
                    //
                    if ((*_pVecSimplifiedArc)[orgPrevArcId].nCriticalNode0 == startNodeId) {
                        objPreEdgePt = (*_offsetPathArcOnMeshPtr)[orgPrevArcId].vecPoints[
                                (*_offsetPathArcOnMeshPtr)[orgPrevArcId].vecPoints.size() - 2];
                        objPreEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[orgPrevArcId][
                                (*_offsetPathArcOnMeshPointTypePtr)[orgPrevArcId].size() - 2];
                    } else {
                        objPreEdgePt = (*_offsetPathArcOnMeshPtr)[orgPrevArcId].vecPoints[1];
                        objPreEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[orgPrevArcId][1];
                    }

                    //
                    if ((*_pVecSimplifiedArc)[prevArcId].nCriticalNode0 == startNodeId) {
                        selfPreEdgePt = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[
                                (*_pathArcOnMeshPtr)[prevArcId].vecPoints.size() - 2];
                        selfPreEdgePtType = (*_pathArcOnMeshPointTypePtr)[prevArcId][
                                (*_pathArcOnMeshPointTypePtr)[prevArcId].size() - 2];
                    } else {
                        selfPreEdgePt = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[1];
                        selfPreEdgePtType = (*_pathArcOnMeshPointTypePtr)[prevArcId][1];
                    }
                    //
                    if (!selfPreEdgePtType.second || !selfEdgePtType.second) {
                        WalkAroundVertexToMakeIntersectionAtFace_edgeCase(dir, startNodeId,
                                                                          selfEdgePtType, selfEdgePt,
                                                                          selfPreEdgePtType, selfPreEdgePt,
                                                                          objEdgePtType, objEdgePt,
                                                                          objPreEdgePtType, objPreEdgePt,
                                                                          outBasisLoop.vecPoints);
                    } else {
                        WalkAroundVertexToMakeIntersectionAtFace(dir, startNodeId,
                                                                 selfEdgePtType, selfEdgePt,
                                                                 selfPreEdgePtType, selfPreEdgePt,
                                                                 objEdgePtType, objEdgePt,
                                                                 objPreEdgePtType, objPreEdgePt,
                                                                 outBasisLoop.vecPoints);
                    }
                } else {
                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints.back());
                }
                //
                //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.back());
                //
                if (!(*_bVecBasisLoopPtr)[i]) {// if the arc is setted to be a straight segment, ignore all other points
                    for (int iarc = (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2;
                         iarc > 0; iarc--) {// ignore the first point
                        outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);

                    }/*
					for ( int iarc = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 2; iarc > 0 ; iarc--)
					{
						outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
					}*/
                }
            } else {// ignore the last point
                sIter = shared_nodes.find(startNodeId);
                if (sIter != shared_nodes.end()) {// need to pull out or in
                    int orgArcId = shared_nodes_pos[startNodeId];
                    int orgPrevArcId = orgArcId - 1;
                    if (orgPrevArcId < 0) {
                        orgPrevArcId = _loopOnMeshPtr->simpArcIndexSet.size() - 1;
                    }
                    orgArcId = _loopOnMeshPtr->simpArcIndexSet[orgArcId];
                    orgPrevArcId = _loopOnMeshPtr->simpArcIndexSet[orgPrevArcId];
                    //
                    selfEdgePt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[1];
                    selfEdgePtType = (*_pathArcOnMeshPointTypePtr)[curArcId][1];
                    //
                    if ((*_pVecSimplifiedArc)[orgArcId].nCriticalNode0 == startNodeId) {
                        objEdgePt = (*_offsetPathArcOnMeshPtr)[orgArcId].vecPoints[
                                (*_offsetPathArcOnMeshPtr)[orgArcId].vecPoints.size() - 2];
                        objEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[orgArcId][
                                (*_offsetPathArcOnMeshPointTypePtr)[orgArcId].size() - 2];
                    } else {
                        objEdgePt = (*_offsetPathArcOnMeshPtr)[orgArcId].vecPoints[1];
                        objEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[orgArcId][1];
                    }
                    if ((*_pVecSimplifiedArc)[orgPrevArcId].nCriticalNode0 == startNodeId) {
                        objPreEdgePt = (*_offsetPathArcOnMeshPtr)[orgPrevArcId].vecPoints[
                                (*_offsetPathArcOnMeshPtr)[orgPrevArcId].vecPoints.size() - 2];
                        objPreEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[orgPrevArcId][
                                (*_offsetPathArcOnMeshPointTypePtr)[orgPrevArcId].size() - 2];
                    } else {
                        objPreEdgePt = (*_offsetPathArcOnMeshPtr)[orgPrevArcId].vecPoints[1];
                        objPreEdgePtType = (*_offsetPathArcOnMeshPointTypePtr)[orgPrevArcId][1];
                    }

                    //
                    if ((*_pVecSimplifiedArc)[prevArcId].nCriticalNode0 == startNodeId) {
                        selfPreEdgePt = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[
                                (*_pathArcOnMeshPtr)[prevArcId].vecPoints.size() - 2];
                        selfPreEdgePtType = (*_pathArcOnMeshPointTypePtr)[prevArcId][
                                (*_pathArcOnMeshPointTypePtr)[prevArcId].size() - 2];
                    } else {
                        selfPreEdgePt = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[1];
                        selfPreEdgePtType = (*_pathArcOnMeshPointTypePtr)[prevArcId][1];
                    }
                    //
                    if (!selfPreEdgePtType.second || !selfEdgePtType.second) {
                        WalkAroundVertexToMakeIntersectionAtFace_edgeCase(dir, startNodeId,
                                                                          selfEdgePtType, selfEdgePt,
                                                                          selfPreEdgePtType, selfPreEdgePt,
                                                                          objEdgePtType, objEdgePt,
                                                                          objPreEdgePtType, objPreEdgePt,
                                                                          outBasisLoop.vecPoints);
                    } else {
                        WalkAroundVertexToMakeIntersectionAtFace(dir, startNodeId,
                                                                 selfEdgePtType, selfEdgePt,
                                                                 selfPreEdgePtType, selfPreEdgePt,
                                                                 objEdgePtType, objEdgePt,
                                                                 objPreEdgePtType, objPreEdgePt,
                                                                 outBasisLoop.vecPoints);
                    }
                } else {
                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints.front());
                }
                //
                if (!(*_bVecBasisLoopPtr)[i]) {// if the arc is setted to be a straight segment, ignore all other points
                    //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.front());
                    for (unsigned int iarc = 1;
                         iarc < (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {// ignore the first point
                        outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);

                    }/*
					for (unsigned int iarc = 1; iarc < (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 1 ; iarc++)
					{
						outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
					}*/
                }
            }
        }
        // push the first point into the loop to close it
        //outLoopOnMesh.vecPoints.push_back(outLoopOnMesh.vecPoints[0]);
        outBasisLoop.vecPoints.push_back(outBasisLoop.vecPoints[0]);

    }

}

void NonOverlappingCycles::ComputeBasisLoopPolygon() {
    // make a mutual decision on _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
    int curArcId = 0;
    int startNodeId = 0;
    int endNodeId = 0;
    Vector3 translatedPt;
    Vector3 prePt, nxtPt;
    if (_basisLoopPtr == _loopOnMeshPtr) {// need to embed every critical node on the path
        for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
            curArcId = _basisLoopPtr->simpArcIndexSet[i];
            startNodeId = _basisLoopPtr->nodeIndexSet[i];
            endNodeId = _basisLoopPtr->nodeIndexSet[i + 1];
            //
            int prevArcId = i - 1;
            if (prevArcId < 0)
                prevArcId = _basisLoopPtr->simpArcIndexSet.size() - 1;
            //
            prevArcId = _basisLoopPtr->simpArcIndexSet[prevArcId];
            //
            if (i == 0) {
                if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet.back()].nCriticalNode0 == startNodeId) {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[
                            (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints.size() - 2];
                } else {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[1];
                }
            } else {
                if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[i - 1]].nCriticalNode0 == startNodeId) {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[
                            (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints.size() - 2];
                } else {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[1];
                }
            }
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction
                double dir = 1.0;
                if (_basisLoopPtr->pathType) {// vertical, so pull inside
                    dir = -1.0;
                }
                nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[(*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2];
                //

                //
                if (unitOutNormalForCritNode.find(startNodeId) == unitOutNormalForCritNode.end()) {
                    std::cout << "can not find normals" << std::endl;
                    exit(0);
                }
                //

                //double smallScale = heightDirection * (*_pathArcOnMeshPtr)[curArcId].vecPoints.back();
                //{
                //	double nxtHeight = nxtPt * heightDirection;
                //	double prevHeight = prePt * heightDirection;
                //	smallScale = abs(smallScale - nxtHeight) < abs(smallScale - prevHeight) ? abs(smallScale - nxtHeight) : abs(smallScale - prevHeight);
                //	if (abs(heightDirection * unitOutNormalForCritNode[startNodeId]) < 1e-10)
                //		smallScale = smallScale * 1e10;
                //	else
                //		smallScale = smallScale / abs(heightDirection * unitOutNormalForCritNode[startNodeId]);
                //}
                //smallScale = smallScale * smallScale;
                //

                double smallScale = 100000000.0;
                for (std::vector<int>::iterator vidIter = inMeshPtr->vecVertex[startNodeId].adjEdges.begin();
                     vidIter != inMeshPtr->vecVertex[startNodeId].adjEdges.end(); vidIter++) {
                    double dist = (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].x -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].x) *
                                  (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].x -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].x) +
                                  (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].y -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].y) *
                                  (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].y -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].y) +
                                  (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].z -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].z) *
                                  (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].z -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].z);
                    if (dist < smallScale)
                        smallScale = dist;
                }
                //
                Vector3 newPullDirection;
                int triangleId = 0;
                int opp_edge_id = 0;
                if ((*_pVecSimplifiedArc)[prevArcId].nCriticalNode0 == startNodeId) {
                    opp_edge_id = (*_pathArcOnMeshPointTypePtr)[prevArcId][
                            (*_pathArcOnMeshPointTypePtr)[prevArcId].size() - 2].first;
                    //newPullDirection = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[(*_pathArcOnMeshPointTypePtr)[prevArcId].size() - 2];
                } else {
                    opp_edge_id = (*_pathArcOnMeshPointTypePtr)[prevArcId][1].first;
                    //newPullDirection = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[1];
                }
                triangleId = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(opp_edge_id, startNodeId);
                newPullDirection = (*(inMeshPtr->meshNormalPtr))[triangleId];
                opp_edge_id = (*_pathArcOnMeshPointTypePtr)[curArcId][(*_pathArcOnMeshPointTypePtr)[curArcId].size() -
                                                                      2].first;
                triangleId = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(opp_edge_id, startNodeId);
                newPullDirection = newPullDirection + (*(inMeshPtr->meshNormalPtr))[triangleId];
                unitize(newPullDirection);
                //
                //newPullDirection = newPullDirection - (*_pathArcOnMeshPtr)[curArcId].vecPoints.back();
                //unitize(newPullDirection);
                //newPullDirection = newPullDirection + unitize((*_pathArcOnMeshPtr)[curArcId].vecPoints[(*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2]
                //-(*_pathArcOnMeshPtr)[curArcId].vecPoints.back());
                //if (newPullDirection * unitOutNormalForCritNode[startNodeId] < 0.0)
                //	newPullDirection = -newPullDirection;
                //newPullDirection = newPullDirection + unitOutNormalForCritNode[startNodeId];
                //unitize(newPullDirection);
                //
                smallScale = sqrt(smallScale);
                TranslatePoint(dir, smallScale,
                               (*_pathArcOnMeshPtr)[curArcId].vecPoints.back(), //newPullDirection, translatedPt);
                               unitOutNormalForCritNode[startNodeId], translatedPt);

                //TranslatePoint(dir, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints.back(),
                //						nxtPt, unitOutNormalForCritNode[startNodeId],
                //						translatedPt);

                outBasisLoop.vecPoints.push_back(translatedPt);
                //
                //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.back());
                //
                for (int iarc = (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2;
                     iarc > 0; iarc--) {// ignore the first point
/*					Vector3 ourNormal;
					WalkEdgeForNormal((*_pathArcOnMeshPointTypePtr)[curArcId][iarc].first, ourNormal);
					TranslatePoint(dir, smallScale, (*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc],
										ourNormal, translatedPt);	*/
                    outBasisLoop.vecPoints.push_back(
                            (*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);//(*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
                /*for ( int iarc = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 2; iarc > 0 ; iarc--)
				{
					outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
				}*/
            } else {// ignore the last point
                double dir = 1.0;
                if (_basisLoopPtr->pathType) {// vertical, so pull inside
                    dir = -1.0;
                }
                nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[1];
                //

                //
                if (unitOutNormalForCritNode.find(startNodeId) == unitOutNormalForCritNode.end()) {
                    std::cout << "can not find normals" << std::endl;
                    exit(0);
                }

                //double smallScale = heightDirection * (*_pathArcOnMeshPtr)[curArcId].vecPoints.front();
                //{
                //	double nxtHeight = nxtPt * heightDirection;
                //	double prevHeight = prePt * heightDirection;
                //	smallScale = abs(smallScale - nxtHeight) < abs(smallScale - prevHeight) ? abs(smallScale - nxtHeight) : abs(smallScale - prevHeight);
                //	if (abs(heightDirection * unitOutNormalForCritNode[startNodeId]) < 1e-10)
                //		smallScale = smallScale * 1e10;
                //	else
                //		smallScale = smallScale / abs(heightDirection * unitOutNormalForCritNode[startNodeId]);
                //}
                //smallScale = smallScale * smallScale;
                double smallScale = 100000000.0;
                for (std::vector<int>::iterator vidIter = inMeshPtr->vecVertex[startNodeId].adjEdges.begin();
                     vidIter != inMeshPtr->vecVertex[startNodeId].adjEdges.end(); vidIter++) {
                    double dist = (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].x -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].x) *
                                  (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].x -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].x) +
                                  (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].y -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].y) *
                                  (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].y -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].y) +
                                  (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].z -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].z) *
                                  (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].z -
                                   inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].z);
                    if (dist < smallScale)
                        smallScale = dist;
                }
                //
                Vector3 newPullDirection;
                int triangleId = 0;
                int opp_edge_id = 0;
                if ((*_pVecSimplifiedArc)[prevArcId].nCriticalNode0 == startNodeId) {
                    opp_edge_id = (*_pathArcOnMeshPointTypePtr)[prevArcId][
                            (*_pathArcOnMeshPointTypePtr)[prevArcId].size() - 2].first;
                    //newPullDirection = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[(*_pathArcOnMeshPointTypePtr)[prevArcId].size() - 2];
                } else {
                    opp_edge_id = (*_pathArcOnMeshPointTypePtr)[prevArcId][1].first;
                    //newPullDirection = (*_pathArcOnMeshPtr)[prevArcId].vecPoints[1];
                }
                triangleId = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(opp_edge_id, startNodeId);
                newPullDirection = (*(inMeshPtr->meshNormalPtr))[triangleId];
                opp_edge_id = (*_pathArcOnMeshPointTypePtr)[curArcId][1].first;
                triangleId = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(opp_edge_id, startNodeId);
                newPullDirection = newPullDirection + (*(inMeshPtr->meshNormalPtr))[triangleId];
                unitize(newPullDirection);
                newPullDirection = newPullDirection + unitOutNormalForCritNode[startNodeId];
                unitize(newPullDirection);
                //
                //
                //newPullDirection = newPullDirection - (*_pathArcOnMeshPtr)[curArcId].vecPoints.front();
                //unitize(newPullDirection);
                //newPullDirection = newPullDirection + unitize((*_pathArcOnMeshPtr)[curArcId].vecPoints[1]
                //-(*_pathArcOnMeshPtr)[curArcId].vecPoints.front());
                //unitize(newPullDirection);
                //if (newPullDirection * unitOutNormalForCritNode[startNodeId] < 0.0)
                //	newPullDirection = -newPullDirection;
                //newPullDirection = newPullDirection + unitOutNormalForCritNode[startNodeId];
                //unitize(newPullDirection);
                //
                //
                smallScale = sqrt(smallScale);
                TranslatePoint(dir, smallScale,
                               (*_pathArcOnMeshPtr)[curArcId].vecPoints.front(), //newPullDirection, translatedPt);
                               unitOutNormalForCritNode[startNodeId], translatedPt);
                //TranslatePoint(dir, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints[0],
                //						nxtPt, unitOutNormalForCritNode[startNodeId],
                //						translatedPt);
                //
                outBasisLoop.vecPoints.push_back(translatedPt);
                //
                //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[0]);
                for (unsigned int iarc = 1;
                     iarc < (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {// ignore the first point
                    //Vector3 ourNormal;
                    //WalkEdgeForNormal((*_pathArcOnMeshPointTypePtr)[curArcId][iarc].first, ourNormal);
                    //TranslatePoint(dir, smallScale, (*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc],
                    //					ourNormal, translatedPt);
                    outBasisLoop.vecPoints.push_back(
                            (*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);//(*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }/*
				for (unsigned int iarc = 1; iarc < (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 1 ; iarc++)
				{
					outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
				}*/

            }
        }
        // push the first point into the loop to close it
        //outLoopOnMesh.vecPoints.push_back(outLoopOnMesh.vecPoints[0]);
        //
        outBasisLoop.vecPoints.push_back(outBasisLoop.vecPoints[0]);
    } else {// only to embed the shared critical node

        std::set<int> shared_nodes;
        std::set<int> basis_nodes_set(_basisLoopPtr->nodeIndexSet.begin(), _basisLoopPtr->nodeIndexSet.end());
        std::set<int>::iterator findIter;
        std::map<int, int> shared_nodes_pos;
        //
        for (unsigned int j = 0; j < _loopOnMeshPtr->nodeIndexSet.size(); j++) {
            findIter = basis_nodes_set.find(_loopOnMeshPtr->nodeIndexSet[j]);
            if (findIter != basis_nodes_set.end()) {
                shared_nodes.insert(_loopOnMeshPtr->nodeIndexSet[j]);
                shared_nodes_pos[_loopOnMeshPtr->nodeIndexSet[j]] = j;
            }
            //if (_basisLoopPtr->nodeIndexSet[i] == _loopOnMeshPtr->nodeIndexSet[j])
            //	shared_nodes.insert(_loopOnMeshPtr->nodeIndexSet[j]);
        }
        //for (unsigned int i = 0; i < _basisLoopPtr->nodeIndexSet.size(); i++)
        //{
        //	for (unsigned int j = 0; j < _loopOnMeshPtr->nodeIndexSet.size(); j++)
        //	{
        //		if (_basisLoopPtr->nodeIndexSet[i] == _loopOnMeshPtr->nodeIndexSet[j])
        //			shared_nodes.insert(_loopOnMeshPtr->nodeIndexSet[j]);
        //	}
        //}
        std::set<int>::iterator sIter;
        //
        for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
            curArcId = _basisLoopPtr->simpArcIndexSet[i];
            startNodeId = _basisLoopPtr->nodeIndexSet[i];
            endNodeId = _basisLoopPtr->nodeIndexSet[i + 1];
            //
            //
            int prevArcId = i - 1;
            if (prevArcId < 0)
                prevArcId = _basisLoopPtr->simpArcIndexSet.size() - 1;
            //
            prevArcId = _basisLoopPtr->simpArcIndexSet[prevArcId];
            //
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction
                sIter = shared_nodes.find(startNodeId);
                if (sIter != shared_nodes.end()) {// need to pull out or in
                    //
                    std::cout << "shered case 1" << std::endl;
                    //
                    int orgArcId = shared_nodes_pos[startNodeId];
                    int orgPrevArcId = orgArcId - 1;
                    if (orgPrevArcId < 0) {
                        orgPrevArcId = _loopOnMeshPtr->simpArcIndexSet.size() - 1;
                    }
                    orgArcId = _loopOnMeshPtr->simpArcIndexSet[orgArcId];
                    orgPrevArcId = _loopOnMeshPtr->simpArcIndexSet[orgPrevArcId];
                    //
                    double dir = 1.0;
                    if (_basisLoopPtr->pathType) {// vertical, so pull inside
                        dir = -1.0;
                    }
                    nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[(*_pathArcOnMeshPtr)[curArcId].vecPoints.size() -
                                                                     2];
                    //
                    if (i == 0) {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet.back()].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[1];
                        }
                    } else {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[i - 1]].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[1];
                        }
                    }
                    //
                    //double smallScale = heightDirection * (*_pathArcOnMeshPtr)[curArcId].vecPoints.back();
                    //{
                    //	double nxtHeight = nxtPt * heightDirection;
                    //	double prevHeight = prePt * heightDirection;
                    //	smallScale = abs(smallScale - nxtHeight) < abs(smallScale - prevHeight) ? abs(smallScale - nxtHeight) : abs(smallScale - prevHeight);
                    //	if (abs(heightDirection * unitOutNormalForCritNode[startNodeId]) < 1e-10)
                    //		smallScale = smallScale * 1e10;
                    //	else
                    //		smallScale = smallScale / abs(heightDirection * unitOutNormalForCritNode[startNodeId]);
                    //}
                    //smallScale = smallScale * smallScale;
                    double smallScale = 100000000.0;
                    for (std::vector<int>::iterator vidIter = inMeshPtr->vecVertex[startNodeId].adjEdges.begin();
                         vidIter != inMeshPtr->vecVertex[startNodeId].adjEdges.end(); vidIter++) {
                        double dist =
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].x -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].x) *
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].x -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].x) +
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].y -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].y) *
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].y -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].y) +
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].z -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].z) *
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].z -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].z);
                        if (dist < smallScale)
                            smallScale = dist;
                    }
                    //
                    smallScale = sqrt(smallScale);
                    //
                    //
                    Vector3 newPullDirection;
                    int triangleId = 0;
                    int opp_edge_id = 0;
                    if ((*_pVecSimplifiedArc)[orgPrevArcId].nCriticalNode0 == startNodeId) {
                        opp_edge_id = (*_pathArcOnMeshPointTypePtr)[orgPrevArcId][
                                (*_pathArcOnMeshPointTypePtr)[orgPrevArcId].size() - 2].first;
                    } else {
                        opp_edge_id = (*_pathArcOnMeshPointTypePtr)[orgPrevArcId][1].first;
                    }
                    triangleId = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(opp_edge_id, startNodeId);
                    newPullDirection = (*(inMeshPtr->meshNormalPtr))[triangleId];
                    //
                    if ((*_pVecSimplifiedArc)[orgArcId].nCriticalNode0 == startNodeId) {
                        opp_edge_id = (*_pathArcOnMeshPointTypePtr)[orgArcId][
                                (*_pathArcOnMeshPointTypePtr)[orgArcId].size() - 2].first;
                    } else {
                        opp_edge_id = (*_pathArcOnMeshPointTypePtr)[orgArcId][1].first;
                    }
                    //opp_edge_id =  (*_pathArcOnMeshPointTypePtr)[orgArcId][(*_pathArcOnMeshPointTypePtr)[orgArcId].size() - 2].first;
                    //
                    triangleId = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(opp_edge_id, startNodeId);
                    newPullDirection = newPullDirection + (*(inMeshPtr->meshNormalPtr))[triangleId];
                    unitize(newPullDirection);
                    newPullDirection = newPullDirection + unitOutNormalForCritNode[startNodeId];
                    unitize(newPullDirection);
                    //
                    //newPullDirection = newPullDirection - (*_pathArcOnMeshPtr)[curArcId].vecPoints.back();
                    //unitize(newPullDirection);
                    //newPullDirection = newPullDirection + unitize((*_pathArcOnMeshPtr)[curArcId].vecPoints[(*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2]
                    //-(*_pathArcOnMeshPtr)[curArcId].vecPoints.back());
                    ////
                    //unitize(newPullDirection);
                    //if (newPullDirection * unitOutNormalForCritNode[startNodeId] < 0.0)
                    //	newPullDirection = -newPullDirection;
                    //newPullDirection = newPullDirection + unitOutNormalForCritNode[startNodeId];
                    //unitize(newPullDirection);
                    //
                    //
                    TranslatePoint(dir, smallScale,
                                   (*_pathArcOnMeshPtr)[curArcId].vecPoints.back(), //newPullDirection, translatedPt);
                                   unitOutNormalForCritNode[startNodeId], translatedPt);
                    //
                    //TranslatePoint(dir, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints.back(),
                    //						nxtPt, unitOutNormalForCritNode[startNodeId],
                    //						translatedPt);
                    //
                    outBasisLoop.vecPoints.push_back(translatedPt);
                } else {
                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints.back());
                }
                //
                //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.back());
                //
                if (!(*_bVecBasisLoopPtr)[i]) {// if the arc is setted to be a straight segment, ignore all other points
                    for (int iarc = (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2;
                         iarc > 0; iarc--) {// ignore the first point
                        outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);

                    }/*
					for ( int iarc = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 2; iarc > 0 ; iarc--)
					{
						outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
					}*/
                }
            } else {// ignore the last point

                sIter = shared_nodes.find(startNodeId);
                if (sIter != shared_nodes.end()) {// need to pull out or in
                    std::cout << "shered case 2" << std::endl;
                    //
                    int orgArcId = shared_nodes_pos[startNodeId];
                    int orgPrevArcId = orgArcId - 1;
                    if (orgPrevArcId < 0) {
                        orgPrevArcId = _loopOnMeshPtr->simpArcIndexSet.size() - 1;
                    }
                    orgArcId = _loopOnMeshPtr->simpArcIndexSet[orgArcId];
                    orgPrevArcId = _loopOnMeshPtr->simpArcIndexSet[orgPrevArcId];
                    //
                    double dir = 1.0;
                    if (_basisLoopPtr->pathType) {// vertical, so pull inside
                        dir = -1.0;
                    }
                    nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[1];
                    //
                    if (i == 0) {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet.back()].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[1];
                        }
                    } else {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[i - 1]].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[1];
                        }
                    }
                    //
                    //double smallScale = heightDirection * (*_pathArcOnMeshPtr)[curArcId].vecPoints.front();
                    //{
                    //	double nxtHeight = nxtPt * heightDirection;
                    //	double prevHeight = prePt * heightDirection;
                    //	smallScale = abs(smallScale - nxtHeight) < abs(smallScale - prevHeight) ? abs(smallScale - nxtHeight) : abs(smallScale - prevHeight);
                    //	if (abs(heightDirection * unitOutNormalForCritNode[startNodeId]) < 1e-10)
                    //		smallScale = smallScale * 1e10;
                    //	else
                    //		smallScale = smallScale / abs(heightDirection * unitOutNormalForCritNode[startNodeId]);
                    //}
                    //smallScale = smallScale * smallScale;
                    std::cout << "before scale" << std::endl;
                    double smallScale = 100000000.0;
                    for (std::vector<int>::iterator vidIter = inMeshPtr->vecVertex[startNodeId].adjEdges.begin();
                         vidIter != inMeshPtr->vecVertex[startNodeId].adjEdges.end(); vidIter++) {
                        double dist =
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].x -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].x) *
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].x -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].x) +
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].y -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].y) *
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].y -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].y) +
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].z -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].z) *
                                (inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v0].z -
                                 inMeshPtr->vecVertex[inMeshPtr->vecEdge[*vidIter].v1].z);
                        if (dist < smallScale)
                            smallScale = dist;
                    }
                    //
                    smallScale = sqrt(smallScale);
                    //
                    //Vector3 newPullDirection;
                    //int triangleId = 0;
                    //int opp_edge_id = 0;
                    //if ((*_pVecSimplifiedArc)[orgPrevArcId].nCriticalNode0 == startNodeId)
                    //{
                    //	opp_edge_id =  (*_pathArcOnMeshPointTypePtr)[orgPrevArcId][(*_pathArcOnMeshPointTypePtr)[orgPrevArcId].size() - 2].first;
                    //}
                    //else
                    //{
                    //	opp_edge_id =  (*_pathArcOnMeshPointTypePtr)[orgPrevArcId][1].first;
                    //}
                    //triangleId = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(opp_edge_id, startNodeId);
                    //newPullDirection = (*(inMeshPtr->meshNormalPtr))[triangleId];
                    ////
                    //if ((*_pVecSimplifiedArc)[orgArcId].nCriticalNode0 == startNodeId)
                    //{
                    //	opp_edge_id =  (*_pathArcOnMeshPointTypePtr)[orgArcId][(*_pathArcOnMeshPointTypePtr)[orgArcId].size() - 2].first;
                    //}
                    //else
                    //{
                    //	opp_edge_id =  (*_pathArcOnMeshPointTypePtr)[orgArcId][1].first;
                    //}
                    ////opp_edge_id =  (*_pathArcOnMeshPointTypePtr)[orgArcId][(*_pathArcOnMeshPointTypePtr)[orgArcId].size() - 2].first;
                    ////
                    //triangleId = inMeshPtr->TriangleAdjacentToOneEdgesOneVertex(opp_edge_id, startNodeId);
                    //newPullDirection = newPullDirection + (*(inMeshPtr->meshNormalPtr))[triangleId];
                    //unitize(newPullDirection);
                    //newPullDirection = newPullDirection + unitOutNormalForCritNode[startNodeId];
                    //unitize(newPullDirection);
                    //
                    //newPullDirection = newPullDirection - (*_pathArcOnMeshPtr)[curArcId].vecPoints.front();
                    //unitize(newPullDirection);
                    //newPullDirection = newPullDirection + unitize((*_pathArcOnMeshPtr)[curArcId].vecPoints[1]
                    //-(*_pathArcOnMeshPtr)[curArcId].vecPoints.front());
                    ////
                    //unitize(newPullDirection);
                    //if (newPullDirection * unitOutNormalForCritNode[startNodeId] < 0.0)
                    //	newPullDirection = -newPullDirection;
                    //newPullDirection = newPullDirection + unitOutNormalForCritNode[startNodeId];
                    //unitize(newPullDirection);
                    //
                    std::cout << "before tran" << std::endl;
                    TranslatePoint(dir, smallScale,
                                   (*_pathArcOnMeshPtr)[curArcId].vecPoints.front(), //newPullDirection,translatedPt);
                                   unitOutNormalForCritNode[startNodeId], translatedPt);
                    //
                    //TranslatePoint(dir, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints.front(),
                    //						nxtPt, unitOutNormalForCritNode[startNodeId],
                    //						translatedPt);
                    //
                    outBasisLoop.vecPoints.push_back(translatedPt);
                } else {
                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints.front());
                }
                //
                if (!(*_bVecBasisLoopPtr)[i]) {// if the arc is setted to be a straight segment, ignore all other points
                    //outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.front());
                    for (unsigned int iarc = 1;
                         iarc < (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {// ignore the first point
                        outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);

                    }/*
					for (unsigned int iarc = 1; iarc < (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 1 ; iarc++)
					{
						outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
					}*/
                }
            }
        }
        // push the first point into the loop to close it
        //outLoopOnMesh.vecPoints.push_back(outLoopOnMesh.vecPoints[0]);
        outBasisLoop.vecPoints.push_back(outBasisLoop.vecPoints[0]);

    }

    return;
}

void NonOverlappingCycles::SharedArcs(std::vector<int> &arcId) {
    std::set<int> loopOnMeshArcSet(_loopOnMeshPtr->simpArcIndexSet.begin(), _loopOnMeshPtr->simpArcIndexSet.end());
    if (_basisLoopPtr == _loopOnMeshPtr) {// all arcs are there
        arcId.assign(_loopOnMeshPtr->simpArcIndexSet.begin(), _loopOnMeshPtr->simpArcIndexSet.end());
    } else {
        std::set<int>::iterator sIter;
        for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
            sIter = loopOnMeshArcSet.find((_basisLoopPtr->simpArcIndexSet)[i]);
            if (sIter != loopOnMeshArcSet.end()) {// find it
                arcId.push_back((_basisLoopPtr->simpArcIndexSet)[i]);
            }
        }
    }
    return;
}

void NonOverlappingCycles::ComputeNewPathesForSharedArcs(
        std::vector<int> sharedArcs,
        std::vector<_Polygon> &sharedArcsNewPoints,
        std::vector<_Polygon> &offsetSharedArcsNewPoints) {
    // prepare the number of points on each edge for computing
    std::map<int, int> edgeId_to_vec_index_mapping;
    std::vector<std::vector<std::pair<Vector3, std::pair<std::pair<int, int>, int> > > > edge_with_all_points;
    //Vector3, std::pair<int, int>
    // point coord, arc_index, loopOnMesh_or_basis (0 on mesh, 1 offset), and point's order on the path
    std::map<int, int>::iterator mIter;
    //
    sharedArcsNewPoints.resize(sharedArcs.size());
    offsetSharedArcsNewPoints.resize(sharedArcs.size());
    //
    int edge_num = 0;
    Vector3 tempPt_a, tempPt_b;
    //
    for (unsigned int i = 0; i < sharedArcs.size(); i++) {
        int arcIdx = sharedArcs[i];
        for (unsigned int j = 0; j < (*_pathArcOnMeshPointTypePtr)[arcIdx].size(); j++) {
            if ((*_pathArcOnMeshPointTypePtr)[arcIdx][j].second) {// this is the edge
                mIter = edgeId_to_vec_index_mapping.find((*_pathArcOnMeshPointTypePtr)[arcIdx][j].first);
                if (mIter == edgeId_to_vec_index_mapping.end()) {// new edge record the index
                    edgeId_to_vec_index_mapping[(*_pathArcOnMeshPointTypePtr)[arcIdx][j].first] = edge_num++;
                    //
                    std::vector<std::pair<Vector3, std::pair<std::pair<int, int>, int> > > temp;
                    //
                    tempPt_a = (*_pathArcOnMeshPtr)[arcIdx].vecPoints[j];
                    //
                    temp.push_back(std::pair<Vector3, std::pair<std::pair<int, int>, int> >(tempPt_a,
                                                                                            std::pair<std::pair<int, int>, int>(
                                                                                                    std::pair<int, int>(
                                                                                                            arcIdx, 0),
                                                                                                    j)));
                    //
                    edge_with_all_points.push_back(temp);
                } else {
                    tempPt_a = (*_pathArcOnMeshPtr)[arcIdx].vecPoints[j];
                    //
                    edge_with_all_points[mIter->second].push_back(
                            std::pair<Vector3, std::pair<std::pair<int, int>, int> >(tempPt_a,
                                                                                     std::pair<std::pair<int, int>, int>(
                                                                                             std::pair<int, int>(arcIdx,
                                                                                                                 0),
                                                                                             j)));
                }
            }
        }
        for (unsigned int j = 0; j < (*_offsetPathArcOnMeshPointTypePtr)[arcIdx].size(); j++) {
            if ((*_offsetPathArcOnMeshPointTypePtr)[arcIdx][j].second) {// this is the edge
                mIter = edgeId_to_vec_index_mapping.find((*_offsetPathArcOnMeshPointTypePtr)[arcIdx][j].first);
                if (mIter == edgeId_to_vec_index_mapping.end()) {
                    edgeId_to_vec_index_mapping[(*_offsetPathArcOnMeshPointTypePtr)[arcIdx][j].first] = edge_num++;
                    //
                    std::vector<std::pair<Vector3, std::pair<std::pair<int, int>, int> > > temp;
                    //
                    tempPt_a = (*_offsetPathArcOnMeshPtr)[arcIdx].vecPoints[j];
                    //
                    temp.push_back(std::pair<Vector3, std::pair<std::pair<int, int>, int> >(tempPt_a,
                                                                                            std::pair<std::pair<int, int>, int>(
                                                                                                    std::pair<int, int>(
                                                                                                            arcIdx, 1),
                                                                                                    j)));
                    //
                    edge_with_all_points.push_back(temp);
                } else {
                    tempPt_a = (*_offsetPathArcOnMeshPtr)[arcIdx].vecPoints[j];
                    edge_with_all_points[mIter->second].push_back(
                            std::pair<Vector3, std::pair<std::pair<int, int>, int> >(tempPt_a,
                                                                                     std::pair<std::pair<int, int>, int>(
                                                                                             std::pair<int, int>(arcIdx,
                                                                                                                 1),
                                                                                             j)));
                }
            }
        }
    }
    // equally distributed all points on the edge by its order to one end point
    for (mIter = edgeId_to_vec_index_mapping.begin();
         mIter != edgeId_to_vec_index_mapping.end();
         mIter++) {
        int edge_id = mIter->first;
        int vecIndex = mIter->second;
        if (edge_with_all_points[vecIndex].size() == 1) {// take the middle point of this edge is fine
            tempPt_a = Vector3(inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v0].x,
                               inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v0].y,
                               inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v0].z);
            //
            tempPt_b = Vector3(inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v1].x,
                               inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v1].y,
                               inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v1].z);
            //
            edge_with_all_points[vecIndex][0].first = 0.5 * (tempPt_a + tempPt_b);
            //
        } else {// sort all points on this edge by their distance to one end point
            std::set<std::pair<double, std::pair<std::pair<int, int>, int> >, doubleIntLessThan> sorted_dist;
            double weight_a = 1.0 / (edge_with_all_points[vecIndex].size() + 1);
            //
            tempPt_a = Vector3(inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v0].x,
                               inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v0].y,
                               inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v0].z);
            //
            tempPt_b = Vector3(inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v1].x,
                               inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v1].y,
                               inMeshPtr->vecVertex[inMeshPtr->vecEdge[edge_id].v1].z);
            for (unsigned int iarc = 0; iarc <
                                        edge_with_all_points[vecIndex].size(); iarc++) {// all point are distinct because it from the previous distinct computation
                sorted_dist.insert(std::pair<double, std::pair<std::pair<int, int>, int> >(
                        norm2(tempPt_a - edge_with_all_points[vecIndex][iarc].first),
                        edge_with_all_points[vecIndex][iarc].second));
            }
            int dist_order = 1;
            for (std::set<std::pair<double, std::pair<std::pair<int, int>, int> >, doubleIntLessThan>::iterator sIter = sorted_dist.begin();
                 sIter != sorted_dist.end();
                 sIter++) {
                edge_with_all_points[vecIndex][dist_order - 1].first =
                        tempPt_a * (1.0 - dist_order * weight_a) + dist_order * weight_a * tempPt_b;
                edge_with_all_points[vecIndex][dist_order - 1].second = sIter->second;
                dist_order++;
            }
        }
    }

    for (unsigned int i = 0; i < sharedArcs.size(); i++) {// assumption each edge cross the edge only once
        int arcIdx = sharedArcs[i];
        int edge_id = 0;
        int vecIndex = 0;
        for (unsigned int j = 0; j < (*_pathArcOnMeshPointTypePtr)[arcIdx].size(); j++) {
            edge_id = (*_pathArcOnMeshPointTypePtr)[arcIdx][j].first;

            //
            if ((*_pathArcOnMeshPointTypePtr)[arcIdx][j].second) {// this is the edge
                vecIndex = edgeId_to_vec_index_mapping[edge_id];
                //
                for (unsigned iarc = 0; iarc < edge_with_all_points[vecIndex].size(); iarc++) {//
                    if (edge_with_all_points[vecIndex][iarc].second ==
                        std::pair<std::pair<int, int>, int>(std::pair<int, int>(arcIdx, 0),
                                                            j)) {// this is the point we want to use
                        sharedArcsNewPoints[i].vecPoints.push_back(edge_with_all_points[vecIndex][iarc].first);
                    }
                }
            } else {
                sharedArcsNewPoints[i].vecPoints.push_back((*_pathArcOnMeshPtr)[arcIdx].vecPoints[j]);
            }
        }
        //
        for (unsigned int j = 0; j < (*_offsetPathArcOnMeshPointTypePtr)[arcIdx].size(); j++) {
            edge_id = (*_offsetPathArcOnMeshPointTypePtr)[arcIdx][j].first;

            //
            if ((*_offsetPathArcOnMeshPointTypePtr)[arcIdx][j].second) {// this is the edge
                vecIndex = edgeId_to_vec_index_mapping[edge_id];
                //
                for (unsigned iarc = 0; iarc < edge_with_all_points[vecIndex].size(); iarc++) {//
                    if (edge_with_all_points[vecIndex][iarc].second ==
                        std::pair<std::pair<int, int>, int>(std::pair<int, int>(arcIdx, 1),
                                                            j)) {// this is the point we want to use
                        offsetSharedArcsNewPoints[i].vecPoints.push_back(edge_with_all_points[vecIndex][iarc].first);
                    }
                }
            } else {
                offsetSharedArcsNewPoints[i].vecPoints.push_back((*_offsetPathArcOnMeshPtr)[arcIdx].vecPoints[j]);
            }
        }
    }
}

void NonOverlappingCycles::ComputeNonOverlappingCycles() {// path is encoded in this way
// v0--e0--v1--e1--....vn--en--v(n+1)[==v0]
// each cycle are passing from high to low
// each edge is taken in this way [a...b)-[b...c)-[c...a)
    int curArcId = 0;
    int startNodeId = 0;
    int endNodeId = 0;
    Vector3 translatedPt;
    Vector3 prePt, nxtPt;
    if (_basisLoopPtr == _loopOnMeshPtr) {// need to embed every critical node on the path
        for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
            curArcId = _basisLoopPtr->simpArcIndexSet[i];
            startNodeId = _basisLoopPtr->nodeIndexSet[i];
            endNodeId = _basisLoopPtr->nodeIndexSet[i + 1];
            //
            if (i == 0) {
                if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet.back()].nCriticalNode0 == startNodeId) {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[
                            (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints.size() - 2];
                } else {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[1];
                }
            } else {
                if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[i - 1]].nCriticalNode0 == startNodeId) {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[
                            (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints.size() - 2];
                } else {
                    prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[1];
                }
            }
            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction
                double dir = 1.0;
                if (_basisLoopPtr->pathType) {// vertical, so pull inside
                    dir = -1.0;
                }
                nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[(*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2];
                //

                //
                TranslatePoint(-1.0, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints.back(),
                               nxtPt, unitOutNormalForCritNode[startNodeId],
                               translatedPt);
                //
                outBasisLoop.vecPoints.push_back(translatedPt);
                //
                outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.back());
                //
                for (int iarc = (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2;
                     iarc > 0; iarc--) {// ignore the first point

                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
                for (int iarc = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 2; iarc > 0; iarc--) {
                    outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
            } else {// ignore the last point
                double dir = 1.0;
                if (_basisLoopPtr->pathType) {// vertical, so pull inside
                    dir = -1.0;
                }
                nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[1];
                //

                //
                TranslatePoint(-1.0, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints[0],
                               nxtPt, unitOutNormalForCritNode[startNodeId],
                               translatedPt);
                //
                outBasisLoop.vecPoints.push_back(translatedPt);
                //
                outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[0]);
                for (unsigned int iarc = 1;
                     iarc < (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {// ignore the first point

                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
                for (unsigned int iarc = 1; iarc < (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {
                    outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }

            }
        }
        // push the first point into the loop to close it
        outLoopOnMesh.vecPoints.push_back(outLoopOnMesh.vecPoints[0]);
        //
        outBasisLoop.vecPoints.push_back(outBasisLoop.vecPoints[0]);
    } else {// only to embed the shared critical node

        std::set<int> shared_nodes;
        for (unsigned int i = 0; i < _basisLoopPtr->nodeIndexSet.size(); i++) {
            for (unsigned int j = 0; j < _loopOnMeshPtr->nodeIndexSet.size(); j++) {
                if (_basisLoopPtr->nodeIndexSet[i] == _loopOnMeshPtr->nodeIndexSet[j])
                    shared_nodes.insert(_loopOnMeshPtr->nodeIndexSet[j]);
            }
        }
        std::set<int>::iterator sIter;
        //
        for (unsigned int i = 0; i < _basisLoopPtr->simpArcIndexSet.size(); i++) {
            curArcId = _basisLoopPtr->simpArcIndexSet[i];
            startNodeId = _basisLoopPtr->nodeIndexSet[i];
            endNodeId = _basisLoopPtr->nodeIndexSet[i + 1];
            //

            if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to trave in the opposite direction
                sIter = shared_nodes.find(startNodeId);
                if (sIter != shared_nodes.end()) {// need to pull out or in
                    double dir = 1.0;
                    if (_basisLoopPtr->pathType) {// vertical, so pull inside
                        dir = -1.0;
                    }
                    nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[(*_pathArcOnMeshPtr)[curArcId].vecPoints.size() -
                                                                     2];
                    //
                    if (i == 0) {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet.back()].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[1];
                        }
                    } else {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[i - 1]].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[1];
                        }
                    }
                    //
                    TranslatePoint(-1.0, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints.back(),
                                   nxtPt, unitOutNormalForCritNode[startNodeId],
                                   translatedPt);
                    //
                    outBasisLoop.vecPoints.push_back(translatedPt);
                } else {
                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints.back());
                }
                //
                outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.back());
                //
                for (int iarc = (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 2;
                     iarc > 0; iarc--) {// ignore the first point
                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);

                }
                for (int iarc = (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 2; iarc > 0; iarc--) {
                    outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
            } else {// ignore the last point
                sIter = shared_nodes.find(startNodeId);
                if (sIter != shared_nodes.end()) {// need to pull out or in
                    double dir = 1.0;
                    if (_basisLoopPtr->pathType) {// vertical, so pull inside
                        dir = -1.0;
                    }
                    nxtPt = (*_pathArcOnMeshPtr)[curArcId].vecPoints[1];
                    //
                    if (i == 0) {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet.back()].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet.back()].vecPoints[1];
                        }
                    } else {
                        if ((*_pVecSimplifiedArc)[_basisLoopPtr->simpArcIndexSet[i - 1]].nCriticalNode0 ==
                            startNodeId) {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[
                                    (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints.size() - 2];
                        } else {
                            prePt = (*_pathArcOnMeshPtr)[_basisLoopPtr->simpArcIndexSet[i - 1]].vecPoints[1];
                        }
                    }
                    //
                    TranslatePoint(-1.0, prePt, (*_pathArcOnMeshPtr)[curArcId].vecPoints.front(),
                                   nxtPt, unitOutNormalForCritNode[startNodeId],
                                   translatedPt);
                    //
                    outBasisLoop.vecPoints.push_back(translatedPt);
                } else {
                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints.front());
                }
                //
                outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.front());
                for (unsigned int iarc = 1;
                     iarc < (*_pathArcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {// ignore the first point
                    outBasisLoop.vecPoints.push_back((*_pathArcOnMeshPtr)[curArcId].vecPoints[iarc]);

                }
                for (unsigned int iarc = 1; iarc < (*_offsetPathArcOnMeshPtr)[curArcId].vecPoints.size() - 1; iarc++) {
                    outLoopOnMesh.vecPoints.push_back((*_offsetPathArcOnMeshPtr)[curArcId].vecPoints[iarc]);
                }
            }
        }
        // push the first point into the loop to close it
        outLoopOnMesh.vecPoints.push_back(outLoopOnMesh.vecPoints[0]);
        outBasisLoop.vecPoints.push_back(outBasisLoop.vecPoints[0]);

    }

    return;
}

