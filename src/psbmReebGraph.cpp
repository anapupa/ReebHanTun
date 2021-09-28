/*
(c) 2012 Fengtao Fan
*/
#include "psbmReebGraph.h"
#include "AdjListGraph.h"
#include "ReebGraphPairing.h"
#include <string>
#include <stack>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "NonOverlappingCycles.h"
#include "LinkNumberComputing.h"
#include "NonOverlappedLevelCycleAndArc.h"
#include "MapLoopsBackToMeshLevelSetAndArc.h"
#include "InverseLinkNumberMatrix.h"

#include "FilesOutputForOptimalCycles.h"

#include "Graph_max_tree.h"

struct myVector3LessThan {
    bool operator()(const Vector3 &lhs, const Vector3 &rhs) {
        bool ret = false;
        if (lhs[0] != rhs[0]) {
            if (lhs[0] < rhs[0])
                ret = true;
            else
                ret = false;
        } else {
            if (lhs[1] != rhs[1]) {
                if (lhs[1] < rhs[1])
                    ret = true;
                else
                    ret = false;
            } else {
                if (lhs[2] != rhs[2]) {
                    if (lhs[2] < rhs[2])
                        ret = true;
                    else
                        ret = false;
                }
            }
        }
        return ret;
    }
};

void SetAllNull(auxMeshEdge &inEdge) {
    inEdge.ptrPos = NULL;
    inEdge.ptrReebArcId = NULL;
}

#define psbmReebGraphSwapVars(type, var1, var2)  \
{\
    type tmp;\
    tmp=(var1);\
    (var1)=(var2);\
    (var2)=tmp;\
}

#define psbmReebGraphAdvanceAlongEdge(repE0, a0) \
{\
    repE0 = repE0->ptrNextPos;\
    if (repE0)\
        a0 = repE0->ptrReebArcId;\
    else\
        a0 = NULL;\
}

/***some comparing function****************/
bool lessThanPairs(const float &lhsF, const int &lhsI,
                   const float &rhsF, const int &rhsI) {
    bool ret = false;
    if (lhsF < rhsF)
        ret = true;
    else if (lhsF == rhsF) {
        if (lhsI < rhsI)
            ret = true;
    }
    return ret;
}

struct myFloatIntPairCompare {
    bool operator()(const std::pair<float, int> &lhs, const std::pair<float, int> &rhs) {
        bool ret = false;
        if (lhs.first < rhs.first)
            ret = true;
        else if (lhs.first == rhs.first) {
            if (lhs.second < rhs.second)
                ret = true;
        }
        return ret;
    }

    bool operator()(const std::pair<float, int> &lhs, const std::pair<float, int> &rhs) const {
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

struct myDoubleIntPairCompare {
    bool operator()(const std::pair<double, int> &lhs, const std::pair<double, int> &rhs) {
        bool ret = false;
        if (lhs.first < rhs.first)
            ret = true;
        else if (lhs.first == rhs.first) {
            if (lhs.second < rhs.second)
                ret = true;
        }
        return ret;
    }

    bool operator()(const std::pair<double, int> &lhs, const std::pair<double, int> &rhs) const {
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

struct myIntDoublePairCompare {
    bool operator()(const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) {
        bool ret = false;
        if (lhs.second < rhs.second)
            ret = true;
        else if (lhs.second == rhs.second) {
            if (lhs.first < rhs.first)
                ret = true;
        }
        return ret;
    }

    bool operator()(const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) const {
        bool ret = false;
        if (lhs.second < rhs.second)
            ret = true;
        else if (lhs.second == rhs.second) {
            if (lhs.first < rhs.first)
                ret = true;
        }
        return ret;
    }
};

struct myIntDoubleIntPairCompare {
    bool
    operator()(const std::pair<std::pair<int, int>, double> &lhs, const std::pair<std::pair<int, int>, double> &rhs) {
        bool ret = false;
        if (lhs.second < rhs.second)
            ret = true;
        else if (lhs.second == rhs.second) {
            if (lhs.first.second < rhs.first.second)
                ret = true;
        }
        return ret;
    }

    bool operator()(const std::pair<std::pair<int, int>, double> &lhs,
                    const std::pair<std::pair<int, int>, double> &rhs) const {
        bool ret = false;
        if (lhs.second < rhs.second)
            ret = true;
        else if (lhs.second == rhs.second) {
            if (lhs.first.second < rhs.first.second)
                ret = true;
        }
        return ret;
    }
};

struct myIntFloatPairCompare {
    bool operator()(const std::pair<int, float> &lhs, const std::pair<int, float> &rhs) {
        bool ret = false;
        if (lhs.second < rhs.second)
            ret = true;
        else if (lhs.second == rhs.second) {
            if (lhs.first < rhs.first)
                ret = true;
        }
        return ret;
    }

    bool operator()(const std::pair<int, float> &lhs, const std::pair<int, float> &rhs) const {
        bool ret = false;
        if (lhs.second < rhs.second)
            ret = true;
        else if (lhs.second == rhs.second) {
            if (lhs.first < rhs.first)
                ret = true;
        }
        return ret;
    }
};

/***************************/
void RotateAxisTobeAlignedWithZ_axis(double u, double v, double w,
                                     double x, double y, double z,
                                     Vector3 &outRes) {
    double uvw = 1.0 / sqrt(u * u + v * v + w * w);
    double uv = 1.0 / sqrt(u * u + v * v);
    outRes[0] = ((u * w * x + w * v * y) * uv - (z) / uv) * uvw;
    outRes[1] = (u * y - v * x) * uv;
    outRes[2] = (u * x + v * y + w * z) * uvw;
}

/******************************************/
void psbmReebGraph::InitReebNodes(const unsigned int nTotSize) {
    if (!pVecReebNode)
        pVecReebNode = new std::vector<class psbmReebNode>;
    pVecReebNode->clear();
    pVecReebNode->resize(nTotSize);

    //
    // add the vertex into processed map
//	if (!ProcessedVertex)
//		ProcessedVertex = new std::map<int, int>;
    for (unsigned int i = 0; i < nTotSize; i++) {
        (*pVecReebNode)[i].nVertexId = i;
        (*pVecReebNode)[i].ptrArcUpId = new std::list<class psbmReebArc *>;
        (*pVecReebNode)[i].ptrArcDownId = new std::list<class psbmReebArc *>;
        (*pVecReebNode)[i].ptrHArc = new std::list<class psbmReebArc *>;
//		(*ProcessedVertex)[i] = i;
    }
    return;
}

int psbmReebGraph::CreateNode(const int nVertexId)//, const float scalarValue)
{// nVertexId is the vertex index in mesh vertices array
    psbmReebNode tmpNode(nVertexId);//, scalarValue);
    if (!pVecReebNode)
        pVecReebNode = new std::vector<class psbmReebNode>;

    pVecReebNode->push_back(tmpNode);

    // add the vertex into processed map
    //if (!ProcessedVertex)
    //	ProcessedVertex = new std::map<int, int>;

    //(*ProcessedVertex)[nVertexId] = pVecReebNode->size() - 1;
    //
    return nVertexId;// (*ProcessedVertex)[nVertexId];
}

void psbmReebGraph::CreateArc(const int N0, const int N1, const int eId) {// pVecAuxMeshEdge is intialized before
    // N0 and N1 are indices in the Reeb graph nodes array
    // N0.value < N1.value
    psbmReebArc *tmpArc = new psbmReebArc(N0, N1, eId);
    //if (!pListReebArc)
    //	pListReebArc = new std::set<class psbmReebArc*>;
    pListReebArc->insert(tmpArc);
    // set mesh edge
    (*pVecAuxMeshEdge)[eId].ptrReebArcId = tmpArc;
    (*pVecAuxMeshEdge)[eId].ptrPos = tmpArc->pArcList->front();
    // set graph node

    (*pVecReebNode)[N0].ptrArcUpId->push_back(tmpArc);
    (*pVecReebNode)[N1].ptrArcDownId->push_back(tmpArc);
    return;
}

void psbmReebGraph::MergePaths(const int N0, const int N1, const int N2,
                               const int e0, const int e1, const int e2) {
    psbmReebArc *a0 = (*pVecAuxMeshEdge)[e0].ptrReebArcId;
    psbmReebArc *a1 = (*pVecAuxMeshEdge)[e1].ptrReebArcId;
    psbmReebArc *a2 = (*pVecAuxMeshEdge)[e2].ptrReebArcId;

    //
    psbmArcListNode *travE0 = (*pVecAuxMeshEdge)[e0].ptrPos;
    psbmArcListNode *travE1 = (*pVecAuxMeshEdge)[e1].ptrPos;
    //
    GlueByMergeSorting(N0, N1, &a0, &a1, &travE0, &travE1);
    //
    travE1 = (*pVecAuxMeshEdge)[e2].ptrPos;
    GlueByMergeSorting(N1, N2, &a0, &a2, &travE0, &travE1);
    //
    return;
}

void psbmReebGraph::DeleteArcListNodeForEdge(const int eid) {
    //std::cout << "deleting edge " << eid << std::endl;
    psbmReebArc *currentReebArc = (*pVecAuxMeshEdge)[eid].ptrReebArcId;
    psbmArcListNode *currentArcListNode = (*pVecAuxMeshEdge)[eid].ptrPos;
    psbmArcListNode *delPtr = NULL;
    while (currentArcListNode) {
        std::list<struct psbmArcListNode *>::iterator listIter =
                currentReebArc->pArcList->begin();
        for (; listIter != currentReebArc->pArcList->end(); listIter++)
            if (*listIter == currentArcListNode)
                break;
        currentReebArc->pArcList->erase(listIter);
        //
        delPtr = currentArcListNode;
        currentArcListNode = delPtr->ptrNextPos;
        if (currentArcListNode)
            currentReebArc = currentArcListNode->ptrReebArcId;
        //
        delete delPtr;
    }
    return;
}

void psbmReebGraph::DeleteArc(psbmReebArc *a_del) {
    // reset arc-node relation
    // delete downArc from a1->nNodeId0
    std::list<class psbmReebArc *>::iterator vIter;
    vIter = std::find((*pVecReebNode)[a_del->nNodeId1].ptrArcDownId->begin(),
                      (*pVecReebNode)[a_del->nNodeId1].ptrArcDownId->end(),
                      a_del);
    (*pVecReebNode)[a_del->nNodeId1].ptrArcDownId->erase(vIter);
    //
    vIter = std::find((*pVecReebNode)[a_del->nNodeId0].ptrArcUpId->begin(),
                      (*pVecReebNode)[a_del->nNodeId0].ptrArcUpId->end(),
                      a_del);
    (*pVecReebNode)[a_del->nNodeId0].ptrArcUpId->erase(vIter);
    //
    std::set<class psbmReebArc *>::iterator iter;
    iter = pListReebArc->find(a_del);//std::find(pListReebArc->begin(),
    //pListReebArc->end(), a_del);
    //(*iter) = NULL;
    a_del->ClearArcList();
    delete a_del;
    pListReebArc->erase(iter);
    return;
}

void psbmReebGraph::GlueByMergeSorting(const int N0, const int N1,
                                       psbmReebArc **a0_in, psbmReebArc **a1_in,
                                       psbmArcListNode **travE0, psbmArcListNode **travE1) {
    psbmReebArc *a0 = *a0_in;
    psbmReebArc *a1 = *a1_in;
    int BottomNode0 = 0;
    int BottomNode1 = 0;

    psbmArcListNode *repE0 = *travE0;
    psbmArcListNode *repE1 = *travE1;


    if (a0 == NULL || a1 == NULL) {// *a0_in = NULL or *a1_in = NULL
        return;
    }
    do {
        BottomNode0 = a0->nNodeId0;
        BottomNode1 = a1->nNodeId0;
        // check if a0 and a1 are the same arcs
        if (a0 == a1) {// just move on, do nothing
            psbmReebGraphAdvanceAlongEdge(repE0, a0);
            psbmReebGraphAdvanceAlongEdge(repE1, a1);
        } else {// check if a0 and a1 have the same downpoint
            if (a0->nNodeId0 == a1->nNodeId0) {// have the same downpoint
                // simply merge them and delete one
                // no edges are duplicated
                std::list<struct psbmArcListNode *>::iterator siter, iterSign;
                iterSign = a0->pArcList->begin();
                a0->pArcList->splice(iterSign, *(a1->pArcList));
                for (siter = a0->pArcList->begin();
                     siter != iterSign; siter++) {// reset the pointer to its host arc
                    (*siter)->ptrReebArcId = a0;
                    // reset parent arc pointer
                    if ((*pVecAuxMeshEdge)[(*siter)->nEdgeId].ptrReebArcId == a1)
                        (*pVecAuxMeshEdge)[(*siter)->nEdgeId].ptrReebArcId = a0;
                }
                DeleteArc(a1);

                //
                psbmReebGraphAdvanceAlongEdge(repE0, a0);
                psbmReebGraphAdvanceAlongEdge(repE1, a1);
                //a1 = NULL;
            } else {// modify the one with long height and merge
                if ((MeshVertexScalarValue[(*pVecReebNode)[BottomNode0].nVertexId] <
                     MeshVertexScalarValue[(*pVecReebNode)[BottomNode1].nVertexId])
                    ||
                    (MeshVertexScalarValue[(*pVecReebNode)[BottomNode0].nVertexId] ==
                     MeshVertexScalarValue[(*pVecReebNode)[BottomNode1].nVertexId]
                     && (*pVecReebNode)[BottomNode0].nVertexId <
                        (*pVecReebNode)[BottomNode1].nVertexId)) //BottomNode0 < BottomNode1
                {//  B1.val > B0.val
                    MergeArcs(a1, a0);
                    psbmReebGraphAdvanceAlongEdge(repE1, a1);
                    repE0 = repE0->ptrNextPos;
                } else {// B0.val > B1.val
                    MergeArcs(a0, a1);
                    psbmReebGraphAdvanceAlongEdge(repE0, a0);
                    repE1 = repE1->ptrNextPos;
                }
            }
        }
        //std::cout << "can not exit " << pVecReebNode->size() << std::endl;
        //std::cout << "arc size " << pListReebArc->size() << std::endl;
    } while (a0 && a1);
    *a0_in = a0;
    *travE0 = repE0;
    return;
}

void psbmReebGraph::MergeArcs(psbmReebArc *a_high, psbmReebArc *a_low) {
    //a_high has bigger f(bottom_node)
    /*		#
		   * *
	(high)*   *
		       *(low)
	*/
    // add all spaces in a_low to a_high
    // modify a_low
    // 1) modify a_low
    // delete downArc from a_low->nNodeId1
    std::list<class psbmReebArc *>::iterator vIter;
    vIter = std::find((*pVecReebNode)[a_low->nNodeId1].ptrArcDownId->begin(),
                      (*pVecReebNode)[a_low->nNodeId1].ptrArcDownId->end(),
                      a_low);
    (*pVecReebNode)[a_low->nNodeId1].ptrArcDownId->erase(vIter);
    //
    //vIter = std::find((*pVecReebNode)[a_low->nNodeId0].ptrArcUpId->begin(),
    //		(*pVecReebNode)[a_low->nNodeId0].ptrArcUpId->end(),
    //		a_low);
    //(*pVecReebNode)[a_low->nNodeId0].ptrArcUpId->erase(vIter);
    // add downArc to a_high->nNodeId1
    (*pVecReebNode)[a_high->nNodeId0].ptrArcDownId->push_back(a_low);
    //(*pVecReebNode)[a_low->nNodeId0].ptrArcUpId->push_back(a_low);
    // if edge in a_low has a_low as highest arc change it to a_high
    a_low->nNodeId1 = a_high->nNodeId0;

    //2) copy every edges in a_low to a_high
    // this makes sure that every edge points to edges in a_low is not changed
    std::list<struct psbmArcListNode *>::iterator siter;
    int a_low_size = (int) a_low->pArcList->size();
    a_high->pArcList->splice(a_high->pArcList->begin(), *(a_low->pArcList));
    //generate new edge nodes in a_low
    struct psbmArcListNode *tmpNode = NULL;
    siter = a_high->pArcList->begin();
    for (int i = 0; i < a_low_size; i++, siter++) {
        tmpNode = new psbmArcListNode;
        tmpNode->nEdgeId = (*siter)->nEdgeId;
        tmpNode->ptrNextPos = (*siter)->ptrNextPos;
        tmpNode->ptrReebArcId = (*siter)->ptrReebArcId; // equal to a_low
        // set new edge node pointer for sister
        (*siter)->ptrNextPos = tmpNode;
        (*siter)->ptrReebArcId = a_high;
        //
        a_low->pArcList->push_back(tmpNode);
        //
        if ((*pVecAuxMeshEdge)[(*siter)->nEdgeId].ptrReebArcId == a_low)
            (*pVecAuxMeshEdge)[(*siter)->nEdgeId].ptrReebArcId = a_high;
    }

}

void psbmReebGraph::AddMeshTriangle(int v0, double f0,
                                    int v1, double f1,
                                    int v2, double f2,
                                    int e01, int e12, int e02) {
    // Add the vertices to the stream
    //std::map<int, int>::iterator sIter;

    // // vertex0
    //sIter = (*ProcessedVertex).find(v0);
    //if(sIter == (*ProcessedVertex).end())
    //{
    //// this vertex hasn't been processed yet, add it
    //	CreateNode(v0);//, f0);
    //}

    // // vertex1
    //sIter = (*ProcessedVertex).find(v1);
    //if(sIter == (*ProcessedVertex).end())
    //{
    //// this vertex hasn't been processed yet, add it
    //	CreateNode(v1);//, f1);
    //}

    // // vertex2
    //sIter = (*ProcessedVertex).find(v2);
    //if(sIter == (*ProcessedVertex).end())
    //{
    //// this vertex hasn't been processed yet, add it
    //	CreateNode(v2);//, f2);
    //}

    // All three vertices are processed before

    int Node0 = v0;//(*ProcessedVertex)[v0],
    int Node1 = v1;//(*ProcessedVertex)[v1],
    int Node2 = v2;//(*ProcessedVertex)[v2];
    // Consistency less check
    if (f2 < f1 || (f2 == f1 && v2 < v1)) // Node2 < Node1
    {// compare index in the mesh so that it can also be used for level set computing
//        psbmReebGraphSwapVars(int, v1, v2);
        std::swap(v1, v2);
        std::swap(Node1, Node2);
        std::swap(f1, f2);
        std::swap(e01, e02);
//        psbmReebGraphSwapVars(int, Node1, Node2);
//        psbmReebGraphSwapVars(double, f1, f2);
//        psbmReebGraphSwapVars(int, e01, e02);


//        std::swap(Node1, Node2);
//        std::swap(f1, f2);
//        std::swap(e01, e02);
    }
    if (f1 < f0 || (f1 == f0 && v1 < v0)) //Node1 < Node0
    {
//        psbmReebGraphSwapVars(int, v0, v1);
//        psbmReebGraphSwapVars(int, Node0, Node1);
//        psbmReebGraphSwapVars(double, f0, f1);
//        psbmReebGraphSwapVars(int, e12, e02);

        std::swap(v1, v0);
        std::swap(Node1, Node0);
        std::swap(f1, f0);
        std::swap(e12, e02);
    }
    if (f2 < f1 || (f2 == f1 && v2 < v1)) //Node2 < Node1
    {
        psbmReebGraphSwapVars(int, v1, v2);
        psbmReebGraphSwapVars(int, Node1, Node2);
        psbmReebGraphSwapVars(double, f1, f2);
        psbmReebGraphSwapVars(int, e01, e02);

//        std::swap(v1, v2);
//        std::swap(Node1, Node2);
//        std::swap(f1, f2);
//        std::swap(e01, e02);
    }
    // f0 < f1 < f2
    // check if edges are processed
    if (!(*pVecAuxMeshEdge)[e01].ptrReebArcId && !(*pVecAuxMeshEdge)[e01].ptrPos) {// both pointers are null
        // this edge is new
        CreateArc(Node0, Node1, e01);
    }
    if (!(*pVecAuxMeshEdge)[e12].ptrReebArcId && !(*pVecAuxMeshEdge)[e12].ptrPos) {// both pointers are null
        // this edge is new
        CreateArc(Node1, Node2, e12);
    }
    if (!(*pVecAuxMeshEdge)[e02].ptrReebArcId && !(*pVecAuxMeshEdge)[e02].ptrPos) {// both pointers are null
        // this edge is new
        CreateArc(Node0, Node2, e02);
    }


    //merge paths
    MergePaths(Node0, Node1, Node2,
               e02, e12, e01);

    //
    if (AddedTriangle[meshData->vecEdge[e01].AdjTri[0]] &&
        AddedTriangle[meshData->vecEdge[e01].AdjTri[1]]) {// remove edge
        DeleteArcListNodeForEdge(e01);
    }
    if (AddedTriangle[meshData->vecEdge[e12].AdjTri[0]] &&
        AddedTriangle[meshData->vecEdge[e12].AdjTri[1]]) {// remove edge
        DeleteArcListNodeForEdge(e12);
    }
    if (AddedTriangle[meshData->vecEdge[e02].AdjTri[0]] &&
        AddedTriangle[meshData->vecEdge[e02].AdjTri[1]]) {// remove edge
        DeleteArcListNodeForEdge(e02);
    }
    return;
}

void psbmReebGraph::AssignData(_SimpleMesh *inMesh, double *scalarField) {
    meshData = inMesh;
    MeshVertexScalarValue = scalarField;
}
//
/*********************************************************************************/
// NEED TO COPY OUT FOR CODE FILE CONSISTENCY
int psbmReebGraph::NonzeroArcValenceInSimplifiedArcGraphForNode(const int nodeId,
                                                                bool UpwardDir) {// Assumption : with phantom nonzero value -1 for arc conncecting to paired node
    std::list<class psbmReebArc *>::iterator listIter;
    int valence = 0;
    if (UpwardDir) {
        for (listIter = (*pVecReebNode)[nodeId].ptrArcUpId->begin();
             listIter != (*pVecReebNode)[nodeId].ptrArcUpId->end();
             listIter++) {
            if ((*listIter)->clusterLabel)
                valence++;
        }
    } else {
        for (listIter = (*pVecReebNode)[nodeId].ptrArcDownId->begin();
             listIter != (*pVecReebNode)[nodeId].ptrArcDownId->end();
             listIter++) {
            if ((*listIter)->clusterLabel)
                valence++;
        }
    }
    return valence;
}

void psbmReebGraph::GoAlong_1_or_2_degree_ReebArc(class psbmReebArc *&travArc,
                                                  const int nodeId,
                                                  const int arc_label,
                                                  bool UpwardDir) {// either v--, or ---v---, or --v
    std::list<class psbmReebArc *>::iterator listIter;
    if (UpwardDir) {
        for (listIter = (*pVecReebNode)[nodeId].ptrArcUpId->begin();
             listIter != (*pVecReebNode)[nodeId].ptrArcUpId->end();
             listIter++) {
            if ((*listIter)->clusterLabel == arc_label) {
                travArc = *listIter;
                break;
            }
        }
    } else {
        for (listIter = (*pVecReebNode)[nodeId].ptrArcDownId->begin();
             listIter != (*pVecReebNode)[nodeId].ptrArcDownId->end();
             listIter++) {
            if ((*listIter)->clusterLabel == arc_label) {
                travArc = *listIter;
                break;
            }
        }
    }
}

int psbmReebGraph::ResetArcLabel(const int nodeId, const int arcId, std::set<int> &pairedCritNodes) {
    bool searchDir_is_upward = false;
    int opposite_node_id = (*pVecSimplifiedArc)[arcId - 1].nCriticalNode0;
    //
    if ((*pVecSimplifiedArc)[arcId - 1].nCriticalNode0 == nodeId) {
        opposite_node_id = (*pVecSimplifiedArc)[arcId - 1].nCriticalNode1;
        searchDir_is_upward = true;
    }
    //
    // check if this node is paired node or not
    //
    int resetValue = arcId; // as phantom nonzero value for arc conncecting to paired node
    std::set<int>::iterator sFindIter;
    sFindIter = pairedCritNodes.find(opposite_node_id);
    if (sFindIter == pairedCritNodes.end()) {// it is NOT paired
        resetValue = 0;
    } else {// paired node , keep the node connected to paired node as arc with value -1
        // don't do anything
        criticalNodeIdSet->insert(nodeId);
        return opposite_node_id;//
    }
    // go over the arc to set all arc label as zero
    int search_node = -1;
    class psbmReebArc *travArc;
    GoAlong_1_or_2_degree_ReebArc(travArc, nodeId, arcId, searchDir_is_upward);
    // set zero arc label
    travArc->clusterLabel = resetValue;
    //
    if (searchDir_is_upward)
        search_node = travArc->nNodeId1;
    else
        search_node = travArc->nNodeId0;
    while (search_node != opposite_node_id) {
        GoAlong_1_or_2_degree_ReebArc(travArc, search_node, arcId, searchDir_is_upward);
        // set zero arc label
        travArc->clusterLabel = resetValue;
        //
        if (searchDir_is_upward)
            search_node = travArc->nNodeId1;
        else
            search_node = travArc->nNodeId0;
    }
    //
    // remove the arc from the simplified arc array
    //
    (*pVecSimplifiedArc)[arcId - 1].edgeLabel = 0;
    //
    return opposite_node_id;
}

void psbmReebGraph::ResetArcLabel(const int nodeId, const int arcId, const int newArcLabel) {
    bool searchDir_is_upward = false;
    int opposite_node_id = (*pVecSimplifiedArc)[arcId - 1].nCriticalNode0;
    //
    if ((*pVecSimplifiedArc)[arcId - 1].nCriticalNode0 == nodeId) {
        opposite_node_id = (*pVecSimplifiedArc)[arcId - 1].nCriticalNode1;
        searchDir_is_upward = true;
    }
    //
    // go over the arc to set all arc label as zero
    int search_node = -1;
    class psbmReebArc *travArc;
    GoAlong_1_or_2_degree_ReebArc(travArc, nodeId, arcId, searchDir_is_upward);
    // set zero arc label
    travArc->clusterLabel = newArcLabel;
    //
    if (searchDir_is_upward)
        search_node = travArc->nNodeId1;
    else
        search_node = travArc->nNodeId0;
    while (search_node != opposite_node_id) {
        GoAlong_1_or_2_degree_ReebArc(travArc, search_node, arcId, searchDir_is_upward);
        // set zero arc label
        travArc->clusterLabel = newArcLabel;
        //
        if (searchDir_is_upward)
            search_node = travArc->nNodeId1;
        else
            search_node = travArc->nNodeId0;
    }
    //
    return;
}

void psbmReebGraph::RemoveSimplifiedArcs() {// use a stack to simulate the depth first walking
    //
    std::stack<int> remove_reeb_node_stack;
    // intializing it with all minimum nodes
    for (std::set<int>::iterator sIter = criticalNodeIdSet->begin();
         sIter != criticalNodeIdSet->end();
         sIter++) {
        if ((*pVecReebNode)[*sIter].crType == MINIMUM_REEB ||
            (*pVecReebNode)[*sIter].crType == MAXIMUM_REEB) {
            remove_reeb_node_stack.push(*sIter);
            //
        }
    }
    //
    std::set<int> pairedCritNodes;
    for (unsigned int i = 0; i < criticalPairing->size(); i++) {
        pairedCritNodes.insert((*criticalPairing)[i].first);
        pairedCritNodes.insert((*criticalPairing)[i].second);
    }
    //
    while (!remove_reeb_node_stack.empty()) {//
        int nodeId = remove_reeb_node_stack.top();
        int opposite_node_id = -1;
        int valence_up = 0;
        int valence_down = 0;
        //
        remove_reeb_node_stack.pop();
        criticalNodeIdSet->erase(nodeId);
        //
        // find the opposite node
        //
        int cur_arc_id = -1;
        // as all arc label > 0;
        std::list<class psbmReebArc *>::iterator listIter;
        for (listIter = (*pVecReebNode)[nodeId].ptrArcDownId->begin();
             listIter != (*pVecReebNode)[nodeId].ptrArcDownId->end();
             listIter++) {
            cur_arc_id = cur_arc_id > (*listIter)->clusterLabel ?
                         cur_arc_id : (*listIter)->clusterLabel;
        }
        for (listIter = (*pVecReebNode)[nodeId].ptrArcUpId->begin();
             listIter != (*pVecReebNode)[nodeId].ptrArcUpId->end();
             listIter++) {
            cur_arc_id = cur_arc_id > (*listIter)->clusterLabel ?
                         cur_arc_id : (*listIter)->clusterLabel;
        }
        //
        opposite_node_id = ResetArcLabel(nodeId, cur_arc_id, pairedCritNodes);
        //

        valence_up = NonzeroArcValenceInSimplifiedArcGraphForNode(opposite_node_id, true);
        valence_down = NonzeroArcValenceInSimplifiedArcGraphForNode(opposite_node_id, false);
        //
        if (valence_up + valence_down == 1) {// remove it as it is a leaf now
            // put it to the stack
            remove_reeb_node_stack.push(opposite_node_id);
        }
    } // end of stack
    return;
}

void psbmReebGraph::MergeTwoSimplifiedArcs() {
    //
    std::vector<int> delete_nodes;
    // iterate all crit nodes to find the node
    for (std::set<int>::iterator sIter = criticalNodeIdSet->begin();
         sIter != criticalNodeIdSet->end();
         sIter++) {//
        if (NonzeroArcValenceInSimplifiedArcGraphForNode(*sIter, true) == 1 &&
            NonzeroArcValenceInSimplifiedArcGraphForNode(*sIter, false) == 1) {
            delete_nodes.push_back(*sIter);
        }
    }
    // merge the two arcs adjacent the 2-degree nodes
    for (unsigned int i = 0; i < delete_nodes.size(); i++) {
        int nodeId = delete_nodes[i];
        int up_arc_id = 0;
        int down_arc_id = 0;
        std::list<class psbmReebArc *>::iterator listIter;
        for (listIter = (*pVecReebNode)[nodeId].ptrArcDownId->begin();
             listIter != (*pVecReebNode)[nodeId].ptrArcDownId->end();
             listIter++) {
            down_arc_id = down_arc_id > (*listIter)->clusterLabel ?
                          down_arc_id : (*listIter)->clusterLabel;
        }
        for (listIter = (*pVecReebNode)[nodeId].ptrArcUpId->begin();
             listIter != (*pVecReebNode)[nodeId].ptrArcUpId->end();
             listIter++) {
            up_arc_id = up_arc_id > (*listIter)->clusterLabel ?
                        up_arc_id : (*listIter)->clusterLabel;
        }
        //
        if (up_arc_id == 0 || down_arc_id == 0) {
            std::cout << "WRONG IN MERGE" << std::endl;
            exit(0);
        }
        // merge up_arc_id to down_arc_id
        ResetArcLabel(nodeId, up_arc_id, down_arc_id);
        // remove node from crit nodes set
        criticalNodeIdSet->erase(nodeId);
        // remove the arc
        (*pVecSimplifiedArc)[up_arc_id - 1].edgeLabel = 0;
        // reset the nodes on the new arclabel
        (*pVecSimplifiedArc)[down_arc_id - 1].nCriticalNode1 =
                (*pVecSimplifiedArc)[up_arc_id - 1].nCriticalNode1;
    }
    //
    delete_nodes.clear();
    // remove the remain leaf adjacent to paired nodes
    //for (std::set<int>::iterator sIter = criticalNodeIdSet->begin();
    //	sIter != criticalNodeIdSet->end();
    //	sIter++)
    //{//
    //	if (NonzeroArcValenceInSimplifiedArcGraphForNode(*sIter, true) +
    //		NonzeroArcValenceInSimplifiedArcGraphForNode(*sIter, false) == 1 )
    //	{
    //		delete_nodes.push_back(*sIter);
    //	}
    //}
    //remove adjacent the 1-degree nodes
    //for (unsigned int i = 0; i < delete_nodes.size(); i++)
    //{
    //	int nodeId = delete_nodes[i];
    //	int up_arc_id = 0;
    //	int down_arc_id = 0;
    //	std::list<class psbmReebArc*>::iterator listIter;
    //	for (listIter = (*pVecReebNode)[nodeId].ptrArcDownId->begin();
    //		listIter != (*pVecReebNode)[nodeId].ptrArcDownId->end();
    //		listIter++)
    //	{
    //		down_arc_id = down_arc_id > (*listIter)->clusterLabel ?
    //					down_arc_id : (*listIter)->clusterLabel ;
    //	}
    //	for (listIter = (*pVecReebNode)[nodeId].ptrArcUpId->begin();
    //		listIter != (*pVecReebNode)[nodeId].ptrArcUpId->end();
    //		listIter++)
    //	{
    //		up_arc_id = up_arc_id > (*listIter)->clusterLabel ?
    //					up_arc_id : (*listIter)->clusterLabel ;
    //	}
    //	// one must be zero
    //	if (down_arc_id > 1)
    //	{
    //		up_arc_id = down_arc_id;
    //	}
    //	// merge up_arc_id to down_arc_id
    //	ResetArcLabel(nodeId, up_arc_id, 0);
    //	// remove node from crit nodes set
    //	criticalNodeIdSet->erase(nodeId);
    //	// remove the arc
    //	(*pVecSimplifiedArc)[up_arc_id - 1].edgeLabel = 0;
    //}


}

void psbmReebGraph::ResetSimplifiedArcIds() {
    std::vector<struct SimplifiedReebGraphArc> newSimplifiedArcSet;
    //
    //
    int sArcCounter = 1;
    for (unsigned int i = 0; i < pVecSimplifiedArc->size(); i++) {
        if ((*pVecSimplifiedArc)[i].edgeLabel) {
            ResetArcLabel((*pVecSimplifiedArc)[i].nCriticalNode0,
                          (*pVecSimplifiedArc)[i].edgeLabel,
                          sArcCounter);
            //
            newSimplifiedArcSet.push_back((*pVecSimplifiedArc)[i]);
            newSimplifiedArcSet.back().edgeLabel = sArcCounter;
            sArcCounter++;
        }
    }
    // update the pVecSimplifiedARc
    delete pVecSimplifiedArc;
    //
    pVecSimplifiedArc = new std::vector<struct SimplifiedReebGraphArc>;
    *pVecSimplifiedArc = newSimplifiedArcSet;
    //
    for (std::set<int>::iterator sIter = criticalNodeIdSet->begin();
         sIter != criticalNodeIdSet->end(); sIter++) {
        int valence_up = NonzeroArcValenceInSimplifiedArcGraphForNode((*sIter), true);
        int valence_down = NonzeroArcValenceInSimplifiedArcGraphForNode((*sIter), false);
        if (valence_up < 1 && valence_down == 1 ||
            valence_down < 1 && valence_up == 1) {
            std::cout << *sIter << std::endl;
            std::cout << "WRONG in after pairing simplifying " << std::endl;
        }

        if (
                (*pVecReebNode)[*sIter].ptrArcUpId->size() > valence_up) {
            std::list<class psbmReebArc *>::iterator sFindIter;
            for (sFindIter = (*pVecReebNode)[*sIter].ptrArcUpId->begin();
                 sFindIter != (*pVecReebNode)[*sIter].ptrArcUpId->end();
                 sFindIter++) {
                if ((*sFindIter)->clusterLabel == 0)
                    break;
            }
            (*pVecReebNode)[*sIter].ptrArcUpId->erase(sFindIter);
        }
        if (
                (*pVecReebNode)[*sIter].ptrArcDownId->size() > valence_down) {
            std::list<class psbmReebArc *>::iterator sFindIter;
            for (sFindIter = (*pVecReebNode)[*sIter].ptrArcDownId->begin();
                 sFindIter != (*pVecReebNode)[*sIter].ptrArcDownId->end();
                 sFindIter++) {
                if ((*sFindIter)->clusterLabel == 0)
                    break;
            }
            (*pVecReebNode)[*sIter].ptrArcDownId->erase(sFindIter);
        }
    }
}

void psbmReebGraph::AfterPairingSimplifyingForSimplifiedArcs() {
    RemoveSimplifiedArcs();
    MergeTwoSimplifiedArcs();
    ResetSimplifiedArcIds();
    return;
}
/*********************************************************************************/
// NEED TO COPY OUT FOR CODE FILE CONSISTENCY
void psbmReebGraph::ComputeReebGraph() {/**/
    std::cout << "Time for computing Reeb graph :" << std::endl;
    //
    AddedTriangle.resize(meshData->vecTriangle.size());
    for (unsigned int i = 0; i < meshData->vecTriangle.size(); i++) {
        AddedTriangle[i] = false;
    }
    //
    //std::vector<int> BFS_triangle_order;
    //BFS_triangle_order.reserve(1.2 * meshData->vecTriangle.size());
    ////
    //std::queue<int> triQ;
    //int edgeIdx[3] = {0};
    //int opposite_tri = 0;
    ////
    //triQ.push(0);
    //AddedTriangle[0] = true;
    //while (!triQ.empty())
    //{
    //	int curTriIdx = triQ.front();
    //	triQ.pop();
    //	BFS_triangle_order.push_back(curTriIdx);
    //	//
    //	edgeIdx[0] = meshData->vecTriangle[curTriIdx].e01;
    //	edgeIdx[1] = meshData->vecTriangle[curTriIdx].e02;
    //	edgeIdx[2] = meshData->vecTriangle[curTriIdx].e12;
    //	for (int i = 0; i < 3; i++)
    //	{
    //		opposite_tri = meshData->vecEdge[edgeIdx[i]].AdjTri[0] + meshData->vecEdge[edgeIdx[i]].AdjTri[1];
    //		opposite_tri = opposite_tri - curTriIdx;
    //		if (!AddedTriangle[opposite_tri])
    //		{
    //			triQ.push(opposite_tri);
    //			AddedTriangle[opposite_tri] = true;
    //		}
    //	}
    //}
    ////
    //for (unsigned int i = 0 ; i < meshData->vecTriangle.size(); i++)
    //{
    //	AddedTriangle[i] = false;
    //}
    //if (BFS_triangle_order.size() != meshData->vecTriangle.size())
    //{
    //	std::cout << "some triangles not VISITED " << std::endl;
    //	exit(0);
    //}
    //memset(&AddedTriangle[0], 0, sizeof(bool) * meshData->vecTriangle.size());
    //
    int vertexId[3] = {0};
    InitReebNodes((int) meshData->vecVertex.size());
    if (!pListReebArc)
        pListReebArc = new std::set<class psbmReebArc *>;
    for (unsigned int i = 0; i < meshData->vecTriangle.size(); i++)
        //for (unsigned int t = 0; t < BFS_triangle_order.size(); t++)
    {
        int t = i;//BFS_triangle_order[t];
        vertexId[0] = (meshData->vecTriangle)[i].v0;
        vertexId[1] = (meshData->vecTriangle)[i].v1;
        vertexId[2] = (meshData->vecTriangle)[i].v2;
        AddedTriangle[i] = true;
        AddMeshTriangle(vertexId[0], MeshVertexScalarValue[vertexId[0]],
                        vertexId[1], MeshVertexScalarValue[vertexId[1]],
                        vertexId[2], MeshVertexScalarValue[vertexId[2]],
                        (meshData->vecTriangle)[i].e01,
                        (meshData->vecTriangle)[i].e12,
                        (meshData->vecTriangle)[i].e02);
        //
        //if (t % 100000 == 1)
        //	std::cout << " t " << t << " vs " << meshData->vecTriangle.size() << std::endl;

        //CheckConsistency(0);
    }

    //WriteReebGraph("filigree_1mt.reeb");
    //exit(8);
    //ReadReebGraph("filigree_1mt.reeb");
    ////CheckConsistency(1);
    //



}

/******maximum spanning tree****************************/
void psbmReebGraph::ComputingCycle_max_tree() {
    vecVertexMappingToSimplifiedReebArc.resize(meshData->vecVertex.size());
    for (unsigned int i = 0; i < vecVertexMappingToSimplifiedReebArc.size(); i++) {
        vecVertexMappingToSimplifiedReebArc[i] = -1;
    }
//    std::cout << "Simplifying " << std::endl;
    if (SimplifyReebGraph()) {// all nodes have degree <= 3
//        std::cout << "Pairing on original graph" << std::endl;
    }
    /******************/
    _simpGraph_max_tree maxGraph;
    std::vector<std::pair<int, std::pair<int, int> > > max_spanning_tree;
    std::vector<std::pair<int, std::pair<int, int> > > non_tree_edges;
    //
    std::map<int, int> critNodeToGraphIndex;
    float cur_edge_weights = 0.0f;
    // visit all simplified arc to construct a graph
    for (unsigned int i = 0; i < pVecSimplifiedArc->size(); i++) {
        std::map<int, int>::iterator mIter;
        int src = (*pVecSimplifiedArc)[i].nCriticalNode0;
        int dst = (*pVecSimplifiedArc)[i].nCriticalNode1;
        int arcId = (*pVecSimplifiedArc)[i].edgeLabel - 1;
        // update the edge weights and maximum edge weights
        cur_edge_weights = MeshVertexScalarValue[src];
        //
        mIter = critNodeToGraphIndex.find(src);
        if (mIter == critNodeToGraphIndex.end()) {
            maxGraph.AddNode(src);
            critNodeToGraphIndex[src] = maxGraph.vecNode.size() - 1;
            src = maxGraph.vecNode.size() - 1;
        } else {
            src = mIter->second;
        }
        //
        mIter = critNodeToGraphIndex.find(dst);
        if (mIter == critNodeToGraphIndex.end()) {
            maxGraph.AddNode(dst);
            critNodeToGraphIndex[dst] = maxGraph.vecNode.size() - 1;
            dst = maxGraph.vecNode.size() - 1;
        } else {
            dst = mIter->second;
        }
        //
        maxGraph.AddEdgeW(src, dst, arcId, cur_edge_weights);
        //
    }
    //updating the weights
    maxGraph.KruskalAlg(max_spanning_tree, non_tree_edges);
    //
    _simpGraph_max_tree span_tree;
    span_tree.InitNodes(maxGraph.vecNode.size());
    for (unsigned int i = 0; i < max_spanning_tree.size(); i++) {
        span_tree.AddEdge(max_spanning_tree[i].second.first, max_spanning_tree[i].second.second,
                          max_spanning_tree[i].first);
    }
    //
    std::vector<int> path;
    std::vector<int> parents;
    //
    ////
    initVerticalLoops = new std::vector<psbmReebPath>;
    criticalPairing = new std::vector<std::pair<int, int> >;
    ////
    //
    std::set<std::pair<int, float>, myIntFloatPairCompare> sortedNonTreeEdges;
    //
    std::set<int> critNodesCounter;
    for (unsigned int i = 0; i < non_tree_edges.size(); i++) {
        int lowNodeId = non_tree_edges[i].second.first;
        if (//maxGraph.vecNode[non_tree_edges[i].first].degree() != 3 ||
                (*pVecReebNode)[maxGraph.vecNode[lowNodeId].nIndexInReebGraph].crType != UP_FORKING_REEB) {
            std::cout << "NOT a upforking node" << std::endl;
            exit(2);
        }
        int reebNodeIndex = maxGraph.vecNode[lowNodeId].nIndexInReebGraph;
        sortedNonTreeEdges.insert(std::pair<int, float>(i, MeshVertexScalarValue[reebNodeIndex]));
        //
        critNodesCounter.insert(lowNodeId);
    }
    if (critNodesCounter.size() != non_tree_edges.size()) {
        std::cout << "common lowest base pt" << std::endl;
        exit(0);
    }
    //
    float highNodeWeight = 0.0f;
    int highNodeIndex = 0;
    for (std::set<std::pair<int, float>, myIntFloatPairCompare>::iterator sIter = sortedNonTreeEdges.begin();
         sIter != sortedNonTreeEdges.end(); sIter++) {
        //for (unsigned int i = 0; i < non_tree_edges.size(); i++)
        //{
        int non_tree_edge_index = sIter->first;
        //
        psbmReebPath cyclePairing;
        //
        path.clear();
        span_tree.BreathFirstSearch(non_tree_edges[non_tree_edge_index].second.second, parents);
        int current = non_tree_edges[non_tree_edge_index].second.first; // lowest point
        while (current != -1) {
            // push this point into path
            path.push_back(current);
            current = parents[current];
        }
        //path.push_back(non_tree_edges[non_tree_edge_index].second.first);
        //
        if (path.size() == 1) {
            std::cout << "isolated points " << std::endl;
            exit(2);
        }
        highNodeIndex = maxGraph.vecNode[path[1]].nIndexInReebGraph;
        highNodeWeight = MeshVertexScalarValue[highNodeIndex];
        for (unsigned int j = 0; j < path.size() - 1; j++) {
            cyclePairing.nodeIndexSet.push_back(maxGraph.vecNode[path[j]].nIndexInReebGraph);
            cyclePairing.simpArcIndexSet.push_back(span_tree.ArcId(path[j], path[j + 1]));

            if (MeshVertexScalarValue[cyclePairing.nodeIndexSet.back()] > highNodeWeight) {
                highNodeIndex = cyclePairing.nodeIndexSet.back();
                highNodeWeight = MeshVertexScalarValue[highNodeIndex];
            }
        }
        //
        cyclePairing.nodeIndexSet.push_back(maxGraph.vecNode[path.back()].nIndexInReebGraph);
        cyclePairing.simpArcIndexSet.push_back(non_tree_edges[non_tree_edge_index].first);
        cyclePairing.nodeIndexSet.push_back(cyclePairing.nodeIndexSet.front());
        //
        // create pairing
        criticalPairing->push_back(std::pair<int, int>(highNodeIndex, cyclePairing.nodeIndexSet.front()));
        //

        std::vector<_Polygon> levelSetPolygons;
        cyclePairing.pathType = DetermineVerticalLoopType(cyclePairing.nodeIndexSet.front(), levelSetPolygons);
        //
        initVerticalLoops->push_back(cyclePairing);
        //
    }
    //
    return;
}

/**********************************/
void psbmReebGraph::ComputeCycleAndPairing() {
    vecVertexMappingToSimplifiedReebArc.resize(meshData->vecVertex.size());
    for (unsigned int i = 0; i < vecVertexMappingToSimplifiedReebArc.size(); i++) {
        vecVertexMappingToSimplifiedReebArc[i] = -1;
    }
    //std::cout << "Simplifying " << std::endl;
    if (SimplifyReebGraph()) {// all nodes have degree <= 3
        //std::cout << "Pairing on original graph" << std::endl;
        ComputingPairing(1);
    } else {
        //std::cout << "Pairing on simplified graph" << std::endl;
        HandlingCritNodeWithDegreeGreaterThan3();
    }
    //for (unsigned int i = 0; i < vecVertexMappingToSimplifiedReebArc.size(); i++)
    //{
    //	if (vecVertexMappingToSimplifiedReebArc[i] == -1)
    //	{
    //		std::cout << "some vertices are not mapped to reeb graph" << std::endl;
    //	}
    //}
    //

    //
    ////
    //ExportMathematicaGraph();
    ////
    //std::cout << "True Pairing " << std::endl;
    //ComputingPairing(1);
    //std::cout << "Own Pairing" << std::endl;
    //ComputingPairing();
    //// compute the tree pairing
    ////ComputingPairing(1);
    //// check the size of the two pairing
    ///*if (trueCriticalPairing.size() != criticalPairing->size())
    //{
    //	std::cout << " WRONG IN PAIRING , ONE HAS TO BE WRONG " << std::endl;
    //}*/

    //std::cout << "Done Pairing " << std::endl;
    ////
    psbmReebPath cyclePairing;
    initVerticalLoops = new std::vector<psbmReebPath>;
    ////
    std::pair<int, int> pairing;
    for (int j = 0; j < int(criticalPairing->size()); j++) {
        pairing.first = (*criticalPairing)[j].second;
        pairing.second = (*criticalPairing)[j].first;
        TraceOutOneCycleForEachPairing(pairing, cyclePairing);
        //
        std::vector<_Polygon> levelSetPolygons;
        cyclePairing.pathType = DetermineVerticalLoopType(pairing.first, levelSetPolygons);
        //

        initVerticalLoops->push_back(cyclePairing);

        cyclePairing.nodeIndexSet.clear();
        cyclePairing.simpArcIndexSet.clear();
    }
    ////
    //SetSimplifiedArcFlagBits();
    ////
    //for (int j = 0; j < int(initVerticalLoops->size()); j++)
    //{
    //	if (j == 29)
    //	{
    //		for (int i = 0; i < int((*initVerticalLoops)[j].simpArcIndexSet.size()); i++)
    //		{

    //				std::cout << (*initVerticalLoops)[j].nodeIndexSet[i] << "--" << (*initVerticalLoops)[j].simpArcIndexSet[i] << "--" ;
    //		}
    //		std::cout << (*initVerticalLoops)[j].nodeIndexSet.back() << std::endl;
    //	}
    //}
}

bool psbmReebGraph::IsUpArcNodesMax(const int nodeIdx) {
    bool ret = false;
    if ((*pVecReebNode)[nodeIdx].ptrArcUpId->size() >= 1 && (*pVecReebNode)[nodeIdx].ptrArcDownId->size() == 1) {
        std::list<psbmReebArc *>::iterator listIter = (*pVecReebNode)[nodeIdx].ptrArcUpId->begin();
        int node = (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode1;
        REEB_GRAPH_CRITICAL_TYPE same_crit_type = (*pVecReebNode)[node].crType;
        bool bSameCriticalType = true;
        for (listIter++; listIter != (*pVecReebNode)[nodeIdx].ptrArcUpId->end(); listIter++) {
            node = (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode1;
            if ((*pVecReebNode)[node].crType != same_crit_type) {
                bSameCriticalType = false;
                break;
            }
        }
        if (bSameCriticalType && same_crit_type == MAXIMUM_REEB) {
            ret = true;
        }
    }
    return ret;
}

bool psbmReebGraph::IsDownArcNodesMin(const int nodeIdx) {
    bool ret = false;
    if ((*pVecReebNode)[nodeIdx].ptrArcDownId->size() >= 1 && (*pVecReebNode)[nodeIdx].ptrArcUpId->size() == 1) {
        std::list<psbmReebArc *>::iterator listIter = (*pVecReebNode)[nodeIdx].ptrArcDownId->begin();
        int node = (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode0;
        REEB_GRAPH_CRITICAL_TYPE same_crit_type = (*pVecReebNode)[node].crType;
        bool bSameCriticalType = true;
        for (listIter++; listIter != (*pVecReebNode)[nodeIdx].ptrArcDownId->end(); listIter++) {
            node = (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode0;
            if ((*pVecReebNode)[node].crType != same_crit_type) {
                bSameCriticalType = false;
                break;
            }
        }
        //
        if (bSameCriticalType && same_crit_type == MINIMUM_REEB) {//
            ret = true;
        }
    }
    return ret;
}

void psbmReebGraph::HandlingCritNodeWithDegreeGreaterThan3() {
    //
    _simpGraph simplifiedReebGraph;

    simplifiedReebGraph.InitNodes(criticalNodeIdSet->size());
    //
    //std::map<int, int> CritNodeToGraphNode;
    boost::unordered_map<int, int> CritNodeToGraphNode;

    std::vector<int> vecSingularCriticalNodes;
    vecSingularCriticalNodes.reserve(criticalNodeIdSet->size());
    //
    int nodeIter = 0;
    for (std::set<int>::iterator sIter = criticalNodeIdSet->begin();
         sIter != criticalNodeIdSet->end(); sIter++) {
        simplifiedReebGraph.vecNode[nodeIter].selfIndex = nodeIter;
        simplifiedReebGraph.vecNode[nodeIter].color = *sIter;
        //
        CritNodeToGraphNode[*sIter] = nodeIter++;
        //
        if ((*pVecReebNode)[*sIter].isCritical && (*pVecReebNode)[*sIter].crType == UNKNOWN_REEB) {
            vecSingularCriticalNodes.push_back(*sIter);
            //
            //for (std::list<psbmReebArc*>::iterator listIter = (*pVecReebNode)[*sIter].ptrArcUpId->begin();
            //	listIter != (*pVecReebNode)[*sIter].ptrArcUpId->end(); listIter++)
            //{
            //	int node = (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode1;
            //	if ((*pVecReebNode)[node].isCritical)
            //		std::cout <<"c";
            //	std::cout << (*pVecReebNode)[node].crType << " ";
            //}
            //std::cout << "\t";
            //for (std::list<psbmReebArc*>::iterator listIter = (*pVecReebNode)[*sIter].ptrArcDownId->begin();
            //	listIter != (*pVecReebNode)[*sIter].ptrArcDownId->end(); listIter++)
            //{
            //	int node = (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode0;
            //	if ((*pVecReebNode)[node].isCritical)
            //		std::cout <<"c";
            //	std::cout << (*pVecReebNode)[node].crType << " ";
            //}
            //std::cout << std::endl;
        }
    }
    //
    for (unsigned int i = 0; i < pVecSimplifiedArc->size(); i++) {
        simplifiedReebGraph.AddEdge(CritNodeToGraphNode[(*pVecSimplifiedArc)[i].nCriticalNode0],
                                    CritNodeToGraphNode[(*pVecSimplifiedArc)[i].nCriticalNode1]);
    }
    //
    for (std::vector<int>::iterator vIter = vecSingularCriticalNodes.begin();
         vIter != vecSingularCriticalNodes.end(); vIter++) {
        if ((*pVecReebNode)[*vIter].isCritical && (*pVecReebNode)[*vIter].crType == UNKNOWN_REEB) {//
            int nCurrentNodeIdx = *vIter;
            if (IsUpArcNodesMax(nCurrentNodeIdx)) {
                do {//
                    for (std::list<psbmReebArc *>::iterator listIter = (*pVecReebNode)[nCurrentNodeIdx].ptrArcUpId->begin();
                         listIter != (*pVecReebNode)[nCurrentNodeIdx].ptrArcUpId->end(); listIter++) {
                        int node = (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode1;
                        simplifiedReebGraph.RemoveEdge(CritNodeToGraphNode[nCurrentNodeIdx], CritNodeToGraphNode[node]);
                    }
                    (*pVecReebNode)[nCurrentNodeIdx].isCritical = false;
                    (*pVecReebNode)[nCurrentNodeIdx].crType = MAXIMUM_REEB;
                    //
                    nCurrentNodeIdx = (*pVecReebNode)[nCurrentNodeIdx].ptrArcDownId->front()->nNodeId0;
                    //
                } while (IsUpArcNodesMax(nCurrentNodeIdx));
            } else {
                while (IsDownArcNodesMin(nCurrentNodeIdx)) {
                    //
                    for (std::list<psbmReebArc *>::iterator listIter = (*pVecReebNode)[nCurrentNodeIdx].ptrArcDownId->begin();
                         listIter != (*pVecReebNode)[nCurrentNodeIdx].ptrArcDownId->end(); listIter++) {
                        int node = (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode0;
                        simplifiedReebGraph.RemoveEdge(CritNodeToGraphNode[nCurrentNodeIdx], CritNodeToGraphNode[node]);
                    }
                    //
                    (*pVecReebNode)[nCurrentNodeIdx].isCritical = false;
                    (*pVecReebNode)[nCurrentNodeIdx].crType = MINIMUM_REEB;
                    //
                    nCurrentNodeIdx = (*pVecReebNode)[nCurrentNodeIdx].ptrArcUpId->front()->nNodeId1;
                }
            }
        }
    }
    //
    for (std::vector<int>::iterator vIter = vecSingularCriticalNodes.begin();
         vIter != vecSingularCriticalNodes.end(); vIter++) {
        if ((*pVecReebNode)[*vIter].isCritical && (*pVecReebNode)[*vIter].crType == UNKNOWN_REEB) {//
            std::vector<int> upMaxNum;
            std::vector<int> downMinNum;
            //
            upMaxNum.reserve(5);
            downMinNum.reserve(5);
            //
            for (std::list<psbmReebArc *>::iterator listIter = (*pVecReebNode)[*vIter].ptrArcUpId->begin();
                 listIter != (*pVecReebNode)[*vIter].ptrArcUpId->end(); listIter++) {
                int node = (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode1;
                if ((*pVecReebNode)[node].crType == MAXIMUM_REEB) {
                    upMaxNum.push_back(node);
                }
            }
            //
            for (std::list<psbmReebArc *>::iterator listIter = (*pVecReebNode)[*vIter].ptrArcDownId->begin();
                 listIter != (*pVecReebNode)[*vIter].ptrArcDownId->end(); listIter++) {
                int node = (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode0;
                if ((*pVecReebNode)[node].crType == MINIMUM_REEB) {
                    downMinNum.push_back(node);
                }
            }
            //
            const int upNodeSize = (*pVecReebNode)[*vIter].ptrArcUpId->size();
            const int downNodeSize = (*pVecReebNode)[*vIter].ptrArcDownId->size();
            const int curNodeOnGraph = CritNodeToGraphNode[*vIter];
            //
            if (upNodeSize - upMaxNum.size() <= 2) {
                if (downNodeSize - downMinNum.size() <= 2) {
                    switch (upNodeSize - upMaxNum.size()) {
                        case 0:
                            if (downNodeSize - downMinNum.size() == 2) { // CASE : 0, 2
                                (*pVecReebNode)[*vIter].crType = DOWN_FORKING_REEB;
                                (*pVecReebNode)[*vIter].isCritical = false;
                                //
                                if (upNodeSize > 0) {
                                    for (unsigned int u = 0; u < upMaxNum.size() - 1; u++)
                                        simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                       CritNodeToGraphNode[upMaxNum[u]]);
                                } else {
                                    std::cout << "up 0, down 2, cannot handle" << std::endl;
                                    exit(6);
                                }
                                //
                                for (unsigned int d = 0; d < downMinNum.size(); d++)
                                    simplifiedReebGraph.RemoveEdge(curNodeOnGraph, CritNodeToGraphNode[downMinNum[d]]);
                            } else {//
                                if (downNodeSize - downMinNum.size() == 1) {//CASE : 0, 1
                                    (*pVecReebNode)[*vIter].crType = DOWN_FORKING_REEB;
                                    (*pVecReebNode)[*vIter].isCritical = false;
                                    //
                                    if (upNodeSize > 0) {
                                        for (unsigned int u = 0; u < upMaxNum.size() - 1; u++)
                                            simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                           CritNodeToGraphNode[upMaxNum[u]]);
                                    } else {
                                        std::cout << "up 0, down 1, cannot handle" << std::endl;
                                        exit(6);
                                    }
                                    //
                                    if (downNodeSize > 1) {
                                        for (unsigned int d = 0; d < downMinNum.size() - 1; d++)
                                            simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                           CritNodeToGraphNode[downMinNum[d]]);
                                    } else {
                                        std::cout << "up 0, down <= 1, cannot handle " << std::endl;
                                        exit(7);
                                    }

                                } else {// (downNodeSize - downMinNum.size() == 0)
                                    std::cout << "Both zeros , CAN NOT handle" << std::endl;
                                    exit(6);
                                }
                            }
                            break;
                        case 1:
                            if (downNodeSize - downMinNum.size() == 2) { //CASE : 1, 2
                                (*pVecReebNode)[*vIter].crType = DOWN_FORKING_REEB;
                                (*pVecReebNode)[*vIter].isCritical = false;
                                //
                                for (unsigned int u = 0; u < upMaxNum.size(); u++)
                                    simplifiedReebGraph.RemoveEdge(curNodeOnGraph, CritNodeToGraphNode[upMaxNum[u]]);
                                //
                                for (unsigned int d = 0; d < downMinNum.size(); d++)
                                    simplifiedReebGraph.RemoveEdge(curNodeOnGraph, CritNodeToGraphNode[downMinNum[d]]);

                            } else {//
                                if (downNodeSize - downMinNum.size() == 1) {// CASE : 1, 1
                                    if (upNodeSize > 1) {//
                                        (*pVecReebNode)[*vIter].crType = UP_FORKING_REEB;
                                        (*pVecReebNode)[*vIter].isCritical = false;
                                        //
                                        for (unsigned int u = 0; u < upMaxNum.size() - 1; u++)
                                            simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                           CritNodeToGraphNode[upMaxNum[u]]);
                                        //
                                        for (unsigned int d = 0; d < downMinNum.size(); d++)
                                            simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                           CritNodeToGraphNode[downMinNum[d]]);
                                    } else {
                                        if (downNodeSize > 1) {
                                            (*pVecReebNode)[*vIter].crType = DOWN_FORKING_REEB;
                                            (*pVecReebNode)[*vIter].isCritical = false;
                                            //
                                            for (unsigned int d = 0; d < downMinNum.size() - 1; d++)
                                                simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                               CritNodeToGraphNode[downMinNum[d]]);
                                        } else {//
                                            std::cout << "middle Node, cannot handle" << std::endl;
                                            exit(5);
                                        }
                                    }
                                } else {// (downNodeSize - downMinNum.size() == 0)
                                    if (upNodeSize > 1) {//CASE : 1, 0
                                        (*pVecReebNode)[*vIter].crType = UP_FORKING_REEB;
                                        (*pVecReebNode)[*vIter].isCritical = false;
                                        //
                                        for (unsigned int u = 0; u < upMaxNum.size() - 1; u++)
                                            simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                           CritNodeToGraphNode[upMaxNum[u]]);
                                        //
                                        if (downNodeSize > 0) {
                                            for (unsigned int d = 0; d < downMinNum.size() - 1; d++)
                                                simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                               CritNodeToGraphNode[downMinNum[d]]);
                                        } else {
                                            std::cout << "up 1 , down 0, cannot handle" << std::endl;
                                            exit(6);
                                        }
                                    } else {
                                        std::cout << "up 1 , down 0, cannot handle" << std::endl;
                                        exit(6);
                                    }
                                }
                            }
                            break;
                        case 2:
                            if (downNodeSize - downMinNum.size() == 2) {//CASE : 2, 2
                                std::cout << "Impossible to handle node with degree >= 4" << std::endl;
                                exit(7);
                            } else {//
                                if (downNodeSize - downMinNum.size() == 1) {//CASE : 2, 1
                                    (*pVecReebNode)[*vIter].crType = UP_FORKING_REEB;
                                    (*pVecReebNode)[*vIter].isCritical = false;
                                    //
                                    for (unsigned int u = 0; u < upMaxNum.size(); u++)
                                        simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                       CritNodeToGraphNode[upMaxNum[u]]);
                                    //
                                    for (unsigned int d = 0; d < downMinNum.size(); d++)
                                        simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                       CritNodeToGraphNode[downMinNum[d]]);
                                } else {// (downNodeSize - downMinNum.size() == 0)
                                    //CASE : 2, 0
                                    (*pVecReebNode)[*vIter].crType = UP_FORKING_REEB;
                                    (*pVecReebNode)[*vIter].isCritical = false;
                                    //
                                    for (unsigned int u = 0; u < upMaxNum.size(); u++)
                                        simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                       CritNodeToGraphNode[upMaxNum[u]]);
                                    //
                                    if (downNodeSize > 0) {
                                        for (unsigned int d = 0; d < downMinNum.size() - 1; d++)
                                            simplifiedReebGraph.RemoveEdge(curNodeOnGraph,
                                                                           CritNodeToGraphNode[downMinNum[d]]);
                                    } else {
                                        std::cout << "case : up 2, down 0; cannot handle " << std::endl;
                                        exit(5);
                                    }
                                }
                            }
                            break;
                    }
                } else {
                    std::cout << "Impossible to handle node with degree >= 4 inner" << std::endl;
                    exit(7);
                }
            } else {
                std::cout << "Impossible to handle node with degree >= 4" << std::endl;
                exit(7);
            }
        }
    }
    //exit(0);
    for (unsigned int i = 0; i < simplifiedReebGraph.vecNode.size(); i++) {
        if (simplifiedReebGraph.vecNode[i].adjList.size() > 3) {
            std::cout << "node with degree 3" << std::endl;
            exit(0);
        }
    }
    //std::cout << " All nodes with degree <= 3" << std::endl;
    //
    std::vector<std::pair<int, int> > grapNodePairing;
    ReebGraphPairing(simplifiedReebGraph, grapNodePairing);
    //
    //exit(0);
    return;
}

bool psbmReebGraph::SimplifyReebGraph() {
    bool ret = true;
    pVecSimplifiedArc = new std::vector<struct SimplifiedReebGraphArc>;
    std::set<class psbmReebArc *>::iterator sIter;

    criticalNodeIdSet = new std::set<int>;

    int curLabelId = 1; // 0 is reserved
    minHeightNodeId = 0;
    sIter = pListReebArc->begin();

    struct SimplifiedReebGraphArc tmpSimpArc;

    for (; sIter != pListReebArc->end(); sIter++) {
        if (!(*sIter)->clusterLabel) {// not assigned yet
            (*sIter)->clusterLabel = curLabelId;
            tmpSimpArc.edgeLabel = curLabelId;
            // down walk
            tmpSimpArc.nCriticalNode0 = DownWalkArcs(*sIter, curLabelId, ret);
            tmpSimpArc.nCriticalNode1 = UpWalkArcs(*sIter, curLabelId, ret);
            //
            //
            vecVertexMappingToSimplifiedReebArc[(*pVecReebNode)[tmpSimpArc.nCriticalNode0].nVertexId] = 0;
            vecVertexMappingToSimplifiedReebArc[(*pVecReebNode)[tmpSimpArc.nCriticalNode1].nVertexId] = 0;
            //
            //
            pVecSimplifiedArc->push_back(tmpSimpArc);
            //
            criticalNodeIdSet->insert(tmpSimpArc.nCriticalNode0);
            criticalNodeIdSet->insert(tmpSimpArc.nCriticalNode1);
            curLabelId++;
        }
    }
    std::cout << "criticalNodeIdSet size : " << criticalNodeIdSet->size() << std::endl;
    std::cout << "ARC size : " << pVecSimplifiedArc->size() << std::endl;
    return ret;

}

void psbmReebGraph::ExportMathematicaGraph() {
    std::map<int, int> graphNode;
    std::set<std::pair<int, int> > edgeSet;
    std::set<std::pair<int, int> >::iterator eIter;
    std::set<int>::iterator sIter;
    std::ofstream ofile;
    unsigned int pos = 0;
    for (sIter = criticalNodeIdSet->begin(); sIter != criticalNodeIdSet->end(); sIter++) {
        graphNode[*sIter] = pos++;
    }
    ofile.open("graph.txt");
    ofile << "Graph[{"; //1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3,
    //3 \[UndirectedEdge] 1}]"
    std::vector<struct SimplifiedReebGraphArc>::iterator vIter;
    pos = 0;
    for (vIter = pVecSimplifiedArc->begin();
         vIter != pVecSimplifiedArc->end();
         vIter++) {
        std::pair<int, int> tmpPair(vIter->nCriticalNode0, vIter->nCriticalNode1);
        eIter = edgeSet.find(tmpPair);
        //if (eIter == edgeSet.end())
        //	edgeSet.insert(tmpPair);
        //else
        //{
        //	tmpPair.first = vIter->nCriticalNode1;
        //	tmpPair.second = vIter->nCriticalNode0;
        //}
        ofile << tmpPair.first;//graphNode[tmpPair.first];
        ofile << " -> ";//"\\[UndirectedEdge]";//
        ofile << tmpPair.second;//graphNode[tmpPair.second];
        if (++pos < pVecSimplifiedArc->size())
            ofile << ",";
        ofile << "\n";
    }
    ofile << "}]";
    ofile << "\n" << criticalNodeIdSet->size() << " " << pVecSimplifiedArc->size() - criticalNodeIdSet->size() + 1 <<
          " " << pVecSimplifiedArc->size() << "\n";
    graphNode.clear();
    edgeSet.clear();
    ofile.close();
}

int psbmReebGraph::DownWalkArcs(class psbmReebArc *inArc, const int label, bool &ret) {
    class psbmReebArc *travArc = inArc;
    int nodeId = travArc->nNodeId0; // smaller index
    unsigned int valenceDown = (*pVecReebNode)[nodeId].ptrArcDownId->size();
    unsigned int valenceUp = (*pVecReebNode)[nodeId].ptrArcUpId->size();
    while (valenceDown == 1 && valenceUp == 1) {
        travArc = (*pVecReebNode)[nodeId].ptrArcDownId->front();
        travArc->clusterLabel = label;
        //
        vecVertexMappingToSimplifiedReebArc[(*pVecReebNode)[nodeId].nVertexId] = label;
        //
        nodeId = travArc->nNodeId0;
        valenceDown = (*pVecReebNode)[nodeId].ptrArcDownId->size();
        valenceUp = (*pVecReebNode)[nodeId].ptrArcUpId->size();
    }
    if ((valenceDown * valenceUp > 0 && valenceDown + valenceUp == 3) || (valenceDown + valenceUp == 1)) {
        if (valenceDown + valenceUp == 3) {
            if (valenceDown == 2) {
                (*pVecReebNode)[nodeId].crType = DOWN_FORKING_REEB;
            } else {
                (*pVecReebNode)[nodeId].crType = UP_FORKING_REEB;
            }
        } else {
            if (valenceDown == 1)
                (*pVecReebNode)[nodeId].crType = MAXIMUM_REEB;
            else {
                if (MeshVertexScalarValue[(*pVecReebNode)[nodeId].nVertexId] <
                    MeshVertexScalarValue[(*pVecReebNode)[minHeightNodeId].nVertexId])
                    minHeightNodeId = nodeId;
                (*pVecReebNode)[nodeId].crType = MINIMUM_REEB;
            }
        }
    } else {
        (*pVecReebNode)[nodeId].isCritical = true;
        ret = false;
        //std::cout << "Node " << nodeId << " has degree " << valenceDown + valenceUp << " is WRONG" << std::endl;
    }
    return nodeId;
}

int psbmReebGraph::UpWalkArcs(class psbmReebArc *inArc, const int label, bool &ret) {
    class psbmReebArc *travArc = inArc;
    int nodeId = travArc->nNodeId1; // lareger index
    int valenceDown = (*pVecReebNode)[nodeId].ptrArcDownId->size();
    int valenceUp = (*pVecReebNode)[nodeId].ptrArcUpId->size();
    while (valenceDown == 1 && valenceUp == 1) {
        travArc = (*pVecReebNode)[nodeId].ptrArcUpId->front();
        travArc->clusterLabel = label;
        //
        vecVertexMappingToSimplifiedReebArc[(*pVecReebNode)[nodeId].nVertexId] = label;
        //
        nodeId = travArc->nNodeId1; // larger index
        valenceDown = (*pVecReebNode)[nodeId].ptrArcDownId->size();
        valenceUp = (*pVecReebNode)[nodeId].ptrArcUpId->size();
    }
    if ((valenceDown * valenceUp > 0 && valenceDown + valenceUp == 3) || (valenceDown + valenceUp == 1)) {
        if (valenceDown + valenceUp == 3) {
            if (valenceDown == 2)
                (*pVecReebNode)[nodeId].crType = DOWN_FORKING_REEB;
            else
                (*pVecReebNode)[nodeId].crType = UP_FORKING_REEB;
        } else {
            if (valenceDown == 1)
                (*pVecReebNode)[nodeId].crType = MAXIMUM_REEB;
            else {
                if (MeshVertexScalarValue[(*pVecReebNode)[nodeId].nVertexId] <
                    MeshVertexScalarValue[(*pVecReebNode)[minHeightNodeId].nVertexId])
                    minHeightNodeId = nodeId;
                (*pVecReebNode)[nodeId].crType = MINIMUM_REEB;
            }
        }
    } else {
        (*pVecReebNode)[nodeId].isCritical = true;
        ret = false;
        //std::cout << "Node " << nodeId << " has degree " << valenceDown + valenceUp << " is WRONG" << std::endl;
        //exit(0);
    }
    return nodeId;
}

void psbmReebGraph::CheckConsistency(int flag) {
    int minNode = 0;
    struct psbmArcListNode *curPtr;
    struct psbmArcListNode *prePtr;
    int maxNode = 0;
    // check whether each edge has pointer to its highest arc
    for (unsigned int i = 0; i < pVecAuxMeshEdge->size(); i++) {
        if ((*pVecAuxMeshEdge)[i].ptrReebArcId) {
            if ((*pVecAuxMeshEdge)[i].ptrReebArcId != (*pVecAuxMeshEdge)[i].ptrPos->ptrReebArcId) {
                std::cout << "Arc identity NOT matched" << std::endl;
            }
            maxNode = meshData->vecEdge[i].v0;//(*ProcessedVertex)[meshData->vecEdge[i].v0];
            minNode = meshData->vecEdge[i].v1;//(*ProcessedVertex)[meshData->vecEdge[i].v1];
            if (maxNode != (*pVecAuxMeshEdge)[i].ptrReebArcId->nNodeId1
                &&
                minNode != (*pVecAuxMeshEdge)[i].ptrReebArcId->nNodeId1) {
                std::cout << "NOT highest arc" << std::endl;
            } else {// check the other endpoint
                if (minNode == (*pVecAuxMeshEdge)[i].ptrReebArcId->nNodeId1)
                    minNode = maxNode;

                // go along next pointer
                prePtr = (*pVecAuxMeshEdge)[i].ptrPos;
                curPtr = (*pVecAuxMeshEdge)[i].ptrPos;
                while (curPtr) {
                    if (curPtr->nEdgeId != i) {
                        std::cout << "Edge index NOT matched" << std::endl;
                    }
                    prePtr = curPtr;
                    curPtr = curPtr->ptrNextPos;
                }
                if (prePtr->ptrReebArcId->nNodeId0 != minNode) {
                    std::cout << "minVertex NOT matched" << std::endl;
                }

            }

        }
    }
    // check node-arc relation
    std::list<class psbmReebArc *>::iterator siter;
    for (unsigned int i = 0; i < pVecReebNode->size(); i++) {
        minNode = (*pVecReebNode)[i].nVertexId;//(*ProcessedVertex)[(*pVecReebNode)[i].nVertexId];
        if (!(*pVecReebNode)[i].ptrArcDownId->empty()) {
            for (siter = (*pVecReebNode)[i].ptrArcDownId->begin();
                 siter != (*pVecReebNode)[i].ptrArcDownId->end();
                 siter++) {
                if ((*siter)->nNodeId1 != minNode) {
                    std::cout << "downARC NOT matched for NODE" << std::endl;
                }
            }
        }
        if (!(*pVecReebNode)[i].ptrArcUpId->empty()) {
            for (siter = (*pVecReebNode)[i].ptrArcUpId->begin();
                 siter != (*pVecReebNode)[i].ptrArcUpId->end();
                 siter++) {
                if ((*siter)->nNodeId0 != minNode) {
                    std::cout << "upARC NOT matched for NODE" << std::endl;
                }
            }
        }
        if (flag) {
            if ((*pVecReebNode)[i].ptrArcUpId->size() +
                (*pVecReebNode)[i].ptrArcDownId->size() +
                (*pVecReebNode)[i].ptrHArc->size() > 3) {
                std::cout << "NONE manifold ??" << std::endl;
                int sin;
                std::cin >> sin;
            }
        }
    }
    return;
}

void psbmReebGraph::WalkFromHighLevelToLowLevel(std::set<int> &highPath, std::set<int> &lowPath,
                                                double &highHeight, int &highestNodeId,
                                                const double &lowHeight, const int &lowestNodeId,
                                                int &pairNode, int &highMINNode) {// no need to take care of termination
// eventually these pathes at least led to minimum node
    int highNodeId;
    int nodeId;
    // used for swaping
    std::set<int> tmpSet;
    // iterator
    std::set<int>::iterator highIter;
    std::set<int>::iterator lowIter;
    // return type for set insert
    std::pair<std::set<int>::iterator, bool> ret;
    std::list<class psbmReebArc *>::iterator arcIter;

    if (highPath.empty())
        return;
    while ((highHeight > lowHeight) ||
           (highHeight == lowHeight && (*pVecReebNode)[highestNodeId].nVertexId >
                                       (*pVecReebNode)[lowestNodeId].nVertexId))//highestNodeId > lowestNodeId
    {
        // initializing  highestHeight
        highNodeId = *(highPath.begin());
        if ((MeshVertexScalarValue[(*pVecReebNode)[highNodeId].nVertexId] < lowHeight) ||
            (MeshVertexScalarValue[(*pVecReebNode)[highNodeId].nVertexId] == lowHeight
             && (*pVecReebNode)[highNodeId].nVertexId <=
                (*pVecReebNode)[lowestNodeId].nVertexId)//highNodeId <= lowestNodeId
                ) {
            highHeight = MeshVertexScalarValue[(*pVecReebNode)[highNodeId].nVertexId];
            highestNodeId = highNodeId;
        } else {
            //if ((*pVecReebNode)[highNodeId].crType != MINIMUM_REEB)
            //{
            arcIter = (*pVecReebNode)[highNodeId].ptrArcDownId->begin();
            //nodeId = (*arcIter)->nNodeId0;
            nodeId = TheOtherEndPoint(highNodeId, (*arcIter)->clusterLabel);
            highHeight = MeshVertexScalarValue[(*pVecReebNode)[nodeId].nVertexId];
            highestNodeId = nodeId;
            //}
        }
        //
        while (!highPath.empty()) {
            highIter = highPath.begin();
            highNodeId = *highIter;
            //
            highPath.erase(highIter);
            //
            if ((MeshVertexScalarValue[(*pVecReebNode)[highNodeId].nVertexId] > lowHeight) ||
                (MeshVertexScalarValue[(*pVecReebNode)[highNodeId].nVertexId] == lowHeight &&
                 (*pVecReebNode)[highNodeId].nVertexId >
                 (*pVecReebNode)[lowestNodeId].nVertexId) //highNodeId > lowestNodeId
                    ) {
                //go down once
                for (arcIter = (*pVecReebNode)[highNodeId].ptrArcDownId->begin();
                     arcIter != (*pVecReebNode)[highNodeId].ptrArcDownId->end();
                     arcIter++) {

                    //nodeId = (*arcIter)->nNodeId0;
                    nodeId = TheOtherEndPoint(highNodeId, (*arcIter)->clusterLabel);
                    //
                    if ((MeshVertexScalarValue[(*pVecReebNode)[nodeId].nVertexId] > highHeight) ||
                        (MeshVertexScalarValue[(*pVecReebNode)[nodeId].nVertexId] == highHeight &&
                         (*pVecReebNode)[nodeId].nVertexId >
                         (*pVecReebNode)[highestNodeId].nVertexId) //nodeId > highestNodeId
                            ) {
                        highHeight = MeshVertexScalarValue[(*pVecReebNode)[nodeId].nVertexId];
                        highestNodeId = nodeId;
                    }
                    //
                    ret = tmpSet.insert(nodeId);
                    //
                    if ((*pVecReebNode)[nodeId].crType == UP_FORKING_REEB) {// check if it is also in right paths
                        lowIter = lowPath.find(nodeId);
                        if (lowIter != lowPath.end()) {// find the a pair , but need to check it is the lowest or not
                            if (pairNode < 0) {// first time
                                pairNode = nodeId;
                                //pairNodeType = UP_FORKING_REEB;
                            } else {
                                if ((MeshVertexScalarValue[(*pVecReebNode)[pairNode].nVertexId] <
                                     MeshVertexScalarValue[(*pVecReebNode)[nodeId].nVertexId]) ||
                                    (MeshVertexScalarValue[(*pVecReebNode)[pairNode].nVertexId] ==
                                     MeshVertexScalarValue[(*pVecReebNode)[nodeId].nVertexId] &&
                                     (*pVecReebNode)[pairNode].nVertexId <
                                     (*pVecReebNode)[nodeId].nVertexId) //pairNode < nodeId
                                        ) {
                                    pairNode = nodeId;
                                    //pairNodeType = UP_FORKING_REEB;
                                }
                            }
                        }
                    } else {// check if it is a minimu or not
                        if ((*pVecReebNode)[nodeId].crType == MINIMUM_REEB) {// check if it is also in right paths
                            lowIter = lowPath.find(nodeId);
                            if (lowIter !=
                                lowPath.end()) {// find the a pair , but need to check it is the lowest or not
                                std::cout << "WRONG minimum has degree more than 1" << std::endl;
                            }
                            if (highMINNode < 0) {// first time
                                highMINNode = nodeId;
                                //pairNodeType = MINIMUM_REEB;
                            } else {
                                if ((MeshVertexScalarValue[(*pVecReebNode)[highMINNode].nVertexId] <
                                     MeshVertexScalarValue[(*pVecReebNode)[nodeId].nVertexId]) ||
                                    (MeshVertexScalarValue[(*pVecReebNode)[highMINNode].nVertexId] ==
                                     MeshVertexScalarValue[(*pVecReebNode)[nodeId].nVertexId]
                                     && (*pVecReebNode)[highMINNode].nVertexId <
                                        (*pVecReebNode)[nodeId].nVertexId) //highMINNode < nodeId
                                        ) {
                                    highMINNode = nodeId;
                                    //pairNodeType = MINIMUM_REEB;
                                }
                            }
                            tmpSet.erase(nodeId);
                        }
                    }
                }
            } else {
                tmpSet.insert(highNodeId);
                if ((MeshVertexScalarValue[(*pVecReebNode)[highNodeId].nVertexId] > highHeight) ||
                    (MeshVertexScalarValue[(*pVecReebNode)[highNodeId].nVertexId] == highHeight
                     && (*pVecReebNode)[highNodeId].nVertexId >
                        (*pVecReebNode)[highestNodeId].nVertexId) //highNodeId > highestNodeId
                        ) {
                    highHeight = MeshVertexScalarValue[(*pVecReebNode)[highNodeId].nVertexId];
                    highestNodeId = highNodeId;
                }
            }
        }
        highPath.swap(tmpSet);

        if (highPath.empty())
            break;
    }
    return;
}

int psbmReebGraph::DownwardWalkingForPairing(const int startingNode, REEB_GRAPH_CRITICAL_TYPE &pairNodeType) {
    int pairNode = -1;
    pairNodeType = UNKNOWN_REEB;
    int leftMinPairNode = -1;
    int rightMinPairNode = -1;
    int upPairNode = -1;
    // each path is represented by its lowest node
    std::set<int> leftPath;
    std::set<int> rightPath;

    bool GoDownOnLeft = false;
    //
    int leftNode = TheOtherEndPoint(startingNode, (*pVecReebNode)[startingNode].ptrArcDownId->front()->clusterLabel);
    int rightNode = TheOtherEndPoint(startingNode, (*pVecReebNode)[startingNode].ptrArcDownId->back()->clusterLabel);

    // highest height in left path or right path
    double rightHighestHeight = 0.0f;
    double leftHighestHeight = 0.0f;
    //
    // this pairing has unique loop or not
    bool uniqueLoop = true;
    //
    if (leftNode == rightNode) {//can not be minimum
        pairNodeType = UP_FORKING_REEB;
        return leftNode;
    } else {
        if ((MeshVertexScalarValue[(*pVecReebNode)[leftNode].nVertexId] >
             MeshVertexScalarValue[(*pVecReebNode)[rightNode].nVertexId]) ||
            ((MeshVertexScalarValue[(*pVecReebNode)[leftNode].nVertexId] ==
              MeshVertexScalarValue[(*pVecReebNode)[rightNode].nVertexId])
             && (*pVecReebNode)[leftNode].nVertexId > (*pVecReebNode)[rightNode].nVertexId) //leftNode > rightNode
                ) {// return higher minimum
            if ((*pVecReebNode)[leftNode].crType == MINIMUM_REEB) {
                pairNodeType = MINIMUM_REEB;
                return leftNode;
            } else {

                if ((*pVecReebNode)[rightNode].crType == MINIMUM_REEB) {
                    rightMinPairNode = rightNode;
                    //pairNodeType = MINIMUM_REEB;
                    //return rightNode;
                }
            }
        } else {
            if ((*pVecReebNode)[rightNode].crType == MINIMUM_REEB) {
                pairNodeType = MINIMUM_REEB;
                return rightNode;
            } else {
                if ((*pVecReebNode)[leftNode].crType == MINIMUM_REEB) {
                    leftMinPairNode = leftNode;
                    //pairNodeType = MINIMUM_REEB;
                    //return leftNode;
                }
            }
        }
    }
    //

    leftPath.insert(leftNode);
    rightPath.insert(rightNode);
    leftHighestHeight = MeshVertexScalarValue[(*pVecReebNode)[leftNode].nVertexId];
    rightHighestHeight = MeshVertexScalarValue[(*pVecReebNode)[rightNode].nVertexId];
    //
    if ((MeshVertexScalarValue[(*pVecReebNode)[leftNode].nVertexId] >
         MeshVertexScalarValue[(*pVecReebNode)[rightNode].nVertexId]) ||
        ((MeshVertexScalarValue[(*pVecReebNode)[leftNode].nVertexId] ==
          MeshVertexScalarValue[(*pVecReebNode)[rightNode].nVertexId])
         && (*pVecReebNode)[leftNode].nVertexId > (*pVecReebNode)[rightNode].nVertexId)//leftNode > rightNode
            ) {

        GoDownOnLeft = true;
    } else {
        GoDownOnLeft = false;
    }

    int leftHighestNodeId = leftNode;
    int rightHighestNodeId = rightNode;
    if (leftMinPairNode > -1)
        leftPath.clear();
    if (rightMinPairNode > -1)
        rightPath.clear();
    while (1) {
        if (GoDownOnLeft) {
            WalkFromHighLevelToLowLevel(leftPath, rightPath,
                                        leftHighestHeight, leftHighestNodeId,
                                        rightHighestHeight, rightHighestNodeId,
                                        pairNode, leftMinPairNode);

            //if (pairNodeType == MINIMUM_REEB)
            //{// check if there is higher minimum in right path
            //	WalkFromHighLevelToLowLevel(rightPath,leftPath,
            //								rightHighestHeight, rightHighestNodeId,
            //								leftHighestHeight, leftHighestNodeId,
            //								pairNode, rightMinPairNode);
            //}
            if (pairNode > 0) {// find the common ancestor
                //if (pairNodeType == MINIMUM_REEB &&
                //	( ( ((*pVecReebNode)[pairNode].value > leftHighestHeight) ||
                //		((*pVecReebNode)[pairNode].value == leftHighestHeight &&
                //		pairNode > leftHighestNodeId) ) &&
                //	   ( ((*pVecReebNode)[pairNode].value > rightHighestHeight) ||
                //		((*pVecReebNode)[pairNode].value == rightHighestHeight &&
                //		pairNode > rightHighestNodeId) )
                //	 )
                //	 )
                //{
                pairNodeType = UP_FORKING_REEB;
                break;
                //}
            }
            if (leftPath.empty())//leftMinPairNode > 0)
            {
                if (((MeshVertexScalarValue[(*pVecReebNode)[leftMinPairNode].nVertexId] > rightHighestHeight) ||
                     (MeshVertexScalarValue[(*pVecReebNode)[leftMinPairNode].nVertexId] == rightHighestHeight &&
                      (*pVecReebNode)[leftMinPairNode].nVertexId >
                      (*pVecReebNode)[rightHighestNodeId].nVertexId)) //leftMinPairNode > rightHighestNodeId
                        ) {
                    break;
                }
            }
            //if (pairNodeType == UP_FORKING_REEB)
            //	break;
            GoDownOnLeft = false;
        } else {
            WalkFromHighLevelToLowLevel(rightPath, leftPath,
                                        rightHighestHeight, rightHighestNodeId,
                                        leftHighestHeight, leftHighestNodeId,
                                        pairNode, rightMinPairNode);
            //if (pairNodeType == MINIMUM_REEB)
            //{
            //	WalkFromHighLevelToLowLevel(leftPath, rightPath,
            //								leftHighestHeight, leftHighestNodeId,
            //								rightHighestHeight, rightHighestNodeId,
            //								pairNode, leftMinPairNode);
            //}
            //if (pairNode > 0)
            //	break;
            if (pairNode > 0) {
                //if (pairNodeType == MINIMUM_REEB &&
                //	( ( ((*pVecReebNode)[pairNode].value > leftHighestHeight) ||
                //		((*pVecReebNode)[pairNode].value == leftHighestHeight &&
                //		pairNode > leftHighestNodeId) ) &&
                //	   ( ((*pVecReebNode)[pairNode].value > rightHighestHeight) ||
                //		((*pVecReebNode)[pairNode].value == rightHighestHeight &&
                //		pairNode > rightHighestNodeId) )
                //	 )
                //	 )
                //{
                pairNodeType = UP_FORKING_REEB;
                break;
                //}
            }
            if (rightPath.empty())//rightMinPairNode > 0)
            {
                if (((MeshVertexScalarValue[(*pVecReebNode)[rightMinPairNode].nVertexId] > leftHighestHeight) ||
                     (MeshVertexScalarValue[(*pVecReebNode)[rightMinPairNode].nVertexId] == leftHighestHeight &&
                      (*pVecReebNode)[rightMinPairNode].nVertexId >
                      (*pVecReebNode)[leftHighestNodeId].nVertexId)) //rightMinPairNode > leftHighestNodeId
                        ) {
                    break;
                }
            }
            //if (pairNodeType == UP_FORKING_REEB)
            //	break;
            GoDownOnLeft = true;
        }


        if (leftPath.empty() && rightPath.empty())
            break;
    }
    if (pairNode < 0) {
        if (rightMinPairNode > 0 && leftMinPairNode > 0) {
            if (MeshVertexScalarValue[(*pVecReebNode)[rightMinPairNode].nVertexId] >
                MeshVertexScalarValue[(*pVecReebNode)[leftMinPairNode].nVertexId] ||
                (MeshVertexScalarValue[(*pVecReebNode)[rightMinPairNode].nVertexId] ==
                 MeshVertexScalarValue[(*pVecReebNode)[leftMinPairNode].nVertexId]
                 && (*pVecReebNode)[rightMinPairNode].nVertexId >
                    (*pVecReebNode)[leftMinPairNode].nVertexId) //rightMinPairNode > leftMinPairNode
                    ) {
                pairNode = rightMinPairNode;
            } else
                pairNode = leftMinPairNode;
        } else {
            if (rightMinPairNode > 0) {
                pairNode = rightMinPairNode;
            } else
                pairNode = leftMinPairNode;

        }
        pairNodeType = MINIMUM_REEB;
    }

    return pairNode;
}

void psbmReebGraph::ComputingPairing() {
    std::set<int>::iterator sIter;
    REEB_GRAPH_CRITICAL_TYPE pairType;
    string a[] = {" non-loop ", " loop "};
    int pairNode;
    std::pair<int, int> tmpPair;
    //
    criticalPairing = new std::vector<std::pair<int, int> >;
    for (sIter = criticalNodeIdSet->begin(); sIter != criticalNodeIdSet->end(); sIter++) {
        if ((*pVecReebNode)[*sIter].crType == DOWN_FORKING_REEB) {
            pairNode = DownwardWalkingForPairing(*sIter, pairType);
            //std::cout << "( " << pairNode ;
            //std::cout	<< " , " << *sIter << " ) "  << a[pairType] << std::endl;
            if (pairType == UP_FORKING_REEB) {
                tmpPair.first = *sIter;
                tmpPair.second = pairNode;//*sIter;
                criticalPairing->push_back(tmpPair);
            }
        }
    }
    for (unsigned int i = 0; i < criticalPairing->size(); i++) {
        std::cout << "(" << (*criticalPairing)[i].first << ", " <<
                  (*criticalPairing)[i].second << ")" << std::endl;

    }
    return;
}

void psbmReebGraph::ReebGraphPairing(_simpGraph &rbGraph, std::vector<std::pair<int, int> > &outPairing) {
    std::pair<int, int> tmpPair;
    //
    std::set<std::pair<int, double>, myIntDoublePairCompare> sortedCritNodes;
    //
    for (unsigned int i = 0; i < rbGraph.vecNode.size(); i++) {
        sortedCritNodes.insert(std::pair<int, double>(i, MeshVertexScalarValue[rbGraph.vecNode[i].color]));
    }
    //
    ReebGraphPairingAlgorithm::ReebGraphPairing rgPair;
    std::set<std::pair<int, double>, myIntDoublePairCompare>::iterator sIter;
    //
    int downEdgeOrVertex_0;
    int downEdgeOrVertex[2];

    for (sIter = sortedCritNodes.begin(); sIter != sortedCritNodes.end(); sIter++) {
        if (!rbGraph.vecNode[sIter->first].adjList.empty()) {
            switch ((*pVecReebNode)[rbGraph.vecNode[sIter->first].color].crType) {
                case MINIMUM_REEB:
                    //
                    rgPair.AddMinimumReebNode(sIter->first, sIter->second);
                    break;
                case UP_FORKING_REEB:
                    //std::cout << "UP_FORKING_REEB " << crit_order++ << std::endl;
                    for (std::set<int>::iterator verIter = rbGraph.vecNode[sIter->first].adjList.begin();
                         verIter != rbGraph.vecNode[sIter->first].adjList.end(); verIter++) {
                        if (MeshVertexScalarValue[rbGraph.vecNode[*verIter].color] <
                            MeshVertexScalarValue[rbGraph.vecNode[sIter->first].color] ||
                            (MeshVertexScalarValue[rbGraph.vecNode[*verIter].color] ==
                             MeshVertexScalarValue[rbGraph.vecNode[sIter->first].color] &&
                             rbGraph.vecNode[*verIter].color < rbGraph.vecNode[sIter->first].color)) {
                            downEdgeOrVertex_0 = *verIter;
                            break;
                        }
                    }
                    //
                    rgPair.AddUpForkingReebNode(sIter->first, sIter->second, downEdgeOrVertex_0);
                    break;
                case DOWN_FORKING_REEB:
                    //std::cout << "DOWN_FORKING_REEB " << crit_order++ << std::endl;
                    downEdgeOrVertex_0 = 0;
                    for (std::set<int>::iterator verIter = rbGraph.vecNode[sIter->first].adjList.begin();
                         verIter != rbGraph.vecNode[sIter->first].adjList.end(); verIter++) {
                        if (MeshVertexScalarValue[rbGraph.vecNode[*verIter].color] <
                            MeshVertexScalarValue[rbGraph.vecNode[sIter->first].color] ||
                            (MeshVertexScalarValue[rbGraph.vecNode[*verIter].color] ==
                             MeshVertexScalarValue[rbGraph.vecNode[sIter->first].color] &&
                             rbGraph.vecNode[*verIter].color < rbGraph.vecNode[sIter->first].color)) {
                            downEdgeOrVertex[downEdgeOrVertex_0] = *verIter;
                            downEdgeOrVertex_0++;
                        }
                    }
                    //
                    rgPair.AddDownForkingReebNode(sIter->first, sIter->second, downEdgeOrVertex[0],
                                                  downEdgeOrVertex[1]);
                    break;
                case MAXIMUM_REEB:
                    //std::cout << "MAXIMUM_REEB " << crit_order++ << std::endl;
                    downEdgeOrVertex_0 = *(rbGraph.vecNode[sIter->first].adjList.begin());
                    //
                    //
                    rgPair.AddMaximumReebNode(sIter->first, sIter->second, downEdgeOrVertex_0);
                    break;
            }
        }
    }
    //
    //std::cout << "components " << rgPair.forest.size() << std::endl;
    //
    std::set<std::pair<std::pair<int, int>, double>, myIntDoubleIntPairCompare> sortedUpforking;
    for (unsigned int i = 0; i < rgPair.resPairing.size(); i++) {
        //std::cout<< "(" << rgPair.resPairing[i].first << ", " <<
        //	rgPair.resPairing[i].second << ")";
        if (((*pVecReebNode)[rbGraph.vecNode[rgPair.resPairing[i].first].color].crType == DOWN_FORKING_REEB ||
             (*pVecReebNode)[rbGraph.vecNode[rgPair.resPairing[i].first].color].crType == UP_FORKING_REEB)
            &&
            ((*pVecReebNode)[rbGraph.vecNode[rgPair.resPairing[i].second].color].crType == DOWN_FORKING_REEB ||
             (*pVecReebNode)[rbGraph.vecNode[rgPair.resPairing[i].second].color].crType == UP_FORKING_REEB)
                ) {
            //std::cout<< "(" << rbGraph.vecNode[rgPair.resPairing[i].first].color << ", " <<
            //	rbGraph.vecNode[rgPair.resPairing[i].second].color << ")";
            std::pair<int, int> tmpPair(rbGraph.vecNode[rgPair.resPairing[i].first].color,
                                        rbGraph.vecNode[rgPair.resPairing[i].second].color);
            sortedUpforking.insert(
                    std::pair<std::pair<int, int>, double>(tmpPair, MeshVertexScalarValue[tmpPair.second]));
            //std::cout << " loop " << std::endl;
        }

    }
    //std::cout << "loop paring size : " << sortedUpforking.size() << std::endl;
    //
    criticalPairing = new std::vector<std::pair<int, int> >;
    for (std::set<std::pair<std::pair<int, int>, double>, myIntDoubleIntPairCompare>::iterator sIter = sortedUpforking.begin();
         sIter != sortedUpforking.end(); sIter++) {
        criticalPairing->push_back(sIter->first);
    }


    return;
}

void psbmReebGraph::ComputingPairing(const int __polymorphism__) {
    std::set<int>::iterator vecIter;
    std::pair<int, int> tmpPair;
    //
    std::set<std::pair<int, double>, myIntDoublePairCompare> sortedCritNodes;
    //
    for (vecIter = criticalNodeIdSet->begin(); vecIter != criticalNodeIdSet->end(); vecIter++) {
        sortedCritNodes.insert(
                std::pair<int, double>(*vecIter, MeshVertexScalarValue[(*pVecReebNode)[*vecIter].nVertexId]));
    }
    //
    ReebGraphPairingAlgorithm::ReebGraphPairing rgPair;
    std::set<std::pair<int, double>, myIntDoublePairCompare>::iterator sIter;
    int downEdgeOrVertex_0, downEdgeOrVertex_1;

    for (sIter = sortedCritNodes.begin(); sIter != sortedCritNodes.end(); sIter++) {
        switch ((*pVecReebNode)[sIter->first].crType) {
            case MINIMUM_REEB:
                //std::cout << "MINIMUM_REEB " << crit_order++ << std::endl;
                //ofile << "rgPair.AddMinimumReebNode(" <<  cirt_map[sIter->first] << "," <<
                //	setprecision(32) << cirt_map[sIter->first] << ");\n";
                //
                rgPair.AddMinimumReebNode(sIter->first, sIter->second);
                break;
            case UP_FORKING_REEB:
                //std::cout << "UP_FORKING_REEB " << crit_order++ << std::endl;
                downEdgeOrVertex_0 = (*pVecReebNode)[sIter->first].ptrArcDownId->front()->clusterLabel - 1;
                downEdgeOrVertex_0 = (*pVecSimplifiedArc)[downEdgeOrVertex_0].nCriticalNode0;
                //
                //ofile << "rgPair.AddUpForkingReebNode(" <<  cirt_map[sIter->first] << "," <<
                //	setprecision(32) << cirt_map[sIter->first] << ","
                //	<< cirt_map[downEdgeOrVertex_0] << ");\n";
                //
                rgPair.AddUpForkingReebNode(sIter->first, sIter->second, downEdgeOrVertex_0);
                break;
            case DOWN_FORKING_REEB:
                //std::cout << "DOWN_FORKING_REEB " << crit_order++ << std::endl;
                downEdgeOrVertex_0 = (*pVecReebNode)[sIter->first].ptrArcDownId->front()->clusterLabel - 1;
                downEdgeOrVertex_0 = (*pVecSimplifiedArc)[downEdgeOrVertex_0].nCriticalNode0;
                downEdgeOrVertex_1 = (*pVecReebNode)[sIter->first].ptrArcDownId->back()->clusterLabel - 1;
                downEdgeOrVertex_1 = (*pVecSimplifiedArc)[downEdgeOrVertex_1].nCriticalNode0;
                //
                //ofile << "rgPair.AddDownForkingReebNode(" <<  cirt_map[sIter->first] << "," <<
                //	setprecision(32) << cirt_map[sIter->first] << ","
                //	<< cirt_map[downEdgeOrVertex_0] << "," << cirt_map[downEdgeOrVertex_1] << ");\n";
                //
                rgPair.AddDownForkingReebNode(sIter->first, sIter->second, downEdgeOrVertex_0, downEdgeOrVertex_1);
                break;
            case MAXIMUM_REEB:
                //std::cout << "MAXIMUM_REEB " << crit_order++ << std::endl;
                downEdgeOrVertex_0 = (*pVecReebNode)[sIter->first].ptrArcDownId->front()->clusterLabel - 1;
                downEdgeOrVertex_0 = (*pVecSimplifiedArc)[downEdgeOrVertex_0].nCriticalNode0;
                //
                //ofile << "rgPair.AddMaximumReebNode(" <<  cirt_map[sIter->first] << "," <<
                //	setprecision(32) << cirt_map[sIter->first] << ","
                //	<< cirt_map[downEdgeOrVertex_0] << ");\n";
                //
                rgPair.AddMaximumReebNode(sIter->first, sIter->second, downEdgeOrVertex_0);
                break;
        }
    }
    //
    std::set<std::pair<std::pair<int, int>, double>, myIntDoubleIntPairCompare> sortedUpforking;
    for (unsigned int i = 0; i < rgPair.resPairing.size(); i++) {
        //std::cout<< "(" << rgPair.resPairing[i].first << ", " <<
        //	rgPair.resPairing[i].second << ")";
        if (((*pVecReebNode)[rgPair.resPairing[i].first].crType == DOWN_FORKING_REEB ||
             (*pVecReebNode)[rgPair.resPairing[i].first].crType == UP_FORKING_REEB)
            &&
            ((*pVecReebNode)[rgPair.resPairing[i].second].crType == DOWN_FORKING_REEB ||
             (*pVecReebNode)[rgPair.resPairing[i].second].crType == UP_FORKING_REEB)
                ) {
            //std::cout<< "(" << rgPair.resPairing[i].first << ", " << rgPair.resPairing[i].second << ")";
            //trueCriticalPairing.push_back(rgPair.resPairing[i]);//std::pair<int, int>(rgPair.resPairing[i].second,rgPair.resPairing[i].first));
            sortedUpforking.insert(std::pair<std::pair<int, int>, double>(rgPair.resPairing[i],
                                                                          MeshVertexScalarValue[rgPair.resPairing[i].second]));
            //std::cout << " loop " << std::endl;
        }
        //else
        //{
        //	std::cout << " non-loop " << std::endl;
        //}

    }
    //std::cout << "loop paring size : " << trueCriticalPairing.size() << std::endl;
    //std::cout << "loop paring size : " << sortedUpforking.size() << std::endl;
    //
    // use the true pairing
    criticalPairing = new std::vector<std::pair<int, int> >;
    for (std::set<std::pair<std::pair<int, int>, double>, myIntDoubleIntPairCompare>::iterator sIter = sortedUpforking.begin();
         sIter != sortedUpforking.end(); sIter++) {
        criticalPairing->push_back(sIter->first);
    }
    //*criticalPairing = trueCriticalPairing;
    return;
}

void psbmReebGraph::CreateReebGraph(int nEdges) {
    for (int i = 0; i < 12; i++)
        CreateNode(i);

    CreateArc(0, 3, 0);
    CreateArc(1, 3, 1);
    CreateArc(2, 4, 2);
    CreateArc(3, 4, 3);
    CreateArc(4, 5, 4);
    CreateArc(5, 6, 5);
    CreateArc(5, 9, 6);
    CreateArc(6, 7, 7);
    CreateArc(6, 10, 8);
    CreateArc(7, 8, 9);
    CreateArc(7, 8, 10);
    CreateArc(8, 9, 11);
    CreateArc(9, 10, 12);
    CreateArc(10, 11, 13);
///////////////////////////////////////////
    //for (int i = 0; i < 14; i++)
    //	CreateNode(i, (float)i);
    //CreateArc(0, 1, 0);
    //CreateArc(1, 5, 1);
    //CreateArc(1, 2, 2);
    //CreateArc(2, 3, 3);
    //CreateArc(2, 7, 4);
    //CreateArc(3, 4, 5);
    //CreateArc(4, 5, 6);
    //CreateArc(4, 6, 7);
    //CreateArc(3, 7, 8);
    //CreateArc(5, 6, 9);
    //CreateArc(6, 8, 10);
    //CreateArc(7, 9, 11);
    //CreateArc(8, 10, 12);
    //CreateArc(8, 10, 13);
    //CreateArc(9, 11, 14);
    //CreateArc(9, 11, 15);
    //CreateArc(11, 12, 16);
    //CreateArc(10, 12, 17);
    //CreateArc(12, 13, 18);
///////////////////////////////////////////
    //for (int i = 0; i < 14; i++)
    //	CreateNode(i, (float)i);

    //CreateArc(0, 2, 0);
    //CreateArc(2, 3, 1);
    //CreateArc(1, 3, 12);
    //CreateArc(2, 8, 2);
    //CreateArc(3, 4, 3);
    //CreateArc(4, 5, 4);
    //CreateArc(4, 7, 5);
    //CreateArc(5, 6, 6);
    //CreateArc(5, 8, 7);
    //CreateArc(6, 10, 8);
    //CreateArc(6, 7, 9);
    //CreateArc(7, 9, 10);
    //CreateArc(8, 11, 11);
    //CreateArc(9, 10, 13);
    //CreateArc(9, 11, 14);
    //CreateArc(10, 12, 15);
    //CreateArc(11, 13, 16);
///////////////////////////////////////////
    //for (int i = 0; i < 10; i++)
    //	CreateNode(i, (float)i);

    //CreateArc(0, 1, 0);
    //CreateArc(1, 2, 1);
    //CreateArc(1, 6, 2);
    //CreateArc(2, 3, 3);
    //CreateArc(2, 8, 4);
    //CreateArc(3, 4, 5);
    //CreateArc(3, 5, 6);
    //CreateArc(4, 6, 7);
    //CreateArc(4, 5, 8);
    //CreateArc(5, 7, 9);
    //CreateArc(6, 8, 10);
    //CreateArc(8, 9, 11);

    std::vector<class psbmReebNode>::iterator vIter;
    criticalNodeIdSet = new std::set<int>;
    for (vIter = pVecReebNode->begin(); vIter != pVecReebNode->end(); vIter++) {
        criticalNodeIdSet->insert(vIter->nVertexId);

        if ((*vIter).ptrArcDownId->size() + (*vIter).ptrArcUpId->size() == 3) {
            if ((*vIter).ptrArcDownId->size() == 2)
                vIter->crType = DOWN_FORKING_REEB;
            else
                vIter->crType = UP_FORKING_REEB;
        } else {
            if ((*vIter).ptrArcDownId->size() + (*vIter).ptrArcUpId->size() == 1) {
                if ((*vIter).ptrArcDownId->size() == 1)
                    vIter->crType = MAXIMUM_REEB;
                else
                    vIter->crType = MINIMUM_REEB;
            } else
                std::cout << "ERROR, NONCRITICAL point" << std::endl;
        }
    }
    //
    pVecSimplifiedArc = new std::vector<struct SimplifiedReebGraphArc>;
    SimplifiedReebGraphArc tmpSimpArc;
    std::set<class psbmReebArc *>::iterator lIter;
    int edgeLabelCounter = 0;
    for (lIter = pListReebArc->begin(); lIter != pListReebArc->end(); lIter++) {
        tmpSimpArc.edgeLabel = edgeLabelCounter;
        (*lIter)->clusterLabel = edgeLabelCounter++;
        tmpSimpArc.nCriticalNode0 = (*lIter)->nNodeId0;
        tmpSimpArc.nCriticalNode1 = (*lIter)->nNodeId1;
        pVecSimplifiedArc->push_back(tmpSimpArc);
    }
    //
    REEB_GRAPH_CRITICAL_TYPE pairType;
    string a[] = {" non-loop ", " loop "};
    for (vIter = pVecReebNode->begin(); vIter != pVecReebNode->end(); vIter++) {
        if ((*vIter).crType == DOWN_FORKING_REEB) {

            std::cout << "( " << DownwardWalkingForPairing(vIter->nVertexId, pairType);
            std::cout << " , " << vIter->nVertexId << " ) " << a[pairType] << std::endl;
        }
    }
    std::cout << std::endl;
    ComputingPairing();

    //
    //ComputeVerticalLoops();
}

void psbmReebGraph::CutLevelSetAlongOneEdgeOnMesh_pair(_Polygon &retPlg,
                                                       std::vector<std::pair<int, int> > &retPair,
                                                       const int edgeIndex,
                                                       const double levelSetValue,
                                                       const int levelSetIndex,
                                                       const int targetArcIdx,
                                                       const bool special_case) {
    // preq : levelSetValue and levelSetIndex are set up

    const double highValue = MeshVertexScalarValue[meshData->vecEdge[edgeIndex].v1];
    //const int   highNodeId = meshData->vecEdge[edgeIndex].v1;
    const double lowValue = MeshVertexScalarValue[meshData->vecEdge[edgeIndex].v0];

    bool equalEndValues = true;
    if (lowValue != highValue)
        equalEndValues = false;
    //
    const int lowNodeIdInMesh = meshData->vecEdge[edgeIndex].v0;
    const int highNodeIdInMesh = meshData->vecEdge[edgeIndex].v1;
    //
    Vector3 tmpPoint0, tmpPoint1;
    //
    int incidentTri = 0;
    int incidentTri_1 = 0;
    int terminationTriId = 0;
    int prevTri = 0;
    int refEdge = edgeIndex;
    //
    std::pair<int, int> point_type;
    //
    if (MeshVertexScalarValue[meshData->vecEdge[refEdge].v0] ==
        MeshVertexScalarValue[meshData->vecEdge[refEdge].v1]) {// take the middle point
        std::cout << "equal edge height" << std::endl;
        tmpPoint0[0] = 0.5 * (meshData->vecVertex[meshData->vecEdge[refEdge].v0].x +
                              meshData->vecVertex[meshData->vecEdge[refEdge].v1].x);
        tmpPoint0[1] = 0.5 * (meshData->vecVertex[meshData->vecEdge[refEdge].v0].y +
                              meshData->vecVertex[meshData->vecEdge[refEdge].v1].y);
        tmpPoint0[2] = 0.5 * (meshData->vecVertex[meshData->vecEdge[refEdge].v0].z +
                              meshData->vecVertex[meshData->vecEdge[refEdge].v1].z);
        // save the point type
        point_type.first = refEdge;
        point_type.second = 1;
    } else {// compute the exact point with the height value
        tmpPoint0.set(meshData->vecVertex[meshData->vecEdge[refEdge].v0].x,
                      meshData->vecVertex[meshData->vecEdge[refEdge].v0].y,
                      meshData->vecVertex[meshData->vecEdge[refEdge].v0].z);
        tmpPoint1.set(meshData->vecVertex[meshData->vecEdge[refEdge].v1].x,
                      meshData->vecVertex[meshData->vecEdge[refEdge].v1].y,
                      meshData->vecVertex[meshData->vecEdge[refEdge].v1].z);
        tmpPoint0 = tmpPoint0 + (levelSetValue - MeshVertexScalarValue[meshData->vecEdge[refEdge].v0]) /
                                (MeshVertexScalarValue[meshData->vecEdge[refEdge].v1] -
                                 MeshVertexScalarValue[meshData->vecEdge[refEdge].v0]) * (tmpPoint1 - tmpPoint0);
        // save the point type
        if (levelSetValue == MeshVertexScalarValue[meshData->vecEdge[refEdge].v0] ||
            levelSetValue == MeshVertexScalarValue[meshData->vecEdge[refEdge].v1]) {
            if (levelSetValue == MeshVertexScalarValue[meshData->vecEdge[refEdge].v0])
                point_type.first = meshData->vecEdge[refEdge].v0;
            else
                point_type.first = meshData->vecEdge[refEdge].v1;
            point_type.second = 0;
        } else {
            point_type.first = refEdge;
            point_type.second = 1;
        }
    }
    // push the first point into the polygon
    retPlg.vecPoints.push_back(tmpPoint0);
    // store the point type
    retPair.push_back(point_type);

    incidentTri = meshData->vecEdge[refEdge].AdjTri[0];
    incidentTri_1 = meshData->vecEdge[refEdge].AdjTri[1];
    terminationTriId = incidentTri;
    prevTri = incidentTri;
    //std::set<int> processedTriangles;
    //std::set<int>::iterator sTriIter;
    boost::unordered_set<int> processedTriangles;
    boost::unordered_set<int>::iterator sTriIter;
    if (!equalEndValues) {
        do {
            processedTriangles.insert(incidentTri);
            int potEdge[2] = {0};
            int potNextTri[2] = {0};
            // search all the other end point across which the level set comes;
            if (meshData->vecTriangle[incidentTri].e01 == refEdge) {
                // compute next cross point at the other two edges
                potEdge[0] = meshData->vecTriangle[incidentTri].e12;
                potEdge[1] = meshData->vecTriangle[incidentTri].e02;
            } else {
                potEdge[0] = meshData->vecTriangle[incidentTri].e01;
                if (meshData->vecTriangle[incidentTri].e02 == refEdge)
                    potEdge[1] = meshData->vecTriangle[incidentTri].e12;
                else
                    potEdge[1] = meshData->vecTriangle[incidentTri].e02;
            }
            // compute the possible next triangle along each edge
            for (int i = 0; i < 2; i++) {
                if (meshData->vecEdge[potEdge[i]].AdjTri[0] == incidentTri)
                    potNextTri[i] = meshData->vecEdge[potEdge[i]].AdjTri[1];
                else
                    potNextTri[i] = meshData->vecEdge[potEdge[i]].AdjTri[0];
            }
            //if (incidentTri != terminationTriId)
            //{
            //	if (potNextTri[0] == terminationTriId ||
            //		potNextTri[1] == terminationTriId)
            //	{
            //		break;
            //	}
            //}
            //pick up the cross edge
            int i = 0;
            int ev0, ev1;
            for (i = 0; i < 2; i++) {
                ev0 = meshData->vecEdge[potEdge[i]].v0;
                ev1 = meshData->vecEdge[potEdge[i]].v1;
                if ((MeshVertexScalarValue[ev0] < levelSetValue &&
                     MeshVertexScalarValue[ev1] > levelSetValue)
                    ||
                    (MeshVertexScalarValue[ev1] < levelSetValue &&
                     MeshVertexScalarValue[ev0] > levelSetValue)
                        ) {
                    sTriIter = processedTriangles.find(potNextTri[i]);
                    if (potNextTri[i] == terminationTriId || incidentTri == terminationTriId ||
                        sTriIter == processedTriangles.end()) {

                        //
                        if (MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0] ==
                            MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v1]) {
                            std::cout << "equal edge height" << std::endl;
                            tmpPoint0[0] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].x +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].x);
                            tmpPoint0[1] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].y +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].y);
                            tmpPoint0[2] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].z +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].z);
                            // save the point type
                            point_type.first = potEdge[i];//refEdge;//
                            point_type.second = 1;
                            //
                        } else {
                            tmpPoint0.set(meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].x,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].y,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].z);
                            tmpPoint1.set(meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].x,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].y,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].z);
                            tmpPoint0 = tmpPoint0 +
                                        (levelSetValue - MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0]) /
                                        (MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v1] -
                                         MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0]) *
                                        (tmpPoint1 - tmpPoint0);
                            // save the point type
                            if (levelSetValue == MeshVertexScalarValue[meshData->vecEdge[refEdge].v0] ||
                                levelSetValue == MeshVertexScalarValue[meshData->vecEdge[refEdge].v1]) {
                                if (levelSetValue == MeshVertexScalarValue[meshData->vecEdge[refEdge].v0])
                                    point_type.first = meshData->vecEdge[refEdge].v0;
                                else
                                    point_type.first = meshData->vecEdge[refEdge].v1;
                                point_type.second = 0;
                            } else {
                                point_type.first = potEdge[i];//refEdge;//
                                point_type.second = 1;
                            }
                            //
                        }
                        refEdge = potEdge[i];
                        break;
                    }
                }
            }
            if (!map_levelset_on_edge_to_this_reeb_arc(targetArcIdx, meshData->vecEdge[refEdge].v0,
                                                       meshData->vecEdge[refEdge].v1, refEdge, special_case)) {
                std::cout << "not level set edge" << std::endl;
            }
            if (meshData->vecEdge[refEdge].AdjTri[0] == incidentTri)
                incidentTri = meshData->vecEdge[refEdge].AdjTri[1];
            else
                incidentTri = meshData->vecEdge[refEdge].AdjTri[0];
            if (incidentTri == terminationTriId)
                break;
            // compute the intersection point
            retPlg.vecPoints.push_back(tmpPoint0);
            // store the point type
            retPair.push_back(point_type);

        } while (incidentTri != terminationTriId);
    } else {
        std::cout << "equal levelset height" << std::endl;
        do {
            processedTriangles.insert(incidentTri);
            int potEdge[2] = {0};
            int potNextTri[2] = {0};
            // search all the other end point across which the level set comes;
            if (meshData->vecTriangle[incidentTri].e01 == refEdge) {
                // compute next cross point at the other two edges
                potEdge[0] = meshData->vecTriangle[incidentTri].e12;
                potEdge[1] = meshData->vecTriangle[incidentTri].e02;
            } else {
                potEdge[0] = meshData->vecTriangle[incidentTri].e01;
                if (meshData->vecTriangle[incidentTri].e02 == refEdge)
                    potEdge[1] = meshData->vecTriangle[incidentTri].e12;
                else
                    potEdge[1] = meshData->vecTriangle[incidentTri].e02;
            }
            // compute the possible next triangle along each edge
            for (int i = 0; i < 2; i++) {
                if (meshData->vecEdge[potEdge[i]].AdjTri[0] == incidentTri)
                    potNextTri[i] = meshData->vecEdge[potEdge[i]].AdjTri[1];
                else
                    potNextTri[i] = meshData->vecEdge[potEdge[i]].AdjTri[0];
            }
            //pick up the cross edge
            int i = 0;
            int ev0, ev1;
            for (i = 0; i < 2; i++) {
                ev0 = meshData->vecEdge[potEdge[i]].v0;
                ev1 = meshData->vecEdge[potEdge[i]].v1;
                if ((ev0 <= levelSetIndex &&
                     ev1 >= levelSetIndex)
                    ||
                    (ev1 <= levelSetIndex &&
                     ev0 >= levelSetIndex)
                        ) {
                    sTriIter = processedTriangles.find(potNextTri[i]);
                    if (potNextTri[i] == terminationTriId || incidentTri == terminationTriId ||
                        sTriIter == processedTriangles.end()) {

                        //
                        if (MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0] ==
                            MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v1]) {
                            std::cout << "equal edge height" << std::endl;
                            tmpPoint0[0] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].x +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].x);
                            tmpPoint0[1] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].y +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].y);
                            tmpPoint0[2] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].z +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].z);
                            // save the point type
                            point_type.first = potEdge[i];//refEdge;//
                            point_type.second = 1;
                            //
                        } else {
                            tmpPoint0.set(meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].x,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].y,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].z);
                            tmpPoint1.set(meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].x,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].y,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].z);
                            tmpPoint0 = tmpPoint0 +
                                        (levelSetValue - MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0]) /
                                        (MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v1] -
                                         MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0]) *
                                        (tmpPoint1 - tmpPoint0);
                            // save the point type
                            if (levelSetValue == MeshVertexScalarValue[meshData->vecEdge[refEdge].v0] ||
                                levelSetValue == MeshVertexScalarValue[meshData->vecEdge[refEdge].v1]) {
                                if (levelSetValue == MeshVertexScalarValue[meshData->vecEdge[refEdge].v0])
                                    point_type.first = meshData->vecEdge[refEdge].v0;
                                else
                                    point_type.first = meshData->vecEdge[refEdge].v1;
                                point_type.second = 0;
                            } else {
                                point_type.first = potEdge[i]; //refEdge;//
                                point_type.second = 1;
                            }
                            //
                        }
                        refEdge = potEdge[i];
                        break;
                    }
                }
            }

            if (meshData->vecEdge[refEdge].AdjTri[0] == incidentTri)
                incidentTri = meshData->vecEdge[refEdge].AdjTri[1];
            else
                incidentTri = meshData->vecEdge[refEdge].AdjTri[0];
            if (incidentTri == terminationTriId)
                break;
            // compute the intersection point
            retPlg.vecPoints.push_back(tmpPoint0);
            // store the point type
            retPair.push_back(point_type);
            //
        } while (incidentTri != terminationTriId);
    }

    retPlg.vecPoints.push_back(retPlg.vecPoints[0]);
    // store the point type
    retPair.push_back(retPair[0]);
    return;
}

void psbmReebGraph::CutLevelSetAlongOneEdgeOnMesh(_Polygon &retPlg,
                                                  const int edgeIndex,
                                                  const double levelSetValue,
                                                  const int levelSetIndex,
                                                  const int targetArcIdx,
                                                  const bool special_case) {
    // preq : levelSetValue and levelSetIndex are set up

    const double highValue = MeshVertexScalarValue[meshData->vecEdge[edgeIndex].v1];
    //const int   highNodeId = meshData->vecEdge[edgeIndex].v1;
    const double lowValue = MeshVertexScalarValue[meshData->vecEdge[edgeIndex].v0];

    bool equalEndValues = true;
    if (lowValue != highValue)
        equalEndValues = false;
    //const int   lowNodeId = meshData->vecEdge[edgeIndex].v0;
    const int lowNodeIdInMesh = meshData->vecEdge[edgeIndex].v0;
    const int highNodeIdInMesh = meshData->vecEdge[edgeIndex].v1;
    //float levelSetValue = highValue;
    //int   levelSetNodeId = highNodeIdInMesh;

    //
    Vector3 tmpPoint0, tmpPoint1;
    //
    int incidentTri = 0;
    int incidentTri_1 = 0;
    int terminationTriId = 0;
    int prevTri = 0;
    int refEdge = edgeIndex;
    //
    //if (meshData->vecEdge[refEdge].v0 == levelSetNodeId ||
    //	meshData->vecEdge[refEdge].v1 == levelSetNodeId )
    //{// cross the edge
    //	levelSetValue = (highValue + lowValue) / 2.0f;
    //	levelSetNodeId = -1;
    //}
    if (MeshVertexScalarValue[meshData->vecEdge[refEdge].v0] ==
        MeshVertexScalarValue[meshData->vecEdge[refEdge].v1]) {// take the middle point
        std::cout << "equal edge height" << std::endl;
        tmpPoint0[0] = 0.5f * (meshData->vecVertex[meshData->vecEdge[refEdge].v0].x +
                               meshData->vecVertex[meshData->vecEdge[refEdge].v1].x);
        tmpPoint0[1] = 0.5f * (meshData->vecVertex[meshData->vecEdge[refEdge].v0].y +
                               meshData->vecVertex[meshData->vecEdge[refEdge].v1].y);
        tmpPoint0[2] = 0.5f * (meshData->vecVertex[meshData->vecEdge[refEdge].v0].z +
                               meshData->vecVertex[meshData->vecEdge[refEdge].v1].z);
    } else {// compute the exact point with the height value
        tmpPoint0.set(meshData->vecVertex[meshData->vecEdge[refEdge].v0].x,
                      meshData->vecVertex[meshData->vecEdge[refEdge].v0].y,
                      meshData->vecVertex[meshData->vecEdge[refEdge].v0].z);
        tmpPoint1.set(meshData->vecVertex[meshData->vecEdge[refEdge].v1].x,
                      meshData->vecVertex[meshData->vecEdge[refEdge].v1].y,
                      meshData->vecVertex[meshData->vecEdge[refEdge].v1].z);
        tmpPoint0 = tmpPoint0 + (levelSetValue - MeshVertexScalarValue[meshData->vecEdge[refEdge].v0]) /
                                (MeshVertexScalarValue[meshData->vecEdge[refEdge].v1] -
                                 MeshVertexScalarValue[meshData->vecEdge[refEdge].v0]) * (tmpPoint1 - tmpPoint0);
    }
    // push the first point into the polygon
    retPlg.vecPoints.push_back(tmpPoint0);

    incidentTri = meshData->vecEdge[refEdge].AdjTri[0];
    incidentTri_1 = meshData->vecEdge[refEdge].AdjTri[1];
    terminationTriId = incidentTri;
    prevTri = incidentTri;
    //std::set<int> processedTriangles;
    //std::set<int>::iterator sTriIter;
    boost::unordered_set<int> processedTriangles;
    boost::unordered_set<int>::iterator sTriIter;
    if (!equalEndValues) {
        do {
            processedTriangles.insert(incidentTri);
            int potEdge[2] = {0};
            int potNextTri[2] = {0};
            // search all the other end point across which the level set comes;
            if (meshData->vecTriangle[incidentTri].e01 == refEdge) {
                // compute next cross point at the other two edges
                potEdge[0] = meshData->vecTriangle[incidentTri].e12;
                potEdge[1] = meshData->vecTriangle[incidentTri].e02;
            } else {
                potEdge[0] = meshData->vecTriangle[incidentTri].e01;
                if (meshData->vecTriangle[incidentTri].e02 == refEdge)
                    potEdge[1] = meshData->vecTriangle[incidentTri].e12;
                else
                    potEdge[1] = meshData->vecTriangle[incidentTri].e02;
            }
            // compute the possible next triangle along each edge
            for (int i = 0; i < 2; i++) {
                if (meshData->vecEdge[potEdge[i]].AdjTri[0] == incidentTri)
                    potNextTri[i] = meshData->vecEdge[potEdge[i]].AdjTri[1];
                else
                    potNextTri[i] = meshData->vecEdge[potEdge[i]].AdjTri[0];
            }
            //if (incidentTri != terminationTriId)
            //{
            //	if (potNextTri[0] == terminationTriId ||
            //		potNextTri[1] == terminationTriId)
            //	{
            //		break;
            //	}
            //}
            //pick up the cross edge
            int i = 0;
            int ev0, ev1;
            for (i = 0; i < 2; i++) {
                ev0 = meshData->vecEdge[potEdge[i]].v0;
                ev1 = meshData->vecEdge[potEdge[i]].v1;
                if ((MeshVertexScalarValue[ev0] < levelSetValue &&
                     MeshVertexScalarValue[ev1] > levelSetValue)
                    ||
                    (MeshVertexScalarValue[ev1] < levelSetValue &&
                     MeshVertexScalarValue[ev0] > levelSetValue)
                        ) {
                    sTriIter = processedTriangles.find(potNextTri[i]);
                    if (potNextTri[i] == terminationTriId || incidentTri == terminationTriId ||
                        sTriIter == processedTriangles.end()) {
                        if (MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0] ==
                            MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v1]) {
                            std::cout << "equal edge height" << std::endl;
                            tmpPoint0[0] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].x +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].x);
                            tmpPoint0[1] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].y +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].y);
                            tmpPoint0[2] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].z +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].z);
                        } else {
                            tmpPoint0.set(meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].x,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].y,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].z);
                            tmpPoint1.set(meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].x,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].y,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].z);
                            tmpPoint0 = tmpPoint0 +
                                        (levelSetValue - MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0]) /
                                        (MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v1] -
                                         MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0]) *
                                        (tmpPoint1 - tmpPoint0);
                        }
                        refEdge = potEdge[i];
                        break;
                    }
                }
            }
            if (!map_levelset_on_edge_to_this_reeb_arc(targetArcIdx, meshData->vecEdge[refEdge].v0,
                                                       meshData->vecEdge[refEdge].v1, refEdge, special_case)) {
                std::cout << "not level set edge" << std::endl;
            }
            if (meshData->vecEdge[refEdge].AdjTri[0] == incidentTri)
                incidentTri = meshData->vecEdge[refEdge].AdjTri[1];
            else
                incidentTri = meshData->vecEdge[refEdge].AdjTri[0];
            if (incidentTri == terminationTriId)
                break;
            // compute the intersection point
            retPlg.vecPoints.push_back(tmpPoint0);


        } while (incidentTri != terminationTriId);
    } else {
        std::cout << "equal level height" << std::endl;
        do {
            processedTriangles.insert(incidentTri);
            int potEdge[2] = {0};
            int potNextTri[2] = {0};
            // search all the other end point across which the level set comes;
            if (meshData->vecTriangle[incidentTri].e01 == refEdge) {
                // compute next cross point at the other two edges
                potEdge[0] = meshData->vecTriangle[incidentTri].e12;
                potEdge[1] = meshData->vecTriangle[incidentTri].e02;
            } else {
                potEdge[0] = meshData->vecTriangle[incidentTri].e01;
                if (meshData->vecTriangle[incidentTri].e02 == refEdge)
                    potEdge[1] = meshData->vecTriangle[incidentTri].e12;
                else
                    potEdge[1] = meshData->vecTriangle[incidentTri].e02;
            }
            // compute the possible next triangle along each edge
            for (int i = 0; i < 2; i++) {
                if (meshData->vecEdge[potEdge[i]].AdjTri[0] == incidentTri)
                    potNextTri[i] = meshData->vecEdge[potEdge[i]].AdjTri[1];
                else
                    potNextTri[i] = meshData->vecEdge[potEdge[i]].AdjTri[0];
            }
            //pick up the cross edge
            int i = 0;
            int ev0, ev1;
            for (i = 0; i < 2; i++) {
                ev0 = meshData->vecEdge[potEdge[i]].v0;
                ev1 = meshData->vecEdge[potEdge[i]].v1;
                if ((ev0 <= levelSetIndex &&
                     ev1 >= levelSetIndex)
                    ||
                    (ev1 <= levelSetIndex &&
                     ev0 >= levelSetIndex)
                        ) {
                    sTriIter = processedTriangles.find(potNextTri[i]);
                    if (potNextTri[i] == terminationTriId || incidentTri == terminationTriId ||
                        sTriIter == processedTriangles.end()) {
                        if (MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0] ==
                            MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v1]) {
                            std::cout << "equal edge height" << std::endl;
                            tmpPoint0[0] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].x +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].x);
                            tmpPoint0[1] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].y +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].y);
                            tmpPoint0[2] = 0.5f * (meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].z +
                                                   meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].z);
                        } else {
                            tmpPoint0.set(meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].x,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].y,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v0].z);
                            tmpPoint1.set(meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].x,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].y,
                                          meshData->vecVertex[meshData->vecEdge[potEdge[i]].v1].z);
                            tmpPoint0 = tmpPoint0 +
                                        (levelSetValue - MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0]) /
                                        (MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v1] -
                                         MeshVertexScalarValue[meshData->vecEdge[potEdge[i]].v0]) *
                                        (tmpPoint1 - tmpPoint0);
                        }
                        refEdge = potEdge[i];
                        break;
                    }
                }
            }

            if (meshData->vecEdge[refEdge].AdjTri[0] == incidentTri)
                incidentTri = meshData->vecEdge[refEdge].AdjTri[1];
            else
                incidentTri = meshData->vecEdge[refEdge].AdjTri[0];
            if (incidentTri == terminationTriId)
                break;
            // compute the intersection point
            retPlg.vecPoints.push_back(tmpPoint0);


        } while (incidentTri != terminationTriId);
    }

    retPlg.vecPoints.push_back(retPlg.vecPoints[0]);
    return;
}

void psbmReebGraph::CutLevelSet_specific_value_pair(_Polygon &retPly,
                                                    std::vector<std::pair<int, int> > &retPtType,
                                                    SimplifiedReebGraphArc &inArc,
                                                    const double levelSetValue) {
    bool bIsEdge = false;
    int nArcAsEdge = -1;
    for (std::vector<int>::iterator vIter = meshData->vecVertex[(*pVecReebNode)[inArc.nCriticalNode0].nVertexId].adjEdges.begin();
         vIter != meshData->vecVertex[(*pVecReebNode)[inArc.nCriticalNode0].nVertexId].adjEdges.end();
         vIter++) {
        int nTheOtherEndPoint = meshData->vecEdge[*vIter].v0 + meshData->vecEdge[*vIter].v1 -
                                (*pVecReebNode)[inArc.nCriticalNode0].nVertexId;
        if (nTheOtherEndPoint == (*pVecReebNode)[inArc.nCriticalNode1].nVertexId) {
            bIsEdge = true;
            nArcAsEdge = *vIter;
            break;
        }
    }
    //
    bool special_case_no_between_pts = false;
    for (std::list<class psbmReebArc *>::iterator listIter = (*pVecReebNode)[inArc.nCriticalNode0].ptrArcUpId->begin();
         listIter != (*pVecReebNode)[inArc.nCriticalNode0].ptrArcUpId->end();
         listIter++) {
        if ((*listIter)->nNodeId1 == inArc.nCriticalNode1) {
            special_case_no_between_pts = true;
            break;
        }
    }
    //
    bool bProcessed = false;
    if (bIsEdge) {
        const int levelSet_index = (*pVecReebNode)[inArc.nCriticalNode0].nVertexId;
        if (special_case_no_between_pts) {
            bProcessed = true;
            CutLevelSetAlongOneEdgeOnMesh_pair(retPly, retPtType, nArcAsEdge, levelSetValue, levelSet_index,
                                               inArc.edgeLabel, special_case_no_between_pts);
        } else {
            if (map_levelset_on_edge_to_this_reeb_arc(inArc.edgeLabel, inArc.nCriticalNode0, inArc.nCriticalNode1,
                                                      nArcAsEdge, special_case_no_between_pts)) {
                bProcessed = true;
                CutLevelSetAlongOneEdgeOnMesh_pair(retPly, retPtType, nArcAsEdge, levelSetValue, levelSet_index,
                                                   inArc.edgeLabel, special_case_no_between_pts);
            }
        }
    }
    if (!bProcessed) {
        double curVerValue = MeshVertexScalarValue[(*pVecReebNode)[inArc.nCriticalNode0].nVertexId];
        // both of them are < levelSetValue
        psbmReebArc *travArcPtr = NULL;
        //
        for (std::list<psbmReebArc *>::iterator listIter = (*pVecReebNode)[inArc.nCriticalNode0].ptrArcUpId->begin();
             listIter != (*pVecReebNode)[inArc.nCriticalNode0].ptrArcUpId->end(); listIter++) {
            if ((*listIter)->clusterLabel == inArc.edgeLabel) {
                travArcPtr = *listIter;
                break;
            }
        }
        //
        curVerValue = MeshVertexScalarValue[(*pVecReebNode)[travArcPtr->nNodeId1].nVertexId];
        //
        while (curVerValue < levelSetValue) {
            travArcPtr = (*pVecReebNode)[travArcPtr->nNodeId1].ptrArcUpId->front();
            curVerValue = MeshVertexScalarValue[(*pVecReebNode)[travArcPtr->nNodeId1].nVertexId];
            //
        }
        // now travel all the edges intersecting this arc to find an edge
        int cross_edge_id = -1;
        int levelSet_index = -1;
        int nVertexOnMeshWithHigherValue = (*pVecReebNode)[travArcPtr->nNodeId1].nVertexId;
        for (std::vector<int>::iterator vIter = meshData->vecVertex[nVertexOnMeshWithHigherValue].adjEdges.begin();
             vIter != meshData->vecVertex[nVertexOnMeshWithHigherValue].adjEdges.end();
             vIter++) {
            int nTheOtherEndPoint =
                    meshData->vecEdge[*vIter].v0 + meshData->vecEdge[*vIter].v1 - nVertexOnMeshWithHigherValue;

            if (MeshVertexScalarValue[nTheOtherEndPoint] <= levelSetValue) {
                if (vecVertexMappingToSimplifiedReebArc[nTheOtherEndPoint] == inArc.edgeLabel ||
                    vecVertexMappingToSimplifiedReebArc[nVertexOnMeshWithHigherValue] == inArc.edgeLabel) {
                    cross_edge_id = *vIter;
                    break;
                } else {
                    if (map_levelset_on_edge_to_this_reeb_arc(inArc.edgeLabel,
                                                              nVertexOnMeshWithHigherValue,//(*ProcessedVertex)[nVertexOnMeshWithHigherValue],
                                                              nTheOtherEndPoint,//(*ProcessedVertex)[nTheOtherEndPoint],
                                                              *vIter, special_case_no_between_pts)) {
                        cross_edge_id = *vIter;
                        break;
                    }
                }

            }
        }
        levelSet_index = (*pVecReebNode)[travArcPtr->nNodeId1].nVertexId;
        if (cross_edge_id == -1) {
            std::cout << "CAN not the edge cross the levelset value" << std::endl;
            exit(0);
        }
        //
        CutLevelSetAlongOneEdgeOnMesh_pair(retPly, retPtType, cross_edge_id, levelSetValue, levelSet_index,
                                           inArc.edgeLabel, special_case_no_between_pts);
    }
}

void psbmReebGraph::CutLevelSet_specific_value(_Polygon &retPly,
                                               SimplifiedReebGraphArc &inArc,
                                               const double levelSetValue) {
    // Preq: levelSetValue is different lowest value and highest value on this arc
    //double prevVerValue = MeshVertexScalarValue[(*pVecReebNode)[inArc.nCriticalNode0].nVertexId];
    bool bIsEdge = false;
    int nArcAsEdge = -1;
    for (std::vector<int>::iterator vIter = meshData->vecVertex[(*pVecReebNode)[inArc.nCriticalNode0].nVertexId].adjEdges.begin();
         vIter != meshData->vecVertex[(*pVecReebNode)[inArc.nCriticalNode0].nVertexId].adjEdges.end();
         vIter++) {
        int nTheOtherEndPoint = meshData->vecEdge[*vIter].v0 + meshData->vecEdge[*vIter].v1 -
                                (*pVecReebNode)[inArc.nCriticalNode0].nVertexId;
        if (nTheOtherEndPoint == (*pVecReebNode)[inArc.nCriticalNode1].nVertexId) {
            bIsEdge = true;
            nArcAsEdge = *vIter;
            break;
        }
    }
    //
    bool special_case_no_between_pts = false;
    for (std::list<class psbmReebArc *>::iterator listIter = (*pVecReebNode)[inArc.nCriticalNode0].ptrArcUpId->begin();
         listIter != (*pVecReebNode)[inArc.nCriticalNode0].ptrArcUpId->end();
         listIter++) {
        if ((*listIter)->nNodeId1 == inArc.nCriticalNode1) {
            special_case_no_between_pts = true;
            break;
        }
    }
    bool bProcessed = false;
    if (bIsEdge) {
        const int levelSet_index = (*pVecReebNode)[inArc.nCriticalNode0].nVertexId;
        if (special_case_no_between_pts) {
            bProcessed = true;
            CutLevelSetAlongOneEdgeOnMesh(retPly, nArcAsEdge, levelSetValue, levelSet_index, inArc.edgeLabel,
                                          special_case_no_between_pts);
        } else {
            if (map_levelset_on_edge_to_this_reeb_arc(inArc.edgeLabel, inArc.nCriticalNode0, inArc.nCriticalNode1,
                                                      nArcAsEdge, special_case_no_between_pts)) {
                bProcessed = true;
                CutLevelSetAlongOneEdgeOnMesh(retPly, nArcAsEdge, levelSetValue, levelSet_index, inArc.edgeLabel,
                                              special_case_no_between_pts);
            }
        }
    }
    if (!bProcessed) {
        double curVerValue = MeshVertexScalarValue[(*pVecReebNode)[inArc.nCriticalNode0].nVertexId];
        // both of them are < levelSetValue
        psbmReebArc *travArcPtr = NULL;
        //
        for (std::list<psbmReebArc *>::iterator listIter = (*pVecReebNode)[inArc.nCriticalNode0].ptrArcUpId->begin();
             listIter != (*pVecReebNode)[inArc.nCriticalNode0].ptrArcUpId->end(); listIter++) {
            if ((*listIter)->clusterLabel == inArc.edgeLabel) {
                travArcPtr = *listIter;
                break;
            }
        }
        //
        curVerValue = MeshVertexScalarValue[(*pVecReebNode)[travArcPtr->nNodeId1].nVertexId];
        while (curVerValue < levelSetValue) {
            //
            travArcPtr = (*pVecReebNode)[travArcPtr->nNodeId1].ptrArcUpId->front();
            curVerValue = MeshVertexScalarValue[(*pVecReebNode)[travArcPtr->nNodeId1].nVertexId];
        }
        // now travel all the edges intersecting this arc to find an edge
        int cross_edge_id = -1;
        int levelSet_index = -1;
        int nVertexOnMeshWithHigherValue = (*pVecReebNode)[travArcPtr->nNodeId1].nVertexId;
        for (std::vector<int>::iterator vIter = meshData->vecVertex[nVertexOnMeshWithHigherValue].adjEdges.begin();
             vIter != meshData->vecVertex[nVertexOnMeshWithHigherValue].adjEdges.end();
             vIter++) {
            int nTheOtherEndPoint =
                    meshData->vecEdge[*vIter].v0 + meshData->vecEdge[*vIter].v1 - nVertexOnMeshWithHigherValue;

            if (MeshVertexScalarValue[nTheOtherEndPoint] <= levelSetValue) {
                if (vecVertexMappingToSimplifiedReebArc[nTheOtherEndPoint] == inArc.edgeLabel ||
                    vecVertexMappingToSimplifiedReebArc[nVertexOnMeshWithHigherValue] == inArc.edgeLabel) {
                    cross_edge_id = *vIter;
                    break;
                } else {
                    if (map_levelset_on_edge_to_this_reeb_arc(inArc.edgeLabel,
                                                              nVertexOnMeshWithHigherValue,//(*ProcessedVertex)[nVertexOnMeshWithHigherValue],
                                                              nTheOtherEndPoint, //(*ProcessedVertex)[nTheOtherEndPoint],
                                                              *vIter, special_case_no_between_pts)) {
                        cross_edge_id = *vIter;
                        break;
                    }
                }

            }
        }
        levelSet_index = (*pVecReebNode)[travArcPtr->nNodeId1].nVertexId;
        if (cross_edge_id == -1) {
            std::cout << "CAN not the edge cross the levelset value" << std::endl;
            exit(0);
        }
        //
        CutLevelSetAlongOneEdgeOnMesh(retPly, cross_edge_id, levelSetValue, levelSet_index, inArc.edgeLabel,
                                      special_case_no_between_pts);
    }
    //
    return;
}

//
void psbmReebGraph::PushIntoQueueForNonVisitedNodes_UP(const int nodeId, const int vecIndex,
                                                       std::vector<std::pair<int, bool> > &processed_node,
                                                       std::vector<bool> &processed_arc_bit,
                                                       std::queue<int> &ready_to_decide_node) {
    if (!processed_node[vecIndex].second) {
        if ((*pVecReebNode)[nodeId].crType == UP_FORKING_REEB) {
            // put it into the queue
            ready_to_decide_node.push(vecIndex);
        } else {
            if ((*pVecReebNode)[nodeId].crType == DOWN_FORKING_REEB) {//
                int left_arc_id = (*pVecReebNode)[nodeId].ptrArcDownId->front()->clusterLabel - 1;
                int right_arc_id = (*pVecReebNode)[nodeId].ptrArcDownId->back()->clusterLabel - 1;
                if (processed_arc_bit[left_arc_id] && processed_arc_bit[right_arc_id]) {
                    ready_to_decide_node.push(vecIndex);
                }
            }
        }

    }
    return;
}

//
void psbmReebGraph::PushIntoQueueForNonVisitedNodes(const int nodeId, const int vecIndex,
                                                    std::vector<std::pair<int, bool> > &processed_node,
                                                    std::vector<bool> &processed_arc_bit,
                                                    std::queue<int> &ready_to_decide_node) {
    if (!processed_node[vecIndex].second) {
        if ((*pVecReebNode)[nodeId].crType == UP_FORKING_REEB) {
            // put it into the queue
            int up_front_lhs = (*pVecReebNode)[nodeId].ptrArcUpId->front()->clusterLabel;
            int up_back_rhs = (*pVecReebNode)[nodeId].ptrArcUpId->front()->clusterLabel;
            int down_front_lhs = (*pVecReebNode)[nodeId].ptrArcDownId->front()->clusterLabel;
            if (processed_arc_bit[up_front_lhs - 1] && processed_arc_bit[up_back_rhs - 1] ||
                processed_arc_bit[down_front_lhs - 1]
                    ) {
                ready_to_decide_node.push(vecIndex);
            }
        } else {
            if ((*pVecReebNode)[nodeId].crType == DOWN_FORKING_REEB) {//
                int up_front_lhs = (*pVecReebNode)[nodeId].ptrArcUpId->front()->clusterLabel;
                int down_front_lhs = (*pVecReebNode)[nodeId].ptrArcDownId->front()->clusterLabel;
                int down_back_rhs = (*pVecReebNode)[nodeId].ptrArcDownId->back()->clusterLabel;
                if (processed_arc_bit[down_front_lhs - 1] && processed_arc_bit[down_back_rhs - 1] ||
                    processed_arc_bit[up_front_lhs - 1]
                        ) {
                    ready_to_decide_node.push(vecIndex);
                }
            }
        }

    }
    return;
}

void psbmReebGraph::PushIntoQueueForNonVisitedNodes_DOWN(const int nodeId, const int vecIndex,
                                                         std::vector<std::pair<int, bool> > &processed_node,
                                                         std::vector<bool> &processed_arc_bit,
                                                         std::queue<int> &ready_to_decide_node) {
    if (!processed_node[vecIndex].second) {
        if ((*pVecReebNode)[nodeId].crType == DOWN_FORKING_REEB) {
            // put it into the queue
            ready_to_decide_node.push(vecIndex);
        } else {
            if ((*pVecReebNode)[nodeId].crType == UP_FORKING_REEB) {//
                int left_arc_id = (*pVecReebNode)[nodeId].ptrArcUpId->front()->clusterLabel - 1;
                int right_arc_id = (*pVecReebNode)[nodeId].ptrArcUpId->back()->clusterLabel - 1;
                if (processed_arc_bit[left_arc_id] && processed_arc_bit[right_arc_id]) {
                    ready_to_decide_node.push(vecIndex);
                }
            }
        }
    }
    return;
}

char psbmReebGraph::LevelSetSeparationProperty(_Polygon &lhs, _Polygon &rhs, const int perpAxis, const int rayAxis) {
    _Polygon projLhs = lhs;
    _Polygon projRhs = rhs;
    //
    //for (unsigned int i = 0; i < projLhs.vecPoints.size(); i++)
    //{//
    //	RotateAxisTobeAlignedWithZ_axis(heightDirection[0], heightDirection[1], heightDirection[2],
    //								projLhs.vecPoints[i][0], projLhs.vecPoints[i][1], projLhs.vecPoints[i][2],
    //								projLhs.vecPoints[i]);
    //}
    //for (unsigned int i = 0; i < projRhs.vecPoints.size(); i++)
    //{
    //	RotateAxisTobeAlignedWithZ_axis(heightDirection[0], heightDirection[1], heightDirection[2],
    //								projRhs.vecPoints[i][0], projRhs.vecPoints[i][1], projRhs.vecPoints[i][2],
    //								projRhs.vecPoints[i]);
    //}
    char ret = 0; // they are separated
    Vector3 basePt = projRhs.vecPoints[0];
    if (InsidePolygon(projLhs, basePt, perpAxis, rayAxis)) { // rhs is inside of lhs, they are nested
        ret = 1;
    } else if (InsidePolygon(projRhs, projLhs.vecPoints[0], perpAxis, rayAxis))
        ret = 2; // lhs is inside of rhs, they are nested
//
    return ret;
}

bool psbmReebGraph::InsideRegion(std::vector<_Polygon> &refPolygons, const Vector3 &basePt, const int perpAxis,
                                 const int rayAxis) {
    double axisIntersect = 0.0;
    int Rcross = 0;
    // shift to make base pt as origin
    for (std::vector<_Polygon>::iterator vIter = refPolygons.begin();
         vIter != refPolygons.end();
         vIter++) {
        for (unsigned int i = 0; i < (*vIter).vecPoints.size(); i++) {
            (*vIter).vecPoints[i][0] -= basePt[0];
            (*vIter).vecPoints[i][1] -= basePt[1];
        }
    }
    std::vector<_Polygon>::iterator vIter = refPolygons.begin();
    // the first polygon is the left side , which we want to know
    for (vIter++;
         vIter != refPolygons.end();
         vIter++) {
        for (unsigned int i = 0; i < (*vIter).vecPoints.size() - 1; i++) // last point is the same as the first point
        {
            int i1 = (i + (*vIter).vecPoints.size() - 2) % ((*vIter).vecPoints.size() - 1);
            if (((*vIter).vecPoints[i][perpAxis] > 0.0f && (*vIter).vecPoints[i1][perpAxis] <= 0.0f) ||
                ((*vIter).vecPoints[i1][perpAxis] > 0.0f && (*vIter).vecPoints[i][perpAxis] <= 0.0f)) {
                axisIntersect = ((*vIter).vecPoints[i][rayAxis] * (*vIter).vecPoints[i1][perpAxis] -
                                 (*vIter).vecPoints[i1][rayAxis] * (*vIter).vecPoints[i][perpAxis]) *
                                ((*vIter).vecPoints[i1][perpAxis] - (*vIter).vecPoints[i][perpAxis]);

                if (axisIntersect > 0.0f)
                    Rcross++;
            }
        }
    }
    //std::cout << "Rcross " << Rcross << std::endl;
    if (Rcross % 2 == 1)
        return true;
    return false;
}

bool psbmReebGraph::InsidePolygon(_Polygon &inRefPly, const Vector3 &basePt, const int perpAxis, const int rayAxis) {
    _Polygon refPly = inRefPly;
    double axisIntersect = 0.0;
    int Rcross = 0;
    // shift to make base pt as origin
    for (unsigned int i = 0; i < refPly.vecPoints.size(); i++) {
        refPly.vecPoints[i][0] -= basePt[0];
        refPly.vecPoints[i][1] -= basePt[1];
    }
    for (unsigned int i = 0; i < refPly.vecPoints.size() - 1; i++) // last point is the same as the first point
    {
        int i1 = (i + refPly.vecPoints.size() - 2) % (refPly.vecPoints.size() - 1);
        if ((refPly.vecPoints[i][perpAxis] > 0.0f && refPly.vecPoints[i1][perpAxis] <= 0.0f) ||
            (refPly.vecPoints[i1][perpAxis] > 0.0f && refPly.vecPoints[i][perpAxis] <= 0.0f)) {
            axisIntersect = (refPly.vecPoints[i][rayAxis] * refPly.vecPoints[i1][perpAxis] -
                             refPly.vecPoints[i1][rayAxis] * refPly.vecPoints[i][perpAxis]) *
                            (refPly.vecPoints[i1][perpAxis] - refPly.vecPoints[i][perpAxis]);

            if (axisIntersect > 0.0f)
                Rcross++;
        }
    }
    if (Rcross % 2 == 1)
        return true;
    return false;
}

//tracing out the vertical path in Reeb graph for each pairing
void psbmReebGraph::ComputeUpwardPathBetweenPair(int const u, int const t, psbmReebPath &finalPath,
                                                 psbmReebPath &currentPath) {
    // u the current visiting node
    // t the end node in the upward tracing
    // finalPath
    if (u == t) {
        // compare the finalPath with the currentPath
        // use the size of the path
        if (finalPath.nodeIndexSet.empty() ||
            currentPath.nodeIndexSet.size() > finalPath.nodeIndexSet.size()) {
            finalPath = currentPath;
        }
        return;
    }
    // record the current node

    currentPath.nodeIndexSet.push_back(u);
    // visit each path recursively
    std::list<class psbmReebArc *>::iterator arcIter;
    for (arcIter = (*pVecReebNode)[u].ptrArcUpId->begin();
         arcIter != (*pVecReebNode)[u].ptrArcUpId->end();
         arcIter++) {
        currentPath.simpArcIndexSet.push_back((*arcIter)->clusterLabel - 1);
        //
        int otherEndPoint = TheOtherEndPoint(u, (*arcIter)->clusterLabel);
        ComputeUpwardPathBetweenPair(otherEndPoint, t, finalPath, currentPath);
        //erasing the simplified arc from the path
        currentPath.simpArcIndexSet.pop_back();
    }
    //erasing the current vertex from the path
    currentPath.nodeIndexSet.pop_back();
    //for each upward connected node w
    //	ComputeUpwardPath( w, finalPath, currentPath, t);

    //
    return;
}

//
void psbmReebGraph::Add_Simplified_Arc_into_Graph(const int node_a,
                                                  const int node_b,
                                                  int &nGraphNodeCounter,
                                                  boost::unordered_map<int, int> &inMapping,
                                                  _simpGraph &outGraph) {
    //std::map<int, int>::iterator mIter;
    boost::unordered_map<int, int>::iterator mIter;
    int nGraphNodeId_a = 0;
    int nGraphNodeId_b = 0;
    mIter = inMapping.find(node_a);
    if (mIter == inMapping.end()) {// create new graph node
        nGraphNodeId_a = nGraphNodeCounter;
        _simpGraphNode tmpNode(nGraphNodeCounter, node_a, MeshVertexScalarValue[(*pVecReebNode)[node_a].nVertexId]);
        inMapping[node_a] = nGraphNodeCounter++;
        outGraph.AddNode(tmpNode);
    } else {
        nGraphNodeId_a = mIter->second;//inMapping[node_a];
    }
    mIter = inMapping.find(node_b);
    if (mIter == inMapping.end()) {// create new graph node
        nGraphNodeId_b = nGraphNodeCounter;
        _simpGraphNode tmpNode(nGraphNodeCounter, node_b, MeshVertexScalarValue[(*pVecReebNode)[node_b].nVertexId]);
        inMapping[node_b] = nGraphNodeCounter++;
        outGraph.AddNode(tmpNode);
    } else {
        nGraphNodeId_b = mIter->second;//inMapping[node_b];
    }
    outGraph.AddEdge(nGraphNodeId_a, nGraphNodeId_b);
    return;
}

//
void psbmReebGraph::PathToSourceViaBFS(const int target,
                                       std::vector<int> &vertexOnPath,
                                       std::vector<int> &parents,
                                       _simpGraph &inGraph) {
    int current = target;
    while (current != -1) {
        // push this point into path
        vertexOnPath.push_back(inGraph.vecNode[current].color);
        current = parents[current];
    }
    return;
}

void psbmReebGraph::FindingCorrespondingSimplifiedArcs(
        std::vector<int> &nodeOnPath,
        std::vector<int> &arcOnPath) {
// expecting to find e_i in the following configuration:
// v0--e0--v1--e1--...v_(n-1)--e_(n-1)--vn
    //// e0 is known and in arcOnPath

    for (unsigned int i = 1; i < nodeOnPath.size() - 1; i++) {// always search upward
        bool bFind = false;
        for (std::list<class psbmReebArc *>::iterator listIter = (*pVecReebNode)[nodeOnPath[i]].ptrArcUpId->begin();
             listIter != (*pVecReebNode)[nodeOnPath[i]].ptrArcUpId->end(); listIter++) {
            if ((*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode0 == nodeOnPath[i + 1] ||
                (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode1 == nodeOnPath[i + 1]) {
                arcOnPath.push_back((*listIter)->clusterLabel - 1);
                bFind = true;
                break;
            }
        }
        if (!bFind) {
            for (std::list<class psbmReebArc *>::iterator listIter = (*pVecReebNode)[nodeOnPath[i]].ptrArcDownId->begin();
                 listIter != (*pVecReebNode)[nodeOnPath[i]].ptrArcDownId->end(); listIter++) {
                if ((*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode0 == nodeOnPath[i + 1] ||
                    (*pVecSimplifiedArc)[(*listIter)->clusterLabel - 1].nCriticalNode1 == nodeOnPath[i + 1]) {
                    arcOnPath.push_back((*listIter)->clusterLabel - 1);
                    bFind = true;
                    break;
                }
            }
        }
        if (!bFind) {
            std::cout << "CAN NOT match arc" << std::endl;
            exit(0);
        }
    }
    nodeOnPath.pop_back();
// pop out vn make the configuration as
// // v0--e0--v1--e1--...v_(n-1)--e_(n-1)
}

//
void psbmReebGraph::OptimizeXYAxisMakeNoParaOrtho(std::vector<_Polygon> &inPolygons) {
//
    // calculating the intersection angle delta
    const double optMy_local_pi = 3.1415927;
    int nSegmentNumber = 0;
    for (unsigned int i = 0; i < inPolygons.size(); i++) {
        nSegmentNumber += (int) (inPolygons[i].vecPoints.size() - 1);
    }
    std::vector<double> delta_angle(nSegmentNumber);
    // process the left polygon
    Vector3 segment;
    double xCoord = 0.0;
    double yCoord = 0.0;
    double invPI = 1.0 / optMy_local_pi;
    double halfPI = 0.5 * optMy_local_pi;
    //
    Vector3 _xAxis(1, 0, 0);
    Vector3 _yAxis(0, 1, 0);
    // process the left polygon
    for (std::vector<_Polygon>::iterator vIter = inPolygons.begin();
         vIter != inPolygons.end();
         vIter++) {
        for (unsigned int i = 0; i < (*vIter).vecPoints.size() - 1; i++) {// project the segment onto the XY-plane
            segment = (*vIter).vecPoints[i + 1] - (*vIter).vecPoints[i];
            xCoord = segment * _xAxis;
            yCoord = segment * _yAxis;
            if (xCoord == 0.0) {
                delta_angle[i] = 0.0;
            } else {
                delta_angle[i] = optMy_local_pi + atan(yCoord / xCoord);
                delta_angle[i] = delta_angle[i] - ((int) (delta_angle[i] * 2.0 * invPI)) * halfPI;
            }
        }
    }
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
    double rotAngle = (maxIntervalBegin + maxIntervalEnd) * (-0.5);
    //
    for (std::vector<_Polygon>::iterator vIter = inPolygons.begin();
         vIter != inPolygons.end();
         vIter++) {
        for (unsigned int i = 0; i < (*vIter).vecPoints.size() - 1; i++) {
            RotateAlongAxis_ReebGraph(0.0, 0.0, 1.0,
                                      (*vIter).vecPoints[i][0], (*vIter).vecPoints[i][1], (*vIter).vecPoints[i][2],
                                      rotAngle,
                                      (*vIter).vecPoints[i]);
        }
        (*vIter).vecPoints.back() = (*vIter).vecPoints.front();
    }
    //

    //
    return;
}

void psbmReebGraph::CutAllLevelsets(const int upforkingReebNodeId, std::vector<_Polygon> &levelSetPolygons) {

    //
    const int leftArcId = (*pVecReebNode)[upforkingReebNodeId].ptrArcUpId->front()->clusterLabel - 1;
    int leftUpwardCritReebNodeId = 0;
    const int rightArcId = (*pVecReebNode)[upforkingReebNodeId].ptrArcUpId->back()->clusterLabel - 1;
    int rightUpwardCritReebNodeId = 0;
    leftUpwardCritReebNodeId = (*pVecSimplifiedArc)[leftArcId].nCriticalNode1;
    rightUpwardCritReebNodeId = (*pVecSimplifiedArc)[rightArcId].nCriticalNode1;
    //
    double levelSetValue = 0.5 * MeshVertexScalarValue[(*pVecReebNode)[upforkingReebNodeId].nVertexId];
    levelSetValue = levelSetValue +
                    0.5 * (MeshVertexScalarValue[(*pVecReebNode)[leftUpwardCritReebNodeId].nVertexId] <
                           MeshVertexScalarValue[(*pVecReebNode)[rightUpwardCritReebNodeId].nVertexId] ?
                           MeshVertexScalarValue[(*pVecReebNode)[leftUpwardCritReebNodeId].nVertexId] :
                           MeshVertexScalarValue[(*pVecReebNode)[rightUpwardCritReebNodeId].nVertexId]);
    _Polygon leftLevelSet, rightLevelSet;
    //compute left level set
    CutLevelSet_specific_value(leftLevelSet, (*pVecSimplifiedArc)[leftArcId], levelSetValue);
    //compute right level set
    CutLevelSet_specific_value(rightLevelSet, (*pVecSimplifiedArc)[rightArcId], levelSetValue);
    //
    //std::vector<_Polygon> levelSetPolygons;
    levelSetPolygons.push_back(leftLevelSet);
    levelSetPolygons.push_back(rightLevelSet);
    //
    //
    for (unsigned int i = 0; i < (*pVecSimplifiedArc).size(); i++) {
        if (i != leftArcId && i != rightArcId) {
            int lowVerOnMeshId = (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode0].nVertexId;
            int highVerOnMeshId = (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode1].nVertexId;
            if (MeshVertexScalarValue[lowVerOnMeshId] <= levelSetValue &&
                MeshVertexScalarValue[highVerOnMeshId] >= levelSetValue) {
                _Polygon tmpLevelSetPlg;
                CutLevelSet_specific_value(tmpLevelSetPlg, (*pVecSimplifiedArc)[i], levelSetValue);
                levelSetPolygons.push_back(tmpLevelSetPlg);
                tmpLevelSetPlg.Clear();
            }
        }
    }
    return;
}

int psbmReebGraph::DetermineVerticalLoopType(const int upforkingReebNodeId, std::vector<_Polygon> &levelSetPolygons) {
    // check the path type: horizontal or vertical

    int point_type = -1;
    int perpAxis = 1;
    int rayAxis = 0;
    //
    const int leftArcId = (*pVecReebNode)[upforkingReebNodeId].ptrArcUpId->front()->clusterLabel - 1;
    int leftUpwardCritReebNodeId = 0;
    const int rightArcId = (*pVecReebNode)[upforkingReebNodeId].ptrArcUpId->back()->clusterLabel - 1;
    int rightUpwardCritReebNodeId = 0;
    leftUpwardCritReebNodeId = (*pVecSimplifiedArc)[leftArcId].nCriticalNode1;
    rightUpwardCritReebNodeId = (*pVecSimplifiedArc)[rightArcId].nCriticalNode1;
    //
    double levelSetValue = 0.5 * MeshVertexScalarValue[(*pVecReebNode)[upforkingReebNodeId].nVertexId];
    levelSetValue = levelSetValue +
                    0.5 * (MeshVertexScalarValue[(*pVecReebNode)[leftUpwardCritReebNodeId].nVertexId] <
                           MeshVertexScalarValue[(*pVecReebNode)[rightUpwardCritReebNodeId].nVertexId] ?
                           MeshVertexScalarValue[(*pVecReebNode)[leftUpwardCritReebNodeId].nVertexId] :
                           MeshVertexScalarValue[(*pVecReebNode)[rightUpwardCritReebNodeId].nVertexId]);
    _Polygon leftLevelSet, rightLevelSet;
    //compute left level set
    CutLevelSet_specific_value(leftLevelSet, (*pVecSimplifiedArc)[leftArcId], levelSetValue);
    //compute right level set
    CutLevelSet_specific_value(rightLevelSet, (*pVecSimplifiedArc)[rightArcId], levelSetValue);
    //
    //std::vector<_Polygon> levelSetPolygons;
    levelSetPolygons.push_back(leftLevelSet);
    levelSetPolygons.push_back(rightLevelSet);
    //
    //
    for (unsigned int i = 0; i < (*pVecSimplifiedArc).size(); i++) {
        if (i != leftArcId && i != rightArcId) {
            int lowVerOnMeshId = (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode0].nVertexId;
            int highVerOnMeshId = (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode1].nVertexId;
            if (MeshVertexScalarValue[lowVerOnMeshId] <= levelSetValue &&
                MeshVertexScalarValue[highVerOnMeshId] >= levelSetValue) {
                _Polygon tmpLevelSetPlg;
                CutLevelSet_specific_value(tmpLevelSetPlg, (*pVecSimplifiedArc)[i], levelSetValue);
                levelSetPolygons.push_back(tmpLevelSetPlg);
                tmpLevelSetPlg.Clear();
            }
        }
    }
    //
    for (std::vector<_Polygon>::iterator vIter = levelSetPolygons.begin();
         vIter != levelSetPolygons.end();
         vIter++) {
        for (unsigned int i = 0; i < (*vIter).vecPoints.size(); i++) {
            RotateAxisTobeAlignedWithZ_axis(heightDirection[0], heightDirection[1], heightDirection[2],
                                            (*vIter).vecPoints[i][0], (*vIter).vecPoints[i][1],
                                            (*vIter).vecPoints[i][2],
                                            (*vIter).vecPoints[i]);
        }
    }
    OptimizeXYAxisMakeNoParaOrtho(levelSetPolygons);
    //
    leftLevelSet = levelSetPolygons[0];
    rightLevelSet = levelSetPolygons[1];

    int sepType = LevelSetSeparationProperty(leftLevelSet, rightLevelSet, perpAxis, rayAxis);
    // because levelSetPolygons[0].vecPoints[0] represents the region outside of the closed curve levelSetPolygons[0]
    // so the interior part is just the opposite
    // true ----  interior of the volume
    // flase ---- exterior of the volume
    Vector3 basePt = levelSetPolygons[0].vecPoints[7];
    bool leftArcFlag = !InsideRegion(levelSetPolygons, basePt, perpAxis, rayAxis);
    //
    //std::cout << "separation " << sepType << std::endl;
    //std::cout << "inside " << leftArcFlag << std::endl;
    if (sepType) {
        // sepType == 1 or sepType == 2
        if (sepType == 1) { // left side is larger
            if (leftArcFlag) {//
                point_type = 0;
            } else {
                point_type = 1;
            }
        } else {// sepType == 2
            // right side is larger
            if (leftArcFlag) {
                point_type = 1;
            } else {
                point_type = 0;
            }
        }
    } else {
        if (leftArcFlag) {
            point_type = 1;
        } else {
            point_type = 0;
        }
    }
    return point_type;
    //
}

void psbmReebGraph::TraceOutOneCycleForEachPairing(std::pair<int, int> &critPair, psbmReebPath &cycleInReeb) {
    /*
	each path is in the followign form
	lowPt--arcId--u--arcId--u1---arcId---(highPt)
	that is, hightPt is not there, but the arcId connecting to it is there
	*/
    psbmReebPath leftPath;
    psbmReebPath rightPath;

    psbmReebPath currentPath;

    int leftUpPt = (*pVecReebNode)[critPair.first].ptrArcUpId->front()->clusterLabel - 1;
    int rightUpPt = (*pVecReebNode)[critPair.first].ptrArcUpId->back()->clusterLabel - 1;
    //
    leftUpPt = (*pVecSimplifiedArc)[leftUpPt].nCriticalNode1;
    rightUpPt = (*pVecSimplifiedArc)[rightUpPt].nCriticalNode1;
    // Initializing pathes
    //if (leftUpPt == critPair.second)
    {
        leftPath.nodeIndexSet.push_back(critPair.first);
        leftPath.simpArcIndexSet.push_back((*pVecReebNode)[critPair.first].ptrArcUpId->front()->clusterLabel - 1);
    }
    //if (rightUpPt == critPair.second)
    {
        rightPath.nodeIndexSet.push_back(critPair.first);
        rightPath.simpArcIndexSet.push_back((*pVecReebNode)[critPair.first].ptrArcUpId->back()->clusterLabel - 1);
    }
    if (leftUpPt != critPair.second ||
        rightUpPt != critPair.second) {// at least one path has more than one arc
        int upforkingPtOnOutputGraph = 0;
        int downforkingPtOnOutputGraph = 0;
        int left_upforking_up_node_in_graph = leftUpPt;
        int right_upforking_up_node_in_graph = rightUpPt;
        _simpGraph simpArcGraph;
        std::vector<int> parents;
        //
        construct_simplified_arc_subgraph_between_two_critical_pionts(upforkingPtOnOutputGraph,
                                                                      downforkingPtOnOutputGraph,
                                                                      left_upforking_up_node_in_graph,
                                                                      right_upforking_up_node_in_graph,
                                                                      critPair,
                                                                      simpArcGraph);
        //

        simpArcGraph.RemoveEdgesAdjacentToVertex(upforkingPtOnOutputGraph);
        simpArcGraph.BreathFirstSearch(downforkingPtOnOutputGraph, parents);
        //
        if (leftUpPt != critPair.second) {
            PathToSourceViaBFS(left_upforking_up_node_in_graph, leftPath.nodeIndexSet, parents, simpArcGraph);
            if (leftPath.nodeIndexSet.back() != critPair.second) {
                std::cout << " left Path wrong " << std::endl;
                exit(0);
            }
            // the path containing the highest downforking node
            FindingCorrespondingSimplifiedArcs(leftPath.nodeIndexSet, leftPath.simpArcIndexSet);
        }
        if (rightUpPt != critPair.second) {
            PathToSourceViaBFS(right_upforking_up_node_in_graph, rightPath.nodeIndexSet, parents, simpArcGraph);
            if (rightPath.nodeIndexSet.back() != critPair.second) {
                std::cout << " right Path wrong " << std::endl;
                exit(0);
            }
            // the path containing the highest downforking node
            FindingCorrespondingSimplifiedArcs(rightPath.nodeIndexSet, rightPath.simpArcIndexSet);
        }
    }
    //
    // get one cycle
    cycleInReeb = leftPath;
    cycleInReeb.nodeIndexSet.push_back(critPair.second);
    ////
    cycleInReeb.pathType = -1;
    //std::cout << cycleInReeb.pathType << " type " << std::endl;
    //
    for (int i = rightPath.nodeIndexSet.size() - 1; i >= 0; i--) {
        cycleInReeb.nodeIndexSet.push_back(rightPath.nodeIndexSet[i]);
        cycleInReeb.simpArcIndexSet.push_back(rightPath.simpArcIndexSet[i]);
    }
    return;
}

//
void psbmReebGraph::construct_simplified_arc_subgraph_between_two_critical_pionts(
        int &upforkingPtOnOutputGraph,
        int &downforkingPtOnOutputGraph,
        int &upforking_up_node_a,
        int &upforking_up_node_b,
        std::pair<int, int> &critPair,
        _simpGraph &outGraph) {

    //
    double lowValue = MeshVertexScalarValue[(*pVecReebNode)[critPair.first].nVertexId];
    double highValue = MeshVertexScalarValue[(*pVecReebNode)[critPair.second].nVertexId];
    //
    //std::map<int, int> reebNodeIndexToGraphNodeIndex;
    boost::unordered_map<int, int> reebNodeIndexToGraphNodeIndex;
    //
    int nGraphNodeCounter = 0;
    //
    for (unsigned int i = 0; i < pVecSimplifiedArc->size(); i++) {
        //look at each simplified arc
        // ignore the arc connecting to minimum or maximum
        //if (((*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode0].crType == UP_FORKING_REEB ||
        //	 (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode0].crType == DOWN_FORKING_REEB ) &&
        //	((*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode1].crType == UP_FORKING_REEB ||
        //	 (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode1].crType == DOWN_FORKING_REEB )
        //	)
        if ((*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode0].crType != MAXIMUM_REEB &&
            (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode0].crType != MINIMUM_REEB &&
            (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode0].crType != UNKNOWN_REEB &&
            (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode1].crType != MAXIMUM_REEB &&
            (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode1].crType != MINIMUM_REEB &&
            (*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode1].crType != UNKNOWN_REEB
                ) {// this is a valid candidate
            // if there is an node out of the value range, this edge is not taken
            if (MeshVertexScalarValue[(*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode0].nVertexId] >= lowValue &&
                MeshVertexScalarValue[(*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode0].nVertexId] <= highValue &&
                MeshVertexScalarValue[(*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode1].nVertexId] >= lowValue &&
                MeshVertexScalarValue[(*pVecReebNode)[(*pVecSimplifiedArc)[i].nCriticalNode1].nVertexId] <=
                highValue) {// need to take this arc into the graph
                Add_Simplified_Arc_into_Graph((*pVecSimplifiedArc)[i].nCriticalNode0,
                                              (*pVecSimplifiedArc)[i].nCriticalNode1,
                                              nGraphNodeCounter,
                                              reebNodeIndexToGraphNodeIndex,
                                              outGraph);

            }
        }
    }
    //
    upforkingPtOnOutputGraph = reebNodeIndexToGraphNodeIndex[critPair.first];
    downforkingPtOnOutputGraph = reebNodeIndexToGraphNodeIndex[critPair.second];
    //
    upforking_up_node_a = reebNodeIndexToGraphNodeIndex[upforking_up_node_a];
    upforking_up_node_b = reebNodeIndexToGraphNodeIndex[upforking_up_node_b];
    //
    return;
}

// find the preimage of each cycle in simplified reeb graph
// that is, embed each cycle back to mesh
void psbmReebGraph::EmbedCyclesOnMesh(const int lowCritPt,
                                      const int highCritPt,
                                      psbmReebPath &arcPath,
                                      _Polygon &retPath) {
    /*
	Reeb graph node index mapping to mesh vertex index
	*/
    //std::map<int, int> NodeIndexToMeshVertexIndex;
    //	// inverse mapping of ProcessedVertex
    //for (int i = 0; i < int(ProcessedVertex->size()); i++)
    //{
    //	NodeIndexToMeshVertexIndex[(*ProcessedVertex)[i]] = i;
    //}
    int currentArcId = 0;// stored arc index in arcPath starting from 0 instead of 1
    int currentVertex = 0; // each time, searching edges attached to this vertex
    int nextVertex = 0; // at each time searching
    bool searchingUpward = true;
// finding the highest node in the cycle
// values are in ascending order
//
    int highestValueIndex = 0;
    for (int i = 0; i < int(arcPath.simpArcIndexSet.size()); i++) {
        currentVertex = arcPath.nodeIndexSet[i];
        nextVertex = arcPath.nodeIndexSet[i + 1];
        // determine which direction to search
        if ((MeshVertexScalarValue[(*pVecReebNode)[nextVertex].nVertexId] >
             MeshVertexScalarValue[(*pVecReebNode)[currentVertex].nVertexId]) ||
            (MeshVertexScalarValue[(*pVecReebNode)[nextVertex].nVertexId] ==
             MeshVertexScalarValue[(*pVecReebNode)[currentVertex].nVertexId] &&
             (*pVecReebNode)[nextVertex].nVertexId > (*pVecReebNode)[currentVertex].nVertexId)
                ) {
            searchingUpward = true;
        } else {
            highestValueIndex = i;
            break;
        }
    }
    // search ascending half of the path
    for (int i = 0; i < highestValueIndex; i++) {// int(arcPath.simpArcIndexSet.size())
        currentVertex = arcPath.nodeIndexSet[i];
        nextVertex = arcPath.nodeIndexSet[i + 1];
        currentArcId = arcPath.simpArcIndexSet[i];
        // determine which direction to search
        if ((MeshVertexScalarValue[(*pVecReebNode)[nextVertex].nVertexId] >
             MeshVertexScalarValue[(*pVecReebNode)[currentVertex].nVertexId]) ||
            (MeshVertexScalarValue[(*pVecReebNode)[nextVertex].nVertexId] ==
             MeshVertexScalarValue[(*pVecReebNode)[currentVertex].nVertexId] &&
             (*pVecReebNode)[nextVertex].nVertexId > (*pVecReebNode)[currentVertex].nVertexId)
                ) {
            searchingUpward = true;
        } else
            searchingUpward = false;
        // get edges
        while (currentVertex != nextVertex) {
            if (searchingUpward) {// searching upward
                class psbmReebArc *arcOnPath;
                for (std::list<class psbmReebArc *>::iterator arcIter = (*pVecReebNode)[currentVertex].ptrArcUpId->begin();
                     arcIter != (*pVecReebNode)[currentVertex].ptrArcUpId->end();
                     arcIter++)
                    if ((*arcIter)->clusterLabel == currentArcId + 1) {
                        arcOnPath = *arcIter;
                        break;
                    }
                // find one edge from the set of edges attached to this arc
                int vIndex = currentVertex;//NodeIndexToMeshVertexIndex[currentVertex];
                int vNextIndex = nextVertex;// NodeIndexToMeshVertexIndex[nextVertex];
                //
                float varInValue = 0.0;
                int varIndices = 0;
                int targetIndex = -1;
                int nEndVertex = 0;
                //
                std::list<struct psbmArcListNode *>::iterator arcNodeIter;
                for (arcNodeIter = arcOnPath->pArcList->begin();
                     arcNodeIter != arcOnPath->pArcList->end();
                     arcNodeIter++) {
                    //1) get the edge with endpoint starting with currentVertex
                    if (meshData->vecEdge[(*arcNodeIter)->nEdgeId].v0 == vIndex ||
                        meshData->vecEdge[(*arcNodeIter)->nEdgeId].v1 == vIndex) {// potential edges
                        if (meshData->vecEdge[(*arcNodeIter)->nEdgeId].v0 == vIndex)
                            nEndVertex = meshData->vecEdge[(*arcNodeIter)->nEdgeId].v1;
                        else
                            nEndVertex = meshData->vecEdge[(*arcNodeIter)->nEdgeId].v0;
                        //2 criteria:	1) height <= nextVertx
                        //				2) arc id matched the currentArcId
                        if (MeshVertexScalarValue[nEndVertex] <=
                            MeshVertexScalarValue[(*pVecReebNode)[nextVertex].nVertexId] &&
                            (*pVecAuxMeshEdge)[(*arcNodeIter)->nEdgeId].ptrReebArcId->clusterLabel ==
                            currentArcId + 1) {
                            // eligible edges
                            if (MeshVertexScalarValue[nEndVertex] > MeshVertexScalarValue[vIndex]
                                && MeshVertexScalarValue[nEndVertex] <=
                                   MeshVertexScalarValue[vNextIndex]) {// there must be at least one endpoint with higher value
                                if ((MeshVertexScalarValue[nEndVertex] == MeshVertexScalarValue[vNextIndex]
                                     && nEndVertex <= vNextIndex)
                                    ||
                                    (MeshVertexScalarValue[nEndVertex] < MeshVertexScalarValue[vNextIndex])
                                        ) {
                                    if (varInValue <
                                        MeshVertexScalarValue[nEndVertex] - MeshVertexScalarValue[vIndex]) {
                                        varInValue = MeshVertexScalarValue[nEndVertex] - MeshVertexScalarValue[vIndex];
                                        varIndices = nEndVertex - vIndex;
                                        targetIndex = nEndVertex;
                                    }
                                }
                            } else {
                                if (varInValue == 0.0 &&
                                    MeshVertexScalarValue[nEndVertex] == MeshVertexScalarValue[vIndex]
                                    && nEndVertex > vIndex
                                        )//&& nEndVertex <= vNextIndex)
                                {
                                    varIndices = nEndVertex - vIndex;
                                    targetIndex = nEndVertex;
                                }
                            }
                        }
                    }
                }
                // record the vertex
                if (targetIndex < 0) {
                    std::cout << "WRONG in EMBEDDING CYCLES ONTO MESH" << std::endl;
                } else {
                    retPath.vecPoints.push_back(Vector3(meshData->vecVertex[targetIndex].x,
                                                        meshData->vecVertex[targetIndex].y,
                                                        meshData->vecVertex[targetIndex].z));
                }
                // update iteration conditions
                currentVertex = targetIndex;//(*ProcessedVertex)[targetIndex];
            } else {// searching downward
                class psbmReebArc *arcOnPath;
                for (std::list<class psbmReebArc *>::iterator arcIter = (*pVecReebNode)[currentVertex].ptrArcDownId->begin();
                     arcIter != (*pVecReebNode)[currentVertex].ptrArcDownId->end();
                     arcIter++)
                    if ((*arcIter)->clusterLabel == currentArcId + 1) {
                        arcOnPath = *arcIter;
                        break;
                    }
                // find one edge from the set of edges attached to this arc
                int vIndex = currentVertex;//NodeIndexToMeshVertexIndex[currentVertex];
                int vNextIndex = nextVertex;// NodeIndexToMeshVertexIndex[nextVertex];
                //
                float varInValue = 0.0;
                int varIndices = 0;
                int targetIndex = -1;
                int nEndVertex = 0;
                //
                std::list<struct psbmArcListNode *>::iterator arcNodeIter;
                for (arcNodeIter = arcOnPath->pArcList->begin();
                     arcNodeIter != arcOnPath->pArcList->end();
                     arcNodeIter++) {
                    //1) get the edge with endpoint starting with currentVertex
                    if (meshData->vecEdge[(*arcNodeIter)->nEdgeId].v0 == vIndex ||
                        meshData->vecEdge[(*arcNodeIter)->nEdgeId].v1 == vIndex) {// potential edges
                        if (meshData->vecEdge[(*arcNodeIter)->nEdgeId].v0 == vIndex)
                            nEndVertex = meshData->vecEdge[(*arcNodeIter)->nEdgeId].v1;
                        else
                            nEndVertex = meshData->vecEdge[(*arcNodeIter)->nEdgeId].v0;
                        //2 criteria:	1) height >= nextVertx
                        //				2) arc id matched the currentArcId
                        if (MeshVertexScalarValue[nEndVertex] >=
                            MeshVertexScalarValue[(*pVecReebNode)[nextVertex].nVertexId] &&
                            (*pVecAuxMeshEdge)[(*arcNodeIter)->nEdgeId].ptrReebArcId->clusterLabel ==
                            currentArcId + 1) {
                            // eligible edges
                            if (MeshVertexScalarValue[nEndVertex] < MeshVertexScalarValue[vIndex]
                                && MeshVertexScalarValue[nEndVertex] >=
                                   MeshVertexScalarValue[vNextIndex]) {// there must be at least one endpoint with higher value
                                if ((MeshVertexScalarValue[nEndVertex] == MeshVertexScalarValue[vNextIndex]
                                     && nEndVertex >= vNextIndex)
                                    ||
                                    (MeshVertexScalarValue[nEndVertex] > MeshVertexScalarValue[vNextIndex])
                                        ) {

                                    if (varInValue <
                                        MeshVertexScalarValue[vIndex] - MeshVertexScalarValue[nEndVertex]) {
                                        varInValue = MeshVertexScalarValue[vIndex] - MeshVertexScalarValue[nEndVertex];
                                        varIndices = nEndVertex - vIndex;
                                        targetIndex = nEndVertex;
                                    }
                                }
                            } else {
                                if (varInValue == 0.0 &&
                                    MeshVertexScalarValue[nEndVertex] == MeshVertexScalarValue[vIndex]
                                    && nEndVertex < vIndex
                                    && nEndVertex >= vNextIndex) {
                                    varIndices = vIndex - nEndVertex;
                                    targetIndex = nEndVertex;
                                }
                            }
                        }
                    }
                }
                // record the vertex
                if (targetIndex < 0) {
                    std::cout << "WRONG in EMBEDDING CYCLES ONTO MESH" << std::endl;
                } else {
                    retPath.vecPoints.push_back(Vector3(meshData->vecVertex[targetIndex].x,
                                                        meshData->vecVertex[targetIndex].y,
                                                        meshData->vecVertex[targetIndex].z));
                }
                // update iteration conditions
                currentVertex = targetIndex;//(*ProcessedVertex)[targetIndex];
            }
/*
			// search one edge along up or down direction
			if( searchingUpward)
			{// searching upward
				int vIndex = NodeIndexToMeshVertexIndex[currentVertex];
				// circulating the edges adjacent to current vertex
				// and pick up one edge with largest variance in height
				float varInValue = 0.0;
				int varIndices = 0;
				int targetIndex = -1;
				for (int j = 0; j < int(meshData->vecVertex[vIndex].adjEdges.size()); j++)
				{
					int nEndVertex = 0;
					if (meshData->vecEdge[meshData->vecVertex[vIndex].adjEdges[j]].v0 == vIndex)
						nEndVertex = meshData->vecEdge[meshData->vecVertex[vIndex].adjEdges[j]].v1;
					else
						nEndVertex = meshData->vecEdge[meshData->vecVertex[vIndex].adjEdges[j]].v0;
					if (MeshVertexScalarValue[nEndVertex] > MeshVertexScalarValue[vIndex] )
					{// there must be at least one endpoint with higher value
						if (varInValue < MeshVertexScalarValue[nEndVertex] - MeshVertexScalarValue[vIndex] )
						{
							varInValue = MeshVertexScalarValue[nEndVertex] - MeshVertexScalarValue[vIndex] ;
							varIndices = nEndVertex - vIndex;
							targetIndex = nEndVertex;
						}
					}
					else
					{
						if (varInValue == 0.0 &&
							MeshVertexScalarValue[nEndVertex] == MeshVertexScalarValue[vIndex]
							&& nEndVertex > vIndex)
						{
							varIndices = nEndVertex - vIndex;
							targetIndex = nEndVertex;
						}
					}
				}
				if (targetIndex < 0)
				{
					std::cout << "WRONG in EMBEDDING CYCLES ONTO MESH" << std::endl;
				}
				else
				{
					retPath.vecPoints.push_back(Vector3(meshData->vecVertex[targetIndex].x,
														meshData->vecVertex[targetIndex].y,
														meshData->vecVertex[targetIndex].z));
				}

			}
			else
			{// searching downward
				int vIndex = NodeIndexToMeshVertexIndex[currentVertex];
				// circulating the edges adjacent to current vertex
				// and pick up one edge with largest variance in height
				float varInValue = 0.0;
				int varIndices = 0;
				int targetIndex = -1;
				for (int j = 0; j < int(meshData->vecVertex[vIndex].adjEdges.size()); j++)
				{
					int nEndVertex = 0;
					if (meshData->vecEdge[meshData->vecVertex[vIndex].adjEdges[j]].v0 == vIndex)
						nEndVertex = meshData->vecEdge[meshData->vecVertex[vIndex].adjEdges[j]].v1;
					else
						nEndVertex = meshData->vecEdge[meshData->vecVertex[vIndex].adjEdges[j]].v0;
					if (MeshVertexScalarValue[nEndVertex] < MeshVertexScalarValue[vIndex])
					{// there must be at least one endpoint with lower value
						if (varInValue < MeshVertexScalarValue[vIndex] - MeshVertexScalarValue[nEndVertex])
						{
							varInValue = MeshVertexScalarValue[nEndVertex] - MeshVertexScalarValue[vIndex] ;
							varIndices = nEndVertex - vIndex;
							targetIndex = nEndVertex;
						}
					}
					else
					{
						if (varInValue == 0.0 &&
							MeshVertexScalarValue[nEndVertex] == MeshVertexScalarValue[vIndex]
							&& nEndVertex < vIndex)
						{
							if (varIndices < vIndex - nEndVertex)
							{
								varIndices = vIndex - nEndVertex;
								targetIndex = nEndVertex;
							}
						}
					}
				}
				if (targetIndex < 0)
				{
					std::cout << "WRONG in EMBEDDING CYCLES ONTO MESH" << std::endl;
				}
				else
				{
					retPath.vecPoints.push_back(Vector3(meshData->vecVertex[targetIndex].x,
														meshData->vecVertex[targetIndex].y,
														meshData->vecVertex[targetIndex].z));
				}
			}
*/
        }// while
    }// for
    // searching descending half of the path
    for (int i = arcPath.nodeIndexSet.size() - 1; i > highestValueIndex; i--) {
        currentVertex = arcPath.nodeIndexSet[i];
        nextVertex = arcPath.nodeIndexSet[i - 1];
        currentArcId = arcPath.simpArcIndexSet[i - 1];

        // get edges
        while (currentVertex != nextVertex) {
            class psbmReebArc *arcOnPath;
            for (std::list<class psbmReebArc *>::iterator arcIter = (*pVecReebNode)[currentVertex].ptrArcUpId->begin();
                 arcIter != (*pVecReebNode)[currentVertex].ptrArcUpId->end();
                 arcIter++)
                if ((*arcIter)->clusterLabel == currentArcId + 1) {
                    arcOnPath = *arcIter;
                    break;
                }
            // find one edge from the set of edges attached to this arc
            int vIndex = currentVertex;//NodeIndexToMeshVertexIndex[currentVertex];
            int vNextIndex = nextVertex;//NodeIndexToMeshVertexIndex[nextVertex];
            //
            float varInValue = 0.0;
            int varIndices = 0;
            int targetIndex = -1;
            int nEndVertex = 0;
            //
            std::list<struct psbmArcListNode *>::iterator arcNodeIter;
            for (arcNodeIter = arcOnPath->pArcList->begin();
                 arcNodeIter != arcOnPath->pArcList->end();
                 arcNodeIter++) {
                //1) get the edge with endpoint starting with currentVertex
                if (meshData->vecEdge[(*arcNodeIter)->nEdgeId].v0 == vIndex ||
                    meshData->vecEdge[(*arcNodeIter)->nEdgeId].v1 == vIndex) {// potential edges
                    if (meshData->vecEdge[(*arcNodeIter)->nEdgeId].v0 == vIndex)
                        nEndVertex = meshData->vecEdge[(*arcNodeIter)->nEdgeId].v1;
                    else
                        nEndVertex = meshData->vecEdge[(*arcNodeIter)->nEdgeId].v0;
                    if (nEndVertex == vNextIndex) {
                        targetIndex = nextVertex;
                        break;
                    }
                    //2 criteria:	1) height <= nextVertx
                    //				2) arc id matched the currentArcId
                    if (MeshVertexScalarValue[nEndVertex] <=
                        MeshVertexScalarValue[(*pVecReebNode)[nextVertex].nVertexId] &&
                        (*pVecAuxMeshEdge)[(*arcNodeIter)->nEdgeId].ptrReebArcId->clusterLabel == currentArcId + 1) {
                        // eligible edges
                        if (MeshVertexScalarValue[nEndVertex] >
                            MeshVertexScalarValue[vIndex]) {// there must be at least one endpoint with higher value
                            if ((MeshVertexScalarValue[nEndVertex] == MeshVertexScalarValue[vNextIndex]
                                 && nEndVertex <= vNextIndex)
                                ||
                                (MeshVertexScalarValue[nEndVertex] < MeshVertexScalarValue[vNextIndex])
                                    ) {
                                if (varInValue <= MeshVertexScalarValue[nEndVertex] - MeshVertexScalarValue[vIndex]) {
                                    varInValue = MeshVertexScalarValue[nEndVertex] - MeshVertexScalarValue[vIndex];
                                    varIndices = nEndVertex - vIndex;
                                    targetIndex = nEndVertex;
                                }
                            }
                        } else {
                            if (varInValue == 0.0 &&
                                MeshVertexScalarValue[nEndVertex] == MeshVertexScalarValue[vIndex]
                                && nEndVertex > vIndex
                                    )//&& nEndVertex <= vNextIndex)
                            {
                                varIndices = nEndVertex - vIndex;
                                targetIndex = nEndVertex;
                            }
                        }
                    }
                }
            }
            // record the vertex
            if (targetIndex < 0) {
                std::cout << "WRONG in EMBEDDING CYCLES ONTO MESH" << std::endl;
            } else {
                retPath.vecPoints.push_back(Vector3(meshData->vecVertex[targetIndex].x,
                                                    meshData->vecVertex[targetIndex].y,
                                                    meshData->vecVertex[targetIndex].z));
            }
            // update iteration conditions
            currentVertex = targetIndex;//(*ProcessedVertex)[targetIndex];
        }// while
    }// for


    return;
}

void psbmReebGraph::add_vertex_one_ring_to_graph(const int vertexId_a,
                                                 const int vertexId_b,
                                                 int &nGraphNodeCounter,
                                                 std::map<int, int> &MeshVertexIndexToGraphNodeIndex,
                                                 _simpGraph &outGraph) {

    std::map<int, int>::iterator mIter;
    int nGraphNodeId_a = 0;
    int nGraphNodeId_b = 0;
    mIter = MeshVertexIndexToGraphNodeIndex.find(vertexId_a);
    if (mIter == MeshVertexIndexToGraphNodeIndex.end()) {// create new graph node
        nGraphNodeId_a = nGraphNodeCounter;
        _simpGraphNode tmpNode(nGraphNodeCounter, vertexId_a);
        MeshVertexIndexToGraphNodeIndex[vertexId_a] = nGraphNodeCounter++;
        outGraph.AddNode(tmpNode);
    } else {
        nGraphNodeId_a = MeshVertexIndexToGraphNodeIndex[vertexId_a];
    }
    mIter = MeshVertexIndexToGraphNodeIndex.find(vertexId_b);
    if (mIter == MeshVertexIndexToGraphNodeIndex.end()) {// create new graph node
        nGraphNodeId_b = nGraphNodeCounter;
        _simpGraphNode tmpNode(nGraphNodeCounter, vertexId_b);
        MeshVertexIndexToGraphNodeIndex[vertexId_b] = nGraphNodeCounter++;
        outGraph.AddNode(tmpNode);
    } else {
        nGraphNodeId_b = MeshVertexIndexToGraphNodeIndex[vertexId_b];
    }
    outGraph.AddEdge(nGraphNodeId_a, nGraphNodeId_b);
    return;
}

void psbmReebGraph::add_vertex_one_ring_to_graph(const int vertexId,
                                                 int &nGraphNodeCounter,
                                                 std::map<int, int> &MeshVertexIndexToGraphNodeIndex,
                                                 _simpGraph &outGraph) {
    std::map<int, int>::iterator mIter;
    int nGraphNodeId_a = 0;
    int nGraphNodeId_b = 0;
    mIter = MeshVertexIndexToGraphNodeIndex.find(vertexId);
    if (mIter == MeshVertexIndexToGraphNodeIndex.end()) {// create new graph node
        nGraphNodeId_a = nGraphNodeCounter;
        _simpGraphNode tmpNode(nGraphNodeCounter, vertexId, MeshVertexScalarValue[vertexId]);
        MeshVertexIndexToGraphNodeIndex[vertexId] = nGraphNodeCounter++;
        //
        outGraph.AddNode(tmpNode);
        // set the value

    } else {
        nGraphNodeId_a = MeshVertexIndexToGraphNodeIndex[vertexId];
    }
    // iterate all edges
    for (unsigned int i = 0; i < meshData->vecVertex[vertexId].adjEdges.size(); i++) {
        int nNodeId = 0;
        nNodeId = meshData->vecEdge[meshData->vecVertex[vertexId].adjEdges[i]].v0;
        if (nNodeId == vertexId)
            nNodeId = meshData->vecEdge[meshData->vecVertex[vertexId].adjEdges[i]].v1;

        mIter = MeshVertexIndexToGraphNodeIndex.find(nNodeId);
        if (mIter == MeshVertexIndexToGraphNodeIndex.end()) {// create new graph node
            nGraphNodeId_b = nGraphNodeCounter;
            _simpGraphNode tmpNode(nGraphNodeCounter, nNodeId, MeshVertexScalarValue[nNodeId]);
            MeshVertexIndexToGraphNodeIndex[nNodeId] = nGraphNodeCounter++;
            outGraph.AddNode(tmpNode);
        } else {
            nGraphNodeId_b = MeshVertexIndexToGraphNodeIndex[nNodeId];
        }
        outGraph.AddEdge(nGraphNodeId_a, nGraphNodeId_b);
    }
    return;
}

void psbmReebGraph::add_vertex_one_ring_to_graph(const int vertexId,
                                                 const double lowValue,
                                                 const double highValue,
                                                 const int lowVertexIdOnMesh,
                                                 const int highVertexIdOnMesh,
                                                 int &nGraphNodeCounter,
                                                 std::map<int, int> &MeshVertexIndexToGraphNodeIndex,
                                                 _simpGraph &outGraph) {

    std::map<int, int>::iterator mIter;
    int nGraphNodeId_a = 0;
    int nGraphNodeId_b = 0;
    mIter = MeshVertexIndexToGraphNodeIndex.find(vertexId);
    if (mIter == MeshVertexIndexToGraphNodeIndex.end()) {// create new graph node
        nGraphNodeId_a = nGraphNodeCounter;
        _simpGraphNode tmpNode(nGraphNodeCounter, vertexId, MeshVertexScalarValue[vertexId]);
        MeshVertexIndexToGraphNodeIndex[vertexId] = nGraphNodeCounter++;
        //
        outGraph.AddNode(tmpNode);
        // set the value

    } else {
        nGraphNodeId_a = MeshVertexIndexToGraphNodeIndex[vertexId];
    }
    // iterate all edges
    for (unsigned int i = 0; i < meshData->vecVertex[vertexId].adjEdges.size(); i++) {
        int nNodeId = 0;
        nNodeId = meshData->vecEdge[meshData->vecVertex[vertexId].adjEdges[i]].v0;
        if (nNodeId == vertexId)
            nNodeId = meshData->vecEdge[meshData->vecVertex[vertexId].adjEdges[i]].v1;
        if ((MeshVertexScalarValue[nNodeId] > lowValue ||
             MeshVertexScalarValue[nNodeId] == lowValue && nNodeId >= lowVertexIdOnMesh)
            &&
            (MeshVertexScalarValue[nNodeId] < highValue ||
             MeshVertexScalarValue[nNodeId] == highValue && nNodeId <= highVertexIdOnMesh)
                ) {
            mIter = MeshVertexIndexToGraphNodeIndex.find(nNodeId);
            if (mIter == MeshVertexIndexToGraphNodeIndex.end()) {// create new graph node
                nGraphNodeId_b = nGraphNodeCounter;
                _simpGraphNode tmpNode(nGraphNodeCounter, nNodeId, MeshVertexScalarValue[nNodeId]);
                MeshVertexIndexToGraphNodeIndex[nNodeId] = nGraphNodeCounter++;
                outGraph.AddNode(tmpNode);
            } else {
                nGraphNodeId_b = MeshVertexIndexToGraphNodeIndex[nNodeId];
            }
            outGraph.AddEdge(nGraphNodeId_a, nGraphNodeId_b);
        }
    }
    return;
}

void psbmReebGraph::compute_representative_point_on_edge(const int lowVertexOnMesh,
                                                         const int highVertexOnMesh,
                                                         const int eid,
                                                         Vector3 &repPt) {
    // Preq : this edge intersects the range by lowVertexOnMesh and highVertexOnMesh
    double inWeight = 0.0;
    double lowHeight = MeshVertexScalarValue[meshData->vecEdge[eid].v0];
    double highHeight = MeshVertexScalarValue[meshData->vecEdge[eid].v1];
    if (lowHeight > highHeight) {
        double tmp = lowHeight;
        lowHeight = highHeight;
        highHeight = tmp;
    }
    //
    lowHeight = lowHeight > MeshVertexScalarValue[lowVertexOnMesh] ? lowHeight : MeshVertexScalarValue[lowVertexOnMesh];
    highHeight =
            highHeight < MeshVertexScalarValue[highVertexOnMesh] ? highHeight : MeshVertexScalarValue[highVertexOnMesh];
    inWeight = (lowHeight + highHeight) * 0.5;
    //
    lowHeight = MeshVertexScalarValue[meshData->vecEdge[eid].v0];
    highHeight = MeshVertexScalarValue[meshData->vecEdge[eid].v1];
    double ratio = (inWeight - lowHeight) / (highHeight - lowHeight);
    //
    repPt[0] = (1 - ratio) * meshData->vecVertex[meshData->vecEdge[eid].v0].x +
               ratio * meshData->vecVertex[meshData->vecEdge[eid].v1].x;
    repPt[1] = (1 - ratio) * meshData->vecVertex[meshData->vecEdge[eid].v0].y +
               ratio * meshData->vecVertex[meshData->vecEdge[eid].v1].y;
    repPt[2] = (1 - ratio) * meshData->vecVertex[meshData->vecEdge[eid].v0].z +
               ratio * meshData->vecVertex[meshData->vecEdge[eid].v1].z;
    return;
}

void psbmReebGraph::construct_path_on_mesh_for_reeb_arc_by_crossing_mesh_edges(
        SimplifiedReebGraphArc &inSimplifiedArc,
        std::vector<std::pair<int, int> > &outPathOnMeshPointType,
        std::vector<int> &parents) {
    /*
	// There is on edge path connceting these two critical points
	// This is the second strategy :
	// Each node except the two critical nodes is on the edge
	// Two edges are connected on the mesh if they belong to the same triangle;
	// Extract the subgraph by walking through the edge-triangle
	// This time, the path is ensured in the output subgraph
	*/
    //
    // data used for constructing graph
    //
    const int lowVertexOnMesh = (*pVecReebNode)[inSimplifiedArc.nCriticalNode0].nVertexId;
    const int highVertexOnMesh = (*pVecReebNode)[inSimplifiedArc.nCriticalNode1].nVertexId;
    //
    const double lowValue = MeshVertexScalarValue[lowVertexOnMesh];
    const double highValue = MeshVertexScalarValue[highVertexOnMesh];
    //
    const int inSimplifiedArcIdx = inSimplifiedArc.edgeLabel;
    //
    std::queue<int> Q;
    bool special_case_no_pts_between_two_crit_pts = false;
    for (std::list<class psbmReebArc *>::iterator listIter = (*pVecReebNode)[inSimplifiedArc.nCriticalNode0].ptrArcUpId->begin();
         listIter != (*pVecReebNode)[inSimplifiedArc.nCriticalNode0].ptrArcUpId->end(); listIter++) {
        if ((*listIter)->nNodeId1 == inSimplifiedArc.nCriticalNode1) {
            special_case_no_pts_between_two_crit_pts = true;
            break;
        }
    }
    //if (special_case_no_pts_between_two_crit_pts)
    //	std::cout << "no interior points in walking" << std::endl;
    //
    //std::vector<int> parents(meshData->vecEdge.size());
    for (unsigned int i = 0; i < parents.size(); i++)
        parents[i] = -2;
    //
    /* Initialize the queue by iterating all edges connecting the lowest vertex */
    for (std::vector<int>::iterator vIter = (*meshData).vecVertex[lowVertexOnMesh].adjEdges.begin();
         vIter != meshData->vecVertex[lowVertexOnMesh].adjEdges.end();
         vIter++) {
        int nTheOtherEndpoint = meshData->vecEdge[*vIter].v0 + meshData->vecEdge[*vIter].v1 - lowVertexOnMesh;
        if ((MeshVertexScalarValue[nTheOtherEndpoint] > lowValue) ||
            (MeshVertexScalarValue[nTheOtherEndpoint] == lowValue && nTheOtherEndpoint > lowVertexOnMesh)
                ) {
            if ((MeshVertexScalarValue[nTheOtherEndpoint] < highValue) ||
                (MeshVertexScalarValue[nTheOtherEndpoint] == highValue && nTheOtherEndpoint < highVertexOnMesh)
                    ) {// this vertex is on the arc, then this edge is in
                if (vecVertexMappingToSimplifiedReebArc[nTheOtherEndpoint] == inSimplifiedArcIdx)
                    Q.push(*vIter);
            } else {// now it needs to consider case it has hihger value
                if (map_levelset_on_edge_to_this_reeb_arc(inSimplifiedArcIdx, inSimplifiedArc.nCriticalNode0,
                                                          nTheOtherEndpoint,
                        //(*ProcessedVertex)[nTheOtherEndpoint],
                                                          *vIter, special_case_no_pts_between_two_crit_pts)) {
                    Q.push(*vIter);
                }
            }
            //
        }
        parents[*vIter] = -1;
    }
    //
    if (Q.empty()) {
        //std::cout << "cases " << cases << std::endl;
        std::cout << " Q is empty , can not continue" << std::endl;
        exit(0);
    }
    /*start searching */
    int opposite_vertices = 0;
    int nEdgeAdjToHighVertex = 0;
    while (!Q.empty()) {
        /* Perform the deapth first search to visit edges until the highest vertex is reached */
        nEdgeAdjToHighVertex = Q.front();
        Q.pop();
        /*Compute a representative point on this edge such that its height falls into the range*/
        //compute_representative_point_on_edge(lowVertexOnMesh, highVertexOnMesh, nEdgeAdjToHighVertex, ptOnEdge);
        //edge_to_pt[nEdgeAdjToHighVertex] = ptOnEdge;
        /*Remember this point and add the node represented by this edge into the graph*/
        /*Check if the hightest vertex is on the two triangles which are adjacent to current edge*/
        int currentTriangleId = meshData->vecEdge[nEdgeAdjToHighVertex].AdjTri[0];
        opposite_vertices = meshData->vecTriangle[currentTriangleId].v0 +
                            meshData->vecTriangle[currentTriangleId].v1 +
                            meshData->vecTriangle[currentTriangleId].v2 -
                            (meshData->vecEdge[nEdgeAdjToHighVertex].v0 + meshData->vecEdge[nEdgeAdjToHighVertex].v1);
        if (opposite_vertices == highVertexOnMesh)
            break;
        currentTriangleId = meshData->vecEdge[nEdgeAdjToHighVertex].AdjTri[1];
        opposite_vertices = meshData->vecTriangle[currentTriangleId].v0 +
                            meshData->vecTriangle[currentTriangleId].v1 +
                            meshData->vecTriangle[currentTriangleId].v2 -
                            (meshData->vecEdge[nEdgeAdjToHighVertex].v0 + meshData->vecEdge[nEdgeAdjToHighVertex].v1);
        if (opposite_vertices == highVertexOnMesh)
            break;
        /*If it is, searching is done*/
        /*If not, push non-visited edges which has point in the range into the queue Q*/
        /*Mark these edge are visited by setting their parents*/
        for (int vid = 0; vid < 2; vid++) {
            currentTriangleId = meshData->vecEdge[nEdgeAdjToHighVertex].AdjTri[vid];
            int triEdgeId[3] = {meshData->vecTriangle[currentTriangleId].e01,
                                meshData->vecTriangle[currentTriangleId].e02,
                                meshData->vecTriangle[currentTriangleId].e12};
            for (int eid = 0; eid < 3; eid++) {
                if (parents[triEdgeId[eid]] == -2) {
                    /*check this edge intersect the range or not*/
                    if (meshData->vecEdge[triEdgeId[eid]].v0 != lowVertexOnMesh &&
                        meshData->vecEdge[triEdgeId[eid]].v1 != lowVertexOnMesh) {
                        int v0 = meshData->vecEdge[triEdgeId[eid]].v0;
                        int v1 = meshData->vecEdge[triEdgeId[eid]].v1;
                        double v0Height = MeshVertexScalarValue[v0];
                        double v1Height = MeshVertexScalarValue[v1];
                        if (v0Height > v1Height ||
                            (v0Height == v1Height && v0 > v1)) {
                            double tmpDouble = v0Height;
                            v0Height = v1Height;
                            v1Height = tmpDouble;
                            int tmpInt = v0;
                            v0 = v1;
                            v1 = tmpInt;
                        }
                        if (!((v1Height < lowValue) ||
                              (v1Height == lowValue && v1 <= lowVertexOnMesh) ||
                              (v0Height > highValue) ||
                              (v0Height == highValue && v0 >= highVertexOnMesh)
                        )
                                ) {
                            if (map_levelset_on_edge_to_this_reeb_arc(inSimplifiedArcIdx, v0,
                                                                      v1, //(*ProcessedVertex)[v0], (*ProcessedVertex)[v1],
                                                                      triEdgeId[eid],
                                                                      special_case_no_pts_between_two_crit_pts)) {
                                parents[triEdgeId[eid]] = nEdgeAdjToHighVertex;
                                Q.push(triEdgeId[eid]);
                            }
                        }
                    }
                }
            }
        }
        //
    }// while
    //
    //_Polygon outPath;
    //outPath.vecPoints.push_back(Vector3(meshData->vecVertex[highVertexOnMesh].x,
    //								meshData->vecVertex[highVertexOnMesh].y,
    //								meshData->vecVertex[highVertexOnMesh].z));
    outPathOnMeshPointType.push_back(std::pair<int, int>(highVertexOnMesh, 0));
    //
    int current_node = nEdgeAdjToHighVertex;
    while (parents[current_node] != -1) {// the point one the edge adjacent to the lowest vertex is ignored
        //outPath.vecPoints.push_back(edge_to_pt[current_node]);
        outPathOnMeshPointType.push_back(std::pair<int, int>(current_node, 1));
        current_node = parents[current_node];
    }
    //
    //outPath.vecPoints.push_back(Vector3(meshData->vecVertex[lowVertexOnMesh].x,
    //								meshData->vecVertex[lowVertexOnMesh].y,
    //								meshData->vecVertex[lowVertexOnMesh].z));
    outPathOnMeshPointType.push_back(std::pair<int, int>(lowVertexOnMesh, 0));
    //

}

void psbmReebGraph::ComputePathToSource(_Polygon &path, std::vector<int> &vertexInMesh, const int target,
                                        std::vector<int> &parents, _simpGraph &inGraph) {
    int current = target;
    while (current != -1) {
        // push this point into path
        path.vecPoints.push_back(Vector3(meshData->vecVertex[inGraph.vecNode[current].color].x,
                                         meshData->vecVertex[inGraph.vecNode[current].color].y,
                                         meshData->vecVertex[inGraph.vecNode[current].color].z));
        vertexInMesh.push_back(inGraph.vecNode[current].color);
        current = parents[current];
    }
    return;
}

void psbmReebGraph::ComputePathToSource_pair(_Polygon &path, std::vector<std::pair<int, int> > &vertexInMesh,
                                             const int target, std::vector<int> &parents, _simpGraph &inGraph) {
    int current = target;
    while (current != -1) {
        // push this point into path
        path.vecPoints.push_back(Vector3(meshData->vecVertex[inGraph.vecNode[current].color].x,
                                         meshData->vecVertex[inGraph.vecNode[current].color].y,
                                         meshData->vecVertex[inGraph.vecNode[current].color].z));
        vertexInMesh.push_back(std::pair<int, int>(inGraph.vecNode[current].color, 0));
        current = parents[current];
    }
    return;
}

void
psbmReebGraph::ComputePathToSource(_Polygon &path, const int target, std::vector<int> &parents, _simpGraph &inGraph) {
    int current = target;
    while (current != -1) {
        // push this point into path
        path.vecPoints.push_back(Vector3(meshData->vecVertex[inGraph.vecNode[current].color].x,
                                         meshData->vecVertex[inGraph.vecNode[current].color].y,
                                         meshData->vecVertex[inGraph.vecNode[current].color].z));
        current = parents[current];
    }
    return;
}

void psbmReebGraph::compute_representative_point_on_edge(const int lowVertexOnMesh,
                                                         const int highVertexOnMesh,
                                                         const int eid,
                                                         const int dirVertrex,
                                                         Vector3 &repPt_org,
                                                         Vector3 &repPt_dir) {
    // Preq : this edge intersects the range by lowVertexOnMesh and highVertexOnMesh
    double inWeight_org = 0.0;
    double inWeight_dir = 0.0;
    double lowHeight = MeshVertexScalarValue[meshData->vecEdge[eid].v0];
    double highHeight = MeshVertexScalarValue[meshData->vecEdge[eid].v1];
    bool bVertex = true; // true ---  highHeight point
    // false --- lowHeight point
    if (meshData->vecEdge[eid].v0 == dirVertrex)
        bVertex = false;
    if (lowHeight > highHeight) {
        double tmp = lowHeight;
        lowHeight = highHeight;
        highHeight = tmp;
        //
        bVertex = !bVertex;
    }
    //
    lowHeight = lowHeight > MeshVertexScalarValue[lowVertexOnMesh] ? lowHeight : MeshVertexScalarValue[lowVertexOnMesh];
    highHeight =
            highHeight < MeshVertexScalarValue[highVertexOnMesh] ? highHeight : MeshVertexScalarValue[highVertexOnMesh];
    if (bVertex) {
        inWeight_org = 0.75 * lowHeight + 0.25 * highHeight;
        inWeight_dir = 0.25 * lowHeight + 0.75 * highHeight;
    } else {
        inWeight_org = 0.25 * lowHeight + 0.75 * highHeight;
        inWeight_dir = 0.75 * lowHeight + 0.25 * highHeight;
    }
    //
    lowHeight = MeshVertexScalarValue[meshData->vecEdge[eid].v0];
    highHeight = MeshVertexScalarValue[meshData->vecEdge[eid].v1];
    double inv_diff_high_low = 1.0 / (highHeight - lowHeight);
    double ratio = (inWeight_org - lowHeight) * inv_diff_high_low;// / (highHeight - lowHeight);
    //
    repPt_org[0] = (1 - ratio) * meshData->vecVertex[meshData->vecEdge[eid].v0].x +
                   ratio * meshData->vecVertex[meshData->vecEdge[eid].v1].x;
    repPt_org[1] = (1 - ratio) * meshData->vecVertex[meshData->vecEdge[eid].v0].y +
                   ratio * meshData->vecVertex[meshData->vecEdge[eid].v1].y;
    repPt_org[2] = (1 - ratio) * meshData->vecVertex[meshData->vecEdge[eid].v0].z +
                   ratio * meshData->vecVertex[meshData->vecEdge[eid].v1].z;
    //
    ratio = (inWeight_dir - lowHeight) * inv_diff_high_low; // / (highHeight - lowHeight);
    //
    repPt_dir[0] = (1 - ratio) * meshData->vecVertex[meshData->vecEdge[eid].v0].x +
                   ratio * meshData->vecVertex[meshData->vecEdge[eid].v1].x;
    repPt_dir[1] = (1 - ratio) * meshData->vecVertex[meshData->vecEdge[eid].v0].y +
                   ratio * meshData->vecVertex[meshData->vecEdge[eid].v1].y;
    repPt_dir[2] = (1 - ratio) * meshData->vecVertex[meshData->vecEdge[eid].v0].z +
                   ratio * meshData->vecVertex[meshData->vecEdge[eid].v1].z;
    //
    return;
}

void psbmReebGraph::ComputeParallelPath_pair_for_path_crossing_mesh_edges(
        std::vector<std::pair<int, int> > &inPathVertexIndexOnMesh,
        std::vector<std::pair<int, int> > &outPathPair,
        _Polygon &inPath, _Polygon &outPath,
        const int lowIndexOnMesh,
        const int highIndexOnMesh) {
//inPathVertexIndexInMesh [index_bit_for_ver_or_edge, vertex(0)_or_edge(1)]
//starting and ending points are vertices on mesh
//All other points are on the edge
    const double lowValue = MeshVertexScalarValue[lowIndexOnMesh];
    const double highValue = MeshVertexScalarValue[highIndexOnMesh];
    outPath.vecPoints.push_back(Vector3(meshData->vecVertex[inPathVertexIndexOnMesh[0].first].x,
                                        meshData->vecVertex[inPathVertexIndexOnMesh[0].first].y,
                                        meshData->vecVertex[inPathVertexIndexOnMesh[0].first].z));
    inPath.vecPoints.push_back(outPath.vecPoints[0]);
    // store the point type
    outPathPair = inPathVertexIndexOnMesh;//.push_back(inPathVertexIndexOnMesh[0]);
    //
    int currentEdgeId = currentEdgeId = inPathVertexIndexOnMesh[1].first;// first edge
    int nextEdgeId = 0;
    int dirVertex = meshData->vecEdge[currentEdgeId].v0;
    int sharedVertexId = 0;
    Vector3 repPt_org;
    Vector3 repPt_dir;
    //
    // invariant : dirVertex is a vertex on the currentEdgeId
    for (unsigned int i = 1; i < inPathVertexIndexOnMesh.size() - 1; i++) {
        //
        //outPathPair.push_back(inPathVertexIndexOnMesh[i]);
        currentEdgeId = inPathVertexIndexOnMesh[i].first;
        //
        /*compute two points on this edge such that the one near dirVertex is on the output path*/
        compute_representative_point_on_edge(lowIndexOnMesh, highIndexOnMesh,
                                             currentEdgeId,
                                             dirVertex,
                                             repPt_org,
                                             repPt_dir);
        //
        inPath.vecPoints.push_back(repPt_org);
        outPath.vecPoints.push_back(repPt_dir);
        // update the dirVertex for next edge
        // No updating for last edge
        if (i < inPathVertexIndexOnMesh.size() - 2) {
            nextEdgeId = inPathVertexIndexOnMesh[i + 1].first;
            sharedVertexId = meshData->vecEdge[nextEdgeId].v0;
            if (meshData->vecEdge[currentEdgeId].v0 != sharedVertexId &&
                meshData->vecEdge[currentEdgeId].v1 != sharedVertexId) {
                sharedVertexId = meshData->vecEdge[nextEdgeId].v1;
            }
            if (dirVertex !=
                sharedVertexId) {// update dirVertex to the vertex on next edge which is not shared by current edge and next edge
                dirVertex = meshData->vecEdge[nextEdgeId].v0 + meshData->vecEdge[nextEdgeId].v1 - sharedVertexId;
            }
        }
    }
    //
    //outPathPair.push_back(inPathVertexIndexOnMesh.back());
    outPath.vecPoints.push_back(Vector3(meshData->vecVertex[inPathVertexIndexOnMesh.back().first].x,
                                        meshData->vecVertex[inPathVertexIndexOnMesh.back().first].y,
                                        meshData->vecVertex[inPathVertexIndexOnMesh.back().first].z));
    inPath.vecPoints.push_back(outPath.vecPoints.back());
    return;
}

void psbmReebGraph::EdgeIntersectArcs(const int edgeId, std::set<int> vecArcId) {
    auxMeshEdge edgePtr = (*pVecAuxMeshEdge)[edgeId];
    struct psbmArcListNode *travEdge = edgePtr.ptrPos;

    while (travEdge) {
        vecArcId.insert(travEdge->ptrReebArcId->clusterLabel);
        //
        travEdge = travEdge->ptrNextPos;
    }
}

void psbmReebGraph::EliminatePathInGraph(const int target, std::vector<int> &parents, _simpGraph &inGraph) {
    int current = target;
    std::vector<int> nodesOnPath;
    // don't erase the target
    current = parents[current];
    if (current == -1)
        return;// the path is only one point
    while (parents[current] != -1)// don't erase the source
    {
        if (!(parents[parents[current]] == -1 && inGraph.vecNode[parents[current]].adjList.size() == 1))
            nodesOnPath.push_back(current);
        current = parents[current];
    }
    for (int i = 0; i < int(nodesOnPath.size()); i++)
        inGraph.RemoveEdgesAdjacentToVertex(nodesOnPath[i]);
}

bool psbmReebGraph::map_levelset_on_edge_to_this_reeb_arc(
        const int target_arc_idx,
        const int nReebNode_0,
        const int nReebNode_1,
        const int inEdgeIdx,
        const bool special_case_no_points_between_two_crit_pts) {
    // preq : 1) height range of the edge intersects the height range of the arc
    // these two conditions should be checked before calling this function
    bool ret = false;
    int ret_arc_idx = -1;
    //
    const int lowVerOnMesh = (*pVecReebNode)[(*pVecSimplifiedArc)[target_arc_idx - 1].nCriticalNode0].nVertexId;
    const int highVerOnMesh = (*pVecReebNode)[(*pVecSimplifiedArc)[target_arc_idx - 1].nCriticalNode1].nVertexId;
    const double lowValue = MeshVertexScalarValue[lowVerOnMesh];
    const double highValue = MeshVertexScalarValue[highVerOnMesh];
    //
    if (vecVertexMappingToSimplifiedReebArc[(*meshData).vecEdge[inEdgeIdx].v0] == target_arc_idx ||
        vecVertexMappingToSimplifiedReebArc[(*meshData).vecEdge[inEdgeIdx].v1] == target_arc_idx) {
        ret = true;
    } else {
        if (((MeshVertexScalarValue[(*meshData).vecEdge[inEdgeIdx].v0] > lowValue ||
              (MeshVertexScalarValue[(*meshData).vecEdge[inEdgeIdx].v0] == lowValue &&
               (*meshData).vecEdge[inEdgeIdx].v0 > lowVerOnMesh)) &&
             (MeshVertexScalarValue[(*meshData).vecEdge[inEdgeIdx].v0] < highValue ||
              (MeshVertexScalarValue[(*meshData).vecEdge[inEdgeIdx].v0] == highValue &&
               (*meshData).vecEdge[inEdgeIdx].v0 < highVerOnMesh))
            ) ||
            ((MeshVertexScalarValue[(*meshData).vecEdge[inEdgeIdx].v1] > lowValue ||
              (MeshVertexScalarValue[(*meshData).vecEdge[inEdgeIdx].v1] == lowValue &&
               (*meshData).vecEdge[inEdgeIdx].v1 > lowVerOnMesh)) &&
             (MeshVertexScalarValue[(*meshData).vecEdge[inEdgeIdx].v1] < highValue ||
              (MeshVertexScalarValue[(*meshData).vecEdge[inEdgeIdx].v1] == highValue &&
               (*meshData).vecEdge[inEdgeIdx].v1 < highVerOnMesh))
            )
                ) {// v0 or v1 falling into the height range of the arc
            ret = false;
        } else {
            // first check if input edge is mapped to both a simplified arc and non-simplified arc
            if (!vecVertexMappingToSimplifiedReebArc[(*meshData).vecEdge[inEdgeIdx].v0] &&
                !vecVertexMappingToSimplifiedReebArc[(*meshData).vecEdge[inEdgeIdx].v1]) {// it is mapped to a simplified arc
                // further check it is mapped to a non-simplified arc
                int lowNode = nReebNode_0;
                if (MeshVertexScalarValue[(*pVecReebNode)[nReebNode_0].nVertexId] >
                    MeshVertexScalarValue[(*pVecReebNode)[nReebNode_1].nVertexId] ||
                    (MeshVertexScalarValue[(*pVecReebNode)[nReebNode_0].nVertexId] ==
                     MeshVertexScalarValue[(*pVecReebNode)[nReebNode_1].nVertexId] &&
                     (*pVecReebNode)[nReebNode_0].nVertexId > (*pVecReebNode)[nReebNode_1].nVertexId)
                        ) {
                    lowNode = nReebNode_1;
                }
                int highNode = nReebNode_0 + nReebNode_1 - lowNode;
                //
                for (std::list<class psbmReebArc *>::iterator listIter = (*pVecReebNode)[lowNode].ptrArcUpId->begin();
                     listIter != (*pVecReebNode)[lowNode].ptrArcUpId->end(); listIter++) {
                    if ((*listIter)->nNodeId1 == highNode) {
                        ret_arc_idx = (*listIter)->clusterLabel;
                        break;
                    }
                }
            }
            if (ret_arc_idx < 0) {// it is not an non-simplified arc
                // Even it is an simplified arc, there could have two simplified arcs connecting these two nodes
                //std::set<int> VisitedEdges;
                boost::unordered_set<int> VisitedEdges;
                std::queue<int> Q;
                bool reach_low_ver = false;
                if (nReebNode_0 == lowVerOnMesh || nReebNode_1 == lowVerOnMesh)
                    reach_low_ver = true;
                bool reach_high_ver = false;
                if (nReebNode_0 == highVerOnMesh || nReebNode_1 == highVerOnMesh)
                    reach_high_ver = true;
                // since this case, lowVerOnMesh and highVerOnMesh is not an edge
                // reach_low_ver and  reach_high_ver can be both true so far
                Q.push(inEdgeIdx);
                VisitedEdges.insert(inEdgeIdx);
                bool bFindBelongings = false;
                while (!Q.empty()) {
                    const int current_edge_idx = Q.front();
                    Q.pop();
                    //
                    int lowEdgeNode = (*meshData).vecEdge[current_edge_idx].v0;
                    if (MeshVertexScalarValue[lowEdgeNode] >
                        MeshVertexScalarValue[(*meshData).vecEdge[current_edge_idx].v1] ||
                        (MeshVertexScalarValue[lowEdgeNode] ==
                         MeshVertexScalarValue[(*meshData).vecEdge[current_edge_idx].v1] &&
                         lowEdgeNode > (*meshData).vecEdge[current_edge_idx].v1)
                            ) {
                        lowEdgeNode = (*meshData).vecEdge[current_edge_idx].v1;
                    }
                    int highEdgeNode =
                            (*meshData).vecEdge[current_edge_idx].v0 + (*meshData).vecEdge[current_edge_idx].v1 -
                            lowEdgeNode;
                    //
                    int nTriangleIdx[2] = {(*meshData).vecEdge[current_edge_idx].AdjTri[0],
                                           (*meshData).vecEdge[current_edge_idx].AdjTri[1]};
                    //
                    for (int i = 0; i < 2; i++) {
                        int nOppositeVertexIdx = (*meshData).vecTriangle[nTriangleIdx[i]].v0 +
                                                 (*meshData).vecTriangle[nTriangleIdx[i]].v1 +
                                                 (*meshData).vecTriangle[nTriangleIdx[i]].v2 -
                                                 ((*meshData).vecEdge[current_edge_idx].v0 +
                                                  (*meshData).vecEdge[current_edge_idx].v1);
                        if (special_case_no_points_between_two_crit_pts && nOppositeVertexIdx == highVerOnMesh) {
                            reach_high_ver = true;
                            if (reach_low_ver) {
                                ret = true;
                                bFindBelongings = true;
                                break;
                            }
                        }
                        if (special_case_no_points_between_two_crit_pts && nOppositeVertexIdx == lowVerOnMesh) {
                            reach_low_ver = true;
                            if (reach_high_ver) {
                                ret = true;
                                bFindBelongings = true;
                                break;
                            }
                        }
                        //
                        if ((MeshVertexScalarValue[nOppositeVertexIdx] > lowValue ||
                             (MeshVertexScalarValue[nOppositeVertexIdx] == lowValue &&
                              nOppositeVertexIdx > lowVerOnMesh)) &&
                            (MeshVertexScalarValue[nOppositeVertexIdx] < highValue ||
                             (MeshVertexScalarValue[nOppositeVertexIdx] == highValue &&
                              nOppositeVertexIdx < highVerOnMesh))
                                ) {
                            if (vecVertexMappingToSimplifiedReebArc[nOppositeVertexIdx] == target_arc_idx) {
                                ret = true;
                            } else {
                                ret = false;
                            }
                            //
                            bFindBelongings = true;
                            break;
                        } else {
                            //
                            int nEdge_0 = 0;
                            if (MeshVertexScalarValue[nOppositeVertexIdx] > highValue ||
                                (MeshVertexScalarValue[nOppositeVertexIdx] == highValue &&
                                 nOppositeVertexIdx > highVerOnMesh) ||
                                nOppositeVertexIdx == highVerOnMesh
                                    ) {
                                meshData->EdgeConnectingTwoVertices(nOppositeVertexIdx, lowEdgeNode, current_edge_idx,
                                                                    nTriangleIdx[i], nEdge_0);
                                if (VisitedEdges.find(nEdge_0) == VisitedEdges.end()) {// not visited before
                                    Q.push(nEdge_0);
                                    VisitedEdges.insert(nEdge_0);
                                }
                            } else {
                                meshData->EdgeConnectingTwoVertices(nOppositeVertexIdx, highEdgeNode, current_edge_idx,
                                                                    nTriangleIdx[i], nEdge_0);
                                if (VisitedEdges.find(nEdge_0) == VisitedEdges.end()) {// not visited before
                                    Q.push(nEdge_0);
                                    VisitedEdges.insert(nEdge_0);
                                }
                            }
                        }
                    }// for 2
                    //
                    if (bFindBelongings)
                        break;
                } // while
            } // if
            else {
                if (ret_arc_idx == target_arc_idx)
                    ret = true;
                else
                    ret = false;
            }
        } // else falling into height range
    }// else == target_arc_idx
    return ret;
}

void psbmReebGraph::compute_path_on_mesh_for_each_simplified_arc() {
    //
    initSimplifiedArcIdForHorizontalLoops.resize(initVerticalLoops->size());
    //
    levelsetOnMesh.resize(initVerticalLoops->size());
    levelsetOnMeshPointType.resize(initVerticalLoops->size());

    //
    pathArcOnMesh.resize(pVecSimplifiedArc->size());
    offsetPathArcOnMesh.resize(pVecSimplifiedArc->size());
    //
    pathArcOnMeshPointType.resize(pVecSimplifiedArc->size());
    offsetPathArcOnMeshPointType.resize(pVecSimplifiedArc->size());
    //
    // collect only necessary arcs
    std::set<int> essentialArcs;
    for (unsigned int i = 0; i < initVerticalLoops->size(); i++) {
        std::copy((*initVerticalLoops)[i].simpArcIndexSet.begin(),
                  (*initVerticalLoops)[i].simpArcIndexSet.end(),
                  std::inserter(essentialArcs, essentialArcs.end()));
    }
    //
    //std::cout << "essential cycle # : " << essentialArcs.size() << std::endl;
    //
    std::vector<int> parents(meshData->vecEdge.size());
    //
    //
    //std::vector<std::vector<std::pair<int, int>>> edgeIntersectArcs(meshData->vecEdge.size());
    edgeIntersectArcs.resize(meshData->vecEdge.size());
    // Only compute pathes for necessary arcs
    //for (unsigned int i = 0; i < pVecSimplifiedArc->size(); i++)
    //{
    //	essentialArcs.insert(i);
    //}
    //
    std::vector<bool> special_flag_for_arcs(pVecSimplifiedArc->size(), false);
    //
    for (std::set<int>::iterator sIter = essentialArcs.begin();
         sIter != essentialArcs.end();
         sIter++) {
        int arcId = *sIter;
        const int lowestVertrexIndexOnMesh = (*pVecReebNode)[(*pVecSimplifiedArc)[*sIter].nCriticalNode0].nVertexId;
        const int highestVertexIndexOnMesh = (*pVecReebNode)[(*pVecSimplifiedArc)[*sIter].nCriticalNode1].nVertexId;
        // clean data structures
        //
        //Check if this simplified arc is just an edge
        // check if the two arcs starting at nCriticalNode0 share the same end points or not
        //
        bool bIsEdge = false;
        bool bProcessedAsAnEdge = false;
        bool special_case_no_interior_pts = false;

        int nArcAsEdge = -1;
        //
        int sepcial_case_counter = 0;
        for (std::list<class psbmReebArc *>::iterator listIter = (*pVecReebNode)[(*pVecSimplifiedArc)[*sIter].nCriticalNode0].ptrArcUpId->begin();
             listIter != (*pVecReebNode)[(*pVecSimplifiedArc)[*sIter].nCriticalNode0].ptrArcUpId->end(); listIter++) {
            if ((*listIter)->nNodeId1 == (*pVecSimplifiedArc)[*sIter].nCriticalNode1) {
                special_case_no_interior_pts = true;
                sepcial_case_counter++;
            }
        }
        special_flag_for_arcs[arcId] = special_case_no_interior_pts;
        //
        if (sepcial_case_counter > 1) {
            std::cout << "Two arcs have the same endpoints but both without interior pts" << std::endl;
            std::cout << "consider another height direction" << std::endl;
            exit(2);
        }
        for (std::vector<int>::iterator vIter = meshData->vecVertex[lowestVertrexIndexOnMesh].adjEdges.begin();
             vIter != meshData->vecVertex[lowestVertrexIndexOnMesh].adjEdges.end();
             vIter++) {
            int nTheOtherEndPoint =
                    meshData->vecEdge[*vIter].v0 + meshData->vecEdge[*vIter].v1 - lowestVertrexIndexOnMesh;
            if (nTheOtherEndPoint == highestVertexIndexOnMesh) {
                bIsEdge = true;
                nArcAsEdge = *vIter;
                break;
            }
        }
        if (bIsEdge) {// this simplified arc is just an edge in the mesh
            //

            if (special_case_no_interior_pts) {
                bProcessedAsAnEdge = true;
                //std::cout << "no interior points " << std::endl;
            } else {
                if (map_levelset_on_edge_to_this_reeb_arc(*sIter + 1, (*pVecSimplifiedArc)[*sIter].nCriticalNode0,
                                                          (*pVecSimplifiedArc)[*sIter].nCriticalNode1, nArcAsEdge))
                    bProcessedAsAnEdge = true;
            }
            if (bProcessedAsAnEdge) {
                //s
                //std::cout << "n0 " << lowestVertrexIndexOnMesh << " n1 " << highestVertexIndexOnMesh << std::endl;
                //std::cout << "arc " << *sIter << std::endl;
                pathArcOnMesh[*sIter].vecPoints.push_back(Vector3(meshData->vecVertex[highestVertexIndexOnMesh].x,
                                                                  meshData->vecVertex[highestVertexIndexOnMesh].y,
                                                                  meshData->vecVertex[highestVertexIndexOnMesh].z));
                pathArcOnMeshPointType[*sIter].push_back(std::pair<int, int>(highestVertexIndexOnMesh, 0));
                //
                pathArcOnMesh[*sIter].vecPoints.push_back(Vector3(meshData->vecVertex[lowestVertrexIndexOnMesh].x,
                                                                  meshData->vecVertex[lowestVertrexIndexOnMesh].y,
                                                                  meshData->vecVertex[lowestVertrexIndexOnMesh].z));
                pathArcOnMeshPointType[*sIter].push_back(std::pair<int, int>(lowestVertrexIndexOnMesh, 0));

                //
                offsetPathArcOnMesh[*sIter] = pathArcOnMesh[*sIter];
                offsetPathArcOnMeshPointType[*sIter] = pathArcOnMeshPointType[*sIter];
                //
                //
                special_flag_for_arcs[arcId] = false;
            }

        }
        if (!bProcessedAsAnEdge) {
            construct_path_on_mesh_for_reeb_arc_by_crossing_mesh_edges((*pVecSimplifiedArc)[*sIter],
                                                                       pathArcOnMeshPointType[*sIter], parents);
            //
            ComputeParallelPath_pair_for_path_crossing_mesh_edges(pathArcOnMeshPointType[*sIter],
                                                                  offsetPathArcOnMeshPointType[*sIter],
                                                                  pathArcOnMesh[*sIter],
                                                                  offsetPathArcOnMesh[*sIter],
                                                                  lowestVertrexIndexOnMesh,
                                                                  highestVertexIndexOnMesh);

            for (int ik = 0; ik < pathArcOnMeshPointType[*sIter].size(); ik++) {
                if (pathArcOnMeshPointType[*sIter][ik].second)
                    edgeIntersectArcs[pathArcOnMeshPointType[*sIter][ik].first].push_back(
                            std::pair<int, int>(arcId, ik));
            }
        }
        //
    }
    //
    // decide which arc to take and what value for level sets
    std::vector<double> levelSetValues(initVerticalLoops->size());
    for (unsigned int i = 0; i < criticalPairing->size(); i++) {//
        double nextUpforkingNodeHeight = 0.0;
        if (i == 36)
            i = 36;
        if (i < criticalPairing->size() - 1) {
            nextUpforkingNodeHeight = MeshVertexScalarValue[(*pVecReebNode)[(*criticalPairing)[i +
                                                                                               1].second].nVertexId];
        } else {
            nextUpforkingNodeHeight = MeshVertexScalarValue[(*pVecReebNode)[(*criticalPairing)[i].first].nVertexId];
        }
        //

        //
        int arcIdx[2] = {(*pVecReebNode)[(*criticalPairing)[i].second].ptrArcUpId->front()->clusterLabel - 1,
                         (*pVecReebNode)[(*criticalPairing)[i].second].ptrArcUpId->back()->clusterLabel - 1};
        //
        std::set<double> distDistruction;
        double maxInterval[2][2];
        for (int j = 0; j < 2; j++) {
            double tempHeight = 0.0;
            double upBound = nextUpforkingNodeHeight;
            //tempHeight = pathArcOnMesh[arcIdx[j]].vecPoints.front() * heightDirection;
            tempHeight = offsetPathArcOnMesh[arcIdx[j]].vecPoints.front() * heightDirection;
            if (upBound > tempHeight)
                upBound = tempHeight;
            //
            distDistruction.insert(upBound);

            //for (std::vector<Vector3>::iterator vIter = pathArcOnMesh[arcIdx[j]].vecPoints.begin();
            //	vIter != pathArcOnMesh[arcIdx[j]].vecPoints.end(); vIter++)
            for (std::vector<Vector3>::iterator vIter = offsetPathArcOnMesh[arcIdx[j]].vecPoints.begin();
                 vIter != offsetPathArcOnMesh[arcIdx[j]].vecPoints.end(); vIter++) {
                tempHeight = (*vIter) * heightDirection;
                if (tempHeight <= upBound)
                    distDistruction.insert(tempHeight);
            }
            //
            maxInterval[j][0] = 1.0; // min
            maxInterval[j][1] = 0.0; // max
            //
            std::set<double>::iterator sIter = distDistruction.begin();
            std::set<double>::iterator prevIter = distDistruction.begin();
            //for (sIter++; sIter != distDistruction.end(); sIter++)
            //{
            //if (*sIter - *prevIter > maxInterval[j][1] - maxInterval[j][0])
            //{
            sIter++;
            maxInterval[j][0] = *prevIter;
            maxInterval[j][1] = *sIter;
            //}
            //prevIter = sIter;
            //}
            //
            distDistruction.clear();
        }// for j == 2
        if (!special_flag_for_arcs[arcIdx[0]] && !special_flag_for_arcs[arcIdx[1]]) {
            if (maxInterval[0][1] - maxInterval[0][0] > maxInterval[1][1] - maxInterval[1][0]) {
                levelSetValues[i] = (maxInterval[0][1] + maxInterval[0][0]) * 0.5;
                initSimplifiedArcIdForHorizontalLoops[i] = arcIdx[0];
            } else {
                levelSetValues[i] = (maxInterval[1][1] + maxInterval[1][0]) * 0.5;
                initSimplifiedArcIdForHorizontalLoops[i] = arcIdx[1];
            }
            if (abs(levelSetValues[i] -
                    MeshVertexScalarValue[(*pVecReebNode)[(*criticalPairing)[i].second].nVertexId]) < 1e-5) {
                std::cout << "level set value is very close to upforking node " << i << std::endl;
            }
        } else {
            if (!special_flag_for_arcs[arcIdx[0]]) {
                levelSetValues[i] = (maxInterval[0][1] + maxInterval[0][0]) * 0.5;
                initSimplifiedArcIdForHorizontalLoops[i] = arcIdx[0];
            } else {
                if (!special_flag_for_arcs[arcIdx[1]]) {
                    levelSetValues[i] = (maxInterval[1][1] + maxInterval[1][0]) * 0.5;
                    initSimplifiedArcIdForHorizontalLoops[i] = arcIdx[1];
                } else {
                    std::cout << "---------Height function is not good-------- " << std::endl;
                    if (maxInterval[0][1] - maxInterval[0][0] > maxInterval[1][1] - maxInterval[1][0]) {
                        levelSetValues[i] = (maxInterval[0][1] + maxInterval[0][0]) * 0.5;
                        initSimplifiedArcIdForHorizontalLoops[i] = arcIdx[0];
                    } else {
                        levelSetValues[i] = (maxInterval[1][1] + maxInterval[1][0]) * 0.5;
                        initSimplifiedArcIdForHorizontalLoops[i] = arcIdx[1];
                    }
                    if (abs(levelSetValues[i] -
                            MeshVertexScalarValue[(*pVecReebNode)[(*criticalPairing)[i].second].nVertexId]) < 1e-5) {
                        std::cout << "level set value is very close to upforking node " << i << std::endl;
                    }
                }
            }
        }
    }
    //compute level sets
    for (unsigned int i = 0; i < initSimplifiedArcIdForHorizontalLoops.size(); i++) {
        const int arcId = initSimplifiedArcIdForHorizontalLoops[i];
        //
        CutLevelSet_specific_value_pair(levelsetOnMesh[i],
                                        levelsetOnMeshPointType[i],
                                        (*pVecSimplifiedArc)[arcId],
                                        levelSetValues[i]);
    }
    //
    sepPathArcOnMesh.resize(pVecSimplifiedArc->size());
    sepOffsetPathArcOnMesh.resize(pVecSimplifiedArc->size());
    //
    for (unsigned int i = 0; i < pathArcOnMesh.size(); i++) {
        sepPathArcOnMesh[i].vecPoints.assign(pathArcOnMesh[i].vecPoints.begin(), pathArcOnMesh[i].vecPoints.end());
        sepOffsetPathArcOnMesh[i].vecPoints.assign(offsetPathArcOnMesh[i].vecPoints.begin(),
                                                   offsetPathArcOnMesh[i].vecPoints.end());
    }
    for (int i = 0; i < edgeIntersectArcs.size(); i++) {
        if (edgeIntersectArcs[i].size() > 0) {
            std::vector<Vector3> pts_on_edge(2 * edgeIntersectArcs[i].size());
            //
            std::set<std::pair<double, int>, myDoubleIntPairCompare> distSorted;
            //
            Vector3 ab_a((*meshData).vecVertex[(*meshData).vecEdge[i].v0].x,
                         (*meshData).vecVertex[(*meshData).vecEdge[i].v0].y,
                         (*meshData).vecVertex[(*meshData).vecEdge[i].v0].z);
            Vector3 ab_b((*meshData).vecVertex[(*meshData).vecEdge[i].v1].x,
                         (*meshData).vecVertex[(*meshData).vecEdge[i].v1].y,
                         (*meshData).vecVertex[(*meshData).vecEdge[i].v1].z);
            //
            int edgeIter = 0;
            //
            for (std::vector<std::pair<int, int> >::iterator vIter = edgeIntersectArcs[i].begin();
                 vIter != edgeIntersectArcs[i].end();
                 vIter++) {
                pts_on_edge[edgeIter] = pathArcOnMesh[vIter->first].vecPoints[vIter->second];
                distSorted.insert(std::pair<double, int>(norm2(pts_on_edge[edgeIter] - ab_a), edgeIter));
                edgeIter++;
                pts_on_edge[edgeIter] = offsetPathArcOnMesh[vIter->first].vecPoints[vIter->second];
                distSorted.insert(std::pair<double, int>(norm2(pts_on_edge[edgeIter] - ab_a), edgeIter));
                edgeIter++;
            }
            //
            const double inv_interval = 1.0 / (pts_on_edge.size() + 1);
            double cur_weight = inv_interval;
            for (std::set<std::pair<double, int>, myDoubleIntPairCompare>::iterator sIter = distSorted.begin();
                 sIter != distSorted.end(); sIter++) {
                pts_on_edge[sIter->second] = (1.0 - cur_weight) * ab_a + cur_weight * ab_b;
                cur_weight = cur_weight + inv_interval;
            }
            //
            edgeIter = 0;
            for (std::vector<std::pair<int, int> >::iterator vIter = edgeIntersectArcs[i].begin();
                 vIter != edgeIntersectArcs[i].end();
                 vIter++) {
                sepPathArcOnMesh[vIter->first].vecPoints[vIter->second] = pts_on_edge[edgeIter++];
                sepOffsetPathArcOnMesh[vIter->first].vecPoints[vIter->second] = pts_on_edge[edgeIter++];
            }
            //
        }
    }
    return;
}

void psbmReebGraph::SetSimplifiedArcFlagBits() {
    // set the flag bits for vertical loops
    arcIntersectWithVerticalLoops.resize(pVecSimplifiedArc->size());
    for (unsigned int i = 0; i < initVerticalLoops->size(); i++) {
        for (unsigned int k = 0; k < (*initVerticalLoops)[i].simpArcIndexSet.size(); k++) {
            arcIntersectWithVerticalLoops[(*initVerticalLoops)[i].simpArcIndexSet[k]].push_back(i);
        }
    }
    // set the flag bits for horizontal loops
    //arcIntersectWithHorizontalLoops.resize(pVecSimplifiedArc->size(), -1);/*
    //for (unsigned int i = 0; i < arcIntersectWithHorizontalLoops.size(); i++)
    //	arcIntersectWithHorizontalLoops.push_back(-1);*/
    //for (unsigned int i = 0; i < initSimplifiedArcIdForHorizontalLoops.size(); i++)
    //	arcIntersectWithHorizontalLoops[initSimplifiedArcIdForHorizontalLoops[i]] = i;
    //
    return;
}

// compute the intersections between two loops
void psbmReebGraph::ComputeIntersections(_Polygon &loop_a,
                                         std::vector<std::pair<int, int> > &point_type_a,
                                         _Polygon &loop_b,
                                         std::vector<std::pair<int, int> > &point_type_b) {
    //
    for (unsigned int i = 0; i < point_type_a.size() - 1; i++) {// [i, i+1] is a segment
        for (unsigned int j = 0;
             j < point_type_b.size() - 1; j++) {// check the intersection between [i, i+1] and [j, j+1]

        }
    }
    return;
}

void psbmReebGraph::LinkNumberMatrixComputing_withLoopsOnFly() {
    std::vector<std::vector<std::pair<int, int> > > myLinkNumbersMatrix;
    LinkNumberPairPolygon linkComputer;
    //
    //
    NonOverlappingCycles ComputePairLines;
    NonOverlappedLevelCycleAndArc levelCycleAndArc;
    //
    // reserve the memory for the sparse matrix
    myLinkNumbersMatrix.resize(initVerticalLoops->size() * 2);
    //
    // compute the link number between the initial vertical loops
    // initialize the common data
    ComputePairLines.InitMeshPtr(meshData);
    //
    ComputePairLines.InitNormalInfo(criticalNodeIdSet,
                                    pVecReebNode);
    //
    ComputePairLines.InitArcPolygons(pVecSimplifiedArc,
                                     &pathArcOnMesh,
                                     &offsetPathArcOnMesh,
                                     &pathArcOnMeshPointType,
                                     &offsetPathArcOnMeshPointType);
    //
    ComputePairLines.InitPairPathes(&(initVerticalLoops->front()),
                                    &((*(initVerticalLoops))[0]));
    //
    ComputePairLines.ComputeOutNormals();
    //
    //

    //
    for (unsigned int j = 0; j < initVerticalLoops->size(); j++) {
        ComputePairLines.ResetLoopOnMesh(&((*(initVerticalLoops))[j]));

        //
        std::set<int> intersect_loops;
        for (unsigned int nArc = 0; nArc < (*(initVerticalLoops))[j].simpArcIndexSet.size(); nArc++) {
            std::copy(arcIntersectWithVerticalLoops[(*(initVerticalLoops))[j].simpArcIndexSet[nArc]].begin(),
                      arcIntersectWithVerticalLoops[(*(initVerticalLoops))[j].simpArcIndexSet[nArc]].end(),
                      std::inserter(intersect_loops, intersect_loops.end()));
        }
        //
        std::cout << " >>>>>>>>>>>>>" << std::endl;
        std::cout << " TYPE " << (*(initVerticalLoops))[j].pathType << std::endl;
        for (std::set<int>::iterator sIter = intersect_loops.begin();
             sIter != intersect_loops.end();
             sIter++) {
            ComputePairLines.ResetLoopOnMesh(&((*(initVerticalLoops))[j]));
            ComputePairLines.ResetBasisLoop(&((*(initVerticalLoops))[*sIter]));
            std::vector<int> sharedArcs;
            std::vector<_Polygon> sharedArcNewPt;
            std::vector<_Polygon> offsetSharedArcNewPt;
            //
            ComputePairLines.SharedArcs(sharedArcs);
            ComputePairLines.ComputeNewPathesForSharedArcs(sharedArcs,
                                                           sharedArcNewPt,
                                                           offsetSharedArcNewPt);
            ComputePairLines.ComputeBasisLoopPolygonOnFly(sharedArcs,
                                                          sharedArcNewPt,
                                                          offsetSharedArcNewPt);
            ComputePairLines.ComputeLoopOnMeshPolygonOnFly(sharedArcs,
                                                           sharedArcNewPt,
                                                           offsetSharedArcNewPt);
            //
            linkComputer.ComputeLinkNumberForTwoPolygons(ComputePairLines.outBasisLoop, ComputePairLines.outLoopOnMesh);
            if (linkComputer._link_number)
                myLinkNumbersMatrix[j].push_back(std::pair<int, int>(*sIter, linkComputer._link_number));
            //std::cout << *sIter << " " << std::endl;
        }
        //std::cout << "................................." << std::endl;
    }
    // compute the link number between the horizontal ones and vertical ones
    // link number between horizontal ones is always zero
    //
    // initialize the data for arc and level
    // common data for all computation
    levelCycleAndArc.InitMeshPtr(meshData);
    levelCycleAndArc.InitHeightDirection(heightDirection);
    // this is needed to be reset twice for each vertical path
    levelCycleAndArc.InitBool(false);
    levelCycleAndArc.InitArcPolygons(pVecSimplifiedArc,
                                     &pathArcOnMesh);//.pathArcOnMesh);

    //
    const int tot_genus = initVerticalLoops->size();
    for (unsigned int i = 0; i < initSimplifiedArcIdForHorizontalLoops.size(); i++) {
        if (i == 5)
            i = 5;
        //std::cout << "LEVEL SET = " << i << std::endl;
        // common data for specific arc and level set
        levelCycleAndArc.InitLevelCycleAndArc(&(levelsetOnMesh[i]),
                                              initSimplifiedArcIdForHorizontalLoops[i]);
        levelCycleAndArc.InitArcLevelPointType(&(levelsetOnMeshPointType[i]),
                                               &(pathArcOnMeshPointType[initSimplifiedArcIdForHorizontalLoops[i]]));
        //&(reebGraph.pathArcOnMeshPointType[reebGraph.initSimplifiedArcIdForHorizontalLoops[pid]]));
        // it is ready to compute the augmented pathes
        levelCycleAndArc.ComputeIntersectionBetweenLevelCycleAndArc();

        for (unsigned int pidx = 0; pidx < initVerticalLoops->size(); pidx++) {
            int vert_path_id = pidx;//arcIntersectWithVerticalLoops[initSimplifiedArcIdForHorizontalLoops[i]][pidx];
            // for each path passing this arc
            //particular vertical path intersects this level set
            levelCycleAndArc.InitVerticalPath(&((*(initVerticalLoops))[vert_path_id]));
            //
            // the vertical loop is on mesh, the level set is pulled
            levelCycleAndArc.ResetBool(true);
            //
            levelCycleAndArc.ComputeOutputCycles();
            //
            linkComputer.ComputeLinkNumberForTwoPolygons(levelCycleAndArc.outLevelCycle,
                                                         levelCycleAndArc.outVerticalCycle);
            //
            if (linkComputer._link_number) {
                myLinkNumbersMatrix[vert_path_id].push_back(
                        std::pair<int, int>(i + tot_genus, linkComputer._link_number));
            }
            //
            // the level set is on the mesh, the vertical loop is pulled
            levelCycleAndArc.ResetBool(false);
            //
            levelCycleAndArc.ComputeOutputCycles();
            //
            linkComputer.ComputeLinkNumberForTwoPolygons(levelCycleAndArc.outLevelCycle,
                                                         levelCycleAndArc.outVerticalCycle);
            //
            if (linkComputer._link_number) {
                myLinkNumbersMatrix[i + tot_genus].push_back(
                        std::pair<int, int>(vert_path_id, linkComputer._link_number));
            }
        }
//
    }
    // report the link number matrix
    std::cout << "M2=zeros(" << tot_genus * 2 << "," << tot_genus * 2 << ");" << std::endl;
    for (int i = 0; i < tot_genus; i++) {
        //std::cout << "Vert Loop " << i << " link against " << std::endl;
        for (unsigned int j = 0; j < myLinkNumbersMatrix[i].size(); j++) {
            std::cout << "M(" << i + 1 << "," << myLinkNumbersMatrix[i][j].first + 1 << ")="
                      << myLinkNumbersMatrix[i][j].second << ";";

        }
        std::cout << std::endl;
    }
    for (int i = 0; i < tot_genus; i++) {
        //std::cout << "Hori Loop " << i + tot_genus << " link against " << std::endl;
        for (unsigned int j = 0; j < myLinkNumbersMatrix[i + tot_genus].size(); j++) {
            std::cout << "M(" << i + tot_genus + 1 << "," << myLinkNumbersMatrix[i + tot_genus][j].first + 1 << ")="
                      << myLinkNumbersMatrix[i + tot_genus][j].second << ";";

        }
        std::cout << std::endl;
    }

    //
    return;
}

void psbmReebGraph::FindPointOnLevelSetIntersectingArc(const std::vector<std::pair<int, int> > &ptType,
                                                       const _Polygon &ptsOnPath,
                                                       const double levelSetValue,
                                                       Vector3 &pt_on_levelset) {
    //preq: levelSetValue is between lowest node and highest node of the arc
    double current_height = 0.0;
    int curVerOnMesh = 0;
    bool bFind = false;
    int ptIter = ptsOnPath.vecPoints.size() - 1;
    curVerOnMesh = ptType.back().first;

    if (MeshVertexScalarValue[curVerOnMesh] == levelSetValue) {
        bFind = true;
        pt_on_levelset[0] = (*meshData).vecVertex[curVerOnMesh].x;
        pt_on_levelset[1] = (*meshData).vecVertex[curVerOnMesh].y;
        pt_on_levelset[2] = (*meshData).vecVertex[curVerOnMesh].z;
    } else {
        curVerOnMesh = ptType.back().first;
        ptIter = ptsOnPath.vecPoints.size() - 2;
        current_height = ptsOnPath.vecPoints[ptIter] * heightDirection;
        if (current_height == levelSetValue) {
            bFind = true;
            pt_on_levelset = ptsOnPath.vecPoints[ptIter];
        } else {
            if (current_height > levelSetValue) {
                bFind = true;
                double weight = (levelSetValue - MeshVertexScalarValue[curVerOnMesh]) /
                                (current_height - MeshVertexScalarValue[curVerOnMesh]);
                pt_on_levelset =
                        (1.0 - weight) * ptsOnPath.vecPoints[ptIter + 1] + weight * ptsOnPath.vecPoints[ptIter];
            }
        }
    }
    if (!bFind) {
        // iterate all edges
        ptIter = ptsOnPath.vecPoints.size() - 2;
        for (; ptIter > 0; ptIter--) {//check the height of this edge containing the levelset value or not
            const int eid = ptType[ptIter].first;
            const int v0 = (*meshData).vecEdge[eid].v0;
            const int v1 = (*meshData).vecEdge[eid].v1;
            //
            if ((MeshVertexScalarValue[v0] >= levelSetValue &&
                 MeshVertexScalarValue[v1] <= levelSetValue) ||
                (MeshVertexScalarValue[v1] >= levelSetValue &&
                 MeshVertexScalarValue[v0] <= levelSetValue)
                    ) {
                Vector3 pt_a((*meshData).vecVertex[v0].x,
                             (*meshData).vecVertex[v0].y,
                             (*meshData).vecVertex[v0].z);
                //
                Vector3 pt_b((*meshData).vecVertex[v1].x,
                             (*meshData).vecVertex[v1].y,
                             (*meshData).vecVertex[v1].z);
                //
                bFind = true;
                //
                if (MeshVertexScalarValue[v0] == levelSetValue || MeshVertexScalarValue[v1] == levelSetValue) {
                    if (MeshVertexScalarValue[v0] == levelSetValue)
                        pt_on_levelset = pt_a;
                    else
                        pt_on_levelset = pt_b;
                } else {
                    double weight = (levelSetValue - MeshVertexScalarValue[v0]) /
                                    (MeshVertexScalarValue[v1] - MeshVertexScalarValue[v0]);
                    pt_on_levelset = (1.0 - weight) * pt_a + weight * pt_b;
                }
                break;
            }
        } // for ptIter

    }
    if (!bFind) {
        curVerOnMesh = ptType.front().first;
        ptIter = 1;
        current_height = ptsOnPath.vecPoints[ptIter] * heightDirection;
        if (MeshVertexScalarValue[curVerOnMesh] == levelSetValue) {
            bFind = true;
            pt_on_levelset[0] = (*meshData).vecVertex[curVerOnMesh].x;
            pt_on_levelset[1] = (*meshData).vecVertex[curVerOnMesh].y;
            pt_on_levelset[2] = (*meshData).vecVertex[curVerOnMesh].z;
        } else {
            if (current_height == levelSetValue) {
                bFind = true;
                pt_on_levelset = ptsOnPath.vecPoints[ptIter];
            } else {
                if (current_height < levelSetValue) {
                    bFind = true;
                    double weight = (levelSetValue - MeshVertexScalarValue[curVerOnMesh]) /
                                    (current_height - MeshVertexScalarValue[curVerOnMesh]);
                    pt_on_levelset =
                            (1.0 - weight) * ptsOnPath.vecPoints[ptIter + 1] + weight * ptsOnPath.vecPoints[ptIter];
                }
            }
        }
    }
    if (!bFind) {
        std::cout << "can not find the point in this arc path" << std::endl;
        exit(0);
    }
    //class psbmReebArc* travArc = NULL;
    //int nodeId = travArc->nNodeId1; // lareger index
    //int valenceDown = (*pVecReebNode)[nodeId].ptrArcDownId->size() ;
    //int valenceUp =  (*pVecReebNode)[nodeId].ptrArcUpId->size();
    ////
    //for (std::list<class psbmReebArc*>::iterator listIter = (*pVecReebNode)[(*pVecSimplifiedArc)[arc_idx].nCriticalNode0].ptrArcUpId->begin();
    //	listIter != (*pVecReebNode)[(*pVecSimplifiedArc)[arc_idx].nCriticalNode0].ptrArcUpId->end(); listIter++)
    //{
    //	if ((*listIter)->clusterLabel == arc_idx)
    //	{
    //		travArc = *listIter;
    //		break;
    //	}
    //}
    ////
    //int curVerOnMesh = (*pVecReebNode)[(*pVecSimplifiedArc)[arc_idx].nCriticalNode0].nVertexId;
    //if (MeshVertexScalarValue[curVerOnMesh] == levelSetValue)
    //{
    //	pt_on_levelset[0] = (*meshData).vecVertex[curVerOnMesh].x;
    //	pt_on_levelset[1] = (*meshData).vecVertex[curVerOnMesh].y;
    //	pt_on_levelset[2] = (*meshData).vecVertex[curVerOnMesh].z;
    //}
    //else
    //{
    //	curVerOnMesh =  (*pVecReebNode)[travArc->nNodeId1].nVertexId;
    //	while (MeshVertexScalarValue[curVerOnMesh] < levelSetValue)
    //	{
    //		travArc = (*pVecReebNode)[travArc->nNodeId1].ptrArcUpId->front();
    //		curVerOnMesh =  (*pVecReebNode)[travArc->nNodeId1].nVertexId;
    //	}
    //	//
    //	if (MeshVertexScalarValue[curVerOnMesh] == levelSetValue)
    //	{
    //		pt_on_levelset[0] = (*meshData).vecVertex[curVerOnMesh].x;
    //		pt_on_levelset[1] = (*meshData).vecVertex[curVerOnMesh].y;
    //		pt_on_levelset[2] = (*meshData).vecVertex[curVerOnMesh].z;
    //	}
    //	else
    //	{// > levelSetValue
    //	 // iterate its edge
    //	 // bug: need to make sure this edge
    //	}
    //}
    //

}

void psbmReebGraph::LinkNumberMatrixComputing() {
    LinkNumberPairPolygon linkComputer;
    //
    //
    NonOverlappingCycles ComputePairLines;
    NonOverlappedLevelCycleAndArc levelCycleAndArc;
    //
    // reserve the memory for the sparse matrix
    linkNumbersMatrix.resize(initVerticalLoops->size() * 2);
    //
    // compute the link number between the initial vertical loops
    // initialize the common data
    ComputePairLines.InitMeshPtr(meshData);
    //
    ComputePairLines.InitNormalInfo(criticalNodeIdSet,
                                    pVecReebNode);
    //
    //ComputePairLines.InitArcPolygons(pVecSimplifiedArc,
    //	&sepOffsetPathArcOnMesh, &sepPathArcOnMesh);
//
    ComputePairLines.SetEdgeIntersectArcsPtr(&edgeIntersectArcs);
    ComputePairLines.InitArcPolygons(pVecSimplifiedArc,
                                     &sepPathArcOnMesh, &sepOffsetPathArcOnMesh,
                                     &pathArcOnMeshPointType, &offsetPathArcOnMeshPointType);
    //&sepPathArcOnMesh, // good for linking number
    //&sepOffsetPathArcOnMesh); // good for linking number
    //
    ComputePairLines.InitPairPathes(&(initVerticalLoops->front()),
                                    &((*(initVerticalLoops))[0]));
    ComputePairLines.SetHeightDiretion(heightDirection);
    //
    //ComputePairLines.ComputeOutNormals();
    //
    std::vector<std::pair<double, double> > arcHeightInSepPath(sepPathArcOnMesh.size());
    std::vector<std::pair<double, double> > arcHeightInOffsetSepPath(sepPathArcOnMesh.size());
    for (unsigned int i = 0; i < sepPathArcOnMesh.size(); i++) {
        double curHeight = 0.0;
        if (!sepPathArcOnMesh[i].vecPoints.empty()) {
            arcHeightInSepPath[i].first = sepPathArcOnMesh[i].vecPoints.front() * heightDirection;
            arcHeightInSepPath[i].second = arcHeightInSepPath[i].first;
            //
            arcHeightInOffsetSepPath[i].first = sepOffsetPathArcOnMesh[i].vecPoints.front() * heightDirection;
            arcHeightInOffsetSepPath[i].second = arcHeightInOffsetSepPath[i].first;
            //
            for (unsigned int j = 0; j < sepPathArcOnMesh[i].vecPoints.size(); j++) {
                curHeight = sepPathArcOnMesh[i].vecPoints[j] * heightDirection;
                if (arcHeightInSepPath[i].first > curHeight)
                    arcHeightInSepPath[i].first = curHeight;
                if (arcHeightInSepPath[i].second < curHeight)
                    arcHeightInSepPath[i].second = curHeight;
                //
                curHeight = sepOffsetPathArcOnMesh[i].vecPoints[j] * heightDirection;
                if (arcHeightInOffsetSepPath[i].first > curHeight)
                    arcHeightInOffsetSepPath[i].first = curHeight;
                if (arcHeightInOffsetSepPath[i].second < curHeight)
                    arcHeightInOffsetSepPath[i].second = curHeight;
            }
        }
    }
    double heightError = 5.0;
    for (unsigned int j = 0; j < initVerticalLoops->size(); j++) {
        std::vector<bool> bArcLoopOnMeshStraight((*initVerticalLoops)[j].simpArcIndexSet.size(), false);
        //
        ComputePairLines.ResetLoopOnMesh(&((*initVerticalLoops)[j]));

        //std::cout << j << " has TYPE " << (*initVerticalLoops)[j].pathType << "\t" ;//std::endl;
        //
        std::pair<double, double> actualMinMax(std::numeric_limits<double>::max(), std::numeric_limits<double>::min());
        for (std::vector<int>::iterator vIter = (*initVerticalLoops)[j].simpArcIndexSet.begin();
             vIter != (*initVerticalLoops)[j].simpArcIndexSet.end();
             vIter++) {
            if (actualMinMax.first > arcHeightInOffsetSepPath[*vIter].first)
                actualMinMax.first = arcHeightInOffsetSepPath[*vIter].first;
            if (actualMinMax.second < arcHeightInOffsetSepPath[*vIter].second)
                actualMinMax.second = arcHeightInOffsetSepPath[*vIter].second;
        }
        //
        for (unsigned int k = 0; k < initVerticalLoops->size(); k++) {
            if (!(MeshVertexScalarValue[(*pVecReebNode)[(*criticalPairing)[k].second].nVertexId] >
                  MeshVertexScalarValue[(*pVecReebNode)[(*criticalPairing)[j].first].nVertexId] ||
                  MeshVertexScalarValue[(*pVecReebNode)[(*criticalPairing)[k].first].nVertexId] <
                  MeshVertexScalarValue[(*pVecReebNode)[(*criticalPairing)[j].second].nVertexId])) {// ignore the case "==", because when the heights are equal, they share some critical point;
                // if this critical point is needed to be moved, the simplified straight line could affect the
                // link number computation;
                for (unsigned int arcIter = 0; arcIter < bArcLoopOnMeshStraight.size(); arcIter++)
                    bArcLoopOnMeshStraight[arcIter] = false;
                ////
                std::vector<bool> bArcBasisLoopStraight((*initVerticalLoops)[k].simpArcIndexSet.size(), false);
                ////check if there is any arc can be used the straight segment connecting two critical nodes
                if (j != k) {// if j==k, both arc sets are the same
                    //
                    int arcIter = 0;
                    //
                    std::pair<int, int> curMinMax(std::numeric_limits<double>::max(),
                                                  std::numeric_limits<double>::min());
                    //
                    for (std::vector<int>::iterator vIter = (*initVerticalLoops)[k].simpArcIndexSet.begin();
                         vIter != (*initVerticalLoops)[k].simpArcIndexSet.end();
                         vIter++) {
                        if (curMinMax.first > arcHeightInSepPath[*vIter].first)
                            curMinMax.first = arcHeightInSepPath[*vIter].first;
                        if (curMinMax.second < arcHeightInSepPath[*vIter].second)
                            curMinMax.second = arcHeightInSepPath[*vIter].second;
                        //
                        if (arcHeightInSepPath[*vIter].first > actualMinMax.second + heightError ||
                            arcHeightInSepPath[*vIter].second + heightError < actualMinMax.first) {
                            bArcBasisLoopStraight[arcIter] = true;
                        }
                        arcIter++;
                    }
                    //
                    arcIter = 0;
                    for (std::vector<int>::iterator vIter = (*initVerticalLoops)[j].simpArcIndexSet.begin();
                         vIter != (*initVerticalLoops)[j].simpArcIndexSet.end();
                         vIter++) {
                        if (arcHeightInOffsetSepPath[*vIter].first > curMinMax.second + heightError ||
                            arcHeightInOffsetSepPath[*vIter].second + heightError < curMinMax.first) {
                            bArcLoopOnMeshStraight[arcIter] = true;
                        }
                        arcIter++;
                    }
                    //

                }
                //
                ComputePairLines.SetBoolVectorPtr(&bArcLoopOnMeshStraight, &bArcBasisLoopStraight);
                //
                ComputePairLines.ResetLoopOnMesh(&((*initVerticalLoops)[j]));
                ComputePairLines.ComputeLoopOnMeshPolygon();
                //
                ComputePairLines.ResetBasisLoop(&((*initVerticalLoops)[k]));
                //std::cout << "in ply_1" << std::endl;
                ComputePairLines.ComputeBasisLoopPolygon_1();
                //

                if (!linkComputer.ComputeLinkNumberForTwoPolygons(ComputePairLines.outLoopOnMesh,
                                                                  ComputePairLines.outBasisLoop)) {// Failure in computing linking number
                    std::cout << j << " and " << k << std::endl;
                    std::cout << (*(initVerticalLoops))[j].nodeIndexSet.front() << std::endl;
                    std::cout << (*(initVerticalLoops))[k].nodeIndexSet.front() << std::endl;
                    std::cout << ComputePairLines.outLoopOnMesh.vecPoints.size() << std::endl;
                    //for (unsigned int iv = 0; iv < ComputePairLines.outLoopOnMesh.vecPoints.size(); iv++)
                    //{
                    //	std::cout << ComputePairLines.outLoopOnMesh.vecPoints[iv][0] << " ";
                    //	std::cout << ComputePairLines.outLoopOnMesh.vecPoints[iv][1] << " ";
                    //	std::cout << ComputePairLines.outLoopOnMesh.vecPoints[iv][2] << std::endl;
                    //}
                    std::cout << std::endl << ComputePairLines.outBasisLoop.vecPoints.size() << std::endl;
                    //for (unsigned int iv = 0; iv < ComputePairLines.outBasisLoop.vecPoints.size(); iv++)
                    //{
                    //	std::cout << ComputePairLines.outBasisLoop.vecPoints[iv][0] <<  " ";
                    //	std::cout << ComputePairLines.outBasisLoop.vecPoints[iv][1] << " ";
                    //	std::cout << ComputePairLines.outBasisLoop.vecPoints[iv][2] << std::endl;
                    //}
                    exit(0);
                }
                if (linkComputer._link_number)
                    linkNumbersMatrix[j].push_back(std::pair<int, int>(k, linkComputer._link_number));
            }
        }

        //std::cout << "................................." << std::endl;
    }

    // compute the link number between the horizontal ones and vertical ones
    // link number between horizontal ones is always zero
    //
    // initialize the data for arc and level
    // common data for all computation
    levelCycleAndArc.InitMeshPtr(meshData);
    levelCycleAndArc.InitHeightDirection(heightDirection);
    // this is needed to be reset twice for each vertical path
    levelCycleAndArc.InitBool(false);
    //levelCycleAndArc.InitArcPolygons(pVecSimplifiedArc,
    //								&pathArcOnMesh);//.pathArcOnMesh);
    levelCycleAndArc.InitArcPolygons(pVecSimplifiedArc,
                                     &offsetPathArcOnMesh);//.pathArcOnMesh);

    //
    const int tot_genus = initVerticalLoops->size();
    for (unsigned int i = 0; i < initSimplifiedArcIdForHorizontalLoops.size(); i++) {
        //std::cout << "LEVEL SET = " << i << "\t" ;//std::endl;

        // common data for specific arc and level set
        levelCycleAndArc.InitLevelCycleAndArc(&(levelsetOnMesh[i]),
                                              initSimplifiedArcIdForHorizontalLoops[i]);
        //levelCycleAndArc.InitArcLevelPointType(&(levelsetOnMeshPointType[i]),
        //	&(pathArcOnMeshPointType[initSimplifiedArcIdForHorizontalLoops[i]]));
        levelCycleAndArc.InitArcLevelPointType(&(levelsetOnMeshPointType[i]),
                                               &(offsetPathArcOnMeshPointType[initSimplifiedArcIdForHorizontalLoops[i]]));
        //&(reebGraph.pathArcOnMeshPointType[reebGraph.initSimplifiedArcIdForHorizontalLoops[pid]]));
        //
        if ((*initVerticalLoops)[i].pathType)
            levelCycleAndArc.SetOwnVerticalType(true);
        else
            levelCycleAndArc.SetOwnVerticalType(false);
        // it is ready to compute the augmented pathes
        levelCycleAndArc.ComputeIntersectionBetweenLevelCycleAndArc();
        //
        const double levelSetValue = levelsetOnMesh[i].vecPoints.front() * heightDirection;
        const int levelArcIdx = initSimplifiedArcIdForHorizontalLoops[i];
        levelCycleAndArc._arcId = levelArcIdx;
        //
        for (unsigned int pidx = 0; pidx < i + 1; pidx++)//initVerticalLoops->size(); pidx++)
        {
            int vert_path_id = pidx;//arcIntersectWithVerticalLoops[initSimplifiedArcIdForHorizontalLoops[i]][pidx];

            if (MeshVertexScalarValue[(*pVecReebNode)[(*criticalPairing)[vert_path_id].second].nVertexId] <=
                levelSetValue &&
                MeshVertexScalarValue[(*pVecReebNode)[(*criticalPairing)[vert_path_id].first].nVertexId] >=
                levelSetValue) {
                bool bHaveLevelArc = false;
                int arc_inside_levelset = 0;
                std::vector<bool> arcHeightRangeIntersectingLevelset(
                        (*initVerticalLoops)[vert_path_id].simpArcIndexSet.size(), false);
                //
                int arcIter_order = 0;
                for (std::vector<int>::iterator vIter = (*initVerticalLoops)[vert_path_id].simpArcIndexSet.begin();
                     vIter != (*initVerticalLoops)[vert_path_id].simpArcIndexSet.end();
                     vIter++) {
                    //arcHeightRangeIntersectingLevelset[arcIter_order] = true;
                    if (*vIter == levelArcIdx) {
                        bHaveLevelArc = true;
                        //
                        arcHeightRangeIntersectingLevelset[arcIter_order] = true;
                    } else {//
                        if (MeshVertexScalarValue[(*pVecReebNode)[(*pVecSimplifiedArc)[*vIter].nCriticalNode0].nVertexId] <=
                            levelSetValue &&
                            MeshVertexScalarValue[(*pVecReebNode)[(*pVecSimplifiedArc)[*vIter].nCriticalNode1].nVertexId] >=
                            levelSetValue) {// this arc range contains the levelSetValue
                            //
                            arcHeightRangeIntersectingLevelset[arcIter_order] = true;
                            //
                            // find a point on this arc which is also on the level set
                            Vector3 pt_on_level_set;
                            FindPointOnLevelSetIntersectingArc(pathArcOnMeshPointType[*vIter], pathArcOnMesh[*vIter],
                                                               levelSetValue, pt_on_level_set);
                            // determine this point is inside or outside of the level set cycle
                            // record it when it is inside
                            _Polygon projPlyg = levelsetOnMesh[i];
                            for (unsigned int vid = 0; vid < projPlyg.vecPoints.size(); vid++) {//
                                RotateAxisTobeAlignedWithZ_axis(heightDirection[0], heightDirection[1],
                                                                heightDirection[2],
                                                                projPlyg.vecPoints[vid][0], projPlyg.vecPoints[vid][1],
                                                                projPlyg.vecPoints[vid][2],
                                                                projPlyg.vecPoints[vid]);
                            }
                            RotateAxisTobeAlignedWithZ_axis(heightDirection[0], heightDirection[1], heightDirection[2],
                                                            pt_on_level_set[0], pt_on_level_set[1], pt_on_level_set[2],
                                                            pt_on_level_set);
                            if (InsidePolygon(projPlyg, pt_on_level_set, 1,
                                              0)) { // pt is inside of plyg, they are nested
                                arc_inside_levelset++;
                            }
                            //arc_inside_levelset++;
                        }
                    }
                    arcIter_order++;
                }
                if (bHaveLevelArc ||
                    arc_inside_levelset) {// ignore the case : 1) it doesn't contain the arc which intersects the level set
                    //					  2) no arc passes through the level set from its inside
                    //                  under these two conditions, level set and the vertical cycle are unlinked;
                    // for each path passing this arc
                    //particular vertical path intersects this level set
                    levelCycleAndArc.InitVerticalPath(&((*(initVerticalLoops))[vert_path_id]));
                    //
                    // the vertical loop is on mesh, the level set is pulled
                    levelCycleAndArc.SetArcHeightRangeIntersectingLevelsetPtr(&arcHeightRangeIntersectingLevelset);
                    levelCycleAndArc.ResetBool(true);
                    //
                    levelCycleAndArc.ComputeOutputCycles();
                    //
                    linkComputer.ComputeLinkNumberForTwoPolygons(levelCycleAndArc.outLevelCycle,
                                                                 levelCycleAndArc.outVerticalCycle);
                    //
                    if (linkComputer._link_number) {
                        linkNumbersMatrix[vert_path_id].push_back(
                                std::pair<int, int>(i + tot_genus, linkComputer._link_number));
                    }
                    //
                    // the level set is on the mesh, the vertical loop is pulled
                    levelCycleAndArc.ResetBool(false);
                    //
                    levelCycleAndArc.ComputeOutputCycles();
                    //
                    linkComputer.ComputeLinkNumberForTwoPolygons(levelCycleAndArc.outVerticalCycle,
                                                                 levelCycleAndArc.outLevelCycle);
                    //
                    if (linkComputer._link_number) {
                        linkNumbersMatrix[i + tot_genus].push_back(
                                std::pair<int, int>(vert_path_id, linkComputer._link_number));
                    }
                } // need to compute it
            }// cycle height range containg level set value
        }

    }
    // report the link number matrix
    InverseLinkNumberMatrix invComputing;
    //invComputing.PrintMatrix("orgMat.txt", linkNumbersMatrix);
    invComputing.SetMatrixHalfDim(initVerticalLoops->size());
    invComputing.SetOrgMatrix(linkNumbersMatrix);
    invComputing.ComputeInverseMatrix();
    //
    invLinkNumbersMatrix = invComputing.inv_mat;
    //invComputing.PrintMatrix("invMat.txt", invLinkNumbersMatrix);

    //std::cout << "M=zeros("<< tot_genus * 2 << "," << tot_genus * 2 << ");" << std::endl;
    //for (int i = 0; i < tot_genus; i++)
    //{
    //	//std::cout << "Vert Loop " << i << " link against " << std::endl;
    //	for (unsigned int j = 0; j < linkNumbersMatrix[i].size(); j++)
    //	{
    //		std::cout << "M(" << i + 1 << "," << linkNumbersMatrix[i][j].first + 1 << ")=" << linkNumbersMatrix[i][j].second << ";";
    //
    //	}
    //	std::cout << std::endl;
    //}
    //for (int i = 0; i < tot_genus; i++)
    //{
    //	//std::cout << "Hori Loop " << i + tot_genus << " link against " << std::endl;
    //	for (unsigned int j = 0; j < linkNumbersMatrix[i + tot_genus].size(); j++)
    //	{
    //		std::cout << "M(" << i + tot_genus + 1 << "," << linkNumbersMatrix[i + tot_genus][j].first + 1<< ")=" << linkNumbersMatrix[i + tot_genus][j].second << ";";
    //
    //	}
    //	std::cout << std::endl;
    //}

    //
    return;
}

void psbmReebGraph::EmbedCycleAsEdgePathOnMesh() {
    // the embedded edge path on mesh for each simplified arc
    std::vector<std::vector<int> > vecArcEdgePathOnMesh_Vertex;
    std::vector<std::vector<int> > vecArcEdgePathOnMesh_Edge;
    //
    MapLoopsBackToMeshLevelSetAndArc cycleMappedOnMesh;
    //
    const int tot_genus = initVerticalLoops->size();
    //
    vecArcEdgePathOnMesh_Vertex.resize(pathArcOnMesh.size());
    vecArcEdgePathOnMesh_Edge.resize(pathArcOnMesh.size());
    // allocate the memory first
    //
    vecVerticalLoopEdgePathOnMesh_Vertex.resize(tot_genus);
    vecVerticalLoopEdgePathOnMesh_Edge.resize(tot_genus);
    //
    vecHorizontalLoopEdgePathOnMesh_Vertex.resize(tot_genus);
    vecHorizontalLoopEdgePathOnMesh_Edge.resize(tot_genus);
    //

    // initialize the data for arc and level
    // common data for all computation
    cycleMappedOnMesh.InitMeshPtr(meshData);
    cycleMappedOnMesh.InitScalarPtr(MeshVertexScalarValue);

    //
    //std::cout << "  initSimplifiedArcIdForHorizontalLoops.size() "
    //	<<  initSimplifiedArcIdForHorizontalLoops.size() << std::endl;
    for (unsigned int i = 0; i < initSimplifiedArcIdForHorizontalLoops.size(); i++) {
        // common data for specific arc and level set
        cycleMappedOnMesh.InitLevelCyclePointType(&(levelsetOnMeshPointType[i]));
        //
        cycleMappedOnMesh.ComputeLevelCycle();
        // store the output
        vecHorizontalLoopEdgePathOnMesh_Vertex[i] = cycleMappedOnMesh.outLevelCycle_VertexOnMesh;
        vecHorizontalLoopEdgePathOnMesh_Edge[i] = cycleMappedOnMesh.outLevelCycle_EdgeOnMesh;
        if (!cycleMappedOnMesh.VerifyingEmbeddedPath(vecHorizontalLoopEdgePathOnMesh_Vertex[i],
                                                     vecHorizontalLoopEdgePathOnMesh_Edge[i]))
            //{
            //	std::cout << "VALID H cycle" <<  i << std::endl;
            //}
            //else
        {
            std::cout << "NON-VALID H cycle" << i << std::endl;
        }
        //else
        //{
        //	std::cout << "size is  " << vecHorizontalLoopEdgePathOnMesh_Edge[i].size() << std::endl;
        //}
    }
    //
    //for (unsigned int i = 0; i < pathArcOnMeshPointType.size(); i++)
    //{
    //	if (!pathArcOnMeshPointType[i].empty() )
    //	{
    //		cycleMappedOnMesh.InitArcCyclePointType( &pathArcOnMeshPointType[i]);
    //		cycleMappedOnMesh.ComputeArcCycle();
    //		// store bakc
    //		vecArcEdgePathOnMesh_Vertex[i] = cycleMappedOnMesh.outArcCycle_VertexOnMesh;
    //		vecArcEdgePathOnMesh_Edge[i] = cycleMappedOnMesh.outArcCycle_EdgeOnMesh;
    //		//
    //		if (!cycleMappedOnMesh.VerifyingEmbeddedPath(vecArcEdgePathOnMesh_Vertex[i] , vecArcEdgePathOnMesh_Edge[i]))
    //		//{
    //		//	std::cout << "VALID arc cycle" <<  i << std::endl;
    //		///}
    //		//else
    //		{
    //			std::cout << "NON-VALID arc cycle" <<  i << std::endl;
    //		}
    //	}
    //}
    // invariant :  current_active_vertex is an end point of the edge refEdge; and both of them on the triangle to decide
    // orientation
    //
    for (unsigned int i = 0; i < initVerticalLoops->size(); i++) {
        //// store the output back
        int curArcId = 0;
        int startNodeId = 0;
        int endNodeId = 0;
        //
        int refEdge = 0;
        int init_side_vertex = 0;
        int nTriangleIdForOrientation = 0;
        int current_active_vertex = 0;
        // initializing the computation
        startNodeId = (*initVerticalLoops)[i].nodeIndexSet[0];
        curArcId = (*initVerticalLoops)[i].simpArcIndexSet[0];
        //
        current_active_vertex = startNodeId;
        //
        vecVerticalLoopEdgePathOnMesh_Vertex[i].push_back(current_active_vertex);
        // simplifying the arcs
        //for (unsigned int j = 0; j < (*initVerticalLoops)[i].simpArcIndexSet.size(); j++)
        //{
        //	curArcId = (*initVerticalLoops)[i].simpArcIndexSet[j];
        //	startNodeId = (*initVerticalLoops)[i].nodeIndexSet[j];
        //	endNodeId = (*initVerticalLoops)[i].nodeIndexSet[j+1];
        //	//
        //	if ((*pVecSimplifiedArc)[curArcId].nCriticalNode0 == startNodeId)
        //	{// since the arc are traversed from high to low,
        //	// need to travel in the opposite direction
        //
        //	}
        //	else
        //	{// ignore the last point
        //	}
        //}
        //
        if (pathArcOnMeshPointType[curArcId].size() == 2) {// it is an edge
            int dst_node = pathArcOnMeshPointType[curArcId].back().first;
            int src_node = pathArcOnMeshPointType[curArcId].front().first;
            int edge_idx = -1;
            //
            for (std::vector<int>::iterator vIter = (*meshData).vecVertex[src_node].adjEdges.begin();
                 vIter != (*meshData).vecVertex[src_node].adjEdges.end();
                 vIter++) {
                if ((*meshData).vecEdge[*vIter].v0 == dst_node ||
                    (*meshData).vecEdge[*vIter].v1 == dst_node) {
                    edge_idx = *vIter;
                    break;
                }
            }
            //
            nTriangleIdForOrientation = meshData->vecEdge[edge_idx].AdjTri[0];
            //
            int third_vertex = meshData->vecTriangle[nTriangleIdForOrientation].v0 +
                               meshData->vecTriangle[nTriangleIdForOrientation].v1 +
                               meshData->vecTriangle[nTriangleIdForOrientation].v2 - (dst_node + src_node);
            //
            cycleMappedOnMesh.FindEdgeConnectingTwoVertices(third_vertex, current_active_vertex, edge_idx, refEdge);
        } else {//
            int edge_idx = 0;
            if ((*pVecSimplifiedArc)[curArcId].nCriticalNode0 == startNodeId) {
                edge_idx = pathArcOnMeshPointType[curArcId][pathArcOnMeshPointType[curArcId].size() - 2].first;
            } else {
                edge_idx = pathArcOnMeshPointType[curArcId][1].first;
            }
            //
            init_side_vertex = meshData->vecEdge[edge_idx].v0;
            //
            nTriangleIdForOrientation = meshData->TriangleAdjacentToOneEdgesOneVertex(edge_idx, current_active_vertex);
            //
            cycleMappedOnMesh.FindEdgeConnectingTwoVertices(current_active_vertex, init_side_vertex, edge_idx, refEdge);

        }
        //
        int startNodes[2] = {0, 0};
        for (unsigned int j = 0; j < (*initVerticalLoops)[i].simpArcIndexSet.size(); j++) {
            curArcId = (*initVerticalLoops)[i].simpArcIndexSet[j];
            startNodeId = (*initVerticalLoops)[i].nodeIndexSet[j];
            endNodeId = (*initVerticalLoops)[i].nodeIndexSet[j + 1];
            //
            if ((*pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to travel in the opposite direction
                if (pathArcOnMeshPointType[curArcId].size() == 2) {// update the start nodes first
                    startNodes[0] = startNodes[1] = pathArcOnMeshPointType[curArcId].front().first;
                }
                    //
                else {
                    int edge_idx = pathArcOnMeshPointType[curArcId][pathArcOnMeshPointType[curArcId].size() - 2].first;
                    startNodes[0] = meshData->vecEdge[edge_idx].v0;
                    startNodes[1] = meshData->vecEdge[edge_idx].v1;
                }

            } else {// ignore the last point
                if (pathArcOnMeshPointType[curArcId].size() == 2) {// update the start nodes first
                    startNodes[0] = startNodes[1] = pathArcOnMeshPointType[curArcId].back().first;
                }
                    //
                else {
                    int edge_idx = pathArcOnMeshPointType[curArcId][1].first;
                    startNodes[0] = meshData->vecEdge[edge_idx].v0;
                    startNodes[1] = meshData->vecEdge[edge_idx].v1;
                }
            }
            //
            // update the orientation
            // invariant :  current_active_vertex is always on refEdge and both of them are on the triange with triangeId
            while (meshData->vecTriangle[nTriangleIdForOrientation].v0 != startNodes[0] &&
                   meshData->vecTriangle[nTriangleIdForOrientation].v0 != startNodes[1] &&
                   meshData->vecTriangle[nTriangleIdForOrientation].v1 != startNodes[0] &&
                   meshData->vecTriangle[nTriangleIdForOrientation].v1 != startNodes[1] &&
                   meshData->vecTriangle[nTriangleIdForOrientation].v2 != startNodes[0] &&
                   meshData->vecTriangle[nTriangleIdForOrientation].v2 !=
                   startNodes[1]) {// go to another triangle along the refEdge
                nTriangleIdForOrientation =
                        meshData->vecEdge[refEdge].AdjTri[0] + meshData->vecEdge[refEdge].AdjTri[1] -
                        nTriangleIdForOrientation;
                //
                int three_edges[3] = {meshData->vecTriangle[nTriangleIdForOrientation].e01,
                                      meshData->vecTriangle[nTriangleIdForOrientation].e02,
                                      meshData->vecTriangle[nTriangleIdForOrientation].e12};
                for (int teid = 0; teid < 3; teid++) {
                    if (three_edges[teid] != refEdge &&
                        (meshData->vecEdge[three_edges[teid]].v0 == current_active_vertex ||
                         meshData->vecEdge[three_edges[teid]].v1 == current_active_vertex)) {
                        refEdge = three_edges[teid];
                        break;
                    }
                }
            }
            //
            if (meshData->vecTriangle[nTriangleIdForOrientation].v0 == startNodes[0] ||
                meshData->vecTriangle[nTriangleIdForOrientation].v1 == startNodes[0] ||
                meshData->vecTriangle[nTriangleIdForOrientation].v2 == startNodes[0]) {
                init_side_vertex = startNodes[0];
            } else {
                init_side_vertex = startNodes[1];
            }
            //
            if ((*pVecSimplifiedArc)[curArcId].nCriticalNode0 ==
                startNodeId) {// since the arc are traversed from high to low,
                // need to travel in the opposite direction
                cycleMappedOnMesh.WalkAlongArcWithInitialSide_InverseArcOrder(init_side_vertex,
                                                                              current_active_vertex,
                                                                              refEdge,
                                                                              nTriangleIdForOrientation,
                                                                              pathArcOnMeshPointType[curArcId],
                                                                              vecVerticalLoopEdgePathOnMesh_Vertex[i],
                                                                              vecVerticalLoopEdgePathOnMesh_Edge[i]);
            } else {// ignore the last point
                cycleMappedOnMesh.WalkAlongArcWithInitialSide_ArcOrder(init_side_vertex,
                                                                       current_active_vertex,
                                                                       refEdge,
                                                                       nTriangleIdForOrientation,
                                                                       pathArcOnMeshPointType[curArcId],
                                                                       vecVerticalLoopEdgePathOnMesh_Vertex[i],
                                                                       vecVerticalLoopEdgePathOnMesh_Edge[i]);
            }

        }
        //

        // push the first point into the loop to close it
        //vecVerticalLoopEdgePathOnMesh_Vertex[i].push_back(vecVerticalLoopEdgePathOnMesh_Vertex[i][0]);

        //vecVerticalLoopEdgePathOnMesh_Vertex[i] = cycleMappedOnMesh.outArcCycle_VertexOnMesh;
        //vecVerticalLoopEdgePathOnMesh_Edge[i] = cycleMappedOnMesh.outArcCycle_EdgeOnMesh;
        ////
        if (!cycleMappedOnMesh.VerifyingEmbeddedPath(vecVerticalLoopEdgePathOnMesh_Vertex[i],
                                                     vecVerticalLoopEdgePathOnMesh_Edge[i]))
            //{
            //	std::cout << "VALID V cycle" <<  i << std::endl;
            //}
            //else
        {
            std::cout << "NON-VALID V cycle " << i << std::endl;
        }
    }

    //////////////////
    //for (unsigned int i = 0; i < initVerticalLoops->size(); i++)
    //{
    //	//// store the output back
    //	int curArcId = 0;
    //	int startNodeId = 0;
    //	int endNodeId  = 0;
    //	// initializing the computation
    //	startNodeId = (*initVerticalLoops)[i].nodeIndexSet[0];
    //	vecVerticalLoopEdgePathOnMesh_Vertex[i].push_back(startNodeId);
    //	if ((
    //	for (unsigned int j = 0; j < (*initVerticalLoops)[i].simpArcIndexSet.size(); j++)
    //	{
    //		curArcId = (*initVerticalLoops)[i].simpArcIndexSet[j];
    //		startNodeId = (*initVerticalLoops)[i].nodeIndexSet[j];
    //		endNodeId = (*initVerticalLoops)[i].nodeIndexSet[j+1];
    //		//
    //		if ((*pVecSimplifiedArc)[curArcId].nCriticalNode0 != startNodeId)
    //		{// since the arc are traversed from high to low,
    //		// need to travel in the opposite direction
    //			// but the order is reversed in cycleMappedOnMesh.ComputeArcCycle()
    //			//
    //			for (unsigned int iarc = 0; iarc < vecArcEdgePathOnMesh_Vertex[curArcId].size() - 1 ; iarc++)
    //			{
    //				vecVerticalLoopEdgePathOnMesh_Vertex[i].push_back(vecArcEdgePathOnMesh_Vertex[curArcId][iarc]);
    //				vecVerticalLoopEdgePathOnMesh_Edge[i].push_back(vecArcEdgePathOnMesh_Edge[curArcId][iarc]);
    //			}
    //
    //		}
    //		else
    //		{// ignore the last point
    //			for ( int iarc = vecArcEdgePathOnMesh_Vertex[curArcId].size() - 1; iarc > 0 ; iarc--)
    //			{// be aware of which arc to take like _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
    //				vecVerticalLoopEdgePathOnMesh_Vertex[i].push_back(vecArcEdgePathOnMesh_Vertex[curArcId][iarc]);
    //				vecVerticalLoopEdgePathOnMesh_Edge[i].push_back(vecArcEdgePathOnMesh_Edge[curArcId][iarc - 1]);
    //			}
    //		}
    //	}
    //	// push the first point into the loop to close it
    //	vecVerticalLoopEdgePathOnMesh_Vertex[i].push_back(vecVerticalLoopEdgePathOnMesh_Vertex[i][0]);

    //	//vecVerticalLoopEdgePathOnMesh_Vertex[i] = cycleMappedOnMesh.outArcCycle_VertexOnMesh;
    //	//vecVerticalLoopEdgePathOnMesh_Edge[i] = cycleMappedOnMesh.outArcCycle_EdgeOnMesh;
    //	////
    //	if (!cycleMappedOnMesh.VerifyingEmbeddedPath(vecVerticalLoopEdgePathOnMesh_Vertex[i], vecVerticalLoopEdgePathOnMesh_Edge[i]))
    //	//{
    //	//	std::cout << "VALID V cycle" <<  i << std::endl;
    //	//}
    //	//else
    //	{
    //		std::cout << "NON-VALID V cycle " <<  i << std::endl;
    //	}
    //}
    return;
}

void psbmReebGraph::WriteSimplifiedReebGraphOBJ(const char *pFileName, std::vector<_SimpleMeshVertex>& verts) {
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    std::ofstream ofile;
    ofile.open(pFileName, std::ifstream::out);
    if (ofile.is_open()) {
        sstr << "# V: " << pVecReebNode->size() << " E: " << pListReebArc->size() << std::endl;
        for (unsigned int i = 0; i < pVecReebNode->size(); i++) {
            assert(i == (*pVecReebNode)[i].nVertexId);
            sstr << "v " << verts[i].x << ' ' << verts[i].y << ' ' << verts[i].z << std::endl;
        }
        for (auto sIter = pVecSimplifiedArc->begin(); sIter != pVecSimplifiedArc->end(); sIter++) {
            sstr << "l " << (*sIter).nCriticalNode0+1 << " " << (*sIter).nCriticalNode1+1 << std::endl;
        }
        //
        ofile << sstr.rdbuf();
        //
        ofile.close();
        ofile.clear();
        //
        sstr.str("");
        sstr.clear();
    } else {
        std::cout << "can not open " << pFileName << std::endl;
        exit(0);
    }
    return;
}

void psbmReebGraph::WriteReebGraphOBJ(const char *pFileName, std::vector<_SimpleMeshVertex>& verts) {
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    std::ofstream ofile;
    ofile.open(pFileName, std::ifstream::out);
    if (ofile.is_open()) {
        sstr << "# V: " << pVecReebNode->size() << " E: " << pListReebArc->size() << std::endl;
        for (unsigned int i = 0; i < pVecReebNode->size(); i++) {
            assert(i == (*pVecReebNode)[i].nVertexId);
            sstr << "v " << verts[i].x << ' ' << verts[i].y << ' ' << verts[i].z << std::endl;
        }
        for (std::set<class psbmReebArc *>::iterator sIter = pListReebArc->begin();
             sIter != pListReebArc->end();
             sIter++) {
            sstr << "l " << (*sIter)->nNodeId0+1 << " " << (*sIter)->nNodeId1+1 << std::endl;
        }
        //
        ofile << sstr.rdbuf();
        //
        ofile.close();
        ofile.clear();
        //
        sstr.str("");
        sstr.clear();
    } else {
        std::cout << "can not open " << pFileName << std::endl;
        exit(0);
    }
    return;
}
void psbmReebGraph::WriteReebGraph(const char *pFileName) {
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    std::ofstream ofile;
    ofile.open(pFileName, std::ifstream::out);
    if (ofile.is_open()) {
//        sstr << "# " << pVecReebNode->size() << " " << pListReebArc->size() << std::endl;
        for (unsigned int i = 0; i < pVecReebNode->size(); i++) {
            assert(i == (*pVecReebNode)[i].nVertexId);
            sstr << (*pVecReebNode)[i].nVertexId << std::endl;
        }
        for (std::set<class psbmReebArc *>::iterator sIter = pListReebArc->begin();
             sIter != pListReebArc->end();
             sIter++) {
            sstr << (*sIter)->nNodeId0 << " " << (*sIter)->nNodeId1 << std::endl;
        }
        //
        ofile << sstr.rdbuf();
        //
        ofile.close();
        ofile.clear();
        //
        sstr.str("");
        sstr.clear();
    } else {
        std::cout << "can not open " << pFileName << std::endl;
        exit(0);
    }
    return;
}

void psbmReebGraph::ReadReebGraph(const char *pFileName) {
    std::cout << "reading... " << std::endl;

    std::ifstream ifile;
    ifile.open(pFileName, std::ifstream::in);
    long long fileSize = 0;
    char *fBuf = NULL;

    std::string sBuf;
    //
    int vertexNum = 0;
    int edgeNum = 0;
    //
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    Vector3 tempVec;
    if (ifile.is_open()) {
        ifile.seekg(0, std::ios::end);
        fileSize = ifile.tellg();
        // move pointer back to beginning
        ifile.seekg(0, std::ios::beg);
        //
        //copy whole file into file buffer
        fBuf = new char[fileSize + 1];
        ifile.read(fBuf, fileSize);
        // add extra symbol
        fBuf[fileSize] = '\n';
        //
        sBuf.assign(fBuf);
        sstr.str(sBuf);

        sBuf.clear();
        // close file
        ifile.close();

        //deallocate memory
        delete[] fBuf;

        //
        int upVerNum = 0;
        int downVerNum = 0;
        int nVertexId = 0;
        int N0 = 0;
        int N1 = 0;
        sstr >> vertexNum >> edgeNum;
        if (!pVecReebNode) {
            pVecReebNode = new std::vector<class psbmReebNode>;
            pVecReebNode->reserve(vertexNum);
        } else {
            pVecReebNode->clear();
            pVecReebNode->reserve(vertexNum);
        }
        if (!pListReebArc) {
            pListReebArc = new std::set<class psbmReebArc *>;
        } else {
            pListReebArc->clear();
        }
        //if (!ProcessedVertex)
        //{
        //	ProcessedVertex = new std::map<int, int>;
        //}
        //else
        //{
        //	ProcessedVertex->clear();
        //}

        for (int i = 0; i < vertexNum; i++) {
            sstr >> nVertexId;
            psbmReebNode tmpNode(nVertexId);
            pVecReebNode->push_back(tmpNode);
            //(*ProcessedVertex)[nVertexId] = i;
            //std::cout<< i << "\t";
        }
        for (int i = 0; i < edgeNum; i++) {
            sstr >> N0 >> N1;
            psbmReebArc *tmpArc = new psbmReebArc(N0, N1, 0);
            pListReebArc->insert(tmpArc);
            // set graph node
            (*pVecReebNode)[N0].ptrArcUpId->push_back(tmpArc);
            (*pVecReebNode)[N1].ptrArcDownId->push_back(tmpArc);
        }

        sstr.clear();
    } else {
        std::cout << "Can NOT open file " << pFileName << std::endl;
        exit(0);
    }
    std::cout << "Done... " << vertexNum << " " << std::endl;
    return;
}
