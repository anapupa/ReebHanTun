/*
(c) 2012 Fengtao Fan
*/
#pragma once

#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "SimpleMesh.h"
#include "psbmReebGraphElements.h"
#include <RenderVector3.h>

void RotateAxisTobeAlignedWithZ_axis(double u, double v, double w,
                                     double x, double y, double z,
                                     Vector3 &outRes);

class __psbmInterimGraph;

void SetAllNull(auxMeshEdge &inEdge);

class psbmReebGraph {

public:
    psbmReebGraph() {
        nMeshEdgeNum = 0;
        //
        //ProcessedVertex = NULL;
        pVecAuxMeshEdge = NULL;
        //
        pVecReebNode = NULL;
        pListReebArc = NULL;
        //
        pVecSimplifiedArc = NULL;
        criticalNodeIdSet = NULL;
        criticalPairing = NULL;
        //
        initVerticalLoops = NULL;
        ReebCycles = NULL;
    }

    psbmReebGraph(const int n) {
        nMeshEdgeNum = n;
        //ProcessedVertex = NULL;
        pVecAuxMeshEdge = NULL;
        //
        pVecReebNode = NULL;
        pListReebArc = NULL;
        //
        pVecSimplifiedArc = NULL;
        criticalNodeIdSet = NULL;
        criticalPairing = NULL;
        //
        initVerticalLoops = NULL;
        ReebCycles = NULL;
        //
        if (nMeshEdgeNum > 0) {// nMeshEdgeNum > 0
            pVecAuxMeshEdge = new std::vector<struct auxMeshEdge>;
            pVecAuxMeshEdge->reserve(nMeshEdgeNum);
            struct auxMeshEdge tmp;
            tmp.ptrPos = NULL;
            tmp.ptrReebArcId = NULL;
            for (int i = 0; i < nMeshEdgeNum; i++) {
                pVecAuxMeshEdge->push_back(tmp);
            }
            std::for_each(pVecAuxMeshEdge->begin(), pVecAuxMeshEdge->end(), SetAllNull);
        }
        //
        //ProcessedVertex = new std::map<int, int>;

    }

    ~psbmReebGraph() {
        //if (ProcessedVertex)
        //{
        //	ProcessedVertex->clear();
        //	delete ProcessedVertex;
        //}
        if (pVecAuxMeshEdge) {
            pVecAuxMeshEdge->clear();
            delete pVecAuxMeshEdge;
        }
        if (pVecReebNode) {
            pVecReebNode->clear();
            delete pVecReebNode;
        }
        if (pVecSimplifiedArc) {
            pVecSimplifiedArc->clear();
            delete pVecSimplifiedArc;
        }
        if (pListReebArc) {
            std::for_each(pListReebArc->begin(), pListReebArc->end(), DeleteObject());
            pListReebArc->clear();
            delete pListReebArc;
        }
        if (criticalNodeIdSet) {
            criticalNodeIdSet->clear();
            delete criticalNodeIdSet;
        }
        if (criticalPairing) {
            criticalPairing->clear();
            delete criticalPairing;
        }
        if (initVerticalLoops) {
            {
                std::vector<psbmReebPath> tmp;
                initVerticalLoops->swap(tmp);
            }
            delete initVerticalLoops;
        }
        if (ReebCycles) {
            {
                std::vector<psbmReebPath> tmp;
                ReebCycles->swap(tmp);
            }
            ReebCycles->clear();
            delete ReebCycles;
        }
    }

    void ReserveSpaceForEdges(const int n) {
        nMeshEdgeNum = n;
        //ProcessedVertex = NULL;
        pVecAuxMeshEdge = NULL;
        //
        pVecReebNode = NULL;
        pListReebArc = NULL;
        //
        pVecSimplifiedArc = NULL;
        criticalNodeIdSet = NULL;
        criticalPairing = NULL;
        //
        initVerticalLoops = NULL;
        ReebCycles = NULL;
        //
        if (nMeshEdgeNum > 0) {// nMeshEdgeNum > 0
            pVecAuxMeshEdge = new std::vector<struct auxMeshEdge>;
            pVecAuxMeshEdge->reserve(nMeshEdgeNum);
            struct auxMeshEdge tmp;
            tmp.ptrPos = NULL;
            tmp.ptrReebArcId = NULL;
            for (int i = 0; i < nMeshEdgeNum; i++) {
                pVecAuxMeshEdge->push_back(tmp);
            }
            std::for_each(pVecAuxMeshEdge->begin(), pVecAuxMeshEdge->end(), SetAllNull);
        }
        //
        //ProcessedVertex = new std::map<int, int>;
    }

    //set the height direction
    void SetHeightDirection(Vector3 &inDir) {
        heightDirection = inDir;
    }

    //
    void CreateReebGraph(int nEdges);

    void AssignData(_SimpleMesh *inMesh, double *scalarField);

    void ComputeReebGraph();

    void AddMeshTriangle(int v0, double f0,
                         int v1, double f1,
                         int v2, double f2,
                         int e01, int e12, int e02);

    //
    void InitReebNodes(const unsigned int nTotSize);

    int CreateNode(const int nVertexId);//, const float scalarValue);
    void CreateArc(const int v0, const int v1, const int eId);

    void MergePaths(const int N0, const int N1, const int N2,
                    const int e0, const int e1, const int e2);

    void GlueByMergeSorting(const int N0, const int N1,
                            psbmReebArc **a0, psbmReebArc **a1,
                            psbmArcListNode **travE0, psbmArcListNode **travE1);

    void DeleteArc(psbmReebArc *a_del);

    void MergeArcs(psbmReebArc *a_high, psbmReebArc *a_low);

    void DeleteArcListNodeForEdge(const int eid);

    //
    void CheckConsistency(int flag);

    //
    bool IsUpArcNodesMax(const int nodeIdx);

    bool IsDownArcNodesMin(const int nodeIdx);

    void HandlingCritNodeWithDegreeGreaterThan3();

    void ReebGraphPairing(_simpGraph &rbGraph, std::vector<std::pair<int, int> > &outPairing);

    //
    bool SimplifyReebGraph();

    int DownWalkArcs(class psbmReebArc *inArc, const int label, bool &ret);

    int UpWalkArcs(class psbmReebArc *inArc, const int label, bool &ret);

    //
    int DownwardWalkingForPairing(const int startingNode, REEB_GRAPH_CRITICAL_TYPE &pairNodeType);

    int TheOtherEndPoint(const int nodeId, const int SimpArcNodeId);

    //
    void WalkFromHighLevelToLowLevel(std::set<int> &highPath, std::set<int> &lowPath,
                                     double &highHeight, int &highNodeId,
                                     const double &lowHeight, const int &lowNodeId,
                                     int &pairNode, int &highMINNode);

    //
    void ComputingPairing();

    void ComputingPairing(const int __polymorphism__);

    //export
    void ExportMathematicaGraph();
    //
    /*********************************************************************************/
    // NEED TO COPY OUT FOR CODE FILE CONSISTENCY
    // Simplifyin the simplified arc again
    // to remove non-necessary loops
    // with phantom nonzero value -1 for arc conncecting to paired node
    int NonzeroArcValenceInSimplifiedArcGraphForNode(const int node_id, bool UpwardDir);

    void
    GoAlong_1_or_2_degree_ReebArc(class psbmReebArc *&travArc, const int node_id, const int arc_label, bool UpwarDir);

    void ResetArcLabel(const int nodeId, const int arcId, const int newArcLabel);

    int ResetArcLabel(const int nodeId, const int arcId, std::set<int> &pairedCritNodes);

    void RemoveSimplifiedArcs();

    void MergeTwoSimplifiedArcs();

    void ResetSimplifiedArcIds();

    void AfterPairingSimplifyingForSimplifiedArcs();
    //
    /***************maximum spanning tree******************************************************************/
    void ComputingCycle_max_tree();
    /*********************************************************************************/
    /*********************************************************************************/

    /*********************************************************************************/
    //cut level set and trace inside / outside properties
    void RotateAlongAxis_ReebGraph(double u, double v, double w,
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

    void OptimizeXYAxisMakeNoParaOrtho(std::vector<_Polygon> &inPolygons);

    int DetermineVerticalLoopType(const int upforkingReebNodeId, std::vector<_Polygon> &levelSetPolygons);

    void CutLevelSet_specific_value(_Polygon &retPly,
                                    SimplifiedReebGraphArc &inArc,
                                    const double levelSetValue);

    void CutLevelSet_specific_value_pair(_Polygon &retPly,
                                         std::vector<std::pair<int, int> > &retPtType,
                                         SimplifiedReebGraphArc &inArc,
                                         const double levelSetValue);

    void CutLevelSetAlongOneEdgeOnMesh(_Polygon &retPlg,
                                       const int edgeIndex,
                                       const double levelSetValue,
                                       const int levelSetIndex,
                                       const int targetArcIdx,
                                       const bool special_case);

    void CutLevelSetAlongOneEdgeOnMesh_pair(_Polygon &retPlg,
                                            std::vector<std::pair<int, int> > &retPair,
                                            const int edgeIndex,
                                            const double levelSetValue,
                                            const int levelSetIndex,
                                            const int targetArcIdx,
                                            const bool special_case);

    /*****************************************************************************************************/
    void PushIntoQueueForNonVisitedNodes_UP(const int nodeId, const int vecIndex,
                                            std::vector<std::pair<int, bool> > &processed_node,
                                            std::vector<bool> &processed_arc_bit,
                                            std::queue<int> &ready_to_decide_node);

    void PushIntoQueueForNonVisitedNodes(const int nodeId, const int vecIndex,
                                         std::vector<std::pair<int, bool> > &processed_node,
                                         std::vector<bool> &processed_arc_bit,
                                         std::queue<int> &ready_to_decide_node);

    void PushIntoQueueForNonVisitedNodes_DOWN(const int nodeId, const int vecIndex,
                                              std::vector<std::pair<int, bool> > &processed_node,
                                              std::vector<bool> &processed_arc_bit,
                                              std::queue<int> &ready_to_decide_node);

    //
    char LevelSetSeparationProperty(_Polygon &lhs, _Polygon &rhs, const int perpAxis, const int rayAxis);

    bool InsidePolygon(_Polygon &refPly, const Vector3 &basePt, const int perpAxis, const int rayAxis);

    bool InsideRegion(std::vector<_Polygon> &refPolygons, const Vector3 &basePt, const int perpAxis, const int rayAxis);
    /*****************************************************************************************************/
    //
    void TraceOutOneCycleForEachPairing(std::pair<int, int> &critPair, psbmReebPath &cycleInReeb);

    //
    void PathToSourceViaBFS(const int target,
                            std::vector<int> &nodeOnPath,
                            std::vector<int> &parents,
                            _simpGraph &inGraph);

    //void FindingCorrespondingSimplifiedArcs(
    //	_simpGraph &graph_map,
    //		std::vector<int> &nodeOnPath,
    //		std::vector<int> &arcOnPath);
    void FindingCorrespondingSimplifiedArcs(
            std::vector<int> &nodeOnPath,
            std::vector<int> &arcOnPath);

    void Add_Simplified_Arc_into_Graph(const int node_a,
                                       const int node_b,
                                       int &nGraphNodeCounter,
                                       boost::unordered_map<int, int> &inMapping, _simpGraph &outGraph);

    void construct_simplified_arc_subgraph_between_two_critical_pionts(int &lowPtInOutputGraph,
                                                                       int &highPtInOutputGraph,
                                                                       int &downforking_down_node_a,
                                                                       int &downforking_down_node_b,
                                                                       std::pair<int, int> &critPair,
                                                                       _simpGraph &outGraph);

    // Compute the vertical cycle in Reeb graph
    void EdgeIntersectArcs(const int edgeId, std::set<int> vecArcId);

    void ComputeUpwardPathBetweenPair(int const u, int const t, psbmReebPath &finalPath,
                                      psbmReebPath &currentPath);

    void TraceOutOneCycleForEachPairing(int const lowPt, int const highPt, psbmReebPath &cycleInReeb);

    // find the preimage of each cycle in simplified reeb graph
    // that is, embed each cycle back to mesh
    void EmbedCyclesOnMesh(const int lowCritPt, const int highCritPt, psbmReebPath &arcPath, _Polygon &retPath);

    void ConstructSubgraph(int &leftLowPtInGraphNode,
                           int &leftHighPtInGraphNode,
                           int &rightLowPtInGraphNode,
                           int &rightHighPtInGraphNode,
                           psbmReebPath &arcPath, _simpGraph &leftPathGraph, _simpGraph &rightPathGraph);

    void ComputePathToSource(_Polygon &path, std::vector<int> &vertexInMesh, int target, std::vector<int> &parents,
                             _simpGraph &inGraph);

    void ComputePathToSource_pair(_Polygon &path, std::vector<std::pair<int, int> > &vertexInMesh, int target,
                                  std::vector<int> &parents, _simpGraph &inGraph);

    void ComputePathToSource(_Polygon &path, int target, std::vector<int> &parents, _simpGraph &inGraph);

    void EliminatePathInGraph(const int target, std::vector<int> &parents, _simpGraph &inGraph);

    void CutAllLevelsets(const int upforkingReebNodeId, std::vector<_Polygon> &levelSetPolygons);

    /*****************************************************************************/
    void ComputeCycleAndPairing();
    /*********************************************************/
    //construct the path for each simplified arc
    void add_vertex_one_ring_to_graph(const int vertexId,
                                      const double lowValue,
                                      const double highValue,
                                      const int lowVertexIdOnMesh,
                                      const int higherVertexIdOnMesh,
                                      int &nGraphNodeCounter,
                                      std::map<int, int> &inMapping, _simpGraph &outGraph);

    void add_vertex_one_ring_to_graph(const int vertexId_a,
                                      const int vertexId_b,
                                      int &nGraphNodeCounter,
                                      std::map<int, int> &inMapping, _simpGraph &outGraph);

    bool map_levelset_on_edge_to_this_reeb_arc(
            const int target_arc_idx,
            const int nReebNode_0,
            const int nReebNode_1,
            const int inEdgeIdx,
            const bool special_case_no_points_between_two_crit_pts = false);

    void add_vertex_one_ring_to_graph(const int vertexId,
                                      int &nGraphNodeCounter,
                                      std::map<int, int> &MeshVertexIndexToGraphNodeIndex,
                                      _simpGraph &outGraph);

    //
    void construct_path_on_mesh_for_reeb_arc_by_crossing_mesh_edges(
            SimplifiedReebGraphArc &inSimplifiedArc,
            //_Polygon &outPath,
            std::vector<std::pair<int, int> > &outPathOnMeshPointType,
            std::vector<int> &parents);

    void ComputeParallelPath_pair_for_path_crossing_mesh_edges(std::vector<std::pair<int, int> > &inPathVer,
                                                               std::vector<std::pair<int, int> > &outPathVer,
                                                               _Polygon &inPath, _Polygon &outPath,
                                                               const int lowIndexInMesh,
                                                               const int highIndexInMesh);

    void compute_representative_point_on_edge(const int lowVertexOnMesh,
                                              const int highVertexOnMesh,
                                              const int eid,
                                              const int dirVertrex,
                                              Vector3 &repPt_org,
                                              Vector3 &repPt_dir);

    void compute_representative_point_on_edge(const int lowVertexOnMesh,
                                              const int highVertexOnMesh,
                                              const int eid,
                                              Vector3 &repPt);


    void compute_path_on_mesh_for_each_simplified_arc();//std::vector<_Polygon> &resPathOnMesh,
    //std::vector<_Polygon> &resOffsetPolygonOnMesh);
    /**********************************************************************************/
    /**********************************************************************************/
    void SetSimplifiedArcFlagBits();
    /**********************************************************************************/
    // compute the intersections between two loops
    void ComputeIntersections(_Polygon &loop_a,
                              std::vector<std::pair<int, int> > &point_type_a,
                              _Polygon &loop_b,
                              std::vector<std::pair<int, int> > &point_type_b);

    void PullInsideOrOutsideForIntersectionPoint();
    /**********************************************************************************/
    // Compute the link between two disjoint loops
    void FindPointOnLevelSetIntersectingArc(const std::vector<std::pair<int, int> > &ptType,
                                            const _Polygon &ptsOnPath,
                                            const double levelSetValue,
                                            Vector3 &pt_on_levelset);

    void LinkNumberMatrixComputing();

    void LinkNumberMatrixComputing_withLoopsOnFly();

    /**********************************************************************************/
    void EmbedCycleAsEdgePathOnMesh();

    /**********************************************************************************/
    void WriteReebGraph(const char *pFileName);

    /**********************************************************************************/
    void WriteReebGraphOBJ(const char *pFileName, std::vector<_SimpleMeshVertex>& verts);
    void WriteSimplifiedReebGraphOBJ(const char *pFileName, std::vector<_SimpleMeshVertex>& verts);

    void ReadReebGraph(const char *pFileName);
    /**********************************************************************************/
public:
    _SimpleMesh *meshData;
    double *MeshVertexScalarValue;
    //
    std::vector<std::vector<std::pair<int, int> > > edgeIntersectArcs;
    //
    int nMeshEdgeNum;
    std::vector<bool> AddedTriangle;
    std::vector<int> vecVertexMappingToSimplifiedReebArc;
    //std::map<int, int> *ProcessedVertex;
    std::vector<struct auxMeshEdge> *pVecAuxMeshEdge;
    std::vector<class psbmReebNode> *pVecReebNode;
    std::set<class psbmReebArc *> *pListReebArc;
    // used for removing memory
    // used for simplified Reeb graph
    std::vector<struct SimplifiedReebGraphArc> *pVecSimplifiedArc;
    // flag bits for simplified arcs
    std::vector<std::vector<int> > arcIntersectWithVerticalLoops; // the size of simplified arc
    std::vector<int> arcIntersectWithHorizontalLoops; // the size of simplified arc
    //
    std::set<int> *criticalNodeIdSet; // index in Reeb graph nodes array
    // pointing out normal direction
    std::map<int, Vector3> unitOutNormalForCritNode; // the same size for critical nodes
    // the first part is lower pt
    // the second part is high pt
    std::vector<std::pair<int, int> > *criticalPairing;
    //std::vector<std::pair<int, int>> trueCriticalPairing;
    // used for traversing
    int minHeightNodeId;
    // pathes
    std::vector<psbmReebPath> *initVerticalLoops;
    std::vector<int> initSimplifiedArcIdForHorizontalLoops; // the same size and order as initVerticalLoops
    std::vector<psbmReebPath> *ReebCycles;
    // direction useda
    int scalarDir;
    // height direction
    Vector3 heightDirection;
    // utility data structure for identifying which vertex is taken for
    // std::vector<bool> vertexTakenForPathOnMesh;
    //std::vector<bool> edgeTakenForPathOnMesh;
    //
    // levelset associated to each simplified arc
    std::vector<_Polygon> levelsetOnMesh; // the same size and order as initVerticalLoops
    /*point type is to indicate this point is vertex or on an edge*/
    std::vector<std::vector<std::pair<int, int> > > levelsetOnMeshPointType;
    ///the embedded path on the mesh for each simplified arc
    std::vector<_Polygon> pathArcOnMesh;
    std::vector<_Polygon> offsetPathArcOnMesh;
    std::vector<_Polygon> sepPathArcOnMesh;
    std::vector<_Polygon> sepOffsetPathArcOnMesh;
    std::vector<std::vector<std::pair<int, int> > > pathArcOnMeshPointType;
    std::vector<std::vector<std::pair<int, int> > > offsetPathArcOnMeshPointType;

    // for visualization purpose
    _simpGraph viewGraph;
    // link numbers for different loops
    // the first half is the order of vertical loops
    // the second half is the corresponding horizontal loops
    //std::vector<std::vector<std::pair<int, int>>> linkNumbersMatrix;
    std::vector<std::list<std::pair<int, int> > > linkNumbersMatrix;
    std::vector<std::list<std::pair<int, int> > > invLinkNumbersMatrix;
    //
    std::vector<std::vector<int> > vecVerticalLoopEdgePathOnMesh_Vertex;
    std::vector<std::vector<int> > vecVerticalLoopEdgePathOnMesh_Edge;
    //
    std::vector<std::vector<int> > vecHorizontalLoopEdgePathOnMesh_Vertex;
    std::vector<std::vector<int> > vecHorizontalLoopEdgePathOnMesh_Edge;
    //


};

inline int psbmReebGraph::TheOtherEndPoint(const int nodeId, const int SimpArcNodeId) {// 0 is reserved
    if ((*pVecSimplifiedArc)[SimpArcNodeId - 1].nCriticalNode0 == nodeId)
        return (*pVecSimplifiedArc)[SimpArcNodeId - 1].nCriticalNode1;
    return (*pVecSimplifiedArc)[SimpArcNodeId - 1].nCriticalNode0;
}
