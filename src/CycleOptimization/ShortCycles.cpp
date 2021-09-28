/*
(c) 2012 Fengtao Fan
*/
#include "SimpleGraph.h"
#include "SimpleMesh.h"
#include "canonical_loops.h"
#include "edge_annotations_gauss.h"
#include "CycleOptimization/Annotation/AnnotationComputation.h"
#include "psbmReebGraph.h"

#include <boost/progress.hpp>

#include <fstream>
#include <sstream>

#include <map>

//_simpGraph_vec meshGraph;
//_SimpleMesh cycMesh;
//std::vector<Vector3>   meshNormal;
void LoadCycleData(psbmReebGraph &reebGraph, std::vector<std::vector<int> > &basis_h_loops,
                   std::vector<std::vector<int> > &basis_v_loops, std::vector<int> &base_pt_vec) {
    //
    std::vector<bool> initVerticalLoopsType;
    // set the flag bit for handle or tunnel loops
    for (unsigned int i = 0; i < reebGraph.initVerticalLoops->size(); i++) {
        if ((*reebGraph.initVerticalLoops)[i].pathType) {// this is a vertical loop
            initVerticalLoopsType.push_back(true);
        } else {
            initVerticalLoopsType.push_back(false);
        }
    }
    // Store the critical points
    base_pt_vec.reserve(2 * reebGraph.criticalPairing->size());
    for (unsigned int j = 0; j < reebGraph.criticalPairing->size(); j++) {
        base_pt_vec.push_back((*reebGraph.criticalPairing)[j].first);
        base_pt_vec.push_back((*reebGraph.criticalPairing)[j].second);
    }
    //
    const int tot_genus = (int) reebGraph.vecVerticalLoopEdgePathOnMesh_Vertex.size();
    basis_v_loops.reserve(tot_genus);
    basis_h_loops.reserve(tot_genus);
    // store v-basis
    for (int i = 0; i < tot_genus; i++) {
        std::vector<int> tempIntVec;
        int edges_num = 0;
        int row_num = i;
        //
        if (!initVerticalLoopsType[i])
            row_num += tot_genus;
        //
        for (std::list<std::pair<int, int> >::iterator listIter = (reebGraph.invLinkNumbersMatrix)[row_num].begin();
             listIter != (reebGraph.invLinkNumbersMatrix)[row_num].end(); listIter++) {
            int col_num = listIter->first;
            int counter = abs(listIter->second);
            if (col_num >= tot_genus) {// it is a horizontal loop
                edges_num +=
                        counter * (int) (reebGraph.vecHorizontalLoopEdgePathOnMesh_Edge)[col_num - tot_genus].size();
            } else {
                edges_num += counter * (int) (reebGraph.vecVerticalLoopEdgePathOnMesh_Edge)[col_num].size();
            }
        }
        //
        //sstr << edges_num << std::endl;
        tempIntVec.reserve(edges_num + 2);
        for (std::list<std::pair<int, int> >::iterator listIter = (reebGraph.invLinkNumbersMatrix)[row_num].begin();
             listIter != (reebGraph.invLinkNumbersMatrix)[row_num].end(); listIter++) {
            int col_num = listIter->first;
            int counter = abs(listIter->second);
            if (col_num >= tot_genus) {// it is a horizontal loop
                for (int j = 0; j < counter; j++) {
                    for (std::vector<int>::iterator vIter = (reebGraph.vecHorizontalLoopEdgePathOnMesh_Edge)[col_num -
                                                                                                             tot_genus].begin();
                         vIter !=
                         (reebGraph.vecHorizontalLoopEdgePathOnMesh_Edge)[col_num - tot_genus].end(); vIter++) {
                        //sstr << *vIter << " ";
                        tempIntVec.push_back(*vIter);
                    }
                }
            } else {
                for (int j = 0; j < counter; j++) {
                    for (std::vector<int>::iterator vIter = (reebGraph.vecVerticalLoopEdgePathOnMesh_Edge)[col_num].begin();
                         vIter != (reebGraph.vecVerticalLoopEdgePathOnMesh_Edge)[col_num].end(); vIter++) {
                        //sstr << *vIter << " ";
                        tempIntVec.push_back(*vIter);
                    }
                }
            }
        }
        //
        basis_v_loops.push_back(tempIntVec);
        //sstr << std::endl;
    }
    //store h-basis
    for (int i = 0; i < tot_genus; i++) {
        std::vector<int> hTempIntVec;
        int edges_num = 0;
        int row_num = i;
        //
        if (initVerticalLoopsType[i])
            row_num += tot_genus;
        //
        for (std::list<std::pair<int, int> >::iterator listIter = (reebGraph.invLinkNumbersMatrix)[row_num].begin();
             listIter != (reebGraph.invLinkNumbersMatrix)[row_num].end(); listIter++) {
            int col_num = listIter->first;
            int counter = abs(listIter->second);
            if (col_num >= tot_genus) {// it is a horizontal loop
                edges_num +=
                        counter * (int) (reebGraph.vecHorizontalLoopEdgePathOnMesh_Edge)[col_num - tot_genus].size();
            } else {
                edges_num += counter * (int) (reebGraph.vecVerticalLoopEdgePathOnMesh_Edge)[col_num].size();
            }
        }
        //
        //sstr << edges_num << std::endl;
        hTempIntVec.reserve(edges_num + 2);
        for (std::list<std::pair<int, int> >::iterator listIter = (reebGraph.invLinkNumbersMatrix)[row_num].begin();
             listIter != (reebGraph.invLinkNumbersMatrix)[row_num].end(); listIter++) {
            int col_num = listIter->first;
            int counter = abs(listIter->second);
            if (col_num >= tot_genus) {// it is a horizontal loop
                for (int j = 0; j < counter; j++) {
                    for (std::vector<int>::iterator vIter = (reebGraph.vecHorizontalLoopEdgePathOnMesh_Edge)[col_num -
                                                                                                             tot_genus].begin();
                         vIter !=
                         (reebGraph.vecHorizontalLoopEdgePathOnMesh_Edge)[col_num - tot_genus].end(); vIter++) {
                        //sstr << *vIter << " ";
                        hTempIntVec.push_back(*vIter);
                    }
                }
            } else {
                for (int j = 0; j < counter; j++) {
                    for (std::vector<int>::iterator vIter = (reebGraph.vecVerticalLoopEdgePathOnMesh_Edge)[col_num].begin();
                         vIter != (reebGraph.vecVerticalLoopEdgePathOnMesh_Edge)[col_num].end(); vIter++) {
                        //sstr << *vIter << " ";
                        hTempIntVec.push_back(*vIter);
                    }
                }
            }
        }
        //sstr << std::endl;
        basis_h_loops.push_back(hTempIntVec);
    }
    return;
}

void AssembleMeshGraph(_simpGraph_vec &meshGraph, _SimpleMesh &mesh, std::vector<float> &sqEdgeLength,
                       std::set<int> &extraVertices, const float fScaleRatio) {
    meshGraph.InitNodes((int) mesh.vecVertex.size());
    //
    std::vector<bool> vertexFlag(mesh.vecVertex.size(), false);
    //
    for (std::set<int>::iterator sIter = extraVertices.begin();
         sIter != extraVertices.end(); sIter++) {
        vertexFlag[*sIter] = true;
    }
    //
    std::vector<int> nonExistEdges;
    //
    float totEdgeLen = 0.0;
    for (unsigned int eid = 0; eid < mesh.vecEdge.size(); eid++) {
        int v_a = mesh.vecEdge[eid].v0;
        int v_b = mesh.vecEdge[eid].v1;

        Vector3 p_a(mesh.vecVertex[v_a].x, mesh.vecVertex[v_a].y, mesh.vecVertex[v_a].z);
        Vector3 p_b(mesh.vecVertex[v_b].x, mesh.vecVertex[v_b].y, mesh.vecVertex[v_b].z);

        //
        sqEdgeLength.push_back(norm(p_a - p_b) * fScaleRatio);
        //

        meshGraph.AddEdge(v_a, v_b, eid);
        //
        totEdgeLen += sqEdgeLength.back();
        //
        if (vertexFlag[v_a] || vertexFlag[v_b]) {
            nonExistEdges.push_back(eid);
        }
    }
    //
    for (unsigned int i = 0; i < nonExistEdges.size(); i++) {
        //
        sqEdgeLength[nonExistEdges[i]] = totEdgeLen;// std::numeric_limits<float>::max(); // totEdgeLen* 2;
    }
    return;
}

void AssembleMeshGraph(_simpGraph_vec &meshGraph, _SimpleMesh &mesh, std::vector<float> &sqEdgeLength,
                       const float fScaleRatio) {
    meshGraph.InitNodes((int) mesh.vecVertex.size());
    //
    for (unsigned int eid = 0; eid < mesh.vecEdge.size(); eid++) {
        int v_a = mesh.vecEdge[eid].v0;
        int v_b = mesh.vecEdge[eid].v1;

        Vector3 p_a(mesh.vecVertex[v_a].x, mesh.vecVertex[v_a].y, mesh.vecVertex[v_a].z);
        Vector3 p_b(mesh.vecVertex[v_b].x, mesh.vecVertex[v_b].y, mesh.vecVertex[v_b].z);

        //
        sqEdgeLength.push_back(norm(p_a - p_b) * fScaleRatio);
        //

        meshGraph.AddEdge(v_a, v_b, eid);
    }
}


float BasisWeight(std::vector<std::set<int> > &basis, std::vector<float> sqEdgeLen) {
    float tot_len = 0.0;
    for (unsigned int i = 0; i < basis.size(); i++) {
        for (std::set<int>::iterator sIter = basis[i].begin();
             sIter != basis[i].end(); sIter++) {
            tot_len += sqEdgeLen[*sIter];
        }
    }
    return tot_len;
}

float BasisWeight(std::vector<std::vector<int> > &basis, std::vector<float> sqEdgeLen) {
    float tot_len = 0.0;
    for (unsigned int i = 0; i < basis.size(); i++) {
        for (std::vector<int>::iterator vIter = basis[i].begin();
             vIter != basis[i].end(); vIter++) {
            tot_len += sqEdgeLen[*vIter];
        }
    }
    return tot_len;
}

void CycleLocalOptimization_bdry(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles,
                                 std::vector<std::set<int> > &out_v_basis_loops,
                                 std::vector<std::set<int> > &out_h_basis_loops,
                                 std::set<int> extraVertices,
                                 const float fScaleRatio) {//
    _simpGraph_vec locMeshGraph;
    //
    edge_annotation_computing edge_anno_proxy;
    //
    //
    ComputeAnnotation(locMesh, edge_anno_proxy.edge_annotations, OrientTriangles);
    edge_anno_proxy.vec_size = (int) edge_anno_proxy.edge_annotations[0].size();
    //
    if (edge_anno_proxy.vec_size == 0) {
        std::cout << "TRIVIAL LOOPS" << std::endl;
        exit(0);
    }
    std::vector<float> sqEdgeLen;
    AssembleMeshGraph(locMeshGraph, locMesh, sqEdgeLen, extraVertices, fScaleRatio);
    //
    edge_anno_proxy.SetMeshPtr(&locMesh);
    //edge_anno_proxy.ReadEdgeAnnotationsFromFile(argv[2]);
//{//
//	//
//	edge_annotation_computing edge_anno_proxy;
//	//
//	//
//	LoadData(argv[1], edge_anno_proxy);
//	//
//	std::set<int> extraVertices;
//	ReadBoundaries(argv[argc - 1], extraVertices);
//	//
//	std::vector<float> sqEdgeLen;
//	AssembleMeshGraph(meshGraph, cycMesh, sqEdgeLen, extraVertices);
//	//
//	edge_anno_proxy.SetMeshPtr(&cycMesh);
//	//edge_anno_proxy.ReadEdgeAnnotationsFromFile(argv[2]);
    //
    std::vector<std::vector<int> > basis_h_loops;
    std::vector<std::vector<int> > basis_v_loops;
    //
    std::vector<Annotation_Type> basis_h_annotations;
    std::vector<Annotation_Type> basis_v_annotations;
    //
    std::vector<int> base_pt_vec;
    //
    LoadCycleData(reebgraph, basis_h_loops, basis_v_loops, base_pt_vec);
    //ReadBasisLoops(argv[4], basis_h_loops);
    //ReadBasisLoops(argv[5], basis_v_loops);

    edge_anno_proxy.CheckTwoGroupVectorOrthogonality(basis_v_loops, basis_v_annotations, basis_h_loops,
                                                     basis_h_annotations);
    //
    canonical_loops_computing cano_loop_proxy;
    //
    //
    //std::set<int> base_pt_set;
    //for (unsigned int i = 0; i < basis_v_loops.size(); i++)
    //{
    //	for (unsigned int j = 0; j < basis_v_loops[i].size(); j++)
    //	{
    //		base_pt_set.insert(mesh.vecEdge[basis_v_loops[i][j]].v0);
    //		base_pt_set.insert(mesh.vecEdge[basis_v_loops[i][j]].v1);
    //	}
    //}
    //std::vector<int> base_pt_vec(base_pt_set.begin(), base_pt_set.end());//
    //
    //std::vector<int> base_pt_vec;
    //ReadBasisCriticalPoints(argv[3], base_pt_vec);
    std::set<int> base_pt_set(base_pt_vec.begin(), base_pt_vec.end());
    //
    //std::cout << "old base_pt_vec " << base_pt_vec.size() << std::endl;
    //
    cano_loop_proxy.SetGraphPtr(&locMeshGraph);
    cano_loop_proxy.SetEdgeWeights(sqEdgeLen);
    cano_loop_proxy.SetBasePointArray(base_pt_vec);
    cano_loop_proxy.SetMeshPtr(&locMesh);
    cano_loop_proxy.SetMeshGenus(edge_anno_proxy.vec_size / 2);
    //
    //cano_loop_proxy.ResetBaseVertices(basis_h_loops, edge_anno_proxy);
    //std::cout << "new base_pt_vec " << cano_loop_proxy.base_pt_index.size() << std::endl;
    //
    //bool bPartialBasis = false;
    std::cout << "Time for shortening handle and tunnel loops : " << std::endl;
    float cur_basis_h_len = BasisWeight(basis_h_loops, sqEdgeLen);
    float cur_basis_v_len = BasisWeight(basis_v_loops, sqEdgeLen);
    std::vector<std::set<int> > tmp_short_h_basis_loop;
    std::vector<std::set<int> > tmp_short_v_basis_loop;
    {
        // starting measuring the time
        boost::progress_timer t;
        //
        std::vector<std::set<int> > cano_loops;
        std::vector<int> non_tree_edges;
        std::vector<Annotation_Type> cano_loop_annotations;
        std::set<std::pair<float, int>, myFloatIntPairLessThan> sorted_cano_loops;
        std::vector<int> non_tree_edges_for_sorted_cano_loops;
        //
        float pre_basis_h_len = 0.0;
        float pre_basis_v_len = 0.0;


        //
        int iterNumber = 0;
        // optimize the all loops together
        do {
            tmp_short_h_basis_loop.clear();
            tmp_short_h_basis_loop.reserve(basis_h_loops.size());
            tmp_short_v_basis_loop.clear();
            tmp_short_v_basis_loop.reserve(basis_v_loops.size());
            base_pt_vec.clear();
            //
            cano_loop_proxy.compute_canonical_loops_fibo(cano_loops, non_tree_edges, cano_loop_annotations,
                                                         sorted_cano_loops, edge_anno_proxy);
            //
            if (edge_anno_proxy.ComputeShortestCanonicalLoops(basis_h_loops, basis_h_annotations, cano_loops,
                                                              non_tree_edges,
                                                              cano_loop_annotations, sorted_cano_loops,
                                                              tmp_short_h_basis_loop,
                                                              non_tree_edges_for_sorted_cano_loops)) {
                std::cout << ".";//ompletly use short loops" << std::endl;
            } else {
                std::cout << "*"; //artially use short loops" << std::endl;
                //bPartialBasis = true;
            }
            /**/
            cano_loop_proxy.base_pt_index.clear();//
            cano_loop_proxy.base_pt_index.reserve(non_tree_edges_for_sorted_cano_loops.size() * 3);
            for (unsigned int i = 0; i < non_tree_edges_for_sorted_cano_loops.size(); i++) {
                int edge_idx = non_tree_edges_for_sorted_cano_loops[i];
                //std::cout << edge_idx << std::endl;
                if ((edge_idx >= 0) && base_pt_set.find(locMesh.vecEdge[edge_idx].v0) == base_pt_set.end()) {
                    base_pt_set.insert(locMesh.vecEdge[edge_idx].v0);
                    cano_loop_proxy.base_pt_index.push_back(locMesh.vecEdge[edge_idx].v0);
                }
                //if (base_pt_set.find(mesh.vecEdge[edge_idx].v1) == base_pt_set.end())
                //{
                //	base_pt_set.insert(mesh.vecEdge[edge_idx].v1);
                //	cano_loop_proxy.base_pt_index.push_back(mesh.vecEdge[edge_idx].v1);
                //}
            }


            if (edge_anno_proxy.ComputeShortestCanonicalLoops(basis_v_loops, basis_v_annotations, cano_loops,
                                                              non_tree_edges,
                                                              cano_loop_annotations, sorted_cano_loops,
                                                              tmp_short_v_basis_loop,
                                                              non_tree_edges_for_sorted_cano_loops)) {
                std::cout << "."; //ompletly use short loops" << std::endl;
            } else {
                std::cout << "*"; //artially use short loops" << std::endl;
                //bPartialBasis = true;
            }
            //
            pre_basis_v_len = cur_basis_v_len;
            pre_basis_h_len = cur_basis_h_len;
            cur_basis_v_len = BasisWeight(tmp_short_v_basis_loop, sqEdgeLen);
            cur_basis_h_len = BasisWeight(tmp_short_h_basis_loop, sqEdgeLen);
            //
            for (unsigned int i = 0; i < non_tree_edges_for_sorted_cano_loops.size(); i++) {
                int edge_idx = non_tree_edges_for_sorted_cano_loops[i];
                //std::cout << edge_idx << std::endl;
                if ((edge_idx >= 0) && base_pt_set.find(locMesh.vecEdge[edge_idx].v0) == base_pt_set.end()) {
                    base_pt_set.insert(locMesh.vecEdge[edge_idx].v0);
                    cano_loop_proxy.base_pt_index.push_back(locMesh.vecEdge[edge_idx].v0);
                }
                //if (base_pt_set.find(mesh.vecEdge[edge_idx].v1) == base_pt_set.end())
                //{
                //	base_pt_set.insert(mesh.vecEdge[edge_idx].v1);
                //	cano_loop_proxy.base_pt_index.push_back(mesh.vecEdge[edge_idx].v1);
                //}
            }
            //std::cout << iterNumber << "\t";//<< " new base pt set size " << base_pt_set.size() << std::endl;
            //std::cout << iterNumber << " new base pt size " << cano_loop_proxy.base_pt_index.size() << std::endl;
            iterNumber++;
            if (cano_loop_proxy.base_pt_index.empty())
                break;
            if (iterNumber > 2000)
                break;

        } while (abs(cur_basis_h_len - pre_basis_h_len) > 0.000001 ||
                 abs(cur_basis_v_len - pre_basis_v_len) > 0.000001);
        //
        std::cout << std::endl;
    }
    //
    //
    out_h_basis_loops = tmp_short_h_basis_loop;
    //
    out_v_basis_loops = tmp_short_v_basis_loop;
    return;
}

void CycleGlobalOptimization(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles,
                             std::vector<std::set<int> > &out_v_basis_loops,
                             std::vector<std::set<int> > &out_h_basis_loops,
                             const float fScaleRatio) {//
    _simpGraph_vec locMeshGraph;
    //
    edge_annotation_computing edge_anno_proxy;
    //
    //
    ComputeAnnotation(locMesh, edge_anno_proxy.edge_annotations, OrientTriangles);
    edge_anno_proxy.vec_size = (int) edge_anno_proxy.edge_annotations[0].size();
    //
    if (edge_anno_proxy.vec_size == 0) {
        std::cout << "TRIVIAL LOOPS" << std::endl;
        exit(0);
    }
    std::vector<float> sqEdgeLen;
    AssembleMeshGraph(locMeshGraph, locMesh, sqEdgeLen, fScaleRatio);
    //
    edge_anno_proxy.SetMeshPtr(&locMesh);
    //edge_anno_proxy.ReadEdgeAnnotationsFromFile(argv[2]);
//{
//	//
//	edge_annotation_computing edge_anno_proxy;
//	//
//	//std::vector<int> meshTriangles;
//	LoadData(argv[1], edge_anno_proxy);
//	//
//	std::vector<float> sqEdgeLen;
//	AssembleMeshGraph(meshGraph, cycMesh, sqEdgeLen);
//	//
//	edge_anno_proxy.SetMeshPtr(&cycMesh);
//	//edge_anno_proxy.ReadEdgeAnnotationsFromFile(argv[2]);
//	//
//	//std::vector<int> fvs;
//	//ComputeFeedbackVertexSet(edge_anno_proxy.edge_annotations, fvs);

    std::vector<std::vector<int> > basis_h_loops;
    std::vector<std::vector<int> > basis_v_loops;
    std::vector<Annotation_Type> basis_h_annotations;
    std::vector<Annotation_Type> basis_v_annotations;
    //
    std::vector<int> base_pt_vec;
    //
    LoadCycleData(reebgraph, basis_h_loops, basis_v_loops, base_pt_vec);
    //
    edge_anno_proxy.CheckTwoGroupVectorOrthogonality(basis_v_loops, basis_v_annotations, basis_h_loops,
                                                     basis_h_annotations);
    //
    canonical_loops_computing cano_loop_proxy;
    //

    //ReadBasisCriticalPoints(argv[3], base_pt_vec);
    std::set<int> base_pt_set(base_pt_vec.begin(), base_pt_vec.end());
    //
    std::cout << "old base_pt_vec " << base_pt_vec.size() << std::endl;
    //
    cano_loop_proxy.SetGraphPtr(&locMeshGraph);
    cano_loop_proxy.SetEdgeWeights(sqEdgeLen);
    cano_loop_proxy.SetBasePointArray(base_pt_vec);
    cano_loop_proxy.SetMeshPtr(&locMesh);
    cano_loop_proxy.SetMeshGenus(edge_anno_proxy.vec_size / 2);
    //
    //cano_loop_proxy.ResetBaseVertices(basis_h_loops, edge_anno_proxy);
    //std::cout << "new base_pt_vec " << cano_loop_proxy.base_pt_index.size() << std::endl;
    //
    std::cout << "in cano loops computation" << std::endl;
    // starting measuring time
    boost::progress_timer t;
    //
    std::vector<std::set<int> > cano_loops;
    std::vector<int> non_tree_edges;
    std::vector<Annotation_Type> cano_loop_annotations;
    std::set<std::pair<float, int>, myFloatIntPairLessThan> sorted_cano_loops;
    std::vector<int> non_tree_edges_for_sorted_cano_loops;
    //
    float pre_basis_h_len = 0.0;
    float pre_basis_v_len = 0.0;
    float cur_basis_h_len = BasisWeight(basis_h_loops, sqEdgeLen);
    float cur_basis_v_len = BasisWeight(basis_v_loops, sqEdgeLen);
    std::vector<std::set<int> > tmp_short_h_basis_loop;
    std::vector<std::set<int> > tmp_short_v_basis_loop;
    //
    int iterNumber = 0;
    // compute the vertical loops first
    do {
        tmp_short_h_basis_loop.clear();
        tmp_short_h_basis_loop.reserve(basis_h_loops.size());
        tmp_short_v_basis_loop.clear();
        tmp_short_v_basis_loop.reserve(basis_v_loops.size());
        base_pt_vec.clear();
        //
        cano_loop_proxy.compute_canonical_loops_fibo(cano_loops, non_tree_edges, cano_loop_annotations,
                                                     sorted_cano_loops, edge_anno_proxy);
        //
        if (edge_anno_proxy.ComputeShortestCanonicalLoops(basis_h_loops, basis_h_annotations, cano_loops,
                                                          non_tree_edges,
                                                          cano_loop_annotations, sorted_cano_loops,
                                                          tmp_short_h_basis_loop,
                                                          non_tree_edges_for_sorted_cano_loops)) {
            std::cout << "c \t ";//ompletly use short loops" << std::endl;
        } else {
            std::cout << "p \t"; //artially use short loops" << std::endl;
        }
        /**/
        cano_loop_proxy.base_pt_index.clear();//
        cano_loop_proxy.base_pt_index.reserve(non_tree_edges_for_sorted_cano_loops.size() * 3);
        for (unsigned int i = 0; i < non_tree_edges_for_sorted_cano_loops.size(); i++) {
            int edge_idx = non_tree_edges_for_sorted_cano_loops[i];
            //std::cout << edge_idx << std::endl;
            if ((edge_idx >= 0) && base_pt_set.find(locMesh.vecEdge[edge_idx].v0) == base_pt_set.end()) {
                base_pt_set.insert(locMesh.vecEdge[edge_idx].v0);
                cano_loop_proxy.base_pt_index.push_back(locMesh.vecEdge[edge_idx].v0);
            }
            //if (base_pt_set.find(mesh.vecEdge[edge_idx].v1) == base_pt_set.end())
            //{
            //	base_pt_set.insert(mesh.vecEdge[edge_idx].v1);
            //	cano_loop_proxy.base_pt_index.push_back(mesh.vecEdge[edge_idx].v1);
            //}
        }


        if (edge_anno_proxy.ComputeShortestCanonicalLoops(basis_v_loops, basis_v_annotations, cano_loops,
                                                          non_tree_edges,
                                                          cano_loop_annotations, sorted_cano_loops,
                                                          tmp_short_v_basis_loop,
                                                          non_tree_edges_for_sorted_cano_loops)) {
            std::cout << "c \t"; //ompletly use short loops" << std::endl;
        } else {
            std::cout << "p \t"; //artially use short loops" << std::endl;
        }
        //
        pre_basis_v_len = cur_basis_v_len;
        pre_basis_h_len = cur_basis_h_len;
        cur_basis_v_len = BasisWeight(tmp_short_v_basis_loop, sqEdgeLen);
        cur_basis_h_len = BasisWeight(tmp_short_h_basis_loop, sqEdgeLen);
        //
        for (unsigned int i = 0; i < non_tree_edges_for_sorted_cano_loops.size(); i++) {
            int edge_idx = non_tree_edges_for_sorted_cano_loops[i];
            //std::cout << edge_idx << std::endl;
            if ((edge_idx >= 0) && base_pt_set.find(locMesh.vecEdge[edge_idx].v0) == base_pt_set.end()) {
                base_pt_set.insert(locMesh.vecEdge[edge_idx].v0);
                cano_loop_proxy.base_pt_index.push_back(locMesh.vecEdge[edge_idx].v0);
            }
            //if (base_pt_set.find(mesh.vecEdge[edge_idx].v1) == base_pt_set.end())
            //{
            //	base_pt_set.insert(mesh.vecEdge[edge_idx].v1);
            //	cano_loop_proxy.base_pt_index.push_back(mesh.vecEdge[edge_idx].v1);
            //}
        }
        std::cout << iterNumber << "\t";//<< " new base pt set size " << base_pt_set.size() << std::endl;
        //std::cout << iterNumber << " new base pt size " << cano_loop_proxy.base_pt_index.size() << std::endl;
        iterNumber++;
        if (cano_loop_proxy.base_pt_index.empty())
            break;
        if (iterNumber > 2000)
            break;

    } while (abs(cur_basis_h_len - pre_basis_h_len) > 0.000001 || abs(cur_basis_v_len - pre_basis_v_len) > 0.000001);
    //
    // use optimzed vertical loops to optimize horizontal loops
    std::set<int> tmpSet;
    cano_loop_proxy.base_pt_index.clear();
    for (unsigned int i = 0; i < tmp_short_v_basis_loop.size(); i++) {
        for (std::set<int>::iterator sIter = tmp_short_v_basis_loop[i].begin();
             sIter != tmp_short_v_basis_loop[i].end(); sIter++) {
            base_pt_set.insert(locMesh.vecEdge[*sIter].v0);
            base_pt_set.insert(locMesh.vecEdge[*sIter].v1);
            tmpSet.insert(locMesh.vecEdge[*sIter].v0);
        }
    }
    //
    for (unsigned int i = 0; i < tmp_short_h_basis_loop.size(); i++) {
        for (std::set<int>::iterator sIter = tmp_short_h_basis_loop[i].begin();
             sIter != tmp_short_h_basis_loop[i].end(); sIter++) {
            base_pt_set.insert(locMesh.vecEdge[*sIter].v0);
            base_pt_set.insert(locMesh.vecEdge[*sIter].v1);
            tmpSet.insert(locMesh.vecEdge[*sIter].v0);
        }
    }
    //
    cano_loop_proxy.base_pt_index.assign(tmpSet.begin(), tmpSet.end());
    //
    tmp_short_h_basis_loop.clear();
    tmp_short_h_basis_loop.reserve(basis_h_loops.size());
    tmp_short_v_basis_loop.clear();
    tmp_short_v_basis_loop.reserve(basis_v_loops.size());
    //
    cano_loop_proxy.compute_canonical_loops_fibo(cano_loops, non_tree_edges, cano_loop_annotations, sorted_cano_loops,
                                                 edge_anno_proxy);
    //

    if (edge_anno_proxy.ComputeShortestCanonicalLoops(basis_v_loops, basis_v_annotations, cano_loops, non_tree_edges,
                                                      cano_loop_annotations, sorted_cano_loops, tmp_short_v_basis_loop,
                                                      non_tree_edges_for_sorted_cano_loops)) {
        std::cout << "c \t"; //ompletly use short loops" << std::endl;
    } else {
        std::cout << "p \t"; //artially use short loops" << std::endl;
    }
    //
    if (edge_anno_proxy.ComputeShortestCanonicalLoops(basis_h_loops, basis_h_annotations, cano_loops, non_tree_edges,
                                                      cano_loop_annotations, sorted_cano_loops, tmp_short_h_basis_loop,
                                                      non_tree_edges_for_sorted_cano_loops)) {
        std::cout << "completly use short loops" << std::endl;
    } else {
        std::cout << "partially use short loops" << std::endl;
    }
    //
    out_h_basis_loops = tmp_short_h_basis_loop;
    //
    out_v_basis_loops = tmp_short_v_basis_loop;
    return;
}
//
//int mainx(int argc, char** argv)
//{
//	int optMethodOption;
//	stringstream sstr(stringstream::in | stringstream::out);
//	if (argc == 8)
//	{
//		sstr.str(argv[7]);
//		sstr >> optMethodOption;
//		sstr.clear();
//		std::cout << optMethodOption << std::endl;
//		if (optMethodOption)
//		{
//			//GlobalOptimization(argc, argv);
//		}
//		else
//		{
//			LocalOptimization(argc, argv);
//		}
//	}
//	else
//	{
//		if (argc == 9)
//			LocalOptimization_bdry(argc, argv);
//		else
//			LocalOptimization(argc, argv);
//	}
//	return 1;
//	//
//}

void CycleLocalOptimization(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles,
                            std::vector<std::set<int> > &out_v_basis_loops,
                            std::vector<std::set<int> > &out_h_basis_loops,
                            const float fScaleRatio) {//
    _simpGraph_vec locMeshGraph;
    //
    edge_annotation_computing edge_anno_proxy;
    //
    //std::cout << "tris size" << OrientTriangles.size() << std::endl;
    //
    ComputeAnnotation(locMesh, edge_anno_proxy.edge_annotations, OrientTriangles);
    edge_anno_proxy.vec_size = (int) edge_anno_proxy.edge_annotations[0].size();
    //
    if (edge_anno_proxy.vec_size == 0) {
        std::cout << "TRIVIAL LOOPS" << std::endl;
        exit(0);
    }
    //std::cout << "before assembe graph " << edge_anno_proxy.vec_size << std::endl;
    std::vector<float> sqEdgeLen;
    AssembleMeshGraph(locMeshGraph, locMesh, sqEdgeLen, fScaleRatio);
    //
    edge_anno_proxy.SetMeshPtr(&locMesh);
    //edge_anno_proxy.ReadEdgeAnnotationsFromFile(argv[2]);
    //
    std::vector<int> basis_h_LowestOnePos;
    std::vector<int> basis_v_LowestOnePos;
    std::vector<Annotation_Type> reduced_basis_h_annotations;
    std::vector<Annotation_Type> reduced_basis_v_annotations;
    //
    std::vector<std::vector<int> > basis_h_loops;
    std::vector<std::vector<int> > basis_v_loops;
    std::vector<Annotation_Type> basis_h_annotations;
    std::vector<Annotation_Type> basis_v_annotations;
    //
    std::vector<int> base_pt_vec;
    //
    LoadCycleData(reebgraph, basis_h_loops, basis_v_loops, base_pt_vec);

    edge_anno_proxy.CheckTwoGroupVectorOrthogonality(basis_v_loops, basis_v_annotations, reduced_basis_v_annotations,
                                                     basis_v_LowestOnePos,
                                                     basis_h_loops, basis_h_annotations, reduced_basis_h_annotations,
                                                     basis_h_LowestOnePos);
    //
    canonical_loops_computing cano_loop_proxy;
    //
    //
    //std::set<int> base_pt_set;
    //for (unsigned int i = 0; i < basis_v_loops.size(); i++)
    //{
    //	for (unsigned int j = 0; j < basis_v_loops[i].size(); j++)
    //	{
    //		base_pt_set.insert(mesh.vecEdge[basis_v_loops[i][j]].v0);
    //		base_pt_set.insert(mesh.vecEdge[basis_v_loops[i][j]].v1);
    //	}
    //}
    //std::vector<int> base_pt_vec(base_pt_set.begin(), base_pt_set.end());//
    //

    //ReadBasisCriticalPoints(argv[3], base_pt_vec);
    std::set<int> base_pt_set(base_pt_vec.begin(), base_pt_vec.end());
    //
    //std::cout << "old base_pt_vec " << base_pt_vec.size() << std::endl;
    //
    cano_loop_proxy.SetGraphPtr(&locMeshGraph);
    cano_loop_proxy.SetEdgeWeights(sqEdgeLen);
    cano_loop_proxy.SetBasePointArray(base_pt_vec);
    cano_loop_proxy.SetMeshPtr(&locMesh);
    cano_loop_proxy.SetMeshGenus(edge_anno_proxy.vec_size / 2);
    //
    //cano_loop_proxy.ResetBaseVertices(basis_h_loops, edge_anno_proxy);
    //std::cout << "new base_pt_vec " << cano_loop_proxy.base_pt_index.size() << std::endl;
    //
    std::cout << "Time for shortening handle and tunnel loops : " << std::endl;
    float cur_basis_h_len = BasisWeight(basis_h_loops, sqEdgeLen);
    float cur_basis_v_len = BasisWeight(basis_v_loops, sqEdgeLen);
    std::vector<std::set<int> > tmp_short_h_basis_loop;
    std::vector<std::set<int> > tmp_short_v_basis_loop;
    {
        // starting measuring the time
        boost::progress_timer t;
        //
        std::vector<std::set<int> > cano_loops;
        std::vector<int> non_tree_edges;
        std::vector<Annotation_Type> cano_loop_annotations;
        std::set<std::pair<float, int>, myFloatIntPairLessThan> sorted_cano_loops;
        std::vector<int> non_tree_edges_for_sorted_cano_loops;
        //
        float pre_basis_h_len = 0.0;
        float pre_basis_v_len = 0.0;


        //
        int iterNumber = 0;
        // optimize the all loops together
        do {
            tmp_short_h_basis_loop.clear();
            tmp_short_h_basis_loop.reserve(basis_h_loops.size());
            tmp_short_v_basis_loop.clear();
            tmp_short_v_basis_loop.reserve(basis_v_loops.size());
            base_pt_vec.clear();
            //

            cano_loop_proxy.compute_canonical_loops_fibo(cano_loops, non_tree_edges, cano_loop_annotations,
                                                         sorted_cano_loops, edge_anno_proxy);

            //
            //std::cout << "in cano 1" << std::endl;
            if (edge_anno_proxy.ComputeShortestCanonicalLoops(basis_h_loops, basis_h_annotations,
                                                              reduced_basis_h_annotations, basis_h_LowestOnePos,
                                                              cano_loops, non_tree_edges,
                                                              cano_loop_annotations, sorted_cano_loops,
                                                              tmp_short_h_basis_loop,
                                                              non_tree_edges_for_sorted_cano_loops)) {
                std::cout << ".";//ompletly use short loops" << std::endl;
            } else {
                std::cout << "*"; //artially use short loops" << std::endl;
            }
            //std::cout << "out cano 1" << std::endl;
            /**/
            cano_loop_proxy.base_pt_index.clear();//
            cano_loop_proxy.base_pt_index.reserve(non_tree_edges_for_sorted_cano_loops.size() * 3);
            for (unsigned int i = 0; i < non_tree_edges_for_sorted_cano_loops.size(); i++) {
                int edge_idx = non_tree_edges_for_sorted_cano_loops[i];
                //std::cout << edge_idx << std::endl;
                if ((edge_idx >= 0) && base_pt_set.find(locMesh.vecEdge[edge_idx].v0) == base_pt_set.end()) {
                    base_pt_set.insert(locMesh.vecEdge[edge_idx].v0);
                    cano_loop_proxy.base_pt_index.push_back(locMesh.vecEdge[edge_idx].v0);
                }
                //if (base_pt_set.find(mesh.vecEdge[edge_idx].v1) == base_pt_set.end())
                //{
                //	base_pt_set.insert(mesh.vecEdge[edge_idx].v1);
                //	cano_loop_proxy.base_pt_index.push_back(mesh.vecEdge[edge_idx].v1);
                //}
            }

            //std::cout << "in cano 2" << std::endl;
            if (edge_anno_proxy.ComputeShortestCanonicalLoops(basis_v_loops, basis_v_annotations,
                                                              reduced_basis_v_annotations, basis_v_LowestOnePos,
                                                              cano_loops, non_tree_edges,
                                                              cano_loop_annotations, sorted_cano_loops,
                                                              tmp_short_v_basis_loop,
                                                              non_tree_edges_for_sorted_cano_loops)) {
                std::cout << "."; //ompletly use short loops" << std::endl;
            } else {
                std::cout << "*"; //artially use short loops" << std::endl;
            }
            //std::cout << "out cano 2" << std::endl;
            //
            pre_basis_v_len = cur_basis_v_len;
            pre_basis_h_len = cur_basis_h_len;
            cur_basis_v_len = BasisWeight(tmp_short_v_basis_loop, sqEdgeLen);
            cur_basis_h_len = BasisWeight(tmp_short_h_basis_loop, sqEdgeLen);
            //
            for (unsigned int i = 0; i < non_tree_edges_for_sorted_cano_loops.size(); i++) {
                int edge_idx = non_tree_edges_for_sorted_cano_loops[i];
                //std::cout << edge_idx << std::endl;
                if ((edge_idx >= 0) && base_pt_set.find(locMesh.vecEdge[edge_idx].v0) == base_pt_set.end()) {
                    base_pt_set.insert(locMesh.vecEdge[edge_idx].v0);
                    cano_loop_proxy.base_pt_index.push_back(locMesh.vecEdge[edge_idx].v0);
                }
                //if (base_pt_set.find(mesh.vecEdge[edge_idx].v1) == base_pt_set.end())
                //{
                //	base_pt_set.insert(mesh.vecEdge[edge_idx].v1);
                //	cano_loop_proxy.base_pt_index.push_back(mesh.vecEdge[edge_idx].v1);
                //}
            }
            //std::cout << iterNumber << "\t";//<< " new base pt set size " << base_pt_set.size() << std::endl;
            //std::cout << iterNumber << " new base pt size " << cano_loop_proxy.base_pt_index.size() << std::endl;
            iterNumber++;
            if (cano_loop_proxy.base_pt_index.empty())
                break;
            if (iterNumber > 2000)
                break;

        } while (fabs(cur_basis_h_len - pre_basis_h_len) > 0.000001 ||
                 fabs(cur_basis_v_len - pre_basis_v_len) > 0.000001);
        //
        std::cout << std::endl;
    }
    //
    out_h_basis_loops = tmp_short_h_basis_loop;
    //
    //
    out_v_basis_loops = tmp_short_v_basis_loop;
    //
    return;
}
