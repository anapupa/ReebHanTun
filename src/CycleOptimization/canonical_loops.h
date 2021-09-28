/*
(c) 2012 Fengtao Fan
*/
#ifndef _MY_CANONICAL_LOOPS_H_
#define _MY_CANONICAL_LOOPS_H_

#include <vector>
#include <iostream>
#include "SimpleGraph.h"
#include "SimpleMesh.h"
#include "edge_annotations_gauss.h"
#include "DijkstraAlgorithm.h"

#ifndef myFloatIntPairLessThan_
#define myFloatIntPairLessThan_
struct myFloatIntPairLessThan {
    bool operator()(const std::pair<float, int> &lhs, const std::pair<float, int> &rhs)
    {
        bool ret = false;
        if (lhs.first < rhs.first)
            ret = true;
        else
            if (lhs.first == rhs.first)
            {
                if (lhs.second < rhs.second)
                    ret = true;
            }
        return ret;
    }
};
#endif

#ifndef My_Annotation_Type
#define My_Annotation_Type
typedef std::vector<int> Annotation_Type;
#endif

//
class canonical_loops_computing {
    //
public:
    // constructor and assign operator
    canonical_loops_computing() {
        Null_all();
    }

    canonical_loops_computing(const canonical_loops_computing &rhs) {
        //
        std::cout << "NEVER USE THIS FUNCTION" << std::endl;
        //
    }

    canonical_loops_computing &operator=(const canonical_loops_computing &rhs) {
        if (this == &rhs)
            return *this;
        //
        std::cout << "NEVER USE THIS FUNCTION" << std::endl;
        //
        return *this;
    }

    // deconstructor and data clean functioin
    void Null_all() {
        _graphPtr = NULL;
    }

    ~canonical_loops_computing() {
        Null_all();
    }

    // actual computation
    void SetGraphPtr(_simpGraph_vec *in_graph_ptr) {
        _graphPtr = in_graph_ptr;
    }

    void SetBasePointArray(std::vector<int> &in_base_pt) {
        base_pt_index = in_base_pt;
    }

    void SetMeshPtr(_SimpleMesh *inMeshPtr) {
        _meshPtr = inMeshPtr;
    }

    void SetMeshGenus(const int inGenus) {
        _genus = inGenus;
    }

    void PreprocessingForEdges(std::vector<std::vector<float> > &adjListWithWeight);

    void SetEdgeWeights(std::vector<float> &in_weight);

    //
    void Dijkstra_shortest_path_tree(const int source,
                                     _simpGraph_vec &graph,
                                     std::vector<std::pair<std::set<int>, float> > &path_to_source_in_tree,
                                     std::vector<float> &weight,
                                     std::set<std::pair<float, int>, myFloatIntPairLessThan> &remaining_vertex,
                                     std::vector<std::pair<int, float> > &parents);

    void Dijkstra_shortest_path_tree(const int source,
                                     _simpGraph_vec &graph,
                                     std::vector<float> &weight,
                                     std::set<std::pair<float, int>, myFloatIntPairLessThan> &remaining_vertex,
                                     std::vector<int> &parents_vertex,
                                     std::vector<int> &parents_edge,
                                     std::vector<float> &distance);

    //
    void compute_canonical_loops(std::vector<std::set<int> > &canonical_loops,
                                 std::vector<Annotation_Type> &canonical_loops_annotation,
                                 std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                 edge_annotation_computing &edge_annotations_proxy);

    //
    void compute_canonical_loops_fibo(std::vector<std::set<int> > &canonical_loops,
                                      std::vector<int> &non_tree_edges,
                                      std::vector<Annotation_Type> &canonical_loops_annotation,
                                      std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                      edge_annotation_computing &edge_annotations_proxy);

    //
    void collect_non_tree_edge_cycles(std::vector<std::pair<std::set<int>, float> > &path_to_source_in_tree,
                                      std::vector<std::set<int> > &canonical_loops,
                                      std::vector<Annotation_Type> &canonical_loops_annotation,
                                      std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                      edge_annotation_computing &edge_annotations_proxy);

    void collect_non_tree_edge_cycles(const int src_node,
                                      std::vector<int> &parents_vertex,
                                      std::vector<int> &parents_edge,
                                      std::vector<float> &distance,
                                      std::vector<unsigned char> &colors,
                                      std::vector<Annotation_Type> &accumulated_parent_annotations,
                                      std::vector<std::set<int> > &canonical_loops,
                                      std::vector<int> &non_tree_edges,
                                      std::vector<Annotation_Type> &canonical_loops_annotation,
                                      std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                      edge_annotation_computing &edge_annotations_proxy);

    bool Find_shortes_cycle_with_same_annotation(const int src_node,
                                                 std::vector<int> &parents_vertex,
                                                 std::vector<int> &parents_edge,
                                                 std::vector<int> &distance,
                                                 std::vector<unsigned char> &colors,
                                                 std::vector<Annotation_Type> &accumulated_parent_annotations,
                                                 const Annotation_Type &org_loop_annotation,
                                                 const int org_loop_weight,
                                                 edge_annotation_computing &edge_annotations_proxy,
                                                 std::set<int> &base_pt_set_unique);

    void ExtractShortestBasisFromCycleSets(const int genus,
                                           std::vector<int> &shortest_basis_indices,
                                           std::vector<Annotation_Type> &canonical_loops_annotation,
                                           std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                           edge_annotation_computing &edge_annotations_proxy);

    void
    ResetBaseVertices(std::vector<std::vector<int> > &basis_loops, edge_annotation_computing &edge_annotations_proxy);
    //
    //bool is_cycle_nontrivial(std::set<int> &cycle, Annotation_Type &cur_edge_annotation);
public:
    _simpGraph_vec *_graphPtr;
    _SimpleMesh *_meshPtr;
    //
    std::vector<int> base_pt_index;
    //
    std::vector<std::vector<std::pair<int, int> > > pair_to_edge_index_mapping;
    //
    std::vector<bool> edge_flag_list;
    //
    std::vector<float> edge_weight_vec;
    //
    std::vector<std::vector<std::vector<int> > > canonical_loops;
    //
    int _genus;
};

#endif //_MY_CANONICAL_LOOPS_H_
