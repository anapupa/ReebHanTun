/*
(c) 2012 Fengtao Fan
*/
#ifndef _My_EDGE_ANNOTATION_COMPUTING_GAUSS_
#define _My_EDGE_ANNOTATION_COMPUTING_GAUSS_

#include "SimpleMesh.h"
#include <vector>
#include <set>

#ifndef My_Annotation_Type
#define My_Annotation_Type
typedef std::vector<char> Annotation_Type;
#endif
//
#ifndef myFloatIntPairLessThan_
#define myFloatIntPairLessThan_

struct myFloatIntPairLessThan {
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
};

#endif

struct Annotation_Type_LessThan {
    bool operator()(const Annotation_Type &lhs, const Annotation_Type &rhs) {
        bool ret = false;
        for (unsigned int i = 0; i < lhs.size(); i++) {
            if (lhs[i] < rhs[i]) {
                ret = true;
                break;
            } else {
                if (lhs[i] > rhs[i]) {
                    ret = false;
                    break;
                }
            }
        }
        return ret;
    }
};

class edge_annotation_computing {
public:
    edge_annotation_computing() {

    }

    ~edge_annotation_computing() {
        meshPtr = NULL;
        //graphPtr = NULL;
    }

    //
    void SetMeshPtr(_SimpleMesh *inMeshPtr) {
        meshPtr = inMeshPtr;
    }

    void SetEdgeSize(const int n) {
        //edge_size = n;
    }

    //
    bool IndependenceCheck(std::vector<Annotation_Type> &bReducedVec,
                           std::vector<int> &lowestOnePtr,
                           Annotation_Type &objVector);

    //
    bool ReduceMatrix(std::vector<Annotation_Type> &annMatrix, std::vector<int> &lowestOnePos);

    bool GaussianElimination(std::vector<Annotation_Type> &vec);

    //
    void ReadEdgeAnnotationsFromFile(const char *file_name);

    void AssignEdgeAnnotations(std::vector<Annotation_Type> &existEdgeVec) {
        if (!edge_annotations.empty()) {
            edge_annotations.clear();
        }
        edge_annotations = existEdgeVec;
    }

    //
    bool CheckZeroAnnotationsAndReturnAnnotation(std::set<int> &cycle, Annotation_Type &res);

    bool CheckZeroAnnotationsAndReturnAnnotation(std::vector<int> &cycle, Annotation_Type &res);

    //
    bool ComputeShortestCanonicalLoops(std::vector<std::vector<int> > &basis_loops,
                                       std::vector<Annotation_Type> &basis_annotations,
                                       std::vector<std::set<int> > &canonical_loops,
                                       std::vector<int> &non_tree_edges,
                                       std::vector<Annotation_Type> &canonical_loops_annotation,
                                       std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                       std::vector<std::set<int> > &short_basis_loops,
                                       std::vector<int> &basis_non_tree_edges);

    //
    bool ComputeShortestCanonicalLoops(std::vector<std::vector<int> > &basis_loops,
                                       std::vector<Annotation_Type> &basis_annotations,
                                       std::vector<Annotation_Type> &reduced_basis_annotations,
                                       std::vector<int> &basis_LowestOnePos,
                                       std::vector<std::set<int> > &canonical_loops,
                                       std::vector<int> &non_tree_edges,
                                       std::vector<Annotation_Type> &canonical_loops_annotation,
                                       std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                       std::vector<std::set<int> > &short_basis_loops,
                                       std::vector<int> &basis_non_tree_edges);

    //
    bool CheckTwoGroupVectorOrthogonality(std::vector<std::vector<int> > &a_loops, std::vector<Annotation_Type> &a_ann,
                                          std::vector<std::vector<int> > &b_loops, std::vector<Annotation_Type> &b_ann);

    bool CheckTwoGroupVectorOrthogonality(std::vector<std::vector<int> > &a_loops, std::vector<Annotation_Type> &a_ann,
                                          std::vector<Annotation_Type> &reduced_a_ann, std::vector<int> &a_lowestOne,
                                          std::vector<std::vector<int> > &b_loops, std::vector<Annotation_Type> &b_ann,
                                          std::vector<Annotation_Type> &reduced_b_ann, std::vector<int> &b_lowestOne);

    void Add(const Annotation_Type &lhs, const Annotation_Type &rhs, Annotation_Type &res) {
        for (unsigned int i = 0; i < lhs.size(); i++) {
            res[i] = lhs[i] ^ rhs[i];
        }
    }

    void Times(const Annotation_Type &lhs, const Annotation_Type &rhs, Annotation_Type &res) {
        for (unsigned int i = 0; i < lhs.size(); i++) {
            res[i] = lhs[i] & rhs[i];
        }
    }

    bool IsZeroAnnotation(const Annotation_Type &ins) {
        int sum = 0;
        for (unsigned int i = 0; i < ins.size(); i++)
            sum += ins[i];
        return !(sum);
    }

    bool IsEqual(const Annotation_Type &lhs, const Annotation_Type &rhs) {
        bool ret = true;
        for (unsigned int i = 0; i < lhs.size(); i++) {
            if (lhs[i] != rhs[i]) {
                ret = false;
                break;
            }
        }
        return ret;
    }

public:
    /*
    */
    _SimpleMesh *meshPtr;
    //int edge_size;
    std::vector<Annotation_Type> edge_annotations;
    int vec_size;
    //
    //std::vector<Annotation_Type>	reducedHandleBasis;
    //std::vector<int>				handleLowestOnePos;
    ////
    //std::vector<Annotation_Type>	reducedTunnelBasis;
    //std::vector<int>				tunnelLowestOnePos;
};

#endif //GAUSS_

