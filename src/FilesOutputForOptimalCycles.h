/*
(c) 2012 Fengtao Fan
*/
#ifndef _MY_FilesOutputForOptimalCycles_H_
#define _MY_FilesOutputForOptimalCycles_H_

#include <vector>
#include <list>
#include <set>
#include <map>

#include "SimpleMesh.h"

class FilesOutputForOptimalCycles {
public:
    // constructor and assign operator
    FilesOutputForOptimalCycles() {
        Null_all();
    }

    FilesOutputForOptimalCycles(const FilesOutputForOptimalCycles &rhs) {
        //
        inMeshPtr = rhs.inMeshPtr;
        //
        vecInitVerticalLoops_Vertex_Ptr = rhs.vecInitVerticalLoops_Vertex_Ptr;
        vecInitVerticalLoops_Edge_Ptr = rhs.vecInitVerticalLoops_Edge_Ptr;
        //
        vecInitHorizontalLoops_Vertex_Ptr = rhs.vecInitHorizontalLoops_Vertex_Ptr;
        vecInitHorizontalLoops_Edge_Ptr = rhs.vecInitHorizontalLoops_Edge_Ptr;
        //
        coeffMatrix = rhs.coeffMatrix;
        //
    }

    FilesOutputForOptimalCycles &operator=(const FilesOutputForOptimalCycles &rhs) {
        if (this == &rhs)
            return *this;
        //
        inMeshPtr = rhs.inMeshPtr;
        //
        vecInitVerticalLoops_Vertex_Ptr = rhs.vecInitVerticalLoops_Vertex_Ptr;
        vecInitVerticalLoops_Edge_Ptr = rhs.vecInitVerticalLoops_Edge_Ptr;
        //
        vecInitHorizontalLoops_Vertex_Ptr = rhs.vecInitHorizontalLoops_Vertex_Ptr;
        vecInitHorizontalLoops_Edge_Ptr = rhs.vecInitHorizontalLoops_Edge_Ptr;
        //
        coeffMatrix = rhs.coeffMatrix;
        //
        return *this;
    }

    // deconstructor and data clean functioin
    void Null_all() {
        inMeshPtr = NULL;
        //
        vecInitVerticalLoops_Vertex_Ptr = NULL;
        vecInitVerticalLoops_Edge_Ptr = NULL;
        //
        vecInitHorizontalLoops_Vertex_Ptr = NULL;
        vecInitHorizontalLoops_Edge_Ptr = NULL;
        coeffMatrix = NULL;
        //
    }

    ~FilesOutputForOptimalCycles() {
        Null_all();
    }

    //
    // actual operations
    void SetLoopTypes(std::vector<bool> &inVertLoopType) {
        //initVerticalLoopsType.clear();
        //initVerticalLoopsType.assign(inVertLoopType.begin(), inVertLoopType.end());
        return;
    }

    void InitMeshPtr(_SimpleMesh *inPtr) {
        inMeshPtr = inPtr;
    }

    void InitLoopsPtr(std::vector<std::vector<int> > *inVec_V_Loop_V_ptr,
                      std::vector<std::vector<int> > *inVec_V_Loop_E_ptr,
                      std::vector<std::vector<int> > *inVec_H_Loop_V_ptr,
                      std::vector<std::vector<int> > *inVec_H_Loop_E_ptr) {
        vecInitVerticalLoops_Vertex_Ptr = inVec_V_Loop_V_ptr;
        vecInitVerticalLoops_Edge_Ptr = inVec_V_Loop_E_ptr;
        //
        vecInitHorizontalLoops_Vertex_Ptr = inVec_H_Loop_V_ptr;
        vecInitHorizontalLoops_Edge_Ptr = inVec_H_Loop_E_ptr;
    }

    //void InitCoeffMatrix(std::vector<std::vector<std::pair<int, int>>> *inCoeffMatrix)
    void InitCoeffMatrix(std::vector<std::list<std::pair<int, int> > > *inCoeffMatrix) {
        coeffMatrix = inCoeffMatrix;
    }

    void WriteOrientedComplex(const char *out_file_name);

    void WriteOrientedCycle(const int cycle_id, const char *out_file_name);

    void WriteCoefficientMatrix(const char *out_file_name);

    void WriteAllCyclesOut();

    void WriteAllCyclesOut(std::vector<bool> &inVertLoopType, const char *out_file_name);

    void WriteAllCyclesOut(std::vector<bool> &inVertLoopType, const char *out_file_name, const int loop_type_num);

    void WritePolygonOut(const char *out_file_name, std::vector<Vector3> &inPly);

    void ReadPolygonIn(const char *in_file_name, std::vector<Vector3> &wantPly);

    void WriteAllCyclesInOneFile(
            std::vector<std::vector<int> > &vertLoops,
            std::vector<std::vector<int> > &horiLoops,
            std::vector<bool> &initVerticalLoopsType,
            const char *out_file_name,
            const bool loop_type);

    void WriteAllCyclesInOneFile(std::vector<std::vector<int> > &vertLoops,
                                 std::vector<std::vector<int> > &horiLoops,
                                 const char *out_file_name);

    void WriteCriticalPoints(const char *out_file_name, std::set<int> &critSet, const bool bVertical);

    void WriteCriticalPoints(const char *out_file_name, std::vector<std::pair<int, int> > &critSet);

    /********************************/
    void OrderCycles(std::set<int> &edge_set, std::map<int, int> &vertex_ordering);

    bool OrderCycles(std::set<int> &edge_set, std::vector<int> &vertex_set);

    void WriteCyclesInformation(const char *out_file_name, std::vector<std::set<int> > &v_basis_loops,
                                std::vector<std::set<int> > &h_basis_loops);

    void WriteGeomviewListFormat(const char *out_file_name, std::vector<std::set<int> > &v_basis_loops,
                                 std::vector<std::set<int> > &h_basis_loops, std::vector<int> &vec_oriented_triangles,
                                 const int orgVertexSize, const int orgTriangleSize, const float fEnlargeFactor);

public:
    _SimpleMesh *inMeshPtr;
    //
    //std::vector<bool> initVerticalLoopsType; // true for vertical, and false for horizontal
    //
    std::vector<std::vector<int> > *vecInitVerticalLoops_Vertex_Ptr;
    std::vector<std::vector<int> > *vecInitVerticalLoops_Edge_Ptr;
    //
    std::vector<std::vector<int> > *vecInitHorizontalLoops_Vertex_Ptr;
    std::vector<std::vector<int> > *vecInitHorizontalLoops_Edge_Ptr;
    //std::vector<std::vector<std::pair<int, int>>> *coeffMatrix;
    std::vector<std::list<std::pair<int, int> > > *coeffMatrix;
    //
};

#endif // _MY_FilesOutputForOptimalCycles_H_
