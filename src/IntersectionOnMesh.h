/*
(c) 2012 Fengtao Fan
*/
#ifndef _MY_INTERSECTION_ON_MESH_H_
#define _MY_INTERSECTION_ON_MESH_H_

#include <vector>
#include "SimpleMesh.h"
#include <RenderVector3.h>

class IntersectionOnMeshComputation {
public:
    // constructor and assign operator
    IntersectionOnMeshComputation() {
        inMeshPtr = NULL;
        //
        A_ptr = NULL;
        B_ptr = NULL;
        C_ptr = NULL;
        D_ptr = NULL;
        //
        bIntersection = false;
        //
        AB_is_intersection = 0;
        CD_is_intersection = 0;
    }

    IntersectionOnMeshComputation(const IntersectionOnMeshComputation &rhs) {
        pt_A_type = rhs.pt_A_type;
        pt_B_type = rhs.pt_B_type;
        pt_C_type = rhs.pt_C_type;
        pt_D_type = rhs.pt_D_type;
        //
        inMeshPtr = rhs.inMeshPtr;
        //
        A_ptr = rhs.A_ptr;
        B_ptr = rhs.B_ptr;
        C_ptr = rhs.C_ptr;
        D_ptr = rhs.D_ptr;
        //
        bIntersection = rhs.bIntersection;
        //
        AB_is_intersection = rhs.AB_is_intersection;
        CD_is_intersection = rhs.CD_is_intersection;
        //
        if (bIntersection)
            resIntersection = rhs.resIntersection;

    }

    IntersectionOnMeshComputation &operator=(const IntersectionOnMeshComputation &rhs) {
        if (this == &rhs)
            return *this;
        pt_A_type = rhs.pt_A_type;
        pt_B_type = rhs.pt_B_type;
        pt_C_type = rhs.pt_C_type;
        pt_D_type = rhs.pt_D_type;
        //
        inMeshPtr = rhs.inMeshPtr;
        //
        A_ptr = rhs.A_ptr;
        B_ptr = rhs.B_ptr;
        C_ptr = rhs.C_ptr;
        D_ptr = rhs.D_ptr;
        //
        bIntersection = rhs.bIntersection;
        //
        AB_is_intersection = rhs.AB_is_intersection;
        CD_is_intersection = rhs.CD_is_intersection;
        //
        if (bIntersection)
            resIntersection = rhs.resIntersection;
        return *this;
    }

    // deconstructor and data clean functioin
    void Null_all() {
        inMeshPtr = NULL;
        //
        A_ptr = NULL;
        B_ptr = NULL;
        C_ptr = NULL;
        D_ptr = NULL;
        //
        bIntersection = false;
        //
        AB_is_intersection = 0;
        CD_is_intersection = 0;
    }

    void Clear() {
        //
        A_ptr = NULL;
        B_ptr = NULL;
        C_ptr = NULL;
        D_ptr = NULL;
        //
        bIntersection = false;
        //
        AB_is_intersection = 0;
        CD_is_intersection = 0;
    }

    ~IntersectionOnMeshComputation() {
        Null_all();
    }

    // real operations functions
    // Initial data
    void InitMeshPtr(_SimpleMesh *inPtr) {
        inMeshPtr = inPtr;
        return;
    }

    void InitFourPoints(std::pair<int, int> &in_A_type,
                        std::pair<int, int> &in_B_type,
                        std::pair<int, int> &in_C_type,
                        std::pair<int, int> &in_D_type,
                        Vector3 *in_A_ptr,
                        Vector3 *in_B_ptr,
                        Vector3 *in_C_ptr,
                        Vector3 *in_D_ptr) {
        pt_A_type = in_A_type;
        pt_B_type = in_B_type;
        pt_C_type = in_C_type;
        pt_D_type = in_D_type;
        //
        A_ptr = in_A_ptr;
        B_ptr = in_B_ptr;
        C_ptr = in_C_ptr;
        D_ptr = in_D_ptr;
        //
        return;
    }

    // report results
    void ReportIntersection() {
        if (bIntersection) {
            switch (intersect_pt_type.second) {
                case 0:
                    std::cout << "VERTEX : ";
                    break;
                case 1:
                    std::cout << "EDGE : ";
                    break;
                case 2:
                    std::cout << "FACE : ";
                    break;
            }
            std::cout << intersect_pt_type.first << std::endl;
            //
            std::cout << "POINT COORD : " << resIntersection[0] << " "
                      << resIntersection[1] << " "
                      << resIntersection[2] << std::endl;

        } else {
            std::cout << "NO INTERSECTION" << std::endl;
        }
    }

    //
    void ComputeIntersection();

    //
    void A_Pair_3D_Line_Intersection(Vector3 *in_A_ptr,
                                     Vector3 *in_B_ptr,
                                     Vector3 *in_C_ptr,
                                     Vector3 *in_D_ptr);

    // 6 cases
    void A_vertex_B_vertex_C_vertex_D_vertex(std::pair<int, int> &in_A_type,
                                             std::pair<int, int> &in_B_type,
                                             std::pair<int, int> &in_C_type,
                                             std::pair<int, int> &in_D_type,
                                             Vector3 *in_A_ptr,
                                             Vector3 *in_B_ptr,
                                             Vector3 *in_C_ptr,
                                             Vector3 *in_D_ptr);

    void A_vertex_B_vertex_C_vertex_D_edge(std::pair<int, int> &in_A_type,
                                           std::pair<int, int> &in_B_type,
                                           std::pair<int, int> &in_C_type,
                                           std::pair<int, int> &in_D_type,
                                           Vector3 *in_A_ptr,
                                           Vector3 *in_B_ptr,
                                           Vector3 *in_C_ptr,
                                           Vector3 *in_D_ptr);

    void A_vertex_B_vertex_C_edge_D_edge(std::pair<int, int> &in_A_type,
                                         std::pair<int, int> &in_B_type,
                                         std::pair<int, int> &in_C_type,
                                         std::pair<int, int> &in_D_type,
                                         Vector3 *in_A_ptr,
                                         Vector3 *in_B_ptr,
                                         Vector3 *in_C_ptr,
                                         Vector3 *in_D_ptr);

    void A_vertex_B_edge_C_vertex_D_edge(std::pair<int, int> &in_A_type,
                                         std::pair<int, int> &in_B_type,
                                         std::pair<int, int> &in_C_type,
                                         std::pair<int, int> &in_D_type,
                                         Vector3 *in_A_ptr,
                                         Vector3 *in_B_ptr,
                                         Vector3 *in_C_ptr,
                                         Vector3 *in_D_ptr);

    void A_vertex_B_edge_C_edge_D_edge(std::pair<int, int> &in_A_type,
                                       std::pair<int, int> &in_B_type,
                                       std::pair<int, int> &in_C_type,
                                       std::pair<int, int> &in_D_type,
                                       Vector3 *in_A_ptr,
                                       Vector3 *in_B_ptr,
                                       Vector3 *in_C_ptr,
                                       Vector3 *in_D_ptr);

    void A_edge_B_edge_C_edge_D_edge(std::pair<int, int> &in_A_type,
                                     std::pair<int, int> &in_B_type,
                                     std::pair<int, int> &in_C_type,
                                     std::pair<int, int> &in_D_type,
                                     Vector3 *in_A_ptr,
                                     Vector3 *in_B_ptr,
                                     Vector3 *in_C_ptr,
                                     Vector3 *in_D_ptr);

public:
    std::pair<int, int> pt_A_type;
    std::pair<int, int> pt_B_type;
    std::pair<int, int> pt_C_type;
    std::pair<int, int> pt_D_type;
    //
    _SimpleMesh *inMeshPtr;
    // coordinates for points
    Vector3 *A_ptr;
    Vector3 *B_ptr;
    Vector3 *C_ptr;
    Vector3 *D_ptr;
    // the actual intersection point
    Vector3 resIntersection;
    // flag to show there is an intersection between AB and CD or not
    // true --- there is an intersection
    // false -- no intersection
    bool bIntersection;
    //
    std::pair<int, int> intersect_pt_type;
    //flag for second bit:	0 -- vertex
    //						1 -- edge
    //						2 -- face
    int AB_is_intersection; // 0 -- no
    // 1 -- A
    // 2 -- B
    int CD_is_intersection; // 0 -- NO
    // 1 -- C
    // 2 -- D
};

#endif //_MY_INTERSECTION_ON_MESH_H_
