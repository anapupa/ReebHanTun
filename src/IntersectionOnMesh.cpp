/*
(c) 2012 Fengtao Fan
*/
#include "IntersectionOnMesh.h"
#include <cstdlib>

void IntersectionOnMeshComputation::A_Pair_3D_Line_Intersection(Vector3 *in_A_ptr,
                                                                Vector3 *in_B_ptr,
                                                                Vector3 *in_C_ptr,
                                                                Vector3 *in_D_ptr) {// preq: there exists an intersection between these two 3d lines
    /*
        L1 = P1 + a V1
        L2 = P2 + b V2
        //
        a V1 = (P2 - P1) + b V2
        //
        a (V1 X V2) = (P2 - P1) X V2
    */
    /*
    t(x_b - x_a) + x_a = s(x_d - x_c) + x_c
    t(y_b - y_a) + y_a = s(y_d - y_c) + y_c;
    common = (x_b - x_a) * (y_d - y_c) - (x_d - x_c) * (y_b - y_a)
    dividend_s =
    */
    Vector3 ab_dir = *in_B_ptr - *in_A_ptr;
    Vector3 cd_dir = *in_D_ptr - *in_C_ptr;
    //
    Vector3 CA = *in_C_ptr - *in_A_ptr;
    Vector3 ab_dir_cross_cd_dir = ab_dir ^ cd_dir;
    Vector3 pt_cross_cd_dir = CA ^ cd_dir;
    if (norm2(ab_dir_cross_cd_dir) == 0.0) {
        std::cout << "ERROR IN 3D LINES INTERSECTION COMPUTATION" << std::endl;
        exit(9);
    } else {
        double ab_ratio = -1.0;
        double max_up = pt_cross_cd_dir[0];
        double max_denom = ab_dir_cross_cd_dir[0];
        for (int i = 1; i < 3; i++) {
            if (fabs((double) (ab_dir_cross_cd_dir[i])) > fabs(max_denom)) {
                max_denom = ab_dir_cross_cd_dir[i];
                max_up = pt_cross_cd_dir[i];
            }
        }
        ab_ratio = max_up / max_denom;
        if (ab_ratio >= 0.0 && ab_ratio <= 1.0) {
            bIntersection = true;
            resIntersection = (*in_A_ptr) + ab_ratio * ab_dir;
            //resIntersection = (*in_C_ptr) + s_cd * ((*in_D_ptr) - (*in_C_ptr));
        } else {
            std::cout << "ERROR IN 3D LINES INTERSECTION COMPUTATION" << std::endl;
            exit(9);
        }
    }
    return;
    //
    //double t_ab = 0.0;
    //double s_cd = 0.0;

    //int Y = 1;
    //int X = 0;
    //if ((*in_A_ptr)[1] == (*in_B_ptr)[1] && (*in_D_ptr)[1] == (*in_C_ptr)[1])
    //{
    //	Y = 2;
    //}
    //else
    //{
    //	if ((*in_A_ptr)[0] == (*in_B_ptr)[0] && (*in_D_ptr)[0] == (*in_C_ptr)[0])
    //	{
    //		X = 2;
    //	}
    //}
    //double dividend_s = ((*in_A_ptr)[Y] - (*in_C_ptr)[Y]) * ((*in_B_ptr)[X] - (*in_A_ptr)[X]) -
    //				((*in_B_ptr)[Y] - (*in_A_ptr)[Y]) * ((*in_A_ptr)[X] - (*in_C_ptr)[X]) ;

    //double divisor = ((*in_D_ptr)[Y] - (*in_C_ptr)[Y]) * ((*in_B_ptr)[X] - (*in_A_ptr)[X]) -
    //				((*in_B_ptr)[Y] - (*in_A_ptr)[Y]) * ((*in_D_ptr)[X] - (*in_C_ptr)[X]) ;

    //double dividend_t = ((*in_C_ptr)[X] - (*in_A_ptr)[X]) * ((*in_D_ptr)[Y] - (*in_C_ptr)[Y]) -
    //				((*in_C_ptr)[Y] - (*in_A_ptr)[Y]) * ((*in_D_ptr)[X] - (*in_C_ptr)[X]) ;

    //s_cd = dividend_s / divisor;
    //t_ab = dividend_t / divisor;

    //if (s_cd >= 0.0 && s_cd <= 1.0)
    //{
    //	bIntersection =  true;
    //	//resIntersection = (*in_A_ptr) + t_ab * ((*in_B_ptr) - (*in_A_ptr));
    //	resIntersection = (*in_C_ptr) + s_cd * ((*in_D_ptr) - (*in_C_ptr));
    //}
    //else
    //{
    //	if (t_ab >= 0.0 && t_ab <= 1.0)
    //	{
    //		bIntersection =  true;
    //		resIntersection = (*in_A_ptr) + t_ab * ((*in_B_ptr) - (*in_A_ptr));
    //
    //	}
    //	else
    //	{
    //		std::cout << "ERROR IN 3D LINES INTERSECTION COMPUTATION" << std::endl;
    //		exit(0);
    //	}
    //}
    return;
}

void IntersectionOnMeshComputation::ComputeIntersection() {
    bIntersection = false;
    //
    AB_is_intersection = 0;
    CD_is_intersection = 0;
    //
    switch (pt_A_type.second + pt_B_type.second + pt_C_type.second + pt_D_type.second) {
        case 0:
            // all are vertices
            A_vertex_B_vertex_C_vertex_D_vertex(pt_A_type, pt_B_type,
                                                pt_C_type, pt_D_type,
                                                A_ptr, B_ptr,
                                                C_ptr, D_ptr);
            break;
        case 1:
            // one is edge, the rest are vertices
            if (pt_A_type.second || pt_B_type.second) {
                if (pt_A_type.second) {
                    A_vertex_B_vertex_C_vertex_D_edge(pt_C_type, pt_D_type,
                                                      pt_B_type, pt_A_type,
                                                      C_ptr, D_ptr,
                                                      B_ptr, A_ptr
                    );
                    int tmp = AB_is_intersection;
                    AB_is_intersection = CD_is_intersection;
                    CD_is_intersection = tmp;
                    if (AB_is_intersection)
                        AB_is_intersection = 3 - AB_is_intersection; // 1-->2, 2-->1
                } else {
                    A_vertex_B_vertex_C_vertex_D_edge(pt_C_type, pt_D_type,
                                                      pt_A_type, pt_B_type,
                                                      C_ptr, D_ptr,
                                                      A_ptr, B_ptr
                    );
                    int tmp = AB_is_intersection;
                    AB_is_intersection = CD_is_intersection;
                    CD_is_intersection = tmp;
                }
            } else {
                if (pt_C_type.second) {
                    A_vertex_B_vertex_C_vertex_D_edge(pt_A_type, pt_B_type,
                                                      pt_D_type, pt_C_type,
                                                      A_ptr, B_ptr,
                                                      D_ptr, C_ptr
                    );
                    if (CD_is_intersection)
                        CD_is_intersection = 3 - CD_is_intersection; // 1-->2, 2-->1
                } else
                    A_vertex_B_vertex_C_vertex_D_edge(pt_A_type, pt_B_type,
                                                      pt_C_type, pt_D_type,
                                                      A_ptr, B_ptr,
                                                      C_ptr, D_ptr
                    );
            }
            break;
        case 2:
            // two vertices and two edges
            if (pt_A_type.second && pt_B_type.second ||
                pt_C_type.second && pt_D_type.second) {
                if (pt_A_type.second && pt_B_type.second) {
                    A_vertex_B_vertex_C_edge_D_edge(pt_C_type, pt_D_type,
                                                    pt_A_type, pt_B_type,
                                                    C_ptr, D_ptr,
                                                    A_ptr, B_ptr
                    );
                    int tmp = AB_is_intersection;
                    AB_is_intersection = CD_is_intersection;
                    CD_is_intersection = tmp;
                } else {
                    A_vertex_B_vertex_C_edge_D_edge(pt_A_type, pt_B_type,
                                                    pt_C_type, pt_D_type,
                                                    A_ptr, B_ptr,
                                                    C_ptr, D_ptr
                    );
                }
            } else {
                if (pt_A_type.second) {
                    if (pt_C_type.second) {
                        A_vertex_B_edge_C_vertex_D_edge(pt_B_type, pt_A_type,
                                                        pt_D_type, pt_C_type,
                                                        B_ptr, A_ptr,
                                                        D_ptr, C_ptr);
                        if (AB_is_intersection)
                            AB_is_intersection = 3 - AB_is_intersection; // 1-->2, 2-->1
                        if (CD_is_intersection)
                            CD_is_intersection = 3 - CD_is_intersection; // 1-->2, 2-->1
                    } else {
                        A_vertex_B_edge_C_vertex_D_edge(pt_B_type, pt_A_type,
                                                        pt_C_type, pt_D_type,
                                                        B_ptr, A_ptr,
                                                        C_ptr, D_ptr);

                        if (AB_is_intersection)
                            AB_is_intersection = 3 - AB_is_intersection; // 1-->2, 2-->1
                    }
                } else {
                    if (pt_C_type.second) {
                        A_vertex_B_edge_C_vertex_D_edge(pt_A_type, pt_B_type,
                                                        pt_D_type, pt_C_type,
                                                        A_ptr, B_ptr,
                                                        D_ptr, C_ptr);
                        if (CD_is_intersection)
                            CD_is_intersection = 3 - CD_is_intersection; // 1-->2, 2-->1
                    } else {
                        A_vertex_B_edge_C_vertex_D_edge(pt_A_type, pt_B_type,
                                                        pt_C_type, pt_D_type,
                                                        A_ptr, B_ptr,
                                                        C_ptr, D_ptr);
                    }
                }

            }
            break;
        case 3:
            // one is vertex, the rest are edges
            if (!pt_A_type.second || !pt_B_type.second) {
                if (!pt_A_type.second) {
                    A_vertex_B_edge_C_edge_D_edge(pt_A_type, pt_B_type,
                                                  pt_C_type, pt_D_type,
                                                  A_ptr, B_ptr,
                                                  C_ptr, D_ptr
                    );
                } else {
                    A_vertex_B_edge_C_edge_D_edge(pt_B_type, pt_A_type,
                                                  pt_C_type, pt_D_type,
                                                  B_ptr, A_ptr,
                                                  C_ptr, D_ptr
                    );
                    if (AB_is_intersection)
                        AB_is_intersection = 3 - AB_is_intersection; // 1-->2, 2-->1
                }
            } else {
                if (!pt_C_type.second) {
                    A_vertex_B_edge_C_edge_D_edge(pt_C_type, pt_D_type,
                                                  pt_A_type, pt_B_type,
                                                  C_ptr, D_ptr,
                                                  A_ptr, B_ptr
                    );
                    int tmp = AB_is_intersection;
                    AB_is_intersection = CD_is_intersection;
                    CD_is_intersection = tmp;
                } else {
                    A_vertex_B_edge_C_edge_D_edge(pt_D_type, pt_C_type,
                                                  pt_A_type, pt_B_type,
                                                  D_ptr, C_ptr,
                                                  A_ptr, B_ptr
                    );
                    int tmp = AB_is_intersection;
                    AB_is_intersection = CD_is_intersection;
                    CD_is_intersection = tmp;
                    if (CD_is_intersection)
                        CD_is_intersection = 3 - CD_is_intersection; // 1-->2, 2-->1
                }
            }
            break;
        case 4:
            // all are edges
            A_edge_B_edge_C_edge_D_edge(pt_A_type, pt_B_type,
                                        pt_C_type, pt_D_type,
                                        A_ptr, B_ptr,
                                        C_ptr, D_ptr);
            break;
    }
    return;
}

void IntersectionOnMeshComputation::A_vertex_B_vertex_C_vertex_D_vertex(std::pair<int, int> &in_A_type,
                                                                        std::pair<int, int> &in_B_type,
                                                                        std::pair<int, int> &in_C_type,
                                                                        std::pair<int, int> &in_D_type,
                                                                        Vector3 *in_A_ptr,
                                                                        Vector3 *in_B_ptr,
                                                                        Vector3 *in_C_ptr,
                                                                        Vector3 *in_D_ptr) {
    if (in_A_type.first == in_C_type.first ||
        in_A_type.first == in_D_type.first ||
        in_B_type.first == in_C_type.first ||
        in_B_type.first == in_D_type.first) {//
        bIntersection = true;
        //
        if (in_A_type.first == in_C_type.first ||
            in_A_type.first == in_D_type.first) {
            resIntersection = *in_A_ptr;
            //
            intersect_pt_type.first = in_A_type.first;
            intersect_pt_type.second = 0; // vertex type
            //
            //
            AB_is_intersection = 1; // for A is the intersection
            if (in_A_type.first == in_C_type.first)
                CD_is_intersection = 1;
            else
                CD_is_intersection = 2;
        } else {
            resIntersection = *in_B_ptr;
            //
            intersect_pt_type.first = in_B_type.first;
            intersect_pt_type.second = 0; // vertex type
            //
            AB_is_intersection = 2; // for A is the intersection
            if (in_B_type.first == in_C_type.first)
                CD_is_intersection = 1;
            else
                CD_is_intersection = 2;
        }
    } else {
        //
        bIntersection = false;
    }

    return;
}

void IntersectionOnMeshComputation::A_vertex_B_vertex_C_vertex_D_edge(std::pair<int, int> &in_A_type,
                                                                      std::pair<int, int> &in_B_type,
                                                                      std::pair<int, int> &in_C_type,
                                                                      std::pair<int, int> &in_D_type,
                                                                      Vector3 *in_A_ptr,
                                                                      Vector3 *in_B_ptr,
                                                                      Vector3 *in_C_ptr,
                                                                      Vector3 *in_D_ptr) {
    if (in_C_type.first == in_A_type.first ||
        in_C_type.first == in_B_type.first) {// C==A || C==B
        bIntersection = true;
        //
        resIntersection = *in_C_ptr;
        //
        intersect_pt_type.first = in_C_type.first;
        intersect_pt_type.second = 0; // vertex type
        //
        CD_is_intersection = 1; // for A is the intersection
        if (in_A_type.first == in_C_type.first)
            AB_is_intersection = 1;
        else
            AB_is_intersection = 2;
    } else {// D is on AB or not
        // get the edge id connecting A and B
        int edge_ab = -1;
        for (unsigned int i = 0; i < inMeshPtr->vecVertex[in_A_type.first].adjEdges.size(); i++) {
            if (inMeshPtr->vecEdge[inMeshPtr->vecVertex[in_A_type.first].adjEdges[i]].v0 == in_B_type.first ||
                inMeshPtr->vecEdge[inMeshPtr->vecVertex[in_A_type.first].adjEdges[i]].v1 == in_B_type.first) {
                edge_ab = inMeshPtr->vecVertex[in_A_type.first].adjEdges[i];
                break;
            }
        }
        if (edge_ab == -1) {
            std::cout << "AB is not an edge " << std::endl;
            exit(0);
        }
        //
        if (edge_ab == in_D_type.first) {// D is on AB
            bIntersection = true;
            resIntersection = *in_D_ptr;
            //
            intersect_pt_type.first = in_D_type.first;
            intersect_pt_type.second = 1; // edge type
            //
            CD_is_intersection = 2; // for D is the intersection
            AB_is_intersection = 0;
        } else
            bIntersection = false;
    }
    return;
}

void IntersectionOnMeshComputation::A_vertex_B_vertex_C_edge_D_edge(std::pair<int, int> &in_A_type,
                                                                    std::pair<int, int> &in_B_type,
                                                                    std::pair<int, int> &in_C_type,
                                                                    std::pair<int, int> &in_D_type,
                                                                    Vector3 *in_A_ptr,
                                                                    Vector3 *in_B_ptr,
                                                                    Vector3 *in_C_ptr,
                                                                    Vector3 *in_D_ptr) {
    int edge_ab = -1;
    for (unsigned int i = 0; i < inMeshPtr->vecVertex[in_A_type.first].adjEdges.size(); i++) {
        if (inMeshPtr->vecEdge[inMeshPtr->vecVertex[in_A_type.first].adjEdges[i]].v0 == in_B_type.first ||
            inMeshPtr->vecEdge[inMeshPtr->vecVertex[in_A_type.first].adjEdges[i]].v1 == in_B_type.first) {
            edge_ab = inMeshPtr->vecVertex[in_A_type.first].adjEdges[i];
            break;
        }
    }
    if (edge_ab == -1) {
        std::cout << "AB is not an edge " << std::endl;
        exit(0);
    }
    if (edge_ab == in_C_type.first ||
        edge_ab == in_D_type.first) {
        bIntersection = true;
        //
        if (edge_ab == in_C_type.first) {
            resIntersection = *in_C_ptr;
            //
            intersect_pt_type.first = in_C_type.first;
            intersect_pt_type.second = 1; // edge type
            //
            AB_is_intersection = 0;
            CD_is_intersection = 1; // for C
        } else {
            resIntersection = *in_D_ptr;
            //
            intersect_pt_type.first = in_D_type.first;
            intersect_pt_type.second = 1; // edge type
            //
            AB_is_intersection = 0;
            CD_is_intersection = 2; // for D
        }
    } else
        bIntersection = false;
    return;
}

void IntersectionOnMeshComputation::A_vertex_B_edge_C_vertex_D_edge(std::pair<int, int> &in_A_type,
                                                                    std::pair<int, int> &in_B_type,
                                                                    std::pair<int, int> &in_C_type,
                                                                    std::pair<int, int> &in_D_type,
                                                                    Vector3 *in_A_ptr,
                                                                    Vector3 *in_B_ptr,
                                                                    Vector3 *in_C_ptr,
                                                                    Vector3 *in_D_ptr) {
    // if there are not on the same triangle, only trivial intersection happens
    int tri_ab = -1;
    int tri_cd = -1;
    if (inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_B_type.first].AdjTri[0]].v0 == in_A_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_B_type.first].AdjTri[0]].v1 == in_A_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_B_type.first].AdjTri[0]].v2 == in_A_type.first) {
        tri_ab = inMeshPtr->vecEdge[in_B_type.first].AdjTri[0];
    } else
        tri_ab = inMeshPtr->vecEdge[in_B_type.first].AdjTri[1];
    //
    if (inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_D_type.first].AdjTri[0]].v0 == in_C_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_D_type.first].AdjTri[0]].v1 == in_C_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_D_type.first].AdjTri[0]].v2 == in_C_type.first) {
        tri_cd = inMeshPtr->vecEdge[in_D_type.first].AdjTri[0];
    } else {
        tri_cd = inMeshPtr->vecEdge[in_D_type.first].AdjTri[1];
    }
    //
    // check trivial intersection first
    if (in_A_type.first == in_C_type.first ||
        in_B_type.first == in_D_type.first) {// this checking is true under this special case
        //always check the vertex first
        if (in_A_type.first == in_C_type.first) {
            bIntersection = true;
            resIntersection = *in_C_ptr;
            //
            intersect_pt_type = in_A_type; // vertex type
            //
            AB_is_intersection = 1;
            CD_is_intersection = 1;
        } else {
            if (*in_B_ptr == *in_D_ptr) {
                bIntersection = true;
                resIntersection = *in_B_ptr;
                //
                intersect_pt_type = in_B_type; // edge type
                //
                AB_is_intersection = 2;
                CD_is_intersection = 2;
            } else {
                bIntersection = false;
            }
        }
    } else {
        if (tri_ab == tri_cd) {// compute non-trivial intersections
            A_Pair_3D_Line_Intersection(in_A_ptr, in_B_ptr, in_C_ptr, in_D_ptr);
            //
            intersect_pt_type.first = tri_ab;
            intersect_pt_type.second = 2;// face type
            bIntersection = true;
            //
            AB_is_intersection = 0;
            CD_is_intersection = 0;
        } else {// no intersection
            bIntersection = false;
        }
    }
}

void IntersectionOnMeshComputation::A_vertex_B_edge_C_edge_D_edge(std::pair<int, int> &in_A_type,
                                                                  std::pair<int, int> &in_B_type,
                                                                  std::pair<int, int> &in_C_type,
                                                                  std::pair<int, int> &in_D_type,
                                                                  Vector3 *in_A_ptr,
                                                                  Vector3 *in_B_ptr,
                                                                  Vector3 *in_C_ptr,
                                                                  Vector3 *in_D_ptr) {
    // if there are not on the same triangle, only trivial intersection happens
    int tri_ab = -1;
    int tri_cd = -1;
    if (inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_B_type.first].AdjTri[0]].v0 == in_A_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_B_type.first].AdjTri[0]].v1 == in_A_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_B_type.first].AdjTri[0]].v2 == in_A_type.first) {
        tri_ab = inMeshPtr->vecEdge[in_B_type.first].AdjTri[0];
    } else
        tri_ab = inMeshPtr->vecEdge[in_B_type.first].AdjTri[1];
    //
    if (inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_D_type.first].AdjTri[0]].e01 == in_C_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_D_type.first].AdjTri[0]].e02 == in_C_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_D_type.first].AdjTri[0]].e12 == in_C_type.first) {
        tri_cd = inMeshPtr->vecEdge[in_D_type.first].AdjTri[0];
    } else {
        tri_cd = inMeshPtr->vecEdge[in_D_type.first].AdjTri[1];
    }
    //
    // check trivial intersection first
    if (*in_B_ptr == *in_C_ptr ||
        *in_B_ptr == *in_D_ptr) {
        bIntersection = true;
        //
        resIntersection = *in_B_ptr;
        //
        intersect_pt_type = in_B_type;
        //
        AB_is_intersection = 2;
        if (*in_B_ptr == *in_C_ptr)
            CD_is_intersection = 1;
        else
            CD_is_intersection = 2;
    } else {
        if (tri_ab == tri_cd) {// if there is a intersection, compute  non-trivial intersections
            // but there are also non-intersection cases
            // find the common intersection point for the edge_edge segment
            int common_pt = -1;
            if (inMeshPtr->vecEdge[in_C_type.first].v0 == inMeshPtr->vecEdge[in_D_type.first].v0 ||
                inMeshPtr->vecEdge[in_C_type.first].v0 == inMeshPtr->vecEdge[in_D_type.first].v1)
                common_pt = inMeshPtr->vecEdge[in_C_type.first].v0;
            else
                common_pt = inMeshPtr->vecEdge[in_C_type.first].v1;
            //
            if (common_pt == in_A_type.first) { // there must exist intersection
                A_Pair_3D_Line_Intersection(in_A_ptr, in_B_ptr, in_C_ptr, in_D_ptr);
                intersect_pt_type.first = tri_ab;
                intersect_pt_type.second = 2;// face type
                bIntersection = true;
                // by default
                // CD_is_intersection = 0;
                // AB_is_intersection = 0;
            } else {// check the distance to the common_pt
                Vector3 common_pt_coord(inMeshPtr->vecVertex[common_pt].x,
                                        inMeshPtr->vecVertex[common_pt].y,
                                        inMeshPtr->vecVertex[common_pt].z);
                if (in_C_type.first == in_B_type.first) {
                    if (norm2((*in_B_ptr) - common_pt_coord) <
                        norm2((*in_C_ptr) - common_pt_coord)) {// there must exist intersection
                        A_Pair_3D_Line_Intersection(in_A_ptr, in_B_ptr, in_C_ptr, in_D_ptr);
                        //
                        intersect_pt_type.first = tri_ab;
                        intersect_pt_type.second = 2;// face type
                        bIntersection = true;
                    } else
                        bIntersection = false;
                } else {
                    if (norm2((*in_B_ptr) - common_pt_coord) <
                        norm2((*in_D_ptr) - common_pt_coord)) {// there must exist intersection
                        A_Pair_3D_Line_Intersection(in_A_ptr, in_B_ptr, in_C_ptr, in_D_ptr);
                        //
                        intersect_pt_type.first = tri_ab;
                        intersect_pt_type.second = 2;// face type
                        bIntersection = true;
                    } else
                        bIntersection = false;
                }
            }
        } else {// no intersection
            bIntersection = false;
        }
    }
}

void IntersectionOnMeshComputation::A_edge_B_edge_C_edge_D_edge(std::pair<int, int> &in_A_type,
                                                                std::pair<int, int> &in_B_type,
                                                                std::pair<int, int> &in_C_type,
                                                                std::pair<int, int> &in_D_type,
                                                                Vector3 *in_A_ptr,
                                                                Vector3 *in_B_ptr,
                                                                Vector3 *in_C_ptr,
                                                                Vector3 *in_D_ptr) {
    // if there are not on the same triangle, only trivial intersection happens
    int tri_ab = -1;
    int tri_cd = -1;
    if (inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_B_type.first].AdjTri[0]].e01 == in_A_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_B_type.first].AdjTri[0]].e02 == in_A_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_B_type.first].AdjTri[0]].e12 == in_A_type.first) {
        tri_ab = inMeshPtr->vecEdge[in_B_type.first].AdjTri[0];
    } else
        tri_ab = inMeshPtr->vecEdge[in_B_type.first].AdjTri[1];
    //
    if (inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_D_type.first].AdjTri[0]].e01 == in_C_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_D_type.first].AdjTri[0]].e02 == in_C_type.first ||
        inMeshPtr->vecTriangle[inMeshPtr->vecEdge[in_D_type.first].AdjTri[0]].e12 == in_C_type.first) {
        tri_cd = inMeshPtr->vecEdge[in_D_type.first].AdjTri[0];
    } else {
        tri_cd = inMeshPtr->vecEdge[in_D_type.first].AdjTri[1];
    }
    //
    // check trivial intersection first
    if (*in_A_ptr == *in_C_ptr ||
        *in_A_ptr == *in_D_ptr ||
        *in_B_ptr == *in_C_ptr ||
        *in_B_ptr == *in_D_ptr
            ) {
        if (*in_A_ptr == *in_C_ptr ||
            *in_A_ptr == *in_D_ptr) {
            bIntersection = true;
            //
            resIntersection = *in_A_ptr;
            //
            intersect_pt_type = in_A_type;// edge type
            //
            AB_is_intersection = 1;
            if (*in_A_ptr == *in_C_ptr)
                CD_is_intersection = 1;
            else
                CD_is_intersection = 2;
        } else {
            bIntersection = true;
            //
            resIntersection = *in_B_ptr;
            //
            intersect_pt_type = in_B_type; //edge type
            //
            AB_is_intersection = 2;
            if (*in_B_ptr == *in_C_ptr)
                CD_is_intersection = 1;
            else
                CD_is_intersection = 2;
        }
    } else {//
        if (tri_ab == tri_cd) {// if there is a intersection, compute  non-trivial intersections
            // but there are also non-intersection cases
            // find the common intersection point for the edge_edge segment
            //std::cout << "IN the same triangle case " << std::endl;
            int ab_common_pt = -1;
            int cd_common_pt = -1;
            //
            if (inMeshPtr->vecEdge[in_A_type.first].v0 == inMeshPtr->vecEdge[in_B_type.first].v0 ||
                inMeshPtr->vecEdge[in_A_type.first].v0 == inMeshPtr->vecEdge[in_B_type.first].v1)
                ab_common_pt = inMeshPtr->vecEdge[in_A_type.first].v0;
            else
                ab_common_pt = inMeshPtr->vecEdge[in_A_type.first].v1;
            //
            if (inMeshPtr->vecEdge[in_C_type.first].v0 == inMeshPtr->vecEdge[in_D_type.first].v0 ||
                inMeshPtr->vecEdge[in_C_type.first].v0 == inMeshPtr->vecEdge[in_D_type.first].v1)
                cd_common_pt = inMeshPtr->vecEdge[in_C_type.first].v0;
            else
                cd_common_pt = inMeshPtr->vecEdge[in_C_type.first].v1;
            //
            if (ab_common_pt == cd_common_pt) {// 4 cases here
                Vector3 common_pt_coord(inMeshPtr->vecVertex[ab_common_pt].x,
                                        inMeshPtr->vecVertex[ab_common_pt].y,
                                        inMeshPtr->vecVertex[ab_common_pt].z);
                double dist_a = norm2((*in_A_ptr) - common_pt_coord);
                double dist_b = norm2((*in_B_ptr) - common_pt_coord);
                double dist_c = norm2((*in_C_ptr) - common_pt_coord);
                double dist_d = norm2((*in_D_ptr) - common_pt_coord);

                if (in_A_type.first == in_C_type.first) {// A and C are on the same line
                    if (dist_a < dist_c && dist_b > dist_d ||
                        dist_a > dist_c && dist_b < dist_d) {// there must exist intersection
                        A_Pair_3D_Line_Intersection(in_A_ptr, in_B_ptr, in_C_ptr, in_D_ptr);
                        //
                        intersect_pt_type.first = tri_ab;
                        intersect_pt_type.second = 2;// face type
                        bIntersection = true;
                    } else
                        bIntersection = false;
                } else {// A and D are on the same line
                    if (dist_a < dist_d && dist_b > dist_c ||
                        dist_a > dist_d && dist_b < dist_c) {// there must exist intersection
                        A_Pair_3D_Line_Intersection(in_A_ptr, in_B_ptr, in_C_ptr, in_D_ptr);
                        //
                        intersect_pt_type.first = tri_ab;
                        intersect_pt_type.second = 2;// face type
                        bIntersection = true;
                    } else
                        bIntersection = false;
                }
            } else {// 4 cases here
                Vector3 common_pt_coord(inMeshPtr->vecVertex[ab_common_pt].x,
                                        inMeshPtr->vecVertex[ab_common_pt].y,
                                        inMeshPtr->vecVertex[ab_common_pt].z);
                // find the common edge
                if (in_A_type.first == in_C_type.first ||
                    in_A_type.first == in_D_type.first) {// common edge A and C | D
                    if (in_A_type.first == in_C_type.first) {// common edge A and C
                        if (norm2((*in_A_ptr) - common_pt_coord) > norm2((*in_C_ptr) - common_pt_coord)) {
                            A_Pair_3D_Line_Intersection(in_A_ptr, in_B_ptr, in_C_ptr, in_D_ptr);
                            //
                            intersect_pt_type.first = tri_ab;
                            intersect_pt_type.second = 2;// face type
                            bIntersection = true;
                        } else
                            bIntersection = false;
                    } else {// common edge A and D
                        if (norm2((*in_A_ptr) - common_pt_coord) > norm2((*in_D_ptr) - common_pt_coord)) {
                            A_Pair_3D_Line_Intersection(in_A_ptr, in_B_ptr, in_C_ptr, in_D_ptr);
                            //
                            intersect_pt_type.first = tri_ab;
                            intersect_pt_type.second = 2;// face type
                            bIntersection = true;
                        } else
                            bIntersection = false;
                    }
                } else {// common edge B and C | D
                    if (in_B_type.first == in_C_type.first) {// common edge B and C
                        if (norm2((*in_B_ptr) - common_pt_coord) > norm2((*in_C_ptr) - common_pt_coord)) {
                            A_Pair_3D_Line_Intersection(in_A_ptr, in_B_ptr, in_C_ptr, in_D_ptr);
                            //
                            intersect_pt_type.first = tri_ab;
                            intersect_pt_type.second = 2;// face type
                            bIntersection = true;
                        } else
                            bIntersection = false;
                    } else {// common edge B and D
                        if (norm2((*in_B_ptr) - common_pt_coord) > norm2((*in_D_ptr) - common_pt_coord)) {
                            A_Pair_3D_Line_Intersection(in_A_ptr, in_B_ptr, in_C_ptr, in_D_ptr);
                            //
                            intersect_pt_type.first = tri_ab;
                            intersect_pt_type.second = 2;// face type
                            bIntersection = true;
                        } else
                            bIntersection = false;
                    }
                }
            }
        } else {// no intersection
            bIntersection = false;
        }
    }
    return;
}
