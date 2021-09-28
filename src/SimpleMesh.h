/*
(c) 2012 Fengtao Fan
*/
#pragma once

#include <vector>
#include <RenderVector3.h>

class _SimpleMeshVertex {
public:
    _SimpleMeshVertex() {
        x = y = z = 0.f;
    }

    _SimpleMeshVertex(const float x0, const float y0, const float z0) {
        x = x0;
        y = y0;
        z = z0;
    }

    _SimpleMeshVertex(const _SimpleMeshVertex &rhs) {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        //
        adjEdges.assign(rhs.adjEdges.begin(), rhs.adjEdges.end());
    }

    _SimpleMeshVertex &operator=(const _SimpleMeshVertex &rhs) {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        //
        adjEdges.assign(rhs.adjEdges.begin(), rhs.adjEdges.end());
        return *this;
    }

    _SimpleMeshVertex &operator-(const _SimpleMeshVertex &rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    float operator*(const _SimpleMeshVertex &rhs) {
        return (x * rhs.x + y * rhs.y + z * rhs.z);
    }

    void set(const float rx, const float ry, const float rz) {
        x = rx;
        y = ry;
        z = rz;
    }

public:
    float x;
    float y;
    float z;
    std::vector<int> adjEdges;
};

class _SimpleMeshEdge {// v0 < v1 , v0 \neq v1 for mesh is manifold
    // v0--e01--e1
public:
    _SimpleMeshEdge() {
        v0 = v1 = AdjTriNum = -1;
        AdjTri[0] = AdjTri[1] = -1;
    }

    _SimpleMeshEdge(const int lowVer, const int highVer, const int flag) {// lowVer < highVer
        v0 = lowVer;
        v1 = highVer;
        AdjTriNum = -1;
        AdjTri[0] = AdjTri[1] = -1;
    }

    _SimpleMeshEdge(const int lowVer, const int highVer) {
        v0 = lowVer < highVer ? lowVer : highVer;
        v1 = lowVer + highVer - v0;
        AdjTriNum = -1;
        AdjTri[0] = AdjTri[1] = -1;
    }

    _SimpleMeshEdge(const _SimpleMeshEdge &rhs) {
        v0 = rhs.v0;
        v1 = rhs.v1;
        AdjTriNum = rhs.AdjTriNum;
        AdjTri[0] = rhs.AdjTri[0];
        AdjTri[1] = rhs.AdjTri[1];
    }

    _SimpleMeshEdge &operator=(const _SimpleMeshEdge &rhs) {
        v0 = rhs.v0;
        v1 = rhs.v1;
        AdjTriNum = rhs.AdjTriNum;
        AdjTri[0] = rhs.AdjTri[0];
        AdjTri[1] = rhs.AdjTri[1];
        return *this;
    }

public:
    int v0;
    int v1;
    int AdjTri[2];
    int AdjTriNum;
};

class _SimpleMeshTriangle {// v0 < v1 < v2 for mesh is manifold
public:
    _SimpleMeshTriangle() {
        v0 = v1 = v2 = -1;
        e01 = e12 = e02 = -1;
    }

    _SimpleMeshTriangle(const int x0, const int y0, const int z0, const int xy, const int yz, const int xz) {
        int tmp = 0;
        v0 = x0;
        v1 = y0;
        v2 = z0;
        e01 = xy;
        e12 = yz;
        e02 = xz;
        if (v0 > v1) {
            tmp = v0;
            v0 = v1;
            v1 = tmp;
            //
            tmp = e12;
            e12 = e02;
            e02 = tmp;
        }
        if (v1 > v2) {
            tmp = v1;
            v1 = v2;
            v2 = tmp;
            //
            tmp = e02;
            e02 = e01;
            e01 = tmp;
        }
        if (v0 > v1) {
            tmp = v0;
            v0 = v1;
            v1 = tmp;
            //
            tmp = e12;
            e12 = e02;
            e02 = tmp;
        }
    }

    _SimpleMeshTriangle(const _SimpleMeshTriangle &rhs) {
        v0 = rhs.v0;
        v1 = rhs.v1;
        v2 = rhs.v2;
        //
        e01 = rhs.e01;
        e12 = rhs.e12;
        e02 = rhs.e02;
    }

    _SimpleMeshTriangle &operator=(const _SimpleMeshTriangle &rhs) {
        v0 = rhs.v0;
        v1 = rhs.v1;
        v2 = rhs.v2;
        //
        e01 = rhs.e01;
        e12 = rhs.e12;
        e02 = rhs.e02;
        return *this;
    }

    void sortVertices() {
        int tmp = 0;
        if (v0 > v1) {
            tmp = v0;
            v0 = v1;
            v1 = tmp;
        }
        if (v1 > v2) {
            tmp = v1;
            v1 = v2;
            v2 = tmp;
        }
        if (v0 > v1) {
            tmp = v0;
            v0 = v1;
            v1 = tmp;
        }
    }

public:
    int v0;
    int v1;
    int v2;
    int e01;
    int e12;
    int e02;
    /*
    v0---e02-----v2
     \          /
      \        /
       e01    e12
        \    /
         \  /
          v1
    */
};

class _SimpleMesh {
public:
    _SimpleMesh() {
    }

    _SimpleMesh(const _SimpleMesh &rhs) {
        vecVertex.assign(rhs.vecVertex.begin(), rhs.vecVertex.end());
        vecEdge.assign(rhs.vecEdge.begin(), rhs.vecEdge.end());
        vecTriangle.assign(rhs.vecTriangle.begin(), rhs.vecTriangle.end());
    }

    _SimpleMesh &operator=(const _SimpleMesh &rhs) {
        vecVertex.assign(rhs.vecVertex.begin(), rhs.vecVertex.end());
        vecEdge.assign(rhs.vecEdge.begin(), rhs.vecEdge.end());
        vecTriangle.assign(rhs.vecTriangle.begin(), rhs.vecTriangle.end());
        return *this;
    }

    ~_SimpleMesh() {
        vecVertex.clear();
        vecEdge.clear();
        vecTriangle.clear();
    }

    void SetMeshNormalPtr(std::vector<Vector3> *inPtr) {
        meshNormalPtr = inPtr;
    }

    inline int TriangleAdjacentToOneEdgesOneVertex(const int e0, const int opp_v) {
        int triangle_id = vecEdge[e0].AdjTri[0];
        if (vecTriangle[triangle_id].v0 != opp_v &&
            vecTriangle[triangle_id].v1 != opp_v &&
            vecTriangle[triangle_id].v2 != opp_v) {//
            triangle_id = vecEdge[e0].AdjTri[1];
        }
        return triangle_id;
    }

    inline bool EdgeOnTriangle(const int eid, const int fid) {
        if (vecTriangle[fid].e01 == eid ||
            vecTriangle[fid].e02 == eid ||
            vecTriangle[fid].e12 == eid) {
            return true;
        } else {
            return false;
        }
    }

    void EdgeConnectingTwoVertices(const int src_x, const int dist_y, const int eid_y, const int nTriangleIdx,
                                   int &outEdge) {
        int edge_u = 0;
        int edge_v = 0;
        //
        edge_u = vecTriangle[nTriangleIdx].e01;
        if (vecTriangle[nTriangleIdx].e01 == eid_y)
            edge_u = vecTriangle[nTriangleIdx].e02;
        //
        edge_v = vecTriangle[nTriangleIdx].e01 +
                 vecTriangle[nTriangleIdx].e02 +
                 vecTriangle[nTriangleIdx].e12 -
                 (eid_y + edge_u);
        //
        if (vecEdge[edge_u].v0 == dist_y ||
            vecEdge[edge_u].v1 == dist_y) {
            outEdge = edge_u;
        } else {
            outEdge = edge_v;
        }
        return;
    }

    std::pair<int, bool> EdgeIndex(const int src, const int dst) {
        std::pair<int, bool> ret = std::pair<int, bool>(0, false);
        for (std::vector<int>::iterator vIter = vecVertex[src].adjEdges.begin();
             vIter != vecVertex[src].adjEdges.end(); vIter++) {
            if (vecEdge[*vIter].v0 == dst ||
                vecEdge[*vIter].v1 == dst) {
                ret.first = *vIter;
                ret.second = true;
                break;
            }
        }
        return ret;
    }

    int EdgeConnectingTwoVertices(const int src_x, const int dist_y, const int eid_y, const int nTriangleIdx) {
        int outEdge = -1;
        int edge_u = 0;
        //
        edge_u = vecTriangle[nTriangleIdx].e01;
        if (vecTriangle[nTriangleIdx].e01 == eid_y)
            edge_u = vecTriangle[nTriangleIdx].e02;
        //
        //
        if (vecEdge[edge_u].v0 == dist_y ||
            vecEdge[edge_u].v1 == dist_y) {
            outEdge = edge_u;
        } else {
            outEdge = vecTriangle[nTriangleIdx].e01 +
                      vecTriangle[nTriangleIdx].e02 +
                      vecTriangle[nTriangleIdx].e12 -
                      (eid_y + edge_u);
        }
        return outEdge;
    }

    void LoadMeshInOFFformat(_SimpleMeshVertex &minBds,
                             _SimpleMeshVertex &maxBds,
                             std::vector<Vector3> &meshNormal,
                             const char *fileName, const float fEnlargeFactor);

    void LoadMeshInOFFformat(_SimpleMeshVertex &minBds,
                             _SimpleMeshVertex &maxBds,
                             std::vector<Vector3> &meshNormal,
                             const char *fileName,
                             std::vector<int> &OrientTriangles, const float fEnlargeFactor);

    void LoadMeshInOFFformat(const char *fileName, std::vector<int> &OrientTriangles, const float fEnlargeFactor);

public:
    std::vector<_SimpleMeshVertex> vecVertex; // vertex array
    std::vector<_SimpleMeshEdge> vecEdge; // edge array
    std::vector<_SimpleMeshTriangle> vecTriangle; // triangle array
    //
    std::vector<Vector3> *meshNormalPtr;
};

struct myPairCompare {
    bool operator()(const std::pair<int, int> &lhs, const std::pair<int, int> &rhs);

    bool operator()(const std::pair<int, int> &lhs, const std::pair<int, int> &rhs) const;
};
