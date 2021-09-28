/*
(c) 2012 Fengtao Fan
*/
#include "SimpleMesh.h"
#include <map>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <cstdlib>

//struct myPairCompare
//{
bool myPairCompare::operator()(const std::pair<int, int> &lhs, const std::pair<int, int> &rhs) {
    bool ret = false;
    if (lhs.first < rhs.first)
        ret = true;
    else if (lhs.first > rhs.first)
        ret = false;
    else if (lhs.second < rhs.second)
        ret = true;
    else
        ret = false;

    return ret;
}

bool myPairCompare::operator()(const std::pair<int, int> &lhs, const std::pair<int, int> &rhs) const {
    bool ret = false;
    if (lhs.first < rhs.first)
        ret = true;
    else if (lhs.first > rhs.first)
        ret = false;
    else if (lhs.second < rhs.second)
        ret = true;
    else
        ret = false;

    return ret;
}

//};
void
_SimpleMesh::LoadMeshInOFFformat(const char *fileName, std::vector<int> &OrientTriangles, const float fEnlargeFactor) {
    std::cout << "reading... " << std::endl;
    std::map<std::pair<int, int>, int, myPairCompare> edgeMapping;

    std::ifstream ifile;
    ifile.open(fileName, std::ifstream::in);
    long long fileSize = 0;
    char *fBuf = NULL;

    std::string sBuf;
    //
    int vertexNum = 0;
    int triangleNum = 0;
    //
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    std::stringstream locSstr(std::stringstream::in | std::stringstream::out);
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
        while (sstr.good()) {
            std::getline(sstr, sBuf);
            if (!sstr.good())
                break;
            //
            if (!sBuf.empty() && sBuf[0] != '#') {
                if (sBuf == "OFF") {
                    std::getline(sstr, sBuf);
                    if (!sstr.good())
                        break;
                    while (sBuf[0] < '0' || sBuf[0] > '9') {
                        std::getline(sstr, sBuf);
                    }
                    //
                    locSstr.str(sBuf);
                    locSstr >> vertexNum >> triangleNum;
                    //sscanf(sBuf.c_str(), "%d %d", &vertexNum, &triangleNum);
                    locSstr.clear();
                    vecVertex.reserve(vertexNum);
                    vecTriangle.reserve(triangleNum);
                    // reserve space for orientation
                    OrientTriangles.reserve(triangleNum * 3);
                } else {
                    if (vertexNum > 0) {
                        _SimpleMeshVertex tmpVer;
                        locSstr.str(sBuf);
                        locSstr >> tmpVer.x;
                        locSstr >> tmpVer.y;
                        locSstr >> tmpVer.z;
                        locSstr.clear();
                        //sscanf(sBuf.c_str(), "%f %f %f", &tmpVer.x, &tmpVer.y, &tmpVer.z);
                        tmpVer.x *= fEnlargeFactor;
                        tmpVer.y *= fEnlargeFactor;
                        tmpVer.z *= fEnlargeFactor;
                        vecVertex.push_back(tmpVer);
                        vertexNum--;
                    } else {
                        if (triangleNum > 0) {
                            _SimpleMeshTriangle tmpTri;
                            _SimpleMeshEdge tmpEdge;
                            locSstr.str(sBuf);
                            locSstr >> tmpTri.v0;
                            locSstr >> tmpTri.v0;
                            locSstr >> tmpTri.v1;
                            locSstr >> tmpTri.v2;
                            locSstr.clear();
                            //sscanf(sBuf.c_str(), "%d %d %d %d", &tmpTri.v0, &tmpTri.v0, &tmpTri.v1, &tmpTri.v2);

                            // store in the oriented triangles
                            OrientTriangles.push_back(tmpTri.v0);
                            OrientTriangles.push_back(tmpTri.v1);
                            OrientTriangles.push_back(tmpTri.v2);
                            // check the existence of tree edges

                            //
                            tmpTri.sortVertices();
                            //
                            std::pair<int, int> tmpEdgePair(tmpTri.v0, tmpTri.v1);
                            std::map<std::pair<int, int>, int, myPairCompare>::iterator mIter;
                            mIter = edgeMapping.find(tmpEdgePair);
                            if (mIter == edgeMapping.end()) {// new edge
                                tmpEdge.v0 = tmpEdgePair.first;
                                tmpEdge.v1 = tmpEdgePair.second;
                                tmpEdge.AdjTri[0] = vecTriangle.size();
                                tmpEdge.AdjTriNum = 1;
                                //
                                vecEdge.push_back(tmpEdge);
                                tmpTri.e01 = vecEdge.size() - 1;
                                edgeMapping[tmpEdgePair] = vecEdge.size() - 1;
                            } else {// existed already
                                vecEdge[mIter->second].AdjTriNum++;
                                vecEdge[mIter->second].AdjTri[1] = vecTriangle.size();
                                //
                                tmpTri.e01 = mIter->second;
                            }
                            //
                            tmpEdgePair.first = tmpTri.v1;
                            tmpEdgePair.second = tmpTri.v2;
                            mIter = edgeMapping.find(tmpEdgePair);
                            if (mIter == edgeMapping.end()) {// new edge
                                tmpEdge.v0 = tmpEdgePair.first;
                                tmpEdge.v1 = tmpEdgePair.second;
                                tmpEdge.AdjTri[0] = vecTriangle.size();
                                tmpEdge.AdjTriNum = 1;
                                //
                                vecEdge.push_back(tmpEdge);
                                tmpTri.e12 = vecEdge.size() - 1;
                                edgeMapping[tmpEdgePair] = vecEdge.size() - 1;
                            } else {// existed already
                                vecEdge[mIter->second].AdjTriNum++;
                                vecEdge[mIter->second].AdjTri[1] = vecTriangle.size();
                                //
                                tmpTri.e12 = mIter->second;
                            }
                            //
                            tmpEdgePair.first = tmpTri.v0;
                            tmpEdgePair.second = tmpTri.v2;
                            mIter = edgeMapping.find(tmpEdgePair);
                            if (mIter == edgeMapping.end()) {// new edge
                                tmpEdge.v0 = tmpEdgePair.first;
                                tmpEdge.v1 = tmpEdgePair.second;
                                tmpEdge.AdjTri[0] = vecTriangle.size();
                                tmpEdge.AdjTriNum = 1;
                                //
                                vecEdge.push_back(tmpEdge);
                                tmpTri.e02 = vecEdge.size() - 1;
                                edgeMapping[tmpEdgePair] = vecEdge.size() - 1;
                            } else {// existed already
                                vecEdge[mIter->second].AdjTriNum++;
                                vecEdge[mIter->second].AdjTri[1] = vecTriangle.size();
                                //
                                tmpTri.e02 = mIter->second;
                            }
                            //
                            vecTriangle.push_back(tmpTri);
                            triangleNum--;
                        }
                    }
                }
            }
        }
        sstr.clear();
    } else {
        std::cout << "Can NOT open file " << fileName << std::endl;
        exit(0);
    }
    // assign incident edges information to vertex
    for (int i = 0; i < int(vecEdge.size()); i++) {
        vecVertex[vecEdge[i].v0].adjEdges.push_back(i);
        vecVertex[vecEdge[i].v1].adjEdges.push_back(i);
    }
    std::cout << "Done... " << vertexNum << " " << triangleNum << std::endl;
    std::cout << "ver... " << vecVertex.size() << " tri " << vecTriangle.size()
              << " edge" << vecEdge.size() << std::endl;
    //
    edgeMapping.clear();
    return;
}

void _SimpleMesh::LoadMeshInOFFformat(_SimpleMeshVertex &minBd,
                                      _SimpleMeshVertex &maxBd,
                                      std::vector<Vector3> &meshNormal,
                                      const char *fileName,
                                      const float fEnlargeFactor) {
    std::cout << "reading... " << std::endl;
    std::map<std::pair<int, int>, int, myPairCompare> edgeMapping;

    std::ifstream ifile;
    ifile.open(fileName, std::ifstream::in);
    long long fileSize = 0;
    char *fBuf = NULL;

    std::string sBuf;
    //
    int vertexNum = 0;
    int triangleNum = 0;
    //
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    std::stringstream locSstr(std::stringstream::in | std::stringstream::out);
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


        while (sstr.good()) {
            std::getline(sstr, sBuf);
            if (!sstr.good())
                break;
            //
            if (!sBuf.empty() && sBuf[0] != '#') {
                if (sBuf == "OFF") {
                    std::getline(sstr, sBuf);
                    if (!sstr.good())
                        break;
                    while (sBuf[0] < '0' || sBuf[0] > '9') {
                        std::getline(sstr, sBuf);
                    }
                    //
                    locSstr.str(sBuf);
                    locSstr >> vertexNum >> triangleNum;
                    //sscanf(sBuf.c_str(), "%d %d", &vertexNum, &triangleNum);
                    locSstr.clear();
                    vecVertex.reserve(vertexNum);
                    vecTriangle.reserve(triangleNum);
                } else {
                    if (vertexNum > 0) {
                        _SimpleMeshVertex tmpVer;
                        locSstr.str(sBuf);
                        locSstr >> tmpVer.x;
                        locSstr >> tmpVer.y;
                        locSstr >> tmpVer.z;
                        locSstr.clear();
                        //sscanf(sBuf.c_str(), "%f %f %f", &tmpVer.x, &tmpVer.y, &tmpVer.z);
                        tmpVer.x *= fEnlargeFactor;
                        tmpVer.y *= fEnlargeFactor;
                        tmpVer.z *= fEnlargeFactor;
                        vecVertex.push_back(tmpVer);
                        if (vecVertex.size() == 1) {
                            minBd = tmpVer;
                            maxBd = tmpVer;
                        } else {
                            minBd.x = tmpVer.x < minBd.x ? tmpVer.x : minBd.x;
                            minBd.y = tmpVer.y < minBd.y ? tmpVer.y : minBd.y;
                            minBd.z = tmpVer.z < minBd.z ? tmpVer.z : minBd.z;
                            //
                            maxBd.x = tmpVer.x > maxBd.x ? tmpVer.x : maxBd.x;
                            maxBd.y = tmpVer.y > maxBd.y ? tmpVer.y : maxBd.y;
                            maxBd.z = tmpVer.z > maxBd.z ? tmpVer.z : maxBd.z;
                        }
                        vertexNum--;
                    } else {
                        if (triangleNum > 0) {
                            _SimpleMeshTriangle tmpTri;
                            _SimpleMeshEdge tmpEdge;
                            locSstr.str(sBuf);
                            locSstr >> tmpTri.v0;
                            locSstr >> tmpTri.v0;
                            locSstr >> tmpTri.v1;
                            locSstr >> tmpTri.v2;
                            locSstr.clear();
                            //sscanf(sBuf.c_str(), "%d %d %d %d", &tmpTri.v0, &tmpTri.v0, &tmpTri.v1, &tmpTri.v2);

                            // check the existence of tree edges
                            //
                            Vector3 leftVec, rightVec;
                            leftVec[0] = vecVertex[tmpTri.v2].x - vecVertex[tmpTri.v1].x;
                            leftVec[1] = vecVertex[tmpTri.v2].y - vecVertex[tmpTri.v1].y;
                            leftVec[2] = vecVertex[tmpTri.v2].z - vecVertex[tmpTri.v1].z;

                            rightVec[0] = vecVertex[tmpTri.v0].x - vecVertex[tmpTri.v1].x;
                            rightVec[1] = vecVertex[tmpTri.v0].y - vecVertex[tmpTri.v1].y;
                            rightVec[2] = vecVertex[tmpTri.v0].z - vecVertex[tmpTri.v1].z;
                            leftVec = leftVec ^ rightVec;
                            unitize(leftVec);
                            //leftVec = leftVec / norm(leftVec);
                            meshNormal.push_back(leftVec);
                            //
                            tmpTri.sortVertices();
                            //
                            std::pair<int, int> tmpEdgePair(tmpTri.v0, tmpTri.v1);
                            std::map<std::pair<int, int>, int, myPairCompare>::iterator mIter;
                            mIter = edgeMapping.find(tmpEdgePair);
                            if (mIter == edgeMapping.end()) {// new edge
                                tmpEdge.v0 = tmpEdgePair.first;
                                tmpEdge.v1 = tmpEdgePair.second;
                                tmpEdge.AdjTri[0] = vecTriangle.size();
                                tmpEdge.AdjTriNum = 1;
                                //
                                vecEdge.push_back(tmpEdge);
                                tmpTri.e01 = vecEdge.size() - 1;
                                edgeMapping[tmpEdgePair] = vecEdge.size() - 1;
                            } else {// existed already
                                vecEdge[mIter->second].AdjTriNum++;
                                vecEdge[mIter->second].AdjTri[1] = vecTriangle.size();
                                //
                                tmpTri.e01 = mIter->second;
                            }
                            //
                            tmpEdgePair.first = tmpTri.v1;
                            tmpEdgePair.second = tmpTri.v2;
                            mIter = edgeMapping.find(tmpEdgePair);
                            if (mIter == edgeMapping.end()) {// new edge
                                tmpEdge.v0 = tmpEdgePair.first;
                                tmpEdge.v1 = tmpEdgePair.second;
                                tmpEdge.AdjTri[0] = vecTriangle.size();
                                tmpEdge.AdjTriNum = 1;
                                //
                                vecEdge.push_back(tmpEdge);
                                tmpTri.e12 = vecEdge.size() - 1;
                                edgeMapping[tmpEdgePair] = vecEdge.size() - 1;
                            } else {// existed already
                                vecEdge[mIter->second].AdjTriNum++;
                                vecEdge[mIter->second].AdjTri[1] = vecTriangle.size();
                                //
                                tmpTri.e12 = mIter->second;
                            }
                            //
                            tmpEdgePair.first = tmpTri.v0;
                            tmpEdgePair.second = tmpTri.v2;
                            mIter = edgeMapping.find(tmpEdgePair);
                            if (mIter == edgeMapping.end()) {// new edge
                                tmpEdge.v0 = tmpEdgePair.first;
                                tmpEdge.v1 = tmpEdgePair.second;
                                tmpEdge.AdjTri[0] = vecTriangle.size();
                                tmpEdge.AdjTriNum = 1;
                                //
                                vecEdge.push_back(tmpEdge);
                                tmpTri.e02 = vecEdge.size() - 1;
                                edgeMapping[tmpEdgePair] = vecEdge.size() - 1;
                            } else {// existed already
                                vecEdge[mIter->second].AdjTriNum++;
                                vecEdge[mIter->second].AdjTri[1] = vecTriangle.size();
                                //
                                tmpTri.e02 = mIter->second;
                            }
                            //
                            vecTriangle.push_back(tmpTri);
                            triangleNum--;
                        }
                    }
                }
            }
        }
        sstr.clear();
    } else {
        std::cout << "Can NOT open file " << fileName << std::endl;
        exit(0);
    }
    // assign incident edges information to vertex
    for (int i = 0; i < int(vecEdge.size()); i++) {
        vecVertex[vecEdge[i].v0].adjEdges.push_back(i);
        vecVertex[vecEdge[i].v1].adjEdges.push_back(i);
    }
    std::cout << "Done... " << vertexNum << " " << triangleNum << std::endl;
    std::cout << "ver... " << vecVertex.size() << " tri " << vecTriangle.size()
              << " edge" << vecEdge.size() << std::endl;
    //
    edgeMapping.clear();
    return;
}

void _SimpleMesh::LoadMeshInOFFformat(_SimpleMeshVertex &minBd,
                                      _SimpleMeshVertex &maxBd,
                                      std::vector<Vector3> &meshNormal,
                                      const char *fileName,
                                      std::vector<int> &OrientTriangles,
                                      const float fEnlargeFactor) {
    std::cout << "Reading input file <" << fileName << ">" << std::endl;
    std::map<std::pair<int, int>, int, myPairCompare> edgeMapping;

    std::ifstream ifile;
    ifile.open(fileName, std::ifstream::in);
    long long fileSize = 0;
    char *fBuf = NULL;

    std::string sBuf;
    //
    int vertexNum = 0;
    int triangleNum = 0;
    //
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    std::stringstream locSstr(std::stringstream::in | std::stringstream::out);
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


        while (sstr.good()) {
            std::getline(sstr, sBuf);
            if (!sstr.good())
                break;
            //
            if (!sBuf.empty() && sBuf[0] != '#' && sBuf[0] != '\n' &&
                sBuf[0] != '\r') // text file with extra character from
            {
                //if (sBuf == "OFF")
                if (sBuf[0] == 'O') {
                    std::getline(sstr, sBuf);
                    if (!sstr.good())
                        break;
                    while (sBuf[0] < '0' || sBuf[0] > '9') {
                        std::getline(sstr, sBuf);
                    }
                    //
                    locSstr.str(sBuf);
                    locSstr >> vertexNum >> triangleNum;
                    //sscanf(sBuf.c_str(), "%d %d", &vertexNum, &triangleNum);
                    locSstr.clear();
                    vecVertex.reserve(vertexNum);
                    vecTriangle.reserve(triangleNum);
                } else {
                    if (vertexNum > 0) {
                        _SimpleMeshVertex tmpVer;
                        locSstr.str(sBuf);
                        locSstr >> tmpVer.x;
                        locSstr >> tmpVer.y;
                        locSstr >> tmpVer.z;
                        locSstr.clear();
                        //sscanf(sBuf.c_str(), "%f %f %f", &tmpVer.x, &tmpVer.y, &tmpVer.z);
                        tmpVer.x *= fEnlargeFactor;
                        tmpVer.y *= fEnlargeFactor;
                        tmpVer.z *= fEnlargeFactor;
                        vecVertex.push_back(tmpVer);
                        if (vecVertex.size() == 1) {
                            minBd = tmpVer;
                            maxBd = tmpVer;
                        } else {
                            minBd.x = tmpVer.x < minBd.x ? tmpVer.x : minBd.x;
                            minBd.y = tmpVer.y < minBd.y ? tmpVer.y : minBd.y;
                            minBd.z = tmpVer.z < minBd.z ? tmpVer.z : minBd.z;
                            //
                            maxBd.x = tmpVer.x > maxBd.x ? tmpVer.x : maxBd.x;
                            maxBd.y = tmpVer.y > maxBd.y ? tmpVer.y : maxBd.y;
                            maxBd.z = tmpVer.z > maxBd.z ? tmpVer.z : maxBd.z;
                        }
                        vertexNum--;
                    } else {
                        if (triangleNum > 0) {
                            _SimpleMeshTriangle tmpTri;
                            _SimpleMeshEdge tmpEdge;
                            locSstr.str(sBuf);
                            locSstr >> tmpTri.v0;
                            locSstr >> tmpTri.v0;
                            locSstr >> tmpTri.v1;
                            locSstr >> tmpTri.v2;
                            locSstr.clear();
//							if (triangleNum > 6100)
//							{
//							    std::cout << sBuf << " \n ";
//							    std::cout << tmpTri.v0 << " " << tmpTri.v1 << " " << tmpTri.v2;
//							}

                            //sscanf(sBuf.c_str(), "%d %d %d %d", &tmpTri.v0, &tmpTri.v0, &tmpTri.v1, &tmpTri.v2);
                            // store in the oriented triangles
                            OrientTriangles.push_back(tmpTri.v0);
                            OrientTriangles.push_back(tmpTri.v1);
                            OrientTriangles.push_back(tmpTri.v2);
                            // check the existence of tree edges
                            //
                            Vector3 leftVec, rightVec;
                            leftVec[0] = vecVertex[tmpTri.v2].x - vecVertex[tmpTri.v1].x;
                            leftVec[1] = vecVertex[tmpTri.v2].y - vecVertex[tmpTri.v1].y;
                            leftVec[2] = vecVertex[tmpTri.v2].z - vecVertex[tmpTri.v1].z;

                            rightVec[0] = vecVertex[tmpTri.v0].x - vecVertex[tmpTri.v1].x;
                            rightVec[1] = vecVertex[tmpTri.v0].y - vecVertex[tmpTri.v1].y;
                            rightVec[2] = vecVertex[tmpTri.v0].z - vecVertex[tmpTri.v1].z;
                            leftVec = leftVec ^ rightVec;
                            unitize(leftVec);
                            //leftVec = leftVec / norm(leftVec);
                            meshNormal.push_back(leftVec);
                            //
                            tmpTri.sortVertices();
                            //
                            std::pair<int, int> tmpEdgePair(tmpTri.v0, tmpTri.v1);
                            std::map<std::pair<int, int>, int, myPairCompare>::iterator mIter;
                            mIter = edgeMapping.find(tmpEdgePair);
                            if (mIter == edgeMapping.end()) {// new edge
                                tmpEdge.v0 = tmpEdgePair.first;
                                tmpEdge.v1 = tmpEdgePair.second;
                                tmpEdge.AdjTri[0] = vecTriangle.size();
                                tmpEdge.AdjTriNum = 1;
                                //
                                vecEdge.push_back(tmpEdge);
                                tmpTri.e01 = vecEdge.size() - 1;
                                edgeMapping[tmpEdgePair] = vecEdge.size() - 1;
                            } else {// existed already
                                vecEdge[mIter->second].AdjTriNum++;
                                vecEdge[mIter->second].AdjTri[1] = vecTriangle.size();
                                //
                                tmpTri.e01 = mIter->second;
                            }
                            //
                            tmpEdgePair.first = tmpTri.v1;
                            tmpEdgePair.second = tmpTri.v2;
                            mIter = edgeMapping.find(tmpEdgePair);
                            if (mIter == edgeMapping.end()) {// new edge
                                tmpEdge.v0 = tmpEdgePair.first;
                                tmpEdge.v1 = tmpEdgePair.second;
                                tmpEdge.AdjTri[0] = vecTriangle.size();
                                tmpEdge.AdjTriNum = 1;
                                //
                                vecEdge.push_back(tmpEdge);
                                tmpTri.e12 = vecEdge.size() - 1;
                                edgeMapping[tmpEdgePair] = vecEdge.size() - 1;
                            } else {// existed already
                                vecEdge[mIter->second].AdjTriNum++;
                                vecEdge[mIter->second].AdjTri[1] = vecTriangle.size();
                                //
                                tmpTri.e12 = mIter->second;
                            }
                            //
                            tmpEdgePair.first = tmpTri.v0;
                            tmpEdgePair.second = tmpTri.v2;
                            mIter = edgeMapping.find(tmpEdgePair);
                            if (mIter == edgeMapping.end()) {// new edge
                                tmpEdge.v0 = tmpEdgePair.first;
                                tmpEdge.v1 = tmpEdgePair.second;
                                tmpEdge.AdjTri[0] = vecTriangle.size();
                                tmpEdge.AdjTriNum = 1;
                                //
                                vecEdge.push_back(tmpEdge);
                                tmpTri.e02 = vecEdge.size() - 1;
                                edgeMapping[tmpEdgePair] = vecEdge.size() - 1;
                            } else {// existed already
                                vecEdge[mIter->second].AdjTriNum++;
                                vecEdge[mIter->second].AdjTri[1] = vecTriangle.size();
                                //
                                tmpTri.e02 = mIter->second;
                            }
                            //
                            vecTriangle.push_back(tmpTri);
                            triangleNum--;
                        }
                    }
                }
            }
        }
        sstr.clear();
    } else {
        std::cout << "ERROR: can NOT open file <" << fileName << ">" << std::endl;
        exit(0);
    }
    // assign incident edges information to vertex
    for (int i = 0; i < int(vecEdge.size()); i++) {
        vecVertex[vecEdge[i].v0].adjEdges.push_back(i);
        vecVertex[vecEdge[i].v1].adjEdges.push_back(i);
    }
    //
    edgeMapping.clear();
    //
    return;
}
