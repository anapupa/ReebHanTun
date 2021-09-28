/*
(c) 2012 Fengtao Fan
*/
#include "FilesOutputForOptimalCycles.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;
using namespace std;

/* Matrix inversion routine.
Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const matrix<T> &input, matrix<T> &inverse) {
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    matrix<T> A(input);

    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A, pm);
    if (res != 0)
        return false;

    // create identity matrix of "inverse"
    inverse.assign(identity_matrix<T>(A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
}

void AssembleMatrix(matrix<double> &A, std::vector<std::vector<std::pair<int, int> > > &sparseLinkMat) {
    for (unsigned int i = 0; i < sparseLinkMat.size(); i++) {
        for (unsigned int j = 0; j < sparseLinkMat.size(); j++) {
            A(i, j) = 0.0;
        }
    }
    for (unsigned int i = 0; i < sparseLinkMat.size(); i++) {
        for (unsigned int j = 0; j < sparseLinkMat[i].size(); j++) {
            A(i, sparseLinkMat[i][j].first) = sparseLinkMat[i][j].second;
        }
    }
    return;
}

void FilesOutputForOptimalCycles::WriteOrientedComplex(const char *out_file_name) {
    // edge orientation is taken as v0-->--v1 in the vecEdge
    std::cout << "writing file ... " << out_file_name << std::endl;

    std::ofstream ofile;
    ofile.open(out_file_name, std::ifstream::out);

    //
    int vertexNum = int(inMeshPtr->vecVertex.size());
    int triangleNum = int(inMeshPtr->vecTriangle.size());
    int edgeNum = int(inMeshPtr->vecEdge.size());
    //
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    //std::stringstream locSstr(std::stringstream::in | std::stringstream::out);
    if (ofile.is_open()) {
        // total number of simplices
        sstr << vertexNum + triangleNum + edgeNum << std::endl;
        // number of skeletons
        sstr << 3 << std::endl;
        // number of vertices
        sstr << vertexNum << std::endl;
        // write all coordinates
        for (int i = 0; i < vertexNum; i++) {
            sstr << inMeshPtr->vecVertex[i].x << " "
                 << inMeshPtr->vecVertex[i].y << " "
                 << inMeshPtr->vecVertex[i].z << std::endl;
        }
        // write the number of edges
        sstr << edgeNum << std::endl;
        // write the edge endpoints in the orientation of v0-->--v1 in the vecEdge
        for (int i = 0; i < edgeNum; i++) {
            sstr << inMeshPtr->vecEdge[i].v0 << " " << inMeshPtr->vecEdge[i].v1 << std::endl;
        }
        // write the number of triangles
        sstr << triangleNum << std::endl;
        // write the triangle endpoints in the orientation of v0--v1--v2 in the vecTriangle
        for (int i = 0; i < triangleNum; i++) {
            sstr << inMeshPtr->vecTriangle[i].e01 << " "
                 << inMeshPtr->vecTriangle[i].e12 << " "
                 << inMeshPtr->vecTriangle[i].e02 << std::endl;
        }

        //
        //ofile.write(sstr.str().c_str(), 10);
        ofile << sstr.rdbuf();
        //
        // delete string stream content
        sstr.str("");
        sstr.clear();
        // close file
        ofile.close();
    } else {
        std::cout << "Can NOT open file " << out_file_name << std::endl;
        exit(0);
    }
    std::cout << "Done... " << vertexNum << " " << edgeNum << " " << triangleNum << std::endl;
    //
    return;
}

void FilesOutputForOptimalCycles::WritePolygonOut(const char *out_file_name, std::vector<Vector3> &inPly) {
    std::cout << "writing file ... " << out_file_name << std::endl;

    std::ofstream ofile;
    ofile.open(out_file_name, std::ifstream::out);

    //
    int vertexNum = int(inPly.size());
    //
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    //std::stringstream locSstr(std::stringstream::in | std::stringstream::out);
    if (ofile.is_open()) {
        // total number of edges
        sstr << vertexNum << " ";
        for (int i = 0; i < vertexNum; i++) {
            sstr << inPly[i][0] << " "
                 << inPly[i][1] << " "
                 << inPly[i][2] << std::endl;
        }

        //
        //ofile.write(sstr.str().c_str(), 10);
        ofile << sstr.rdbuf();
        //
        // delete string stream content
        sstr.str("");
        sstr.clear();
        // close file
        ofile.close();
    } else {
        std::cout << "Can NOT open file " << out_file_name << std::endl;
        exit(0);
    }
    std::cout << "Done... " << vertexNum << " " << std::endl;
}

void FilesOutputForOptimalCycles::ReadPolygonIn(const char *in_file_name, std::vector<Vector3> &wantPly) {
    std::cout << "reading... " << std::endl;

    std::ifstream ifile;
    ifile.open(in_file_name, std::ifstream::in);
    long long fileSize = 0;
    char *fBuf = NULL;

    std::string sBuf;
    //
    int vertexNum = 0;
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

        sstr >> vertexNum;
        for (int i = 0; i < vertexNum; i++) {
            sstr >> tempVec[0] >> tempVec[1] >> tempVec[2];
        }

        sstr.clear();
    } else {
        std::cout << "Can NOT open file " << in_file_name << std::endl;
        exit(0);
    }
    std::cout << "Done... " << vertexNum << " " << std::endl;
    return;

}

void FilesOutputForOptimalCycles::WriteOrientedCycle(const int in_cycle_id, const char *out_file_name) {
    // edge orientation is taken as v0-->--v1 in the vecEdge
    int cycle_id = in_cycle_id;
    std::vector<std::vector<int> > *vecVertex_Ptr = vecInitVerticalLoops_Vertex_Ptr;
    std::vector<std::vector<int> > *vecEdge_Ptr = vecInitVerticalLoops_Edge_Ptr;

    if (cycle_id >= vecInitVerticalLoops_Vertex_Ptr->size()) {
        cycle_id -= int(vecInitVerticalLoops_Vertex_Ptr->size());
        vecVertex_Ptr = vecInitHorizontalLoops_Vertex_Ptr;
        vecEdge_Ptr = vecInitHorizontalLoops_Edge_Ptr;
    }

    //if (in_cycle_id == 1)
    //{
    //	int own_orientation = -1;
    //	std::cout << "writing file ... " << out_file_name << std::endl;

    //	std::ofstream ofile;
    //	ofile.open(out_file_name, std::ifstream::out);

    //	//
    //	int vertexNum = int((*vecVertex_Ptr)[cycle_id].size());
    //	int edgeNum = int ((*vecEdge_Ptr)[cycle_id].size());
    //	//
    //	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //	//
    //	//std::stringstream locSstr(std::stringstream::in | std::stringstream::out);
    //	if (ofile.is_open())
    //	{
    //		// total number of edges
    //		sstr << edgeNum + (*vecInitHorizontalLoops_Edge_Ptr)[1].size() << " ";
    //		for (int i = 0; i < edgeNum; i++)
    //		{
    //			// vi--ei--v(i+1)
    //			// check the orientation
    //			if( (inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][i]].v0 != (*vecVertex_Ptr)[cycle_id][i] &&
    //				inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][i]].v0 != (*vecVertex_Ptr)[cycle_id][i+1] ) ||
    //				(inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][i]].v1 != (*vecVertex_Ptr)[cycle_id][i] &&
    //				inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][i]].v1 != (*vecVertex_Ptr)[cycle_id][i+1] ) )
    //			{
    //				std::cout << "EDGE doesn't match vertices " << std::endl;
    //				exit(0);
    //			}
    //			if (inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][i]].v0 == (*vecVertex_Ptr)[cycle_id][i])
    //			{// vi-->---v(i+1) has the same orientation as ei
    //				sstr << 1 * own_orientation << " " ;
    //			}
    //			else
    //			{// opposite orientation
    //				sstr << -1 * own_orientation << " ";
    //			}
    //			// write the edge
    //			sstr << (*vecEdge_Ptr)[cycle_id][i] << " ";
    //		}

    //		for (int i = 0; i < (*vecInitHorizontalLoops_Edge_Ptr)[1].size(); i++)
    //		{
    //			// vi--ei--v(i+1)
    //			// check the orientation

    //			if (inMeshPtr->vecEdge[(*vecInitHorizontalLoops_Edge_Ptr)[1][i]].v0 == (*vecInitHorizontalLoops_Vertex_Ptr)[1][i])
    //			{// vi-->---v(i+1) has the same orientation as ei
    //				sstr << 1  << " " ;
    //			}
    //			else
    //			{// opposite orientation
    //				sstr << -1  << " ";
    //			}
    //			// write the edge
    //			sstr << (*vecInitHorizontalLoops_Edge_Ptr)[1][i] << " ";
    //		}
    //
    //		//
    //		//ofile.write(sstr.str().c_str(), 10);
    //		ofile << sstr.rdbuf();
    //		//
    //		// delete string stream content
    //		sstr.str("");
    //		sstr.clear();
    //		// close file
    //		ofile.close();
    //	}
    //	else
    //	{
    //		std::cout << "Can NOT open file " << out_file_name << std::endl;
    //		exit(0);
    //	}
    //	std::cout << "Done... " << vertexNum << " "  << edgeNum << " " << std::endl;
    //	return;

    //}
    if (in_cycle_id >= 0 && in_cycle_id < 2 * vecInitVerticalLoops_Vertex_Ptr->size()) {
        std::cout << "writing file ... " << out_file_name << std::endl;

        std::ofstream ofile;
        ofile.open(out_file_name, std::ifstream::out);

        //
        int vertexNum = int((*vecVertex_Ptr)[cycle_id].size());
        int edgeNum = int((*vecEdge_Ptr)[cycle_id].size());
        //
        std::stringstream sstr(std::stringstream::in | std::stringstream::out);
        //
        //std::stringstream locSstr(std::stringstream::in | std::stringstream::out);
        if (ofile.is_open()) {
            // total number of edges
            sstr << edgeNum << " ";
            for (int i = 0; i < edgeNum; i++) {
                // vi--ei--v(i+1)
                // check the orientation
                if ((inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][i]].v0 != (*vecVertex_Ptr)[cycle_id][i] &&
                     inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][i]].v0 != (*vecVertex_Ptr)[cycle_id][i + 1]) ||
                    (inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][i]].v1 != (*vecVertex_Ptr)[cycle_id][i] &&
                     inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][i]].v1 != (*vecVertex_Ptr)[cycle_id][i + 1])) {
                    std::cout << "EDGE doesn't match vertices " << std::endl;
                    exit(0);
                }
                if (inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][i]].v0 ==
                    (*vecVertex_Ptr)[cycle_id][i]) {// vi-->---v(i+1) has the same orientation as ei
                    sstr << 1 << " ";
                } else {// opposite orientation
                    sstr << -1 << " ";
                }
                // write the edge
                sstr << (*vecEdge_Ptr)[cycle_id][i] << " ";
            }

            //
            //ofile.write(sstr.str().c_str(), 10);
            ofile << sstr.rdbuf();
            //
            // delete string stream content
            sstr.str("");
            sstr.clear();
            // close file
            ofile.close();
        } else {
            std::cout << "Can NOT open file " << out_file_name << std::endl;
            exit(0);
        }
        std::cout << "Done... " << vertexNum << " " << edgeNum << " " << std::endl;
    } else {
        std::cout << "cycle index " << cycle_id << " out of RANGE " << std::endl;
    }
    //
    return;
}

void FilesOutputForOptimalCycles::WriteCoefficientMatrix(const char *out_file_name) {
    //int tot_genus = coeffMatrix->size() / 2 ;

    //std::cout << "M=zeros("<< tot_genus * 2 << "," << tot_genus * 2 << ");" << std::endl;
    //for (int i = 0; i < tot_genus; i++)
    //{
    //	//std::cout << "Vert Loop " << i << " link against " << std::endl;
    //	for (unsigned int j = 0; j < (*coeffMatrix)[i].size(); j++)
    //	{
    //		std::cout << "M(" << i + 1 << "," << (*coeffMatrix)[i][j].first + 1 << ")=" << (*coeffMatrix)[i][j].second << ";";
    //
    //	}
    //	std::cout << std::endl;
    //}
    //for (int i = 0; i < tot_genus; i++)
    //{
    //	//std::cout << "Hori Loop " << i + tot_genus << " link against " << std::endl;
    //	for (unsigned int j = 0; j < (*coeffMatrix)[i + tot_genus].size(); j++)
    //	{
    //		std::cout << "M(" << i + tot_genus + 1 << "," << (*coeffMatrix)[i + tot_genus][j].first + 1<< ")=" << (*coeffMatrix)[i + tot_genus][j].second << ";";
    //
    //	}
    //	std::cout << std::endl;
    //}
}

void FilesOutputForOptimalCycles::OrderCycles(std::set<int> &edge_set, std::map<int, int> &vertex_ordering) {
    vertex_ordering.clear();
    //
    int verCounter = 0;
    std::map<int, int>::iterator mFindIter;
    for (std::set<int>::iterator sIter = edge_set.begin();
         sIter != edge_set.end(); sIter++) {
        int node_a = inMeshPtr->vecEdge[*sIter].v0;
        int node_b = inMeshPtr->vecEdge[*sIter].v1;
        //
        mFindIter = vertex_ordering.find(node_a);
        if (mFindIter == vertex_ordering.end()) {
            vertex_ordering[node_a] = verCounter++;
        }
        //
        mFindIter = vertex_ordering.find(node_b);
        if (mFindIter == vertex_ordering.end()) {
            vertex_ordering[node_b] = verCounter++;
        }
    }
    return;
}

bool FilesOutputForOptimalCycles::OrderCycles(std::set<int> &edge_set, std::vector<int> &vertex_set) {
    int bSuccess = true;
    vertex_set.clear();
    vertex_set.reserve(edge_set.size() + 1);
    //
    std::map<int, std::pair<int, int> > vertex_edge_conn;
    std::map<int, std::pair<int, int> >::iterator mFindIter;
    for (std::set<int>::iterator sIter = edge_set.begin();
         sIter != edge_set.end(); sIter++) {
        int node_a = inMeshPtr->vecEdge[*sIter].v0;
        int node_b = inMeshPtr->vecEdge[*sIter].v1;
        //
        mFindIter = vertex_edge_conn.find(node_a);
        if (mFindIter == vertex_edge_conn.end()) {
            vertex_edge_conn[node_a] = std::pair<int, int>(*sIter, -1);
        } else {
            mFindIter->second.second = *sIter;
        }
        //
        mFindIter = vertex_edge_conn.find(node_b);
        if (mFindIter == vertex_edge_conn.end()) {
            vertex_edge_conn[node_b] = std::pair<int, int>(*sIter, -1);
        } else {
            mFindIter->second.second = *sIter;
        }
    }
    //
    std::set<int>::iterator sFindIter;
    int current_edge = *edge_set.begin();
    int end_vertex = inMeshPtr->vecEdge[current_edge].v0;
    int current_vertex = inMeshPtr->vecEdge[current_edge].v1;
    //
    while (current_vertex != end_vertex) {
        vertex_set.push_back(current_vertex);

        mFindIter = vertex_edge_conn.find(current_vertex);
        if (mFindIter->second.second < 0) {
            //std::cout << "NOT a closed loop" << std::endl;
            bSuccess = false;
            vertex_set.clear();
            break;
        }
        if (mFindIter->second.first == current_edge) {
            current_edge = mFindIter->second.second;
        } else {
            current_edge = mFindIter->second.first;
        }
        //
        if (current_vertex == inMeshPtr->vecEdge[current_edge].v0)
            current_vertex = inMeshPtr->vecEdge[current_edge].v1;
        else
            current_vertex = inMeshPtr->vecEdge[current_edge].v0;
    }
    if (bSuccess)
        vertex_set.push_back(current_vertex);
    return bSuccess;
}

void FilesOutputForOptimalCycles::WriteCyclesInformation(const char *out_file_prefix,
                                                         std::vector<std::set<int> > &v_basis_loops,
                                                         std::vector<std::set<int> > &h_basis_loops) {
//
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    string file_prefix(out_file_prefix);
    string file_name;
    // assemble the matrix
    const int tot_genus = v_basis_loops.size();
    //
    file_name = "loops_" + file_prefix + ".lop";
    //
    std::cout << file_name << " \n";
    //
    std::ofstream ofile;
    ofile.open(file_name.c_str(), std::ifstream::out);
    //ofile.open(file_name.c_str());
    if (ofile.is_open()) {
        sstr << "# there are " << tot_genus << " handle loops" << std::endl << std::endl;
        std::vector<int> cycleVertex;
        for (unsigned int i = 0; i < tot_genus; i++) {

            //if (!OrderCycles(h_basis_loops[i], cycleVertex))
            {
                //Vertices are not ordered
                //Only the set of vertices are there
                std::set<int> verSet;
                for (std::set<int>::iterator sIter = h_basis_loops[i].begin();
                     sIter != h_basis_loops[i].end(); sIter++) {
                    verSet.insert(inMeshPtr->vecEdge[*sIter].v0);
                    verSet.insert(inMeshPtr->vecEdge[*sIter].v1);
                }
                cycleVertex.assign(verSet.begin(), verSet.end());
            }
            sstr << "handle loop " << i << "size(" << cycleVertex.size() << "): ";
            for (unsigned int vid = 0; vid < cycleVertex.size(); vid++) {
                sstr << cycleVertex[vid] << ",";
            }
            sstr << std::endl;
        }
        //
        sstr << std::endl << "# there are " << tot_genus << " tunnel loops" << std::endl << std::endl;
        for (unsigned int i = 0; i < tot_genus; i++) {

            //if (!OrderCycles(v_basis_loops[i], cycleVertex))
            {
                //Vertices are not ordered
                //Only the set of vertices are there
                std::set<int> verSet;
                for (std::set<int>::iterator sIter = v_basis_loops[i].begin();
                     sIter != v_basis_loops[i].end(); sIter++) {
                    verSet.insert(inMeshPtr->vecEdge[*sIter].v0);
                    verSet.insert(inMeshPtr->vecEdge[*sIter].v1);
                }
                cycleVertex.assign(verSet.begin(), verSet.end());
            }
            sstr << "tunnel loop " << i << "size(" << cycleVertex.size() << "): ";
            for (unsigned int vid = 0; vid < cycleVertex.size(); vid++) {
                sstr << cycleVertex[vid] << ",";
            }
            sstr << std::endl;
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
        std::cout << "can not open " << file_name << std::endl;
        exit(0);
    }
    return;
}

void FilesOutputForOptimalCycles::WriteGeomviewListFormat(const char *out_file_prefix,
                                                          std::vector<std::set<int> > &v_basis_loops,
                                                          std::vector<std::set<int> > &h_basis_loops,
                                                          std::vector<int> &vec_oriented_triangles,
                                                          const int orgVertexSize, const int orgTriangleSize,
                                                          const float fEnlargeFactor) {
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    string file_prefix(out_file_prefix);
    string file_name;
    // assemble the matrix
    const int tot_genus = v_basis_loops.size();
    //
    file_name = "loops_" + file_prefix + ".list";
    std::cout << file_name << std::endl;
    std::ofstream ofile;
    bool full_mesh_flag = false;
    ofile.open(file_name.c_str(), std::ifstream::out);
    if (ofile.is_open()) {
        sstr << "LIST" << std::endl;
        sstr << "{" << std::endl;
        sstr << "# there are " << tot_genus << " handle loops" << std::endl << std::endl;
        std::map<int, int> cycleVertex;
        for (unsigned int i = 0; i < tot_genus; i++) {

            OrderCycles(h_basis_loops[i], cycleVertex);
            sstr << "# handle loop " << i << std::endl;
            sstr << "appearance {linewidth 6}" << std::endl;
            sstr << "{" << std::endl;
            sstr << "OFF" << std::endl;
            sstr << cycleVertex.size() << " " << h_basis_loops[i].size() << " 0" << std::endl;
            //
            std::vector<int> vecCycVertex(cycleVertex.size());
            for (std::map<int, int>::iterator mIter = cycleVertex.begin();
                 mIter != cycleVertex.end(); mIter++) {
                vecCycVertex[mIter->second] = mIter->first;
                if (!full_mesh_flag && mIter->first >= orgVertexSize)
                    full_mesh_flag = true;
            }
            for (std::vector<int>::iterator vIter = vecCycVertex.begin();
                 vIter != vecCycVertex.end(); vIter++) {
                sstr << inMeshPtr->vecVertex[*vIter].x * fEnlargeFactor << " ";
                sstr << inMeshPtr->vecVertex[*vIter].y * fEnlargeFactor << " ";
                sstr << inMeshPtr->vecVertex[*vIter].z * fEnlargeFactor << std::endl;
            }
            for (std::set<int>::iterator sIter = h_basis_loops[i].begin();
                 sIter != h_basis_loops[i].end(); sIter++) {
                sstr << "2 ";
                sstr << cycleVertex[inMeshPtr->vecEdge[*sIter].v0] << " ";
                sstr << cycleVertex[inMeshPtr->vecEdge[*sIter].v1] << " ";
                sstr << "0.000000 0.800000 0.000000 1.0" << std::endl;
            }
            sstr << "}" << std::endl;
        }
        //
        sstr << std::endl << "# there are " << tot_genus << " tunnel loops" << std::endl << std::endl;
        for (unsigned int i = 0; i < tot_genus; i++) {

            OrderCycles(v_basis_loops[i], cycleVertex);
            sstr << "# tunnel loop " << i << std::endl;
            sstr << "appearance {linewidth 6}" << std::endl;
            sstr << "{" << std::endl;
            sstr << "OFF" << std::endl;
            sstr << cycleVertex.size() << " " << v_basis_loops[i].size() << " 0" << std::endl;
            //
            std::vector<int> vecCycVertex(cycleVertex.size());
            for (std::map<int, int>::iterator mIter = cycleVertex.begin();
                 mIter != cycleVertex.end(); mIter++) {
                vecCycVertex[mIter->second] = mIter->first;
                if (!full_mesh_flag && mIter->first >= orgVertexSize)
                    full_mesh_flag = true;
            }
            for (std::vector<int>::iterator vIter = vecCycVertex.begin();
                 vIter != vecCycVertex.end(); vIter++) {
                sstr << inMeshPtr->vecVertex[*vIter].x * fEnlargeFactor << " ";
                sstr << inMeshPtr->vecVertex[*vIter].y * fEnlargeFactor << " ";
                sstr << inMeshPtr->vecVertex[*vIter].z * fEnlargeFactor << std::endl;
            }
            for (std::set<int>::iterator sIter = v_basis_loops[i].begin();
                 sIter != v_basis_loops[i].end(); sIter++) {
                sstr << "2 ";
                sstr << cycleVertex[inMeshPtr->vecEdge[*sIter].v0] << " ";
                sstr << cycleVertex[inMeshPtr->vecEdge[*sIter].v1] << " ";
                sstr << "0.800000 0.000000 0.000000 1.0" << std::endl;
            }
            //
            sstr << "}" << std::endl;
        }
        //
        // write the mesh back to file
        //
        sstr << "# the mesh" << std::endl;
        sstr << "appearance {+transparent linewidth 1}" << std::endl;
        sstr << "{" << std::endl;
        // vertices
        int outVertexSize = orgVertexSize;
        int outTriangleSize = orgTriangleSize;
        if (full_mesh_flag) {
            outVertexSize = inMeshPtr->vecVertex.size();
            outTriangleSize = inMeshPtr->vecTriangle.size();
        }
        sstr << "OFF" << std::endl;
        sstr << outVertexSize << " " << outTriangleSize << " 0" << std::endl;

        for (unsigned int i = 0; i < outVertexSize; i++) {
            sstr << inMeshPtr->vecVertex[i].x * fEnlargeFactor << " ";
            sstr << inMeshPtr->vecVertex[i].y * fEnlargeFactor << " ";
            sstr << inMeshPtr->vecVertex[i].z * fEnlargeFactor << std::endl;
        }
        // triangles
        for (unsigned int i = 0; i < outTriangleSize; i++) {
            sstr << "3 ";
            sstr << vec_oriented_triangles[3 * i] << " ";
            sstr << vec_oriented_triangles[3 * i + 1] << " ";
            sstr << vec_oriented_triangles[3 * i + 2] << " ";
            sstr << "0.9 0.9 0.9 0.4" << std::endl;
        }
        sstr << "}" << std::endl;
        sstr << "}" << std::endl; // pair with LIST
        ofile << sstr.rdbuf();
        //
        ofile.close();
        ofile.clear();
        //
        sstr.str("");
        sstr.clear();
    } else {
        std::cout << "can not open " << file_name << std::endl;
        exit(0);
    }
}

void FilesOutputForOptimalCycles::WriteAllCyclesOut(std::vector<bool> &inVertLoopType, const char *out_file_prefix) {
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    string file_prefix(out_file_prefix);
    string file_name;
    // assemble the matrix
    const int tot_genus = vecInitVerticalLoops_Vertex_Ptr->size();
    //
    file_name = file_prefix + "_all_v.txt";
    std::ofstream ofile;
    ofile.open(file_name.c_str(), std::ifstream::out);
    if (ofile.is_open()) {
        sstr << tot_genus << std::endl;
        for (unsigned int i = 0; i < tot_genus; i++) {
            int edges_num = 0;
            int row_num = i;
            //
            if (!inVertLoopType[i])
                row_num += tot_genus;
            //
            for (std::list<std::pair<int, int> >::iterator listIter = (*coeffMatrix)[row_num].begin();
                 listIter != (*coeffMatrix)[row_num].end(); listIter++) {
                int col_num = listIter->first;
                int counter = abs(listIter->second);
                if (col_num >= tot_genus) {// it is a horizontal loop
                    edges_num += counter * (*vecInitHorizontalLoops_Edge_Ptr)[col_num - tot_genus].size();
                } else {
                    edges_num += counter * (*vecInitVerticalLoops_Edge_Ptr)[col_num].size();
                }
            }
            //
            sstr << edges_num << std::endl;
            for (std::list<std::pair<int, int> >::iterator listIter = (*coeffMatrix)[row_num].begin();
                 listIter != (*coeffMatrix)[row_num].end(); listIter++) {
                int col_num = listIter->first;
                int counter = abs(listIter->second);
                if (col_num >= tot_genus) {// it is a horizontal loop
                    for (int j = 0; j < counter; j++) {
                        for (std::vector<int>::iterator vIter = (*vecInitHorizontalLoops_Edge_Ptr)[col_num -
                                                                                                   tot_genus].begin();
                             vIter != (*vecInitHorizontalLoops_Edge_Ptr)[col_num - tot_genus].end(); vIter++) {
                            sstr << *vIter << " ";
                        }
                    }
                } else {
                    for (int j = 0; j < counter; j++) {
                        for (std::vector<int>::iterator vIter = (*vecInitVerticalLoops_Edge_Ptr)[col_num].begin();
                             vIter != (*vecInitVerticalLoops_Edge_Ptr)[col_num].end(); vIter++) {
                            sstr << *vIter << " ";
                        }
                    }
                }
            }
            sstr << std::endl;
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
        std::cout << "can not open " << file_name << std::endl;
        exit(0);
    }
    //
    file_name = file_prefix + "_all_h.txt";
    ofile.open(file_name.c_str(), std::ifstream::out);
    if (ofile.is_open()) {
        sstr << tot_genus << std::endl;
        for (unsigned int i = 0; i < tot_genus; i++) {
            int edges_num = 0;
            int row_num = i;
            //
            if (inVertLoopType[i])
                row_num += tot_genus;
            //
            for (std::list<std::pair<int, int> >::iterator listIter = (*coeffMatrix)[row_num].begin();
                 listIter != (*coeffMatrix)[row_num].end(); listIter++) {
                int col_num = listIter->first;
                int counter = abs(listIter->second);
                if (col_num >= tot_genus) {// it is a horizontal loop
                    edges_num += counter * (*vecInitHorizontalLoops_Edge_Ptr)[col_num - tot_genus].size();
                } else {
                    edges_num += counter * (*vecInitVerticalLoops_Edge_Ptr)[col_num].size();
                }
            }
            //
            sstr << edges_num << std::endl;
            for (std::list<std::pair<int, int> >::iterator listIter = (*coeffMatrix)[row_num].begin();
                 listIter != (*coeffMatrix)[row_num].end(); listIter++) {
                int col_num = listIter->first;
                int counter = abs(listIter->second);
                if (col_num >= tot_genus) {// it is a horizontal loop
                    for (int j = 0; j < counter; j++) {
                        for (std::vector<int>::iterator vIter = (*vecInitHorizontalLoops_Edge_Ptr)[col_num -
                                                                                                   tot_genus].begin();
                             vIter != (*vecInitHorizontalLoops_Edge_Ptr)[col_num - tot_genus].end(); vIter++) {
                            sstr << *vIter << " ";
                        }
                    }
                } else {
                    for (int j = 0; j < counter; j++) {
                        for (std::vector<int>::iterator vIter = (*vecInitVerticalLoops_Edge_Ptr)[col_num].begin();
                             vIter != (*vecInitVerticalLoops_Edge_Ptr)[col_num].end(); vIter++) {
                            sstr << *vIter << " ";
                        }
                    }
                }
            }
            sstr << std::endl;
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
        std::cout << "can not open " << file_name << std::endl;
        exit(0);
    }
    //
    return;
}

void FilesOutputForOptimalCycles::WriteAllCyclesOut(std::vector<bool> &inVertLoopType, const char *out_file_name,
                                                    const int loop_type_num) {
    std::vector<std::vector<int> > vertLoops_Edge;
    std::vector<std::vector<int> > horiLoops_Edge;
    std::vector<int> *activeResLoopPtr;
    //

    std::stringstream locSstr(std::stringstream::in | std::stringstream::out);
    //
    string v_file_name = "v_loop_";
    string h_file_name = "h_loop_";
    string file_name;
    // assemble the matrix
    const int tot_genus = vecInitVerticalLoops_Vertex_Ptr->size();
    matrix<double> linkNumMat(2 * tot_genus, 2 * tot_genus);
    matrix<double> invLinkNumMat(2 * tot_genus, 2 * tot_genus);
    //
    //AssembleMatrix(linkNumMat, *coeffMatrix);
    //InvertMatrix(linkNumMat, invLinkNumMat);
    //
    //for (unsigned int i = 0; i < 2 * tot_genus; i++)
    //{		for (unsigned int j = 0; j < 2 * tot_genus; j++)
    //		{
    //			std::cout << invLinkNumMat(i , j ) << " " ;
    //		}
    //		std::cout << std::endl;
    //}
    //std::cout << std::endl;
    //
    // reserve memory for loops
    vertLoops_Edge.resize(tot_genus);
    horiLoops_Edge.resize(tot_genus);
    //
    for (unsigned int i = 0; i < 2 * tot_genus; i++) {
        //
        if (i < tot_genus) {
            locSstr << i << ".txt";
            file_name = v_file_name + locSstr.str();
            locSstr.str("");
            locSstr.clear();

            // set the active loop
            activeResLoopPtr = &vertLoops_Edge[i];
        } else {
            locSstr << i - tot_genus << ".txt";
            file_name = h_file_name + locSstr.str();
            locSstr.str("");
            locSstr.clear();

            // set the active loop
            activeResLoopPtr = &horiLoops_Edge[i - tot_genus];
        }
        //
        std::cout << "writing file ... " << file_name << std::endl;

        std::ofstream ofile;
        ofile.open(file_name.c_str(), std::ifstream::out);
        if (ofile.is_open()) {
            int totEdgeNum = 0;
            for (unsigned int col = 0; col < 2 * tot_genus; col++) {
                if (invLinkNumMat(i, col) < -0.5 || invLinkNumMat(i, col) > 0.5) {
                    //std::cout << invLinkNumMat(i, col)   << " " ;
                    if (col < tot_genus)
                        totEdgeNum += (*vecInitVerticalLoops_Edge_Ptr)[col].size();
                    else
                        totEdgeNum += (*vecInitHorizontalLoops_Edge_Ptr)[col - tot_genus].size();
                }
            }
            //std::cout << std::endl;
            std::stringstream sstr(std::stringstream::in | std::stringstream::out);

            // total number of edges
            sstr << totEdgeNum << " ";
            for (unsigned int col = 0; col < 2 * tot_genus; col++) {
                int coeff = int(float(invLinkNumMat(i, col)));
                std::cout << coeff << " ";
                if (coeff != 0) {// add this into the out
                    int cycle_id = col;
                    std::vector<std::vector<int> > *vecVertex_Ptr = vecInitVerticalLoops_Vertex_Ptr;
                    std::vector<std::vector<int> > *vecEdge_Ptr = vecInitVerticalLoops_Edge_Ptr;


                    if (cycle_id >= vecInitVerticalLoops_Vertex_Ptr->size()) {
                        cycle_id -= int(vecInitVerticalLoops_Vertex_Ptr->size());
                        vecVertex_Ptr = vecInitHorizontalLoops_Vertex_Ptr;
                        vecEdge_Ptr = vecInitHorizontalLoops_Edge_Ptr;
                    }
                    //
                    int vertexNum = int((*vecVertex_Ptr)[cycle_id].size());
                    int edgeNum = int((*vecEdge_Ptr)[cycle_id].size());
                    //

                    for (int eid = 0; eid < edgeNum; eid++) {
                        // vi--ei--v(i+1)
                        // check the orientation
                        if ((inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][eid]].v0 != (*vecVertex_Ptr)[cycle_id][eid] &&
                             inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][eid]].v0 !=
                             (*vecVertex_Ptr)[cycle_id][eid + 1]) ||
                            (inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][eid]].v1 != (*vecVertex_Ptr)[cycle_id][eid] &&
                             inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][eid]].v1 !=
                             (*vecVertex_Ptr)[cycle_id][eid + 1])) {
                            std::cout << "EDGE doesn't match vertices " << std::endl;
                            exit(0);
                        }
                        if (inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][eid]].v0 ==
                            (*vecVertex_Ptr)[cycle_id][eid]) {// vi-->---v(i+1) has the same orientation as ei
                            sstr << 1 * coeff << " ";
                        } else {// opposite orientation
                            sstr << -1 * coeff << " ";
                        }
                        // write the edge
                        sstr << (*vecEdge_Ptr)[cycle_id][eid] << " ";

                        // save the edge to the vector
                        for (unsigned int nCounter = 0; nCounter < abs(coeff); nCounter++)
                            activeResLoopPtr->push_back((*vecEdge_Ptr)[cycle_id][eid]);
                    }

                } //coeff != 0
            } // col
            std::cout << std::endl;
            //
            //ofile.write(sstr.str().c_str(), 10);
            ofile << sstr.rdbuf();
            //
            // delete string stream content
            sstr.str("");
            sstr.clear();
            // close file
            ofile.close();
        }//
        else {
            std::cout << "Can NOT open file " << file_name << std::endl;
            exit(0);
        }

    }// i
    //
    if (loop_type_num == 1)
        WriteAllCyclesInOneFile(vertLoops_Edge, horiLoops_Edge, inVertLoopType, out_file_name, true);
    if (loop_type_num == 2) {
        WriteAllCyclesInOneFile(vertLoops_Edge, horiLoops_Edge, inVertLoopType, "all_h_loops.txt", false);
        WriteAllCyclesInOneFile(vertLoops_Edge, horiLoops_Edge, inVertLoopType, "all_v_loops.txt", true);
    }
    if (loop_type_num == 3)
        WriteAllCyclesInOneFile(vertLoops_Edge, horiLoops_Edge, inVertLoopType, out_file_name, false);
    //
    return;
}

void FilesOutputForOptimalCycles::WriteAllCyclesOut() {
    std::vector<std::vector<int> > vertLoops_Edge;
    std::vector<std::vector<int> > horiLoops_Edge;
    //

    std::stringstream locSstr(std::stringstream::in | std::stringstream::out);
    string v_file_name = "v_loop_";
    string h_file_name = "h_loop_";
    string file_name;
    // assemble the matrix
    const int tot_genus = vecInitVerticalLoops_Vertex_Ptr->size();
    matrix<double> linkNumMat(2 * tot_genus, 2 * tot_genus);
    matrix<double> invLinkNumMat(2 * tot_genus, 2 * tot_genus);
    //
    //AssembleMatrix(linkNumMat, *coeffMatrix);
    //InvertMatrix(linkNumMat, invLinkNumMat);
    ////////////////////////////////////////////
    //
    //for (unsigned int i = 0; i < 2 * tot_genus; i++)
    //{		for (unsigned int j = 0; j < 2 * tot_genus; j++)
    //		{
    //			std::cout << invLinkNumMat(i , j ) << " " ;
    //		}
    //		std::cout << std::endl;
    //}
    //std::cout << std::endl;
    //
    for (unsigned int i = 0; i < 2 * tot_genus; i++) {
        //
        if (i < tot_genus) {
            locSstr << i << ".txt";
            file_name = v_file_name + locSstr.str();
            locSstr.str("");
            locSstr.clear();

        } else {
            locSstr << i - tot_genus << ".txt";
            file_name = h_file_name + locSstr.str();
            locSstr.str("");
            locSstr.clear();
        }
        //
        std::cout << "writing file ... " << file_name << std::endl;

        std::ofstream ofile;
        ofile.open(file_name.c_str(), std::ifstream::out);
        if (ofile.is_open()) {
            int totEdgeNum = 0;
            for (unsigned int col = 0; col < 2 * tot_genus; col++) {
                if (invLinkNumMat(i, col) < -0.5 || invLinkNumMat(i, col) > 0.5) {
                    //std::cout << invLinkNumMat(i, col)   << " " ;
                    if (col < tot_genus)
                        totEdgeNum += (*vecInitVerticalLoops_Edge_Ptr)[col].size();
                    else
                        totEdgeNum += (*vecInitHorizontalLoops_Edge_Ptr)[col - tot_genus].size();
                }
            }
            //std::cout << std::endl;
            std::stringstream sstr(std::stringstream::in | std::stringstream::out);

            // total number of edges
            sstr << totEdgeNum << " ";
            for (unsigned int col = 0; col < 2 * tot_genus; col++) {
                int coeff = int(float(invLinkNumMat(i, col)));
                //std::cout << coeff << " " ;
                if (coeff != 0) {// add this into the out
                    int cycle_id = col;
                    std::vector<std::vector<int> > *vecVertex_Ptr = vecInitVerticalLoops_Vertex_Ptr;
                    std::vector<std::vector<int> > *vecEdge_Ptr = vecInitVerticalLoops_Edge_Ptr;


                    if (cycle_id >= vecInitVerticalLoops_Vertex_Ptr->size()) {
                        cycle_id -= int(vecInitVerticalLoops_Vertex_Ptr->size());
                        vecVertex_Ptr = vecInitHorizontalLoops_Vertex_Ptr;
                        vecEdge_Ptr = vecInitHorizontalLoops_Edge_Ptr;
                    }
                    //
                    int vertexNum = int((*vecVertex_Ptr)[cycle_id].size());
                    int edgeNum = int((*vecEdge_Ptr)[cycle_id].size());
                    //

                    for (int eid = 0; eid < edgeNum; eid++) {
                        // vi--ei--v(i+1)
                        // check the orientation
                        if ((inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][eid]].v0 != (*vecVertex_Ptr)[cycle_id][eid] &&
                             inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][eid]].v0 !=
                             (*vecVertex_Ptr)[cycle_id][eid + 1]) ||
                            (inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][eid]].v1 != (*vecVertex_Ptr)[cycle_id][eid] &&
                             inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][eid]].v1 !=
                             (*vecVertex_Ptr)[cycle_id][eid + 1])) {
                            std::cout << "EDGE doesn't match vertices " << std::endl;
                            exit(0);
                        }
                        if (inMeshPtr->vecEdge[(*vecEdge_Ptr)[cycle_id][eid]].v0 ==
                            (*vecVertex_Ptr)[cycle_id][eid]) {// vi-->---v(i+1) has the same orientation as ei
                            sstr << 1 * coeff << " ";
                        } else {// opposite orientation
                            sstr << -1 * coeff << " ";
                        }
                        // write the edge
                        sstr << (*vecEdge_Ptr)[cycle_id][eid] << " ";
                    }

                } //coeff != 0
            } // col
            //std::cout << std::endl;
            //
            //ofile.write(sstr.str().c_str(), 10);
            ofile << sstr.rdbuf();
            //
            // delete string stream content
            sstr.str("");
            sstr.clear();
            // close file
            ofile.close();
        }//
        else {
            std::cout << "Can NOT open file " << file_name << std::endl;
            exit(0);
        }

    }// i
    return;
}

void FilesOutputForOptimalCycles::WriteAllCyclesInOneFile(std::vector<std::vector<int> > &vertLoops,
                                                          std::vector<std::vector<int> > &horiLoops,
                                                          const char *out_file_name) {// loop_type ---- true for vertical loop
    // loop_type ---- false for horizontal loop
    std::cout << "writing file ... " << out_file_name << std::endl;

    std::ofstream ofile;
    ofile.open(out_file_name, std::ifstream::out);

    //
    int tot_num = vertLoops.size();
    //
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    if (ofile.is_open()) {
        // total number of edges
        sstr << "h " << tot_num << std::endl;;
        for (int i = 0; i < tot_num; i++) {
            sstr << horiLoops[i].size() << std::endl;
            for (std::vector<int>::iterator sIter = horiLoops[i].begin();
                 sIter != horiLoops[i].end();
                 sIter++) {
                sstr << *sIter << " ";
            }
            sstr << std::endl;
        }
        sstr << std::endl << "v " << tot_num << std::endl;
        for (int i = 0; i < tot_num; i++) {
            sstr << vertLoops[i].size() << std::endl;
            for (std::vector<int>::iterator sIter = vertLoops[i].begin();
                 sIter != vertLoops[i].end();
                 sIter++) {
                sstr << *sIter << " ";
            }
            sstr << std::endl;
        }
        //
        //ofile.write(sstr.str().c_str(), 10);
        ofile << sstr.rdbuf();
        //
        // delete string stream content
        sstr.str("");
        sstr.clear();
        // close file
        ofile.close();
    } else {
        std::cout << "Can NOT open file " << out_file_name << std::endl;
        exit(0);
    }
    std::cout << "Done... " << tot_num << std::endl;
    return;

}

void FilesOutputForOptimalCycles::WriteAllCyclesInOneFile(std::vector<std::vector<int> > &vertLoops,
                                                          std::vector<std::vector<int> > &horiLoops,
                                                          std::vector<bool> &initVerticalLoopsType,
                                                          const char *out_file_name,
                                                          const bool loop_type) {// loop_type ---- true for vertical loop
    // loop_type ---- false for horizontal loop
    std::cout << "writing file ... " << out_file_name << std::endl;

    std::ofstream ofile;
    ofile.open(out_file_name, std::ifstream::out);

    //
    int tot_num = initVerticalLoopsType.size();
    int edgeNum = 0;
    //
    std::vector<std::vector<int> > *posLoops;
    std::vector<std::vector<int> > *negLoops;
    if (loop_type) {
        posLoops = &vertLoops;
        negLoops = &horiLoops;
    } else {
        posLoops = &horiLoops;
        negLoops = &vertLoops;
    }
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    if (ofile.is_open()) {
        // total number of edges
        sstr << tot_num << std::endl;;
        for (int i = 0; i < tot_num; i++) {
            // for each loop
            std::vector<int> *curEdgeLoopPtr = NULL;
            if (loop_type) {
                if (initVerticalLoopsType[i] == true) { //vertical
                    curEdgeLoopPtr = &((vertLoops)[i]);
                } else {// horizontal
                    curEdgeLoopPtr = &((horiLoops)[i]);
                }
            } else {
                if (initVerticalLoopsType[i] == true) { //vertical
                    curEdgeLoopPtr = &((horiLoops)[i]);
                } else {// horizontal
                    curEdgeLoopPtr = &((vertLoops)[i]);
                }
            }
            sstr << curEdgeLoopPtr->size() << std::endl;
            for (std::vector<int>::iterator sIter = curEdgeLoopPtr->begin();
                 sIter != curEdgeLoopPtr->end();
                 sIter++) {
                sstr << *sIter << " ";
            }
            sstr << std::endl;
        }

        //
        //ofile.write(sstr.str().c_str(), 10);
        ofile << sstr.rdbuf();
        //
        // delete string stream content
        sstr.str("");
        sstr.clear();
        // close file
        ofile.close();
    } else {
        std::cout << "Can NOT open file " << out_file_name << std::endl;
        exit(0);
    }
    std::cout << "Done... " << tot_num << std::endl;
    return;

}

void FilesOutputForOptimalCycles::WriteCriticalPoints(const char *out_file_name,
                                                      std::vector<std::pair<int, int> > &critSet) {
    std::cout << "writing file ... " << out_file_name << std::endl;

    std::ofstream ofile;
    ofile.open(out_file_name, std::ifstream::out);
    //
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    if (ofile.is_open()) {
        sstr << critSet.size() << std::endl;

        for (unsigned int i = 0; i < critSet.size(); i++) {
            sstr << critSet[i].first << " " << critSet[i].second << std::endl;
        }
        //ofile.write(sstr.str().c_str(), 10);
        ofile << sstr.rdbuf();
        //
        // delete string stream content
        sstr.str("");
        sstr.clear();
        // close file
        ofile.close();
    } else {
        std::cout << "Can NOT open file " << out_file_name << std::endl;
        exit(0);
    }
    std::cout << "Done... " << critSet.size() << std::endl;
    return;
}

void FilesOutputForOptimalCycles::WriteCriticalPoints(const char *out_file_name, std::set<int> &critSet,
                                                      const bool bVertical) {
    std::cout << "writing file ... " << out_file_name << std::endl;

    std::ofstream ofile;
    ofile.open(out_file_name, std::ifstream::out);

    std::vector<std::vector<int> > *activeLoops;
    if (bVertical)
        activeLoops = vecInitHorizontalLoops_Vertex_Ptr;
    else
        activeLoops = vecInitVerticalLoops_Vertex_Ptr;
    //
//	int tot_num = critSet.size();
    int tot_num = 0;
    for (unsigned int i = 0; i < vecInitHorizontalLoops_Edge_Ptr->size(); i++) {
        tot_num += (*activeLoops)[i].size();
    }
    //
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
    //
    if (ofile.is_open()) {
        // total number of edges
        //sstr << tot_num << std::endl;
        //for (std::set<int>::iterator sIter = critSet.begin();
        //	sIter != critSet.end();
        //	sIter++)
        //{
        //	sstr << *sIter << " " ;
        //}
        //sstr << std::endl;

        sstr << tot_num << std::endl;
        std::set<int> uniqueHoriVertex;
        for (unsigned int i = 0; i < activeLoops->size(); i++) {
            for (unsigned int k = 0; k < (*activeLoops)[i].size(); k++) {
                uniqueHoriVertex.insert((*activeLoops)[i][k]);
            }
        }
        for (std::set<int>::iterator sIter = uniqueHoriVertex.begin();
             sIter != uniqueHoriVertex.end();
             sIter++) {
            sstr << *sIter << " ";
        }
        sstr << std::endl;
        //
        //ofile.write(sstr.str().c_str(), 10);
        ofile << sstr.rdbuf();
        //
        // delete string stream content
        sstr.str("");
        sstr.clear();
        // close file
        ofile.close();
    } else {
        std::cout << "Can NOT open file " << out_file_name << std::endl;
        exit(0);
    }
    std::cout << "Done... " << tot_num << std::endl;
}
