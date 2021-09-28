/*
(c) 2012 Fengtao Fan
*/
#include "psbmReebGraph.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include "SimpleMesh.h"
#include "FilesOutputForOptimalCycles.h"


#include <time.h>
#include <sstream>
#include <boost/progress.hpp>
#include <boost/program_options.hpp>

//
bool CheckBoundaries(_SimpleMesh &inMesh, std::vector<std::vector<std::pair<int, int> > > &boundries) {
    bool ret = false;
    std::vector<bool> edge_flag(inMesh.vecEdge.size(), false);
    for (unsigned int i = 0; i < inMesh.vecEdge.size(); i++) {
        if (!edge_flag[i] && inMesh.vecEdge[i].AdjTriNum == 1) {// this is a boundary edge
            ret = true;
            std::vector<std::pair<int, int> > bdrEdges;
            //
            int terminate_vertex = inMesh.vecEdge[i].v0;
            int current_vertex = inMesh.vecEdge[i].v1;
            int current_edge = i;
            int loop_counter = 0;
            while (current_vertex != terminate_vertex) {
                bdrEdges.push_back(std::pair<int, int>(current_edge, inMesh.vecEdge[current_edge].AdjTri[0]));
                edge_flag[current_edge] = true;
                //
                unsigned int evid = 0;
                for (; evid < inMesh.vecVertex[current_vertex].adjEdges.size(); evid++) {
                    int rot_edge_id = inMesh.vecVertex[current_vertex].adjEdges[evid];
                    if (rot_edge_id != current_edge && inMesh.vecEdge[rot_edge_id].AdjTriNum == 1)
                        break;
                }
                if (evid < inMesh.vecVertex[current_vertex].adjEdges.size()) {
                    current_edge = inMesh.vecVertex[current_vertex].adjEdges[evid];
                    current_vertex = inMesh.vecEdge[current_edge].v0 + inMesh.vecEdge[current_edge].v1 - current_vertex;
                } else {
                    std::cout << "READING MESH MODEL ERROR : " << std::endl;
                    std::cout << "--- CAN NOT FIND ADJACENT BOUDNARY EDGE !" << std::endl;
                    std::cout << "--- PLEASE CHECK THE MESH MODEL" << std::endl;
                    exit(0);
                }
                loop_counter++;
                if (loop_counter > inMesh.vecEdge.size()) {
                    std::cout << "READING MESH MODEL ERROR : " << std::endl;
                    std::cout << "--- PLEASE CHECK THE MESH MODEL" << std::endl;
                    exit(0);
                }
            }
            edge_flag[current_edge] = true;
            bdrEdges.push_back(std::pair<int, int>(current_edge, inMesh.vecEdge[current_edge].AdjTri[0]));
            //
            boundries.push_back(bdrEdges);
        }
        edge_flag[i] = true;
    }
    return ret;
}

void FindEdgeOrientation(_SimpleMesh &inMesh, const int fid, std::vector<int> &OrientTriangles,
                         std::pair<int, int> &EdgeOrientation) {
    int start_index = 3 * fid;
    std::vector<int> TriOrientation(OrientTriangles.begin() + start_index, OrientTriangles.begin() + start_index + 3);
    TriOrientation.push_back(TriOrientation.front());
    //
    for (unsigned int i = 0; i < 3; i++) {
        if ((EdgeOrientation.first == TriOrientation[i] && EdgeOrientation.second == TriOrientation[i + 1]) ||
            (EdgeOrientation.second == TriOrientation[i] && EdgeOrientation.first == TriOrientation[i + 1])) {
            EdgeOrientation.first = TriOrientation[i + 1];
            EdgeOrientation.second = TriOrientation[i];
        }
    }
    return;
}

void CloseHoles(_SimpleMesh &inMesh, std::vector<std::vector<std::pair<int, int> > > &meshBoundaries,
                std::vector<Vector3> &inMeshNormal,
                std::vector<int> &OrientTriangles, std::set<int> &extraVertices) {
    for (unsigned int b = 0; b < meshBoundaries.size(); b++) {
        Vector3 centroid(0.0, 0.0, 0.0);
        int curEdge = 0;
        int nextEdge = 0;
        int curVertex = 0;
        for (unsigned int i = 0; i < meshBoundaries[b].size(); i++) {// the edges on the boundary is ordered
            curEdge = meshBoundaries[b][i].first;
            nextEdge = meshBoundaries[b][0].first;
            if (i < meshBoundaries[b].size() - 1)
                nextEdge = meshBoundaries[b][i + 1].first;
            // shared vertex
            curVertex = inMesh.vecEdge[curEdge].v0;
            if (curVertex != inMesh.vecEdge[nextEdge].v0 &&
                curVertex != inMesh.vecEdge[nextEdge].v1)
                curVertex = inMesh.vecEdge[curEdge].v1;
            //
            centroid = centroid + Vector3(inMesh.vecVertex[curVertex].x,
                                          inMesh.vecVertex[curVertex].y,
                                          inMesh.vecVertex[curVertex].z);
            //
        }
        centroid = centroid / meshBoundaries[b].size();
        //
        // add vertex
        int centroid_index = inMesh.vecVertex.size();
        extraVertices.insert(centroid_index);
        _SimpleMeshVertex tmpVer;
        tmpVer.x = centroid[0];
        tmpVer.y = centroid[1];
        tmpVer.z = centroid[2];
        //
        inMesh.vecVertex.push_back(tmpVer);
        //
        // add edges and triangles
        std::map<int, int> edge_mapping;
        for (unsigned int i = 0; i < meshBoundaries[b].size(); i++) {
            // Find orientation of this edge
            curEdge = meshBoundaries[b][i].first;
            std::pair<int, int> EdgeOrientation(inMesh.vecEdge[curEdge].v0, inMesh.vecEdge[curEdge].v1);
            FindEdgeOrientation(inMesh, meshBoundaries[b][i].second, OrientTriangles, EdgeOrientation);
            //
            OrientTriangles.push_back(EdgeOrientation.first);
            OrientTriangles.push_back(EdgeOrientation.second);
            OrientTriangles.push_back(centroid_index);
            //
            _SimpleMeshTriangle tmpTri;
            _SimpleMeshEdge tmpEdge;
            //
            tmpTri.v0 = EdgeOrientation.first;
            tmpTri.v1 = EdgeOrientation.second;
            tmpTri.v2 = centroid_index;
            //
            //
            Vector3 leftVec, rightVec;
            leftVec[0] = inMesh.vecVertex[tmpTri.v2].x - inMesh.vecVertex[tmpTri.v1].x;
            leftVec[1] = inMesh.vecVertex[tmpTri.v2].y - inMesh.vecVertex[tmpTri.v1].y;
            leftVec[2] = inMesh.vecVertex[tmpTri.v2].z - inMesh.vecVertex[tmpTri.v1].z;

            rightVec[0] = inMesh.vecVertex[tmpTri.v0].x - inMesh.vecVertex[tmpTri.v1].x;
            rightVec[1] = inMesh.vecVertex[tmpTri.v0].y - inMesh.vecVertex[tmpTri.v1].y;
            rightVec[2] = inMesh.vecVertex[tmpTri.v0].z - inMesh.vecVertex[tmpTri.v1].z;
            leftVec = leftVec ^ rightVec;
            unitize(leftVec);
            //leftVec = leftVec / norm(leftVec);
            inMeshNormal.push_back(leftVec);
            //
            tmpTri.sortVertices();
            // [v0 v1] existing edge
            inMesh.vecEdge[curEdge].AdjTriNum++;
            inMesh.vecEdge[curEdge].AdjTri[1] = inMesh.vecTriangle.size();
            //
            tmpTri.e01 = curEdge;
            //
            std::map<int, int>::iterator mIter;
            mIter = edge_mapping.find(tmpTri.v1);
            if (mIter == edge_mapping.end()) {// new edge
                tmpEdge.v0 = tmpTri.v1;
                tmpEdge.v1 = tmpTri.v2;
                tmpEdge.AdjTri[0] = inMesh.vecTriangle.size();
                tmpEdge.AdjTriNum = 1;
                //
                inMesh.vecEdge.push_back(tmpEdge);
                tmpTri.e12 = inMesh.vecEdge.size() - 1;
                edge_mapping[tmpTri.v1] = tmpTri.e12;
                //
                inMesh.vecVertex[tmpEdge.v0].adjEdges.push_back(tmpTri.e12);
                inMesh.vecVertex[tmpEdge.v1].adjEdges.push_back(tmpTri.e12);
            } else {// existed edge
                inMesh.vecEdge[mIter->second].AdjTriNum++;
                inMesh.vecEdge[mIter->second].AdjTri[1] = inMesh.vecTriangle.size();
                //
                tmpTri.e12 = mIter->second;
            }
            mIter = edge_mapping.find(tmpTri.v0);
            if (mIter == edge_mapping.end()) {
                // new edge
                tmpEdge.v0 = tmpTri.v0;
                tmpEdge.v1 = tmpTri.v2;
                tmpEdge.AdjTri[0] = inMesh.vecTriangle.size();
                tmpEdge.AdjTriNum = 1;
                //
                inMesh.vecEdge.push_back(tmpEdge);
                tmpTri.e02 = inMesh.vecEdge.size() - 1;
                edge_mapping[tmpTri.v0] = tmpTri.e02;
                //
                inMesh.vecVertex[tmpEdge.v0].adjEdges.push_back(tmpTri.e02);
                inMesh.vecVertex[tmpEdge.v1].adjEdges.push_back(tmpTri.e02);
            } else {// existed already
                inMesh.vecEdge[mIter->second].AdjTriNum++;
                inMesh.vecEdge[mIter->second].AdjTri[1] = inMesh.vecTriangle.size();
                //
                tmpTri.e02 = mIter->second;
            }
            //
            inMesh.vecTriangle.push_back(tmpTri);
        }
    }
    return;
}

int LoadData(_SimpleMesh &mesh, std::vector<Vector3> &meshNormal,
             const char *mesh_file_name, double &BoundingBoxRadius, std::vector<int> &tris,
             std::set<int> &extraVertices, int &orgTriangleSize,
             const float fEnlargeFactor) {
    int genus = 0;
    _SimpleMeshVertex minBd;
    _SimpleMeshVertex maxBd;
    mesh.LoadMeshInOFFformat(minBd, maxBd, meshNormal, mesh_file_name, tris,
                             fEnlargeFactor);//"E:\\RProgramming\\Models\\torus.off");// "D:\\MeshModels\\OFF-models\\HighGenusCubeHC1.off");//
    mesh.SetMeshNormalPtr(&meshNormal);
    // triangles in TRIS are in the same order as triangles in mesh.vecTriangle;
    //
    orgTriangleSize = mesh.vecTriangle.size();
    //
    std::vector<std::vector<std::pair<int, int> > > meshBoundaries;
    if (CheckBoundaries(mesh, meshBoundaries)) {// it is a mesh with bondary
        int EulerCharacteristic = mesh.vecVertex.size() + mesh.vecTriangle.size() - mesh.vecEdge.size();
        genus = 1 - (EulerCharacteristic + meshBoundaries.size()) / 2;
        if (genus) {
            CloseHoles(mesh, meshBoundaries, meshNormal, tris, extraVertices);
        }
    } else {// it is a closed mesh
        int EulerCharacteristic = mesh.vecVertex.size() + mesh.vecTriangle.size() - mesh.vecEdge.size();
        genus = 1 - EulerCharacteristic / 2;
    }
    if (!genus) {
        std::cout << "NOTHING IS COMPUTED : " << std::endl;
        std::cout << " ---- MESH HAS GENUS 0!" << std::endl;
        exit(1);
    }
    std::cout << "Mesh has genus : " << genus << std::endl;
    //

// compute the bounding box

    return genus;
}

void RandomUniqueDirection(_SimpleMesh &mesh, Vector3 &uniDirection) {
    //
    std::vector<Vector3> vecEdgeDirections(2 * mesh.vecEdge.size());
    // all vectors are pointed to positive Z direction
    for (unsigned int i = 0; i < mesh.vecEdge.size(); i++) {
        vecEdgeDirections[i][0] = mesh.vecVertex[mesh.vecEdge[i].v0].x - mesh.vecVertex[mesh.vecEdge[i].v1].x;
        vecEdgeDirections[i][1] = mesh.vecVertex[mesh.vecEdge[i].v0].y - mesh.vecVertex[mesh.vecEdge[i].v1].y;
        vecEdgeDirections[i][2] = mesh.vecVertex[mesh.vecEdge[i].v0].z - mesh.vecVertex[mesh.vecEdge[i].v1].z;
        //
        if (vecEdgeDirections[i][2] < 0) {
            for (int j = 0; j < 3; j++)
                vecEdgeDirections[i][j] = -vecEdgeDirections[i][j];
        }
        //
        unitize(vecEdgeDirections[i]);
    }
    int mesh_edge_size = mesh.vecEdge.size();
    int opp_v_0;
    int opp_v_1;
    for (unsigned int i = 0; i < mesh.vecEdge.size(); i++) {
        if (mesh.vecEdge[i].AdjTriNum == 2) {
            opp_v_0 = mesh.vecEdge[i].AdjTri[0];
            if (mesh.vecTriangle[opp_v_0].v0 != mesh.vecEdge[i].v0 &&
                mesh.vecTriangle[opp_v_0].v0 != mesh.vecEdge[i].v1) {
                opp_v_0 = mesh.vecTriangle[opp_v_0].v0;
            } else {
                if (mesh.vecTriangle[opp_v_0].v1 != mesh.vecEdge[i].v0 &&
                    mesh.vecTriangle[opp_v_0].v1 != mesh.vecEdge[i].v1) {
                    opp_v_0 = mesh.vecTriangle[opp_v_0].v1;
                } else {
                    opp_v_0 = mesh.vecTriangle[opp_v_0].v2;
                }
            }
            //
            opp_v_1 = mesh.vecEdge[i].AdjTri[1];
            if (mesh.vecTriangle[opp_v_1].v0 != mesh.vecEdge[i].v0 &&
                mesh.vecTriangle[opp_v_1].v0 != mesh.vecEdge[i].v1) {
                opp_v_1 = mesh.vecTriangle[opp_v_1].v0;
            } else {
                if (mesh.vecTriangle[opp_v_1].v1 != mesh.vecEdge[i].v0 &&
                    mesh.vecTriangle[opp_v_1].v1 != mesh.vecEdge[i].v1) {
                    opp_v_1 = mesh.vecTriangle[opp_v_1].v1;
                } else {
                    opp_v_1 = mesh.vecTriangle[opp_v_1].v2;
                }
            }
            //
            vecEdgeDirections[i + mesh_edge_size][0] = mesh.vecVertex[opp_v_0].x - mesh.vecVertex[opp_v_1].x;
            vecEdgeDirections[i + mesh_edge_size][1] = mesh.vecVertex[opp_v_0].y - mesh.vecVertex[opp_v_1].y;
            vecEdgeDirections[i + mesh_edge_size][2] = mesh.vecVertex[opp_v_0].z - mesh.vecVertex[opp_v_1].z;
            //
            if (vecEdgeDirections[i + mesh_edge_size][2] < 0) {
                for (int j = 0; j < 3; j++)
                    vecEdgeDirections[i + mesh_edge_size][j] = -vecEdgeDirections[i + mesh_edge_size][j];
            }
            //
            unitize(vecEdgeDirections[i + mesh_edge_size]);
        }
    }
    //
    srand(time(NULL));
    //
    int runTimesCounter = 0;
    int halfRandMax = RAND_MAX >> 2;
    bool bFindDirection = false;
    double minAangleValue = 1.0;
    double error_precision = 1e-7;
    double max_running_counter = 2000;
    while (!bFindDirection) {
        runTimesCounter++;
        for (int i = 0; i < 3; i++) {
            uniDirection[i] = rand() - halfRandMax;
        }
        if (uniDirection[2] < 0) {
            for (int i = 0; i < 3; i++)
                uniDirection[i] = -uniDirection[i];
        }
        //
        unitize(uniDirection);
        //
        bFindDirection = true;
        minAangleValue = 1.0;
        for (unsigned int i = 0; i < vecEdgeDirections.size(); i++) {
            double angleValue = uniDirection * vecEdgeDirections[i];
            if (angleValue < error_precision && angleValue > -error_precision) {
                bFindDirection = false;
                break;
            }
            if (std::abs(minAangleValue) > std::abs(angleValue))
                minAangleValue = angleValue;
        }
        if (runTimesCounter > max_running_counter)
            break;
    }
    //
    //std::cout << "run times : " << runTimesCounter << std::endl;
    //std::cout << "min angle : " << minAangleValue << std::endl;
    //
    return;
}

//
void CycleLocalOptimization(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles,
                            std::vector<std::set<int> > &v_basis,
                            std::vector<std::set<int> > &h_basis, const float fEnlargeFactor);

void CycleLocalOptimization_bdry(_SimpleMesh &locMesh, psbmReebGraph &reebgraph, std::vector<int> &OrientTriangles,
                                 std::vector<std::set<int> > &out_v_basis_loops,
                                 std::vector<std::set<int> > &out_h_basis_loops,
                                 std::set<int> extraVertices, const float fEnlargeFactor);

////
char *strLicense = "THIS SOFTWARE IS PROVIDED \"AS-IS\". THERE IS NO WARRANTY OF ANY KIND. "
                   "NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE FOR "
                   "ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY. \n"
                   "\n"
                   "This software was developed (and is copyrighted by) the Jyamiti group at "
                   "The Ohio State University. Please do not redistribute this software. "
                   "This program is for academic research use only. This software uses the "
                   "CGAL library (www.cgal.org), Boost library (www.boost.org) and Ann library "
                   "(www.cs.umd.edu/~mount/ANN/) which are covered under their own licenses.\n"
                   "\n"
                   "The CGAL library's license "
                   "(which applies to the CGAL library ONLY and NOT to this program itself) is "
                   "as follows:\n"
                   "\n"
                   "LICENSE\n"
                   "---------------------------------------------------------------------------\n"
                   "\n"
                   "The CGAL software consists of several parts, each of which is licensed under "
                   "an open source license. It is also possible to obtain commercial licenses "
                   "from GeometryFactory (www.geometryfactory.com) for all or parts of CGAL. \n"
                   "\n"
                   "The source code of the CGAL library can be found in the directories "
                   "\"src/CGAL\", \"src/CGALQt\" and \"include/CGAL\". It is specified in each file of "
                   "the CGAL library which license applies to it. This is either the GNU Lesser "
                   "General Public License (as published by the Free Software Foundation; "
                   "version 2.1 of the License) or the Q Public License (version 1.0). The texts "
                   "of both licenses can be found in the files LICENSE.LGPL and LICENSE.QPL. \n"
                   "\n"
                   "Distributed along with CGAL (for the users' convenience), but not part of "
                   "CGAL, are the following third-party libraries, available under their own "
                   "licenses: \n"
                   "\n"
                   "- CORE, in the directories \"include/CORE\" and \"src/Core\", is licensed under"
                   " the QPL (see LICENSE.QPL). \n"
                   "- OpenNL, in the directory \"include/OpenNL\", is licensed under the LGPL"
                   " (see include/OpenNL/LICENSE.OPENNL). \n"
                   "- ImageIO, in the directory \"examples/Surface_mesher/ImageIO\", is licensed"
                   " under the LGPL (see LICENSE.LGPL). \n"
                   "\n"
                   "All other files that do not have an explicit copyright notice (e.g., all "
                   "examples and some demos) are licensed under a very permissive license. The "
                   "exact license text can be found in the file LICENSE.FREE_USE. Note that some "
                   "subdirectories have their own copy of LICENSE.FREE_USE. These copies have "
                   "the same license text and differ only in the copyright holder.\n"
                   "---------------------------------------------------------------------------\n"
                   "\n"
                   "The Boost library's license "
                   "(which applies to the Boost library ONLY and NOT to this program itself) is "
                   "as follows:\n"
                   "\n"
                   "LICENSE\n"
                   "---------------------------------------------------------------------------\n"
                   "Boost Software License - Version 1.0 - August 17th, 2003\n"
                   "\n"
                   "Permission is hereby granted, free of charge, to any person or organization "
                   "obtaining a copy of the software and accompanying documentation covered by "
                   "this license (the \"Software\") to use, reproduce, display, distribute, "
                   "execute, and transmit the Software, and to prepare derivative works of the "
                   "Software, and to permit third-parties to whom the Software is furnished to "
                   "do so, all subject to the following: \n"
                   "\n"
                   "The copyright notices in the Software and this entire statement, including "
                   "the above license grant, this restriction and the following disclaimer, "
                   "must be included in all copies of the Software, in whole or in part, and "
                   "all derivative works of the Software, unless such copies or derivative "
                   "works are solely in the form of machine-executable object code generated by "
                   "a source language processor. \n"
                   "\n"
                   "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR "
                   "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, "
                   "FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT "
                   "SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE "
                   "FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, "
                   "ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER "
                   "DEALINGS IN THE SOFTWARE. \n"
                   "---------------------------------------------------------------------------\n"
                   "\n"
                   "ANN library's license "
                   "(which applies to the ANN library ONLY and NOT to this program itself) is "
                   "as follows: \n"
                   "\n"
                   "LICENSE\n"
                   "---------------------------------------------------------------------------\n"
                   "The ANN Library (all versions) is provided under the terms and "
                   "conditions of the GNU Lesser General Public Library, which is stated "
                   "below.  It can also be found at: \n"
                   "\n"
                   "   http:\//www.gnu.org/copyleft/lesser.html \n"
                   "---------------------------------------------------------------------------\n";

//char * strLicense = "THIS SOFTWARE IS PROVIDED \"AS-IS\". http:\//www.gnu.org/copyleft/lesser.html THERE IS NO WARRANTY OF ANY KIND. \n";
bool ParseCommand(int argc, char **argv, std::string &InputFile, std::string &OutputFile) {
    try {
        /* Define the program options description
        */
        namespace po = boost::program_options;
        po::options_description desc("ReebHanTun Usage");
        desc.add_options()
                (",h", "Help information")
                (",l", "License information")
                (",I", po::value<std::string>(&InputFile)->required(), "Input file name")
                (",O", po::value<std::string>(&OutputFile)->required(), "Output file name prefix");

        // Parser map
        po::variables_map vm;
        try {
            po::store(po::parse_command_line(argc, argv, desc), vm);

            //
            if (vm.count("-h")) {
                std::cout << desc << std::endl;
            }
            //
            if (vm.count("-l")) {
                std::cout << strLicense << std::endl;
            }
            //
            po::notify(vm);
        }
        catch (boost::program_options::required_option &e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            return false;
        }
        catch (boost::program_options::error &e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            return false;
        }
        //if (vm.count("-I"))
        //{
        //	std::cout << vm["-I"].as<std::string>() << std::endl;
        //}
        //if (vm.count("-O"))
        //{
        //	std::cout << vm["-O"].as<std::string>() << std::endl;
        //}
    }
    catch (std::exception &e) {
        std::cerr << "Unhandled Exception reached the top of main: "
                  << e.what() << ", application will now exit" << std::endl;
        return false;

    }
    return true;
}

int main(int argc, char **argv) {
    std::string InputFileName;
    std::string OutputFileName;
    Vector3 distinctDirection;
    //
    _SimpleMesh mesh;
    const float fEnlargeFactor = 10000.f;
    std::vector<std::set<int> > v_basis_loops;
    std::vector<std::set<int> > h_basis_loops;
    //
    psbmReebGraph reebGraph;
    std::vector<Vector3> meshNormal;
    //
    if (ParseCommand(argc, argv, InputFileName, OutputFileName)) {
        //
        std::vector<int> OrientTriangles;
        std::set<int> extraVertices;
        int nOrgTriangleSize = 0;
        double BoundingBoxRadius;
        LoadData(mesh, meshNormal, InputFileName.c_str(), BoundingBoxRadius, OrientTriangles, extraVertices,
                 nOrgTriangleSize, fEnlargeFactor);
        //////

        RandomUniqueDirection(mesh, distinctDirection);
        //
        //std::cout << distinctDirection[0] << " " <<
        //			 distinctDirection[1] << " " <<
        //			 distinctDirection[2] << std::endl;

        //

        reebGraph.ReserveSpaceForEdges(mesh.vecEdge.size());
        //
        double *scalarField = new double[mesh.vecVertex.size()];
        int perDir = 0;
        int rayDir = 0;
        for (unsigned int i = 0; i < mesh.vecVertex.size(); i++) {
            scalarField[i] = mesh.vecVertex[i].x * distinctDirection[0] +
                             mesh.vecVertex[i].y * distinctDirection[1] +
                             mesh.vecVertex[i].z * distinctDirection[2];
        }
        reebGraph.SetHeightDirection(distinctDirection);
        reebGraph.AssignData(&mesh, scalarField);
        reebGraph.scalarDir = 1;// x=0, y=1, z=2
        //
        //minimumDiffBetweenVertices(scalarField);
        //
        std::cout << std::endl;
        {
            boost::progress_timer t;
            //
            reebGraph.ComputeReebGraph();
            //
        }
        {
        	reebGraph.WriteReebGraphOBJ("casting.obj", mesh.vecVertex);
        	std::cout << "Vertices in RG : " << reebGraph.pVecReebNode->size() << std::endl;
        	std::cout << "Edges in RG : " << reebGraph.pListReebArc->size() << std::endl;
        	std::cout << "Genus is : " << reebGraph.pListReebArc->size() - reebGraph.pVecReebNode->size() + 1  << std::endl;
        }
        {
            boost::progress_timer t;
            ////
            std::cout << "Time for mapping and linking :" << std::endl;
            //reebGraph.ComputeCycleAndPairing();
            reebGraph.ComputingCycle_max_tree();

//            reebGraph.WriteSimplifiedReebGraphOBJ("casting-sim.obj", mesh.vecVertex);

            //std::cout << "mapping" << std::endl;
            //// computing the cycle on surface
            reebGraph.compute_path_on_mesh_for_each_simplified_arc();//pathArcOnMesh, offsetPathArcOnMesh);

            //std::cout << "embed" << std::endl;
            reebGraph.EmbedCycleAsEdgePathOnMesh();
            //std::cout << "linking" << std::endl;
            reebGraph.LinkNumberMatrixComputing();

        }
        //

        //
        if (extraVertices.empty())
            CycleLocalOptimization(mesh, reebGraph, OrientTriangles, v_basis_loops, h_basis_loops,
                                   1.f / fEnlargeFactor);
        else
            CycleLocalOptimization_bdry(mesh, reebGraph, OrientTriangles, v_basis_loops, h_basis_loops, extraVertices,
                                        1.f / fEnlargeFactor);
        //
        std::cout << "Handle and tunnel loops written in files :  \n";
        FilesOutputForOptimalCycles files_out_op;
        files_out_op.InitMeshPtr(&mesh);
        files_out_op.WriteCyclesInformation(OutputFileName.c_str(), v_basis_loops, h_basis_loops);
        int orgVertexSize = mesh.vecVertex.size();
        if (!extraVertices.empty())
            orgVertexSize = *extraVertices.begin();
        files_out_op.WriteGeomviewListFormat(OutputFileName.c_str(), v_basis_loops, h_basis_loops, OrientTriangles,
                                             orgVertexSize, nOrgTriangleSize, 1.f / fEnlargeFactor);

        ///////////////////////////////////////////////////
    }
    return 1;
}

