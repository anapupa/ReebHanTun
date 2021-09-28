/*
(c) 2012 Fengtao Fan
*/
#include "edge_annotations_gauss.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

bool edge_annotation_computing::IndependenceCheck(std::vector<Annotation_Type> &bReducedVec,
                                                  std::vector<int> &lowestOnePtr,
                                                  Annotation_Type &objVector) {
    bool dependence = true;
    Annotation_Type testVec(objVector);
    for (int row = vec_size - 1; row >= 0; row--) {
        if (testVec[row]) {
            if (lowestOnePtr[row] >= 0) {
                int vecIndex = lowestOnePtr[row];
                for (int j = row; j >= 0; j--) {
                    testVec[j] = testVec[j] ^ bReducedVec[vecIndex][j];
                }
            } else {
                dependence = false;
                break;
            }
        }
    }
    return dependence;
}

bool edge_annotation_computing::ReduceMatrix(std::vector<Annotation_Type> &annMatrix, std::vector<int> &lowestOnePos) {
    bool dependence = false;
    int nVecSize = (int) annMatrix.front().size();
    //
    lowestOnePos.resize(nVecSize);
    for (unsigned int i = 0; i < lowestOnePos.size(); i++)
        lowestOnePos[i] = -1;
    //
    for (unsigned int row = 0; row < annMatrix.size(); row++) {
        bool bZeroAnnotation = true;
        for (int col = nVecSize - 1; col >= 0; col--) {
            if (annMatrix[row][col]) {// check the bit is unique or not
                if (lowestOnePos[col] >= 0) {// this bit is occupied
                    int curBit = lowestOnePos[col];
                    for (int jcol = col; jcol >= 0; jcol--)
                        annMatrix[row][jcol] = annMatrix[row][jcol] ^ annMatrix[curBit][jcol];
                } else {// it is not occupied
                    lowestOnePos[col] = row;
                    bZeroAnnotation = false;
                    break;
                }
            }
        }
        if (bZeroAnnotation) {
            dependence = true;
            break;
        }
    }
    return dependence;
}

bool edge_annotation_computing::CheckTwoGroupVectorOrthogonality(std::vector<std::vector<int> > &a_loops,
                                                                 std::vector<Annotation_Type> &a_ann,
                                                                 std::vector<Annotation_Type> &reduced_a_ann,
                                                                 std::vector<int> &a_lowestOne,
                                                                 std::vector<std::vector<int> > &b_loops,
                                                                 std::vector<Annotation_Type> &b_ann,
                                                                 std::vector<Annotation_Type> &reduced_b_ann,
                                                                 std::vector<int> &b_lowestOne) {
    //std::vector<Annotation_Type> a_annotations;
    //std::vector<Annotation_Type> b_annotations;

    // get annotation representation for a loops
    Annotation_Type tempAnnotation;
    for (unsigned int i = 0; i < a_loops.size(); i++) {
        if (!CheckZeroAnnotationsAndReturnAnnotation(a_loops[i], tempAnnotation)) {
            a_ann.push_back(tempAnnotation);
        } else {
            std::cout << "ERROR: containing zero annotations !" << std::endl;
            exit(0);
        }
    }
    //
    reduced_a_ann = a_ann;
    //
    if (ReduceMatrix(reduced_a_ann, a_lowestOne)) {
        std::cout << "ERROR: the loops are NOT independent !" << std::endl;
        exit(9);
    }
    //else
    //{
    //	std::cout << "a_loops are independent" << std::endl;
    //}
    // get annotation representation for b loops
    for (unsigned int i = 0; i < b_loops.size(); i++) {
        if (!CheckZeroAnnotationsAndReturnAnnotation(b_loops[i], tempAnnotation)) {
            b_ann.push_back(tempAnnotation);

        } else {
            std::cout << "ERROR: containing zero annotations !" << std::endl;
            exit(0);
        }
    }
    //
    reduced_b_ann = b_ann;
    //
    if (ReduceMatrix(reduced_b_ann, b_lowestOne)) {
        std::cout << "ERROR: the loops are NOT independent !" << std::endl;
        exit(9);
    }
    //else
    //{
    //	std::cout << "b_loops are independent" << std::endl;
    //}
    //
    bool indepBasis = true;
    string echo_msg = "These two basis vectors are independent";
    for (unsigned int i = 0; i < a_ann.size(); i++) {
        b_ann.push_back(a_ann[i]);
        if (GaussianElimination(b_ann)) {
            b_ann.pop_back();
            echo_msg = "These two vectors are NOT independent";
            indepBasis = false;
            break;
        }
        b_ann.pop_back();
    }
    //
    if (indepBasis) {
        for (unsigned int i = 0; i < b_ann.size(); i++) {
            a_ann.push_back(b_ann[i]);
            if (GaussianElimination(a_ann)) {
                a_ann.pop_back();
                echo_msg = "These two vectors are NOT independent";
                indepBasis = false;
                break;
            }
            a_ann.pop_back();
        }
    }
    if (!indepBasis) {
        std::cout << "ERROR: " << echo_msg << std::endl;
        exit(9);
    }
    //
    return false;
}

bool edge_annotation_computing::CheckTwoGroupVectorOrthogonality(std::vector<std::vector<int> > &a_loops,
                                                                 std::vector<Annotation_Type> &a_anno,
                                                                 std::vector<std::vector<int> > &b_loops,
                                                                 std::vector<Annotation_Type> &b_anno) {
    std::vector<Annotation_Type> a_annotations;
    std::vector<Annotation_Type> b_annotations;

    // get annotation representation for a loops
    Annotation_Type tempAnnotation;
    for (unsigned int i = 0; i < a_loops.size(); i++) {
        if (!CheckZeroAnnotationsAndReturnAnnotation(a_loops[i], tempAnnotation)) {
            a_annotations.push_back(tempAnnotation);
        } else {
            std::cout << "Containing zero annotations" << std::endl;
            exit(0);
        }
    }
    //
    a_anno = a_annotations;
    //
    if (GaussianElimination(a_annotations)) {
        std::cout << "ERROR : the loops are NOT independent !" << std::endl;
        exit(9);
    }
    //else
    //{
    //	std::cout << "ERROR : the loops are independent !" << std::endl;
    //}
    // get annotation representation for b loops
    for (unsigned int i = 0; i < b_loops.size(); i++) {
        if (!CheckZeroAnnotationsAndReturnAnnotation(b_loops[i], tempAnnotation)) {
            b_annotations.push_back(tempAnnotation);

        } else {
            std::cout << "ERROR: containing zero annotations !" << std::endl;
            exit(9);
        }
    }
    //
    b_anno = b_annotations;
    //
    if (GaussianElimination(b_annotations)) {
        std::cout << "ERROR : the loops are NOT independent !" << std::endl;
        exit(9);
    }
    //else
    //{
    //	std::cout << "b_loops are independent" << std::endl;
    //}
    //
    bool indepBasis = true;
    string echo_msg = "These two basis vectors are independent";
    for (unsigned int i = 0; i < a_annotations.size(); i++) {
        b_annotations.push_back(a_annotations[i]);
        if (GaussianElimination(b_annotations)) {
            b_annotations.pop_back();
            echo_msg = "These two vectors are NOT independent";
            indepBasis = false;
            break;
        }
        b_annotations.pop_back();
    }
    //
    if (indepBasis) {
        for (unsigned int i = 0; i < b_annotations.size(); i++) {
            a_annotations.push_back(b_annotations[i]);
            if (GaussianElimination(a_annotations)) {
                a_annotations.pop_back();
                echo_msg = "These two vectors are NOT independent";
                indepBasis = false;
                break;
            }
            a_annotations.pop_back();
        }
    }
    if (!indepBasis) {
        std::cout << "ERROR: " << echo_msg << std::endl;
        exit(9);
    }
    //

    return false;
}

bool edge_annotation_computing::GaussianElimination(std::vector<Annotation_Type> &in_vec) {
    std::vector<Annotation_Type> vec(in_vec);
    //std::set<int> processed_row;
    ////
    //for (unsigned int i = 0; i < vec.size(); i++)
    //{
    //	processed_row.insert(i);
    //}
    int row = 0;
    int col = 1;
    bool linearly_dependent = false;
    int vec_size = 0;
    if (!in_vec.empty())
        vec_size = (int) in_vec[0].size();
    //while(!processed_row.empty())
    //{
    //row = *(processed_row.begin());
    //processed_row.erase(processed_row.begin());
    //
    for (unsigned int k = 0; k < vec.size(); k++) {
        row = k;
        col = -1;
        for (int j = 0; j < vec_size; j++) {
            if (vec[row][j]) {
                col = j;
                break;
            }
        }
        //
        if (col == -1) {// linearly dependent
            linearly_dependent = true;
            break;
        }
        for (unsigned int i = 0; i < vec.size(); i++) {
            if (row != i && vec[i][col]) {// sum this row with the candidate row to make this bit zero
                for (int j = 0; j < vec_size; j++)
                    vec[i][j] = (vec[i][j] ^ vec[row][j]);//% 2;
            }
        }
    }
    return linearly_dependent;
}

bool edge_annotation_computing::CheckZeroAnnotationsAndReturnAnnotation(std::set<int> &cycle, Annotation_Type &res) {
    // set res as zero vector
    res.clear();
    res.resize(vec_size, 0);
    int sum_res = 0;
    for (std::set<int>::iterator sIter = cycle.begin();
         sIter != cycle.end();
         sIter++) {
        for (int i = 0; i < vec_size; i++) {
            res[i] = (res[i] + edge_annotations[*sIter][i]) % 2;
        }
    }
    for (int i = 0; i < vec_size; i++) {
        sum_res += res[i];
    }
    if (!sum_res)
        return true;

    return false;
}

bool edge_annotation_computing::CheckZeroAnnotationsAndReturnAnnotation(std::vector<int> &cycle, Annotation_Type &res) {
    // set res as zero vector
    res.clear();
    res.resize(vec_size, 0);
    int sum_res = 0;

    for (std::vector<int>::iterator sIter = cycle.begin();
         sIter != cycle.end();
         sIter++) {
        for (int i = 0; i < vec_size; i++) {
            res[i] = (res[i] + edge_annotations[*sIter][i]) % 2;
        }
    }
    for (int i = 0; i < vec_size; i++) {
        sum_res += res[i];
    }
    if (!sum_res)
        return true;

    return false;
}

void edge_annotation_computing::ReadEdgeAnnotationsFromFile(const char *file_name) {
    std::cout << "reading... " << file_name << std::endl;

    std::ifstream ifile;
    ifile.open(file_name, std::ifstream::in);
    long long fileSize = 0;
    char *fBuf = NULL;

    std::string sBuf;
    //
    int edgeNum = (int) meshPtr->vecEdge.size();
    /*
    reserve space for edge_annotations
    */
    //edge_annotations.resize(edgeNum);
    //
    std::stringstream sstr(std::stringstream::in | std::stringstream::out);
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

        sstr >> vec_size;
        //std::vector<int> tmpInt(vec_size, 0);
        Annotation_Type tmpInt(vec_size, 0);
        /*
        reserve space for edge_annotations
        */
        Annotation_Type temp_z2(vec_size, 0);//(vec_size, tmpInt.begin(), tmpInt.begin() + vec_size);
        edge_annotations.resize(edgeNum);

        /*	for (unsigned int i = 0; i < edgeNum; i++)
            {
                edge_annotations.push_back(temp_z2);
            }*/
        //

        int vertex_a = 0;
        int vertex_b = 0;
        int edge_id = 0;
        for (int i = 0; i < edgeNum; i++) {
            sstr >> vertex_a >> vertex_b;
            //
            // find the edge id here
            edge_id = -1;
            for (std::vector<int>::iterator eidIter = meshPtr->vecVertex[vertex_a].adjEdges.begin();
                 eidIter != meshPtr->vecVertex[vertex_a].adjEdges.end(); eidIter++) {
                if (meshPtr->vecEdge[*eidIter].v0 == vertex_b ||
                    meshPtr->vecEdge[*eidIter].v1 == vertex_b) {
                    edge_id = *eidIter;
                    break;
                }
            }
            //
            if (edge_id < 0) {
                std::cout << "can not find edge (" << vertex_a << "," << vertex_b << ")" << std::endl;
                exit(0);
            }
            sBuf.clear();
            sBuf = "";
            sstr >> sBuf;
            for (int j = 0; j < vec_size; j++)
                tmpInt[j] = sBuf[j] - '0';
            //
            edge_annotations[edge_id] = tmpInt;
        }

        sstr.clear();
    } else {
        std::cout << "Can NOT open file " << file_name << std::endl;
        exit(0);
    }
    //
    std::cout << "Done... " << edgeNum << std::endl;
    //
}

bool edge_annotation_computing::ComputeShortestCanonicalLoops(std::vector<std::vector<int> > &basis_loops,
                                                              std::vector<Annotation_Type> &basis_annotations,
                                                              std::vector<std::set<int> > &canonical_loops,
                                                              std::vector<int> &non_tree_edges,
                                                              std::vector<Annotation_Type> &canonical_loops_annotation,
                                                              std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                                              std::vector<std::set<int> > &short_basis_loops,
                                                              std::vector<int> &basis_non_tree_edges) {
    bool success = true;
    basis_non_tree_edges.clear();
    basis_non_tree_edges.reserve(basis_loops.size());
    //
    //std::vector<Annotation_Type> basis_annotations(in_basis_annotations);//basis_loops.size());
    // get annotations for each edge
    //Annotation_Type sumAnnotation;
    //for (unsigned int i = 0; i < basis_loops.size(); i++)
    //{
    //	if (CheckZeroAnnotationsAndReturnAnnotation(basis_loops[i], sumAnnotation))
    //	{// conati
    //		std::cout << "Contains zero annotation in basis " << std::endl;
    //		exit(0);
    //	}
    //	else
    //	{
    //		basis_annotations[i] = sumAnnotation;
    //	}
    //}
    //
    //
    std::vector<Annotation_Type> picked_vec_annotation;
    //
//	for (unsigned int i = 0; i < tst_vec.size(); i++)
//	{
//		basis.push_back(tst_vec[i]);
//		if (!CGAL::linearly_independent(basis.begin(), basis.end()))
//		{
//			picked_vec.push_back(tst_vec[i]);
//			if (!CGAL::linearly_independent(picked_vec.begin(), picked_vec.end()))
//			{//
//				picked_vec.pop_back();
//			}
//			if (picked_vec.size() == want_basis_vec_num)
//				break;
//		}
//		basis.pop_back();
//	}
    for (std::set<std::pair<float, int>, myFloatIntPairLessThan>::iterator sIter = sorted_loops.begin();
         sIter != sorted_loops.end();
         sIter++) {
        basis_annotations.push_back(canonical_loops_annotation[sIter->second]);
        //
        if (GaussianElimination(basis_annotations)) {// this vector is in the space spanned by the basis
            picked_vec_annotation.push_back(canonical_loops_annotation[sIter->second]);
            if (GaussianElimination(picked_vec_annotation)) {// this vector is in the space spanned by the picked basis
                picked_vec_annotation.pop_back();
            } else {// pick it as a new vector
                short_basis_loops.push_back(canonical_loops[sIter->second]);
                basis_non_tree_edges.push_back(non_tree_edges[sIter->second]);
                if (short_basis_loops.size() == basis_loops.size()) {
                    basis_annotations.pop_back();
                    break;
                }
            }
        }
        basis_annotations.pop_back();
    }
    if (short_basis_loops.size() < basis_loops.size()) {// need to include the extra ones
        success = false;
        for (unsigned int i = 0; i < basis_loops.size(); i++) {
            picked_vec_annotation.push_back(basis_annotations[i]);
            if (GaussianElimination(picked_vec_annotation)) {//this vector is in the space spanned by the picked basis
                picked_vec_annotation.pop_back();
            } else {
                //exit(0);
                std::set<int> tmp(basis_loops[i].begin(), basis_loops[i].end());
                short_basis_loops.push_back(tmp);//basis_loops[i]);
                //
                basis_non_tree_edges.push_back(-1);
            }
        }
    }
    if (short_basis_loops.size() != basis_loops.size()) {
        std::cout << "larger or small size " << std::endl;
        exit(9);
    }
    return success;
}

bool edge_annotation_computing::ComputeShortestCanonicalLoops(std::vector<std::vector<int> > &basis_loops,
                                                              std::vector<Annotation_Type> &basis_annotations,
                                                              std::vector<Annotation_Type> &reduced_basis_annotations,
                                                              std::vector<int> &basis_LowestOnePos,
                                                              std::vector<std::set<int> > &canonical_loops,
                                                              std::vector<int> &non_tree_edges,
                                                              std::vector<Annotation_Type> &canonical_loops_annotation,
                                                              std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                                              std::vector<std::set<int> > &short_basis_loops,
                                                              std::vector<int> &basis_non_tree_edges) {
    bool success = true;
    basis_non_tree_edges.clear();
    basis_non_tree_edges.reserve(basis_loops.size());
    //
    //std::vector<Annotation_Type> basis_annotations(in_basis_annotations);//basis_loops.size());
    // get annotations for each edge
    //Annotation_Type sumAnnotation;
    //for (unsigned int i = 0; i < basis_loops.size(); i++)
    //{
    //	if (CheckZeroAnnotationsAndReturnAnnotation(basis_loops[i], sumAnnotation))
    //	{// conati
    //		std::cout << "Contains zero annotation in basis " << std::endl;
    //		exit(0);
    //	}
    //	else
    //	{
    //		basis_annotations[i] = sumAnnotation;
    //	}
    //}
    //
    //
    std::vector<Annotation_Type> picked_vec_annotation;
    //
//	for (unsigned int i = 0; i < tst_vec.size(); i++)
//	{
//		basis.push_back(tst_vec[i]);
//		if (!CGAL::linearly_independent(basis.begin(), basis.end()))
//		{
//			picked_vec.push_back(tst_vec[i]);
//			if (!CGAL::linearly_independent(picked_vec.begin(), picked_vec.end()))
//			{//
//				picked_vec.pop_back();
//			}
//			if (picked_vec.size() == want_basis_vec_num)
//				break;
//		}
//		basis.pop_back();
//	}
    for (std::set<std::pair<float, int>, myFloatIntPairLessThan>::iterator sIter = sorted_loops.begin();
         sIter != sorted_loops.end();
         sIter++) {
        //basis_annotations.push_back(canonical_loops_annotation[sIter->second]);
        //
        if (IndependenceCheck(reduced_basis_annotations, basis_LowestOnePos,
                              canonical_loops_annotation[sIter->second])) {
            //if (GaussianElimination(basis_annotations))
            //{// this vector is in the space spanned by the basis
            picked_vec_annotation.push_back(canonical_loops_annotation[sIter->second]);
            if (GaussianElimination(picked_vec_annotation)) {// this vector is in the space spanned by the picked basis
                picked_vec_annotation.pop_back();
            } else {// pick it as a new vector
                short_basis_loops.push_back(canonical_loops[sIter->second]);
                basis_non_tree_edges.push_back(non_tree_edges[sIter->second]);
                if (short_basis_loops.size() == basis_loops.size()) {
                    //basis_annotations.pop_back();
                    break;
                }
            }
        }
        //basis_annotations.pop_back();
    }
    if (short_basis_loops.size() < basis_loops.size()) {// need to include the extra ones
        success = false;
        for (unsigned int i = 0; i < basis_loops.size(); i++) {
            picked_vec_annotation.push_back(basis_annotations[i]);
            if (GaussianElimination(picked_vec_annotation)) {//this vector is in the space spanned by the picked basis
                picked_vec_annotation.pop_back();
            } else {
                //exit(0);
                std::set<int> tmp(basis_loops[i].begin(), basis_loops[i].end());
                short_basis_loops.push_back(tmp);//basis_loops[i]);
                //
                basis_non_tree_edges.push_back(-1);
            }
        }
    }
    if (short_basis_loops.size() != basis_loops.size()) {
        std::cout << "larger or small size " << std::endl;
        exit(9);
    }
    return success;
}
