/*
(c) 2012 Fengtao Fan
*/
#include "InverseLinkNumberMatrix.h"
#include <fstream>

void InverseLinkNumberMatrix::add_rows(const int dst_row, const int src_row, const int multiplier,
                                       std::vector<std::list<std::pair<int, int> > > &mat) {
    // elements in each row list are ordered by the colomn number
    std::list<std::pair<int, int> >::iterator srcIter = mat[src_row].begin();
    std::list<std::pair<int, int> >::iterator dstIter = mat[dst_row].begin();
    // scan both list by the increasing order of columns
    while (srcIter != mat[src_row].end() && dstIter != mat[dst_row].end()) {
        if (srcIter->first < dstIter->first) {
            // insert a new element with value of srcIter->second * muliplier
            mat[dst_row].insert(dstIter, std::pair<int, int>(srcIter->first, srcIter->second * multiplier));
            srcIter++;
        } else {
            if (srcIter->first > dstIter->first) {
                dstIter++;
            } else {// both equal
                dstIter->second = dstIter->second + srcIter->second * multiplier;
                if (dstIter->second == 0) {// remove it from dst row
                    dstIter = mat[dst_row].erase(
                            dstIter); // dstIter now points to the element after the erasted element
                }
                srcIter++;
            }
        }
    }
    if (dstIter == mat[dst_row].end() &&
        srcIter != mat[src_row].end()) {// add all elements starting srcIter to dst rows
        do {
            mat[dst_row].push_back(std::pair<int, int>(srcIter->first, srcIter->second * multiplier));
            srcIter++;
        } while (srcIter != mat[src_row].end());
    }
    //
    return;
}

void InverseLinkNumberMatrix::ComputeInverseMatrix() {
    ComputeInverseMatrix(n, org_mat);
}

void InverseLinkNumberMatrix::ComputeInverseMatrix(const int halfDim,
                                                   std::vector<std::list<std::pair<int, int> > > &inOrgMat) {
    /*
      * * * 1 * *
      * * * 0 1 *
      * * * 0 0 1
      1 0 0 0 0 0
      * 1 0 0 0 0
      * * 1 0 0 0
      The submatrix M(n:2n, 1:n) and M(1:n, n:2n) are triangular matrix
    */
    // initialize inverse matrix as identity matrix
    inv_mat.resize(2 * halfDim);
    for (int row = 0; row < 2 * halfDim; row++) {
        inv_mat[row].push_back(std::pair<int, int>(row, 1));
    }
    //
    for (int i = halfDim; i < 2 * halfDim; i++) {
        const int pivot_value = inOrgMat[i].front().second; // pivot_value == 1 or -1
        for (int row = 0; row < 2 * halfDim; row++) {
            if (row != i) {
                std::list<std::pair<int, int> >::iterator listIter = inOrgMat[row].begin();
                int multiplier = 0;
                for (; listIter != inOrgMat[row].end(); listIter++) {
                    if (listIter->first == inOrgMat[i].front().first) {
                        break;
                    }
                }
                //
                if (listIter != inOrgMat[row].end()) {
                    multiplier = -(pivot_value * listIter->second);
                    //
                    add_rows(row, i, multiplier, inOrgMat);
                    add_rows(row, i, multiplier, inv_mat);
                }
            }
        }
    }
    //
    for (int i = halfDim - 1; i >= 0; i--) {
        const int pivot_value = inOrgMat[i].front().second;
        for (int row = 0; row < halfDim; row++) {
            if (row != i) {
                std::list<std::pair<int, int> >::iterator listIter = inOrgMat[row].begin();
                int multiplier = 0;
                for (; listIter != inOrgMat[row].end(); listIter++) {
                    if (listIter->first == inOrgMat[i].front().first) {
                        break;
                    }
                }
                //
                if (listIter != inOrgMat[row].end()) {
                    multiplier = -(pivot_value * listIter->second);
                    //
                    add_rows(row, i, multiplier, inOrgMat);
                    add_rows(row, i, multiplier, inv_mat);
                }
            }
        }
    }
}

void InverseLinkNumberMatrix::PrintMatrix(std::vector<std::list<std::pair<int, int> > > &mat) {

    const int dim = (int) mat.size();
    std::cout << "M=zeros(" << dim << "," << dim << ");" << std::endl;
    for (int i = 0; i < dim; i++) {
        //std::cout << "Vert Loop " << i << " link against " << std::endl;
        for (std::list<std::pair<int, int> >::iterator listIter = mat[i].begin();
             listIter != mat[i].end(); listIter++) {
            std::cout << "M(" << i + 1 << "," << listIter->first + 1 << ")=" << listIter->second << ";";

        }
        std::cout << std::endl;
    }
    return;
}

void InverseLinkNumberMatrix::PrintMatrix(const char *pFileName, std::vector<std::list<std::pair<int, int> > > &mat) {
    std::ofstream ofile;
    ofile.open(pFileName, std::ifstream::out);
    if (ofile.is_open()) {
        const int dim = (int) mat.size();
        ofile << "M=zeros(" << dim << "," << dim << ");" << std::endl;
        for (int i = 0; i < dim; i++) {
            //std::cout << "Vert Loop " << i << " link against " << std::endl;
            for (std::list<std::pair<int, int> >::iterator listIter = mat[i].begin();
                 listIter != mat[i].end(); listIter++) {
                ofile << "M(" << i + 1 << "," << listIter->first + 1 << ")=" << listIter->second << ";";

            }
            ofile << std::endl;
        }
        //
        ofile.close();
    }
    return;
}
