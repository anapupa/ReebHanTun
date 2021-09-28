/*
(c) 2012 Fengtao Fan
*/
/*compute the inverse matrix by using the speciality of the linking number matrix*/
/*
  * * * 1 * *
  * * * 0 1 *
  * * * 0 0 1
  1 0 0 0 0 0
  * 1 0 0 0 0
  * * 1 0 0 0
  The submatrix M(n:2n, 1:n) and M(1:n, n:2n) are triangular matrix
*/
#ifndef _INVERSE_LINK_NUMBER_MATRIX_H_
#define _INVERSE_LINK_NUMBER_MATRIX_H_

#include <iostream>
#include <vector>
#include <list>

class InverseLinkNumberMatrix {
public:
    InverseLinkNumberMatrix() {
        n = 0;
    }

    InverseLinkNumberMatrix(const InverseLinkNumberMatrix &lhs) {
        org_mat = lhs.org_mat;
        inv_mat = lhs.inv_mat;
        n = lhs.n;
    }

    ~InverseLinkNumberMatrix() {
        Clear();
    }

    void SetOrgMatrix(const std::vector<std::list<std::pair<int, int> > > &inOrgMat) {
        org_mat = inOrgMat;
    }

    void SetMatrixHalfDim(const int inHalfDim) {
        n = inHalfDim;
    }

    void Clear() {
        org_mat.clear();
        inv_mat.clear();
    }

    //
    void add_rows(const int dst_row, const int src_row, const int multiplier,
                  std::vector<std::list<std::pair<int, int> > > &mat);

    void ComputeInverseMatrix();

    void ComputeInverseMatrix(const int n, std::vector<std::list<std::pair<int, int> > > &inOrgMat);

    void PrintMatrix(std::vector<std::list<std::pair<int, int> > > &mat);

    void PrintMatrix(const char *pFileName, std::vector<std::list<std::pair<int, int> > > &mat);

public:
    std::vector<std::list<std::pair<int, int> > > inv_mat;
    std::vector<std::list<std::pair<int, int> > > org_mat;
    int n; // matrix is 2nx2n
};

#endif //_INVERSE_LINK_NUMBER_MATRIX_H_
