/*
(c) 2012 Fengtao Fan
*/
#ifndef _My_ANNOTATION_COMP_H_
#define _My_ANNOTATION_COMP_H_

#include <vector>
#include "SimpleMesh.h"

extern int ComputeAnnotation(_SimpleMesh &inMesh, std::vector<std::vector<char> > &vecEdgeAnno, std::vector<int> &tris);

#endif
