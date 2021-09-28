/*
(c) 2012 Fengtao Fan
*/
#ifndef _DIJKSTRA_ALGORITHM_H_
#define _DIJKSTRA_ALGORITHM_H_

#include <iostream>
#include <vector>
#include "FibonacciHeap.h"

typedef float DIJ_real_type;
namespace FIBO = FibonacciHeap;
namespace DijkstraAlgorithm {

    void DijkstraShortestPath(const int srcIndex,
                              std::vector<FIBO::Node *> &vecNode,
                              enum FIBO::State unTouchedState,
                              enum FIBO::State scannedState,
                              std::vector<int> &parents_vertex,
                              std::vector<int> &parents_edge,
                              std::vector<float> &distance);

    //
    void DijkstraShortestPath(const int srcIndex,
                              std::vector<FIBO::Node *> &vecNode,
                              enum FIBO::State unTouchedState,
                              enum FIBO::State scannedState,
                              std::vector<int> &parents_vertex,
                              std::vector<int> &parents_edge,
                              std::vector<int> &distance);

}
#endif //_DIJKSTRA_ALGORITHM_H_
