/*
(c) 2012 Fengtao Fan
*/
#include "DijkstraAlgorithm.h"
#include <algorithm>
#include <limits>

namespace FIBO = FibonacciHeap;
namespace DijkstraAlgorithm {

    void DijkstraShortestPath(const int srcIndex,
                              std::vector<FIBO::Node *> &vecNode,
                              enum FIBO::State unTouchedState,
                              enum FIBO::State scannedState,
                              std::vector<int> &parents_vertex,
                              std::vector<int> &parents_edge,
                              std::vector<float> &distance) {// only for connected graph

        FIBO::FibonacciHeap *heap = new FIBO::FibonacciHeap();

        vecNode[srcIndex]->state = FIBO::LABELED;
        vecNode[srcIndex]->key = 0.f;

        heap->insertVertex(vecNode[srcIndex]);

        //
        parents_vertex[srcIndex] = -1;
        parents_edge[srcIndex] = -1;
        distance[srcIndex] = 0.f;
        // Scan
        do {
            // Delete minimum path
            FIBO::Node *v = heap->deleteMin();

            if (v->key == std::numeric_limits<float>::max())
                break; // it is not connected

            v->state = scannedState;//SCANNED;

            for (unsigned int j = 0; j < v->incomingEdges.size(); j++) {
                FIBO::Edge *currentEdge = v->incomingEdges[j];
                FIBO::Node *headOfCurrentEdge = currentEdge->tail;

                if (headOfCurrentEdge->state != scannedState)//SCANNED)
                {
                    if (headOfCurrentEdge->state == unTouchedState)//UNLABELED)
                    {
                        // Insert a vertex with infinite key
                        headOfCurrentEdge->state = FIBO::LABELED;
                        headOfCurrentEdge->pred = v;
                        headOfCurrentEdge->key = v->key + currentEdge->length;
                        heap->insertVertex(headOfCurrentEdge);
                        // record the parents info
                        parents_vertex[headOfCurrentEdge->data] = v->data;
                        parents_edge[headOfCurrentEdge->data] = currentEdge->edge_idx;
                    } else {
                        if (headOfCurrentEdge->key > v->key + currentEdge->length) {
                            // decrease the key of a vertex with finite key
                            headOfCurrentEdge->pred = v;
                            heap->decreaseKey(v->key + currentEdge->length, headOfCurrentEdge);
                            //
                            // record the parents info
                            parents_vertex[headOfCurrentEdge->data] = v->data;
                            parents_edge[headOfCurrentEdge->data] = currentEdge->edge_idx;
                        }
                    }
                }
            }
        } while (!heap->isEmpty());

        while (!heap->isEmpty())
            heap->deleteMin();
        // record the distance
        for (unsigned int vid = 0; vid < vecNode.size(); vid++) {
            distance[vid] = vecNode[vid]->key;
        }
    }

    void DijkstraShortestPath(const int srcIndex,
                              std::vector<FIBO::Node *> &vecNode,
                              enum FIBO::State unTouchedState,
                              enum FIBO::State scannedState,
                              std::vector<int> &parents_vertex,
                              std::vector<int> &parents_edge,
                              std::vector<int> &distance) {// only for connected graph

        FIBO::FibonacciHeap *heap = new FIBO::FibonacciHeap();

        vecNode[srcIndex]->state = FIBO::LABELED;
        vecNode[srcIndex]->key = 0;

        heap->insertVertex(vecNode[srcIndex]);

        //
        parents_vertex[srcIndex] = -1;
        parents_edge[srcIndex] = -1;
        distance[srcIndex] = 0;
        // Scan
        do {
            // Delete minimum path
            FIBO::Node *v = heap->deleteMin();

            if (v->key == std::numeric_limits<int>::max())
                break; // it is not connected

            v->state = scannedState;//SCANNED;

            for (unsigned int j = 0; j < v->incomingEdges.size(); j++) {
                FIBO::Edge *currentEdge = v->incomingEdges[j];
                FIBO::Node *headOfCurrentEdge = currentEdge->tail;

                if (headOfCurrentEdge->state != scannedState)//SCANNED)
                {
                    if (headOfCurrentEdge->state == unTouchedState)//UNLABELED)
                    {
                        // Insert a vertex with infinite key
                        headOfCurrentEdge->state = FIBO::LABELED;
                        headOfCurrentEdge->pred = v;
                        headOfCurrentEdge->key = v->key + currentEdge->length;
                        heap->insertVertex(headOfCurrentEdge);
                        // record the parents info
                        parents_vertex[headOfCurrentEdge->data] = v->data;
                        parents_edge[headOfCurrentEdge->data] = currentEdge->edge_idx;
                    } else {
                        if (headOfCurrentEdge->key > v->key + currentEdge->length) {
                            // decrease the key of a vertex with finite key
                            headOfCurrentEdge->pred = v;
                            heap->decreaseKey(v->key + currentEdge->length, headOfCurrentEdge);
                            //
                            // record the parents info
                            parents_vertex[headOfCurrentEdge->data] = v->data;
                            parents_edge[headOfCurrentEdge->data] = currentEdge->edge_idx;
                        }
                    }
                }
            }
        } while (!heap->isEmpty());

        while (!heap->isEmpty())
            heap->deleteMin();
        // record the distance
        for (unsigned int vid = 0; vid < vecNode.size(); vid++) {
            distance[vid] = vecNode[vid]->key;
        }
    }

}
