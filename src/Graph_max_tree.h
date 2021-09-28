/*
(c) 2012 Fengtao Fan
*/
#ifndef _MY_GRAPH_MAX_TREE_H_
#define _MY_GRAPH_MAX_TREE_H_

#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
#include <iterator>

#include <boost/pending/disjoint_sets.hpp>

class _simpGraphNode_max_tree {
public:
    _simpGraphNode_max_tree() : nIndexInReebGraph(-1) {
        value = 0.0;
        // set the value
    }

    _simpGraphNode_max_tree(const int index) : nIndexInReebGraph(index) {
        value = 0.0;
    }

    _simpGraphNode_max_tree(const int index, const float inVal) : nIndexInReebGraph(index), value(inVal) {

    }

    _simpGraphNode_max_tree(const _simpGraphNode_max_tree &rhs) {
        nIndexInReebGraph = rhs.nIndexInReebGraph;
        //
        value = rhs.value;

        adjVertexIdVec.assign(rhs.adjVertexIdVec.begin(), rhs.adjVertexIdVec.end());
        adjArcIdVec.assign(rhs.adjArcIdVec.begin(), rhs.adjArcIdVec.end());
        //
        /*	std::copy (	rhs.adjList.begin(),
                        rhs.adjList.end(),
                        std::inserter(adjList, adjList.begin()));*/
    }

    void CopyElementsVector(std::vector<int> &dst, const std::vector<int> &inVec) {
        dst.clear();
        dst.assign(inVec.begin(), inVec.end());
    }

    _simpGraphNode_max_tree &operator=(const _simpGraphNode_max_tree &rhs) {
        if (this == &rhs)
            return *this;
        CopyElementsVector(adjVertexIdVec, rhs.adjVertexIdVec);
        //
        CopyElementsVector(adjArcIdVec, rhs.adjArcIdVec);
        //
        nIndexInReebGraph = rhs.nIndexInReebGraph;
        //
        value = rhs.value;
        //
        return *this;
    }

    ~_simpGraphNode_max_tree() {
        adjVertexIdVec.clear();
        adjArcIdVec.clear();
    }

    /************************************/
    bool is_neighbor(const int u) {
        std::set<int> sorted_ver_id(adjVertexIdVec.begin(), adjVertexIdVec.end());
        std::set<int>::const_iterator sIter;
        sIter = sorted_ver_id.find(u);
        if (sIter == sorted_ver_id.end())
            return false;
        return true;
    }

    bool is_neighbor(const int u) const {
        std::set<int> sorted_ver_id(adjVertexIdVec.begin(), adjVertexIdVec.end());
        std::set<int>::const_iterator sIter;
        sIter = sorted_ver_id.find(u);
        if (sIter == sorted_ver_id.end())
            return false;
        return true;
    }

    int degree() {
        return int(adjVertexIdVec.size());
    }

public:
    int nIndexInReebGraph;
    std::vector<int> adjVertexIdVec;
    std::vector<int> adjArcIdVec;
    //
    float value; // the height of this point
};

class _simpGraph_max_tree {
private:
    struct FourCompLessThan {
        bool operator()(std::pair<float, std::pair<int, std::pair<int, int> > > &lhs,
                        std::pair<float, std::pair<int, std::pair<int, int> > > &rhs) {
            return lhs.first < rhs.first;
        }
    };

public:
    //
    typedef std::pair<int, std::pair<int, int> > TripleInt;
    typedef std::pair<float, TripleInt> FourComp;

    //
    _simpGraph_max_tree() : edge_size(0) {
    }

    _simpGraph_max_tree(const int n) : edge_size(0) {
        _simpGraphNode_max_tree tmpNode;
        vecNode.reserve(n);
        for (int i = 0; i < n; i++) {
            tmpNode.nIndexInReebGraph = i;
            vecNode.push_back(tmpNode);
        }
    }

    _simpGraph_max_tree(const _simpGraph_max_tree &rhs) {
        vecNode.assign(rhs.vecNode.begin(),
                       rhs.vecNode.end());
        edge_size = rhs.edge_size;
        vecEdgeWeight.assign(rhs.vecEdgeWeight.begin(),
                             rhs.vecEdgeWeight.end());
    }

    _simpGraph_max_tree &operator=(const _simpGraph_max_tree &rhs) {
        if (this == &rhs)
            return *this;
        vecNode.clear();
        vecNode.assign(rhs.vecNode.begin(),
                       rhs.vecNode.end());
        edge_size = rhs.edge_size;
        vecEdgeWeight.assign(rhs.vecEdgeWeight.begin(),
                             rhs.vecEdgeWeight.end());
        return *this;
    }

    ~_simpGraph_max_tree() {
        vecNode.clear();
    }

    /********************************/
    void AddEdge(const int i, const int j) {
        vecNode[i].adjVertexIdVec.push_back(j);
        vecNode[j].adjVertexIdVec.push_back(i);
        // use the default order for edge
        vecNode[i].adjArcIdVec.push_back(edge_size);
        vecNode[j].adjArcIdVec.push_back(edge_size);
        //
        edge_size++;
        //
        return;
    }

    void AddEdge(const int i, const int j, const int eid) {
        vecNode[i].adjVertexIdVec.push_back(j);
        vecNode[j].adjVertexIdVec.push_back(i);
        // use the input order for edge
        vecNode[i].adjArcIdVec.push_back(eid);
        vecNode[j].adjArcIdVec.push_back(eid);
        //
        edge_size++;
        return;
    }

    /********************************/
    void AddEdgeW(const int i, const int j, const int eid, float w) {
        vecNode[i].adjVertexIdVec.push_back(j);
        vecNode[j].adjVertexIdVec.push_back(i);
        // use the input order for edge
        vecNode[i].adjArcIdVec.push_back(eid);
        vecNode[j].adjArcIdVec.push_back(eid);
        //
        edge_size++;
        //
        vecEdgeWeight.push_back(FourComp(w, TripleInt(eid, std::pair<int, int>(i, j))));
        return;
    }

    /************used only in the tree****************/
    int ArcId(const int i, const int j) {
        int src = i;
        int dst = j;
        if (vecNode[i].adjArcIdVec.size() > vecNode[j].adjArcIdVec.size()) {
            src = j;
            dst = i;
        }
        unsigned int aid = 0;
        for (; aid < vecNode[src].adjVertexIdVec.size(); aid++)
            if (vecNode[src].adjVertexIdVec[aid] == dst)
                break;
        if (aid > 0 && aid == vecNode[src].adjVertexIdVec.size()) {
            std::cout << "No edge connect this two vertices" << std::endl;
            exit(0);
        }
        return vecNode[src].adjArcIdVec[aid];

    }

    /********************************/
    void SortEdges() {
        std::make_heap(vecEdgeWeight.begin(), vecEdgeWeight.end(), FourCompLessThan());
        std::sort_heap(vecEdgeWeight.begin(), vecEdgeWeight.end(), FourCompLessThan());
        return;
    }

    void AddNode(const int nIndexInReebGraph) {
        _simpGraphNode_max_tree tmpNode(nIndexInReebGraph);
        vecNode.push_back(tmpNode);
    }

    void AddNode(const int nIndexInReebGraph, const float hval) {
        _simpGraphNode_max_tree tmpNode(nIndexInReebGraph, hval);
        vecNode.push_back(tmpNode);
    }

    void InitNodes(const int n) {
        for (int nIndexInReebGraph = 0; nIndexInReebGraph < n; nIndexInReebGraph++) {
            _simpGraphNode_max_tree tmpNode(nIndexInReebGraph, nIndexInReebGraph);
            vecNode.push_back(tmpNode);
        }
    }

    void Clear() {
        std::vector<_simpGraphNode_max_tree> tmp;
        vecNode.swap(tmp);
        vecEdgeWeight.clear();
    }

    //
    void KruskalAlg(std::vector<std::pair<int, std::pair<int, int> > > &max_span_tree,
                    std::vector<std::pair<int, std::pair<int, int> > > &non_tree_edges) {
        //
        max_span_tree.clear();
        max_span_tree.reserve(vecNode.size() + 1);
        //
        SortEdges();
        //
        std::vector<int> rank(vecNode.size() + 2);
        std::vector<int> parent(vecNode.size() + 2);
        //
        boost::disjoint_sets<int *, int *> ds(&rank[0], &parent[0]);

        for (int i = 0; i < vecNode.size(); i++)
            ds.make_set(i);
        //
        std::vector<FourComp>::reverse_iterator re_iter = vecEdgeWeight.rbegin();
        for (; re_iter != vecEdgeWeight.rend(); re_iter++) {
            //
            std::pair<int, int> node_pair = re_iter->second.second;
            if (max_span_tree.size() < vecNode.size() - 1) {
                int src = ds.find_set(node_pair.first);
                int dst = ds.find_set(node_pair.second);
                if (src != dst) {
                    max_span_tree.push_back(re_iter->second);
                    ds.link(src, dst);
                } else {
                    non_tree_edges.push_back(re_iter->second);
                }
            } else {
                non_tree_edges.push_back(re_iter->second);
            }
            //
            //if (max_span_tree.size() == vecNode.size() - 1)
            //	break;
        }
        return;
    }

    //
    void BreathFirstSearch(const int sourceVertex, std::vector<int> &parents) {// color : 0 --- untouched
        //			1 --- grey (touched)
        //			2 --- black (finished)
        std::vector<int> color;
        color.resize(vecNode.size());
        parents.resize(vecNode.size());
        // initializing color and parents
        for (int i = 0; i < int(color.size()); i++) {
            color[i] = 0; // untouched
            parents[i] = -1;// no parents
        }
        // set up a queue
        std::queue<int> Q;
        // push up source into queue
        // set its color as touched
        color[sourceVertex] = 1;
        Q.push(sourceVertex);
        while (!Q.empty()) {
            int curVertex = Q.front();
            Q.pop();
            for (std::vector<int>::iterator sIter = vecNode[curVertex].adjVertexIdVec.begin();
                 sIter != vecNode[curVertex].adjVertexIdVec.end();
                 ++sIter) {
                if (color[*sIter] == 0) {// current vertex is not touched
                    color[*sIter] = 1;
                    parents[*sIter] = curVertex;
                    Q.push(*sIter);
                }
            }
            color[curVertex] = 2;
        }
    }
    //
public:
    std::vector<_simpGraphNode_max_tree> vecNode;
    int edge_size;
    std::vector<FourComp> vecEdgeWeight;
};

#endif //_MY_GRAPH_MAX_TREE_H_
