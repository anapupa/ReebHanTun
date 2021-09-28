/*
(c) 2012 Fengtao Fan
*/
#pragma once

#include "psbmReebGraph.h"

class __psbmGraphNode {
public:
    __psbmGraphNode();

    void Clear();

    __psbmGraphNode(const __psbmGraphNode &rhs);

    __psbmGraphNode &operator=(const __psbmGraphNode &rhs);

    /////////
    int selfIdxInVecNode;    // index in the graph node array
    int mappedIdx;            // index mapped into another array
    float value;            // vertex value
    //int color;
    //int parent;
    //int DiscoverTime;
    //int FinishTime;

    int upNode[2];            // adjacent vertices with higher value
    int downNode[2];        // adjacent vertices with lower value
    int upEdgeLabel[2];        // edge label connecting to higher vertices
    int downEdgeLabel[2];    // edge label connecting to lower vertices
    int upValence;            // # of higher vertices
    int downValence;        // # of lower vertices
};

class __psbmInterimGraph {
public:
    __psbmInterimGraph();

    ~__psbmInterimGraph();

    __psbmInterimGraph(const __psbmInterimGraph &rhs);

    __psbmInterimGraph &operator=(const __psbmInterimGraph &rhs);

    __psbmInterimGraph(psbmReebGraph &rhs);

    __psbmInterimGraph(psbmReebGraph &rhs,
                       const int highNode,
                       const int lowNode);

/////////////////////////////////
    void Clear();

    void DFS_VISIT_STACK(const int SourceNodeIdx,
                         char *Color,
                         int *Parent,
                         int *DiscoverTime,
                         int *FinishTime);

    void DFS_Reachable_Subgraph(
            __psbmInterimGraph &retGraph,
            const int SourceNodeIdx,
            const float HighValue,
            const float LowValue);

    bool DFS_Reachable_STACK(psbmReebGraph &rhs,
                             const int SourceNodeIdx,
                             const int TargetNodeIdx,
                             char *Color,
                             const float HighValue,
                             const float LowValue);

    bool DFS_UP_Reachable_STACK(//psbmReebGraph &rhs,
            const int SourceNodeIdx,
            const int TargetNodeIdx,
            // char *Color,
            const float HighValue,
            const float LowValue);


public:
    std::vector<__psbmGraphNode> vecNode; // vertex array

};
//
//class __psbmGraphNode
//{
//public:
//	__psbmGraphNode();
//	void Clear();
//	__psbmGraphNode(const __psbmGraphNode& rhs);
//	__psbmGraphNode& operator=(const __psbmGraphNode& rhs);
//	/////////
//	int selfIdxInVecNode;
//	int mappedIdx;
//	float value;
//	//int color;
//	//int parent;
//	//int DiscoverTime;
//	//int FinishTime;
//	int upNode[2];
//	int downNode[2];
//	int upEdgeLabel[2];
//	int downEdgeLabel[2];
//};
//
//class __psbmInterimGraph
//{
//public:
//	__psbmInterimGraph()
//	{
//	}
//	~__psbmInterimGraph()
//	{
//		vecNode.clear();
//	}
//	__psbmInterimGraph(const __psbmInterimGraph& rhs)
//	{
//		vecNode.clear();
//		vecNode.assign(rhs.vecNode.begin(), rhs.vecNode.end());
//
//	}
//	__psbmInterimGraph& operator=(const __psbmInterimGraph& rhs)
//	{
//		vecNode.clear();
//		vecNode.assign(rhs.vecNode.begin(), rhs.vecNode.end());
//
//	}
//	__psbmInterimGraph(psbmReebGraph &rhs, 
//						const int highNode,
//						const int lowNode);
///////////////////////////////////
//	void DFS_VISIT_STACK(const int SourceNodeIdx,
//						 char *Color,
//						 int *Parent,
//						 int *DiscoverTime,
//						 int *FinishTime);
//	__psbmInterimGraph& DFS_Reachable_Subgraph(
//							const int SourceNodeIdx,
//							char* Color,
//							const float HighValue,
//							const float LowValue);
//	bool DFS_Reachable_STACK(psbmReebGraph &rhs, 
//							const int SourceNodeIdx,
//							 const int TargetNodeIdx,
//							 char *Color,
//							 const float HighValue,
//							 const float LowValue);
//	bool DFS_UP_Reachable_STACK(psbmReebGraph &rhs, 
//							const int SourceNodeIdx,
//							 const int TargetNodeIdx,
//							 char *Color,
//							 const float HighValue,
//							 const float LowValue);
//
//
//public:
//	std::vector<__psbmGraphNode> vecNode;
//
//};
/*
//#pragma once
//
//class __psbmGraphNode
//{
//public:
//	__psbmGraphNode();
//	void Clear();
//	__psbmGraphNode(const __psbmGraphNode& rhs);
//	__psbmGraphNode& operator=(const __psbmGraphNode& rhs);
//	/////////
//	int selfIdxInVecNode;
//	int mappedIdx;
//	float value;
//	//int color;
//	//int parent;
//	//int DiscoverTime;
//	//int FinishTime;
//	int upNode[2];
//	int downNode[2];
//	int upEdgeLabel[2];
//	int downEdgeLabel[2];
//};
//
//class __psbmInterimGraph
//{
//public:
//	__psbmInterimGraph()
//	{
//	}
//	~__psbmInterimGraph()
//	{
//		vecNode.clear();
//	}
//	__psbmInterimGraph(const __psbmInterimGraph& rhs)
//	{
//		vecNode.clear();
//		vecNode.assign(rhs.vecNode.begin(), rhs.vecNode.end());
//
//	}
//	__psbmInterimGraph& operator=(const __psbmInterimGraph& rhs)
//	{
//		vecNode.clear();
//		vecNode.assign(rhs.vecNode.begin(), rhs.vecNode.end());
//
//	}
//	__psbmInterimGraph(psbmReebGraph &rhs, 
//						const int highNode,
//						const int lowNode);
///////////////////////////////////
//	void DFS_VISIT_STACK(const int SourceNodeIdx,
//						 char *Color,
//						 int *Parent,
//						 int *DiscoverTime,
//						 int *FinishTime);
//	__psbmInterimGraph& DFS_Reachable_Subgraph(
//							const int SourceNodeIdx,
//							char* Color,
//							const float HighValue,
//							const float LowValue);
//	bool DFS_Reachable_STACK(psbmReebGraph &rhs, 
//							const int SourceNodeIdx,
//							 const int TargetNodeIdx,
//							 char *Color,
//							 const float HighValue,
//							 const float LowValue);
//	bool DFS_UP_Reachable_STACK(psbmReebGraph &rhs, 
//							const int SourceNodeIdx,
//							 const int TargetNodeIdx,
//							 char *Color,
//							 const float HighValue,
//							 const float LowValue);
//
//
//public:
//	std::vector<__psbmGraphNode> vecNode;
//
//};
*/