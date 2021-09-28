/*
(c) 2012 Fengtao Fan
*/
//#include "psbmReebGraph.h"
#include "AdjListGraph.h"
#include <string>
#include <stack>
#include <cstring>

//__psbmGraphNode
__psbmGraphNode::__psbmGraphNode() {
    value = 0.f;
    selfIdxInVecNode = -1;
    mappedIdx = -1;
    upNode[0] = upNode[1] = -1;
    //
    downNode[0] = downNode[1] = -1;
    //
    upEdgeLabel[0] = upEdgeLabel[1] = -1;
    //
    downEdgeLabel[0] = downEdgeLabel[1] = -1;
    //
    upValence = 0;
    downValence = 0;

}

void __psbmGraphNode::Clear() {
    value = 0.f;
    selfIdxInVecNode = -1;
    mappedIdx = -1;
    upNode[0] = upNode[1] = -1;
    //
    downNode[0] = downNode[1] = -1;
    //
    upEdgeLabel[0] = upEdgeLabel[1] = -1;
    //
    downEdgeLabel[0] = downEdgeLabel[1] = -1;
    //
    upValence = 0;
    downValence = 0;

}

__psbmGraphNode::__psbmGraphNode(const __psbmGraphNode &rhs) {
    value = rhs.value;
    selfIdxInVecNode = rhs.selfIdxInVecNode;
    mappedIdx = rhs.mappedIdx;
    //
    upNode[0] = rhs.upNode[0];
    upNode[1] = rhs.upNode[1];
    //
    downNode[0] = rhs.downNode[0];
    downNode[1] = rhs.downNode[1];
    //
    upEdgeLabel[0] = rhs.upEdgeLabel[0];
    upEdgeLabel[1] = rhs.upEdgeLabel[1];
    //
    downEdgeLabel[0] = rhs.downEdgeLabel[0];
    downEdgeLabel[1] = rhs.downEdgeLabel[1];
    //
    upValence = rhs.upValence;
    downValence = rhs.downValence;
}

__psbmGraphNode &__psbmGraphNode::operator=(const __psbmGraphNode &rhs) {
    value = rhs.value;
    selfIdxInVecNode = rhs.selfIdxInVecNode;
    mappedIdx = rhs.mappedIdx;
    //
    upNode[0] = rhs.upNode[0];
    upNode[1] = rhs.upNode[1];
    //
    downNode[0] = rhs.downNode[0];
    downNode[1] = rhs.downNode[1];
    //
    upEdgeLabel[0] = rhs.upEdgeLabel[0];
    upEdgeLabel[1] = rhs.upEdgeLabel[1];
    //
    downEdgeLabel[0] = rhs.downEdgeLabel[0];
    downEdgeLabel[1] = rhs.downEdgeLabel[1];
    //
    upValence = rhs.upValence;
    downValence = rhs.downValence;
    return *this;
}

///////__psbmInterimGraph//////////////
__psbmInterimGraph::__psbmInterimGraph() {
}

__psbmInterimGraph::~__psbmInterimGraph() {
    vecNode.clear();
}

__psbmInterimGraph::__psbmInterimGraph(const __psbmInterimGraph &rhs) {
    vecNode.clear();
    vecNode.assign(rhs.vecNode.begin(), rhs.vecNode.end());

}

__psbmInterimGraph &__psbmInterimGraph::operator=(const __psbmInterimGraph &rhs) {
    vecNode.clear();
    vecNode.assign(rhs.vecNode.begin(), rhs.vecNode.end());
    return *this;
}

//////////////////
__psbmInterimGraph::__psbmInterimGraph(psbmReebGraph &rhs,
                                       const int highNode,
                                       const int lowNode) {
    __psbmGraphNode tmpNode;
    int orderCounter = 0;
    unsigned int valence = 0;
    const float highValue = rhs.MeshVertexScalarValue[(*rhs.pVecReebNode)[highNode].nVertexId];
    const float lowValue = rhs.MeshVertexScalarValue[(*rhs.pVecReebNode)[lowNode].nVertexId];
    //
    std::set<int>::iterator sIter;
    std::list<class psbmReebArc *>::iterator hIter;
    //
    std::map<int, int> ReebToInterimNodeMap;
    //
    for (sIter = rhs.criticalNodeIdSet->begin();
         sIter != rhs.criticalNodeIdSet->end();
         sIter++) {
        //tmpNode.Clear();
        tmpNode.value = rhs.MeshVertexScalarValue[(*rhs.pVecReebNode)[*sIter].nVertexId];
        if (tmpNode.value >= lowValue &&
            tmpNode.value <= lowValue) {
            tmpNode.selfIdxInVecNode = orderCounter++;
            tmpNode.mappedIdx = *sIter; // record index in reeb node array
            //tmpNode.value =
            //
            ReebToInterimNodeMap[*sIter] = tmpNode.selfIdxInVecNode;
            //
            vecNode.push_back(tmpNode);
        }

    }
    int ReebNodeId = 0;
    int graphNodeId = 0;
    for (unsigned ivec = 0; ivec < vecNode.size(); ivec++) {
        valence = 0;
        //
        ReebNodeId = vecNode[ivec].mappedIdx;
        //
        valence = (*rhs.pVecReebNode)[ReebNodeId].ptrArcDownId->size();
        vecNode[ivec].downValence = valence;
        if (valence > 0) {
            hIter = (*rhs.pVecReebNode)[ReebNodeId].ptrArcDownId->begin();
            for (unsigned int i = 0; i < valence; hIter++, i++) {
                graphNodeId = rhs.TheOtherEndPoint(ReebNodeId, (*hIter)->clusterLabel);
                if (rhs.MeshVertexScalarValue[(*rhs.pVecReebNode)[graphNodeId].nVertexId] >= lowValue &&
                    rhs.MeshVertexScalarValue[(*rhs.pVecReebNode)[graphNodeId].nVertexId] <= lowValue) {
                    vecNode[ivec].downEdgeLabel[i] = (*hIter)->clusterLabel;
                    vecNode[ivec].downNode[i] = ReebToInterimNodeMap[graphNodeId];
                }
            }
            valence = 0;
        }
        //
        valence = (*rhs.pVecReebNode)[ReebNodeId].ptrArcUpId->size();
        vecNode[ivec].upValence = valence;
        if (valence > 0) {
            hIter = (*rhs.pVecReebNode)[ReebNodeId].ptrArcUpId->begin();
            for (unsigned int i = 0; i < valence; hIter++, i++) {
                graphNodeId = rhs.TheOtherEndPoint(ReebNodeId, (*hIter)->clusterLabel);
                if (rhs.MeshVertexScalarValue[(*rhs.pVecReebNode)[graphNodeId].nVertexId] >= lowValue &&
                    rhs.MeshVertexScalarValue[(*rhs.pVecReebNode)[graphNodeId].nVertexId] <= lowValue) {
                    vecNode[ivec].upEdgeLabel[i] = (*hIter)->clusterLabel;
                    vecNode[ivec].upNode[i] = ReebToInterimNodeMap[graphNodeId];
                }
            }
            valence = 0;
        }
    }
    ReebToInterimNodeMap.clear();
}

__psbmInterimGraph::__psbmInterimGraph(psbmReebGraph &rhs) {
    __psbmGraphNode tmpNode;
    int orderCounter = 0;
    unsigned int valence = 0;
    //
    std::set<int>::iterator sIter;
    std::list<class psbmReebArc *>::iterator hIter;
    //
    std::map<int, int> ReebToInterimNodeMap;
    //Add all critical nodes into vertex set
    for (sIter = rhs.criticalNodeIdSet->begin();
         sIter != rhs.criticalNodeIdSet->end();
         sIter++) {
        //tmpNode.Clear();
        tmpNode.value = rhs.MeshVertexScalarValue[(*rhs.pVecReebNode)[*sIter].nVertexId];
        tmpNode.selfIdxInVecNode = orderCounter++;
        tmpNode.mappedIdx = *sIter; // record index in reeb graph node array
        //
        ReebToInterimNodeMap[*sIter] = tmpNode.selfIdxInVecNode;
        //
        vecNode.push_back(tmpNode);
    }
    int ReebNodeId = 0;
    int graphNodeId = 0;
    // Add all edges into adjacent list for each vertex
    for (unsigned ivec = 0; ivec < vecNode.size(); ivec++) {
        valence = 0;
        //
        ReebNodeId = vecNode[ivec].mappedIdx;
        //
        valence = (*rhs.pVecReebNode)[ReebNodeId].ptrArcDownId->size();
        vecNode[ivec].downValence = valence;
        if (valence > 0) {
            hIter = (*rhs.pVecReebNode)[ReebNodeId].ptrArcDownId->begin();
            for (unsigned int i = 0; i < valence; hIter++, i++) {
                graphNodeId = rhs.TheOtherEndPoint(ReebNodeId, (*hIter)->clusterLabel);
                vecNode[ivec].downEdgeLabel[i] = (*hIter)->clusterLabel;
                vecNode[ivec].downNode[i] = ReebToInterimNodeMap[graphNodeId];
            }
            valence = 0;
        }
        //
        valence = (*rhs.pVecReebNode)[ReebNodeId].ptrArcUpId->size();
        vecNode[ivec].upValence = valence;
        if (valence > 0) {
            hIter = (*rhs.pVecReebNode)[ReebNodeId].ptrArcUpId->begin();
            for (unsigned int i = 0; i < valence; hIter++, i++) {
                graphNodeId = rhs.TheOtherEndPoint(ReebNodeId, (*hIter)->clusterLabel);
                vecNode[ivec].upEdgeLabel[i] = (*hIter)->clusterLabel;
                vecNode[ivec].upNode[i] = ReebToInterimNodeMap[graphNodeId];
            }
            valence = 0;
        }
    }
    ReebToInterimNodeMap.clear();
}

void __psbmInterimGraph::DFS_VISIT_STACK(const int idx,
                                         char *Color,
                                         int *Parent,
                                         int *DiscoverTime,
                                         int *FinishTime) {
    std::stack<int> nodeStack;
    int _globalTime = 0;

    DiscoverTime[idx] = ++_globalTime;
    Color[idx] = 1;

    nodeStack.push(idx);

    bool nodeFinished = true; // to add the possibility that this node is isolated

    std::list<int>::iterator sIter;

    int curNode = -1;

    while (!nodeStack.empty()) {
        curNode = nodeStack.top();
        nodeFinished = true;

        for (int d = 0; d < 2; d++) {
            if (vecNode[curNode].downNode[d] >= 0) {
                int childNode = vecNode[curNode].downNode[d];
                if (Color[childNode] == 0) {
                    Parent[childNode] = idx;
                    DiscoverTime[childNode] = ++_globalTime;
                    Color[childNode] = 1;

                    nodeStack.push(childNode);

                    nodeFinished = false;

                    break;
                }
            }
        }
        for (int u = 0; u < 2; u++) {
            if (vecNode[curNode].upNode[u] >= 0) {
                int childNode = vecNode[curNode].upNode[u];
                if (Color[childNode] == 0) {
                    Parent[childNode] = idx;
                    DiscoverTime[childNode] = ++_globalTime;
                    Color[childNode] = 1;

                    nodeStack.push(childNode);

                    nodeFinished = false;

                    break;
                }
            }
        }
        if (nodeFinished) {// all incident nodes are discovered
            // this node is finised
            Color[curNode] = 2;
            FinishTime[curNode] = ++_globalTime;

            //
            nodeStack.pop();
        }

    }

}

void __psbmInterimGraph::Clear() {
    vecNode.clear();
    return;
}

void __psbmInterimGraph::DFS_Reachable_Subgraph(
        __psbmInterimGraph &retGraph,
        const int SourceNodeIdx,
        const float HighValue,
        const float LowValue) {// SourceNodeIdx is an index of AdjList graph
    __psbmInterimGraph tmpSubgraph;
    // clear return graph
    retGraph.vecNode.clear();
    // initialize the color to zero (White)
    char *insideColor = new char[vecNode.size()];
    memset(insideColor, 0, vecNode.size() * sizeof(char));
    //
    __psbmGraphNode tmpNode;
    //
    std::map<int, int> oldMappedToNew; // mapping from old to new
    std::map<int, int>::iterator mIter;
    //
    std::stack<int> nodeStack; // used for DFS travel
    int _globalTime = 0;
    // add source vertex into subgraph
    tmpNode.selfIdxInVecNode = _globalTime;
    tmpNode.mappedIdx = SourceNodeIdx; // record index in Adjlist graph
    tmpNode.value = vecNode[SourceNodeIdx].value;
    //
    tmpSubgraph.vecNode.push_back(tmpNode);
    //
    oldMappedToNew[SourceNodeIdx] = _globalTime++;

    // set color to gray
    insideColor[SourceNodeIdx] = 1;
    // push it into the stack
    nodeStack.push(SourceNodeIdx);

    bool nodeFinished = true; // to add the possibility that this node is isolated

    std::list<int>::iterator sIter;

    int curNode = -1;

    while (!nodeStack.empty()) {
        curNode = nodeStack.top();
        nodeFinished = true;

        for (int d = 0; d < vecNode[curNode].downValence; d++) {
            int childNode = vecNode[curNode].downNode[d];
            if (vecNode[childNode].value >= LowValue &&
                vecNode[childNode].value <= HighValue) {

                if (insideColor[childNode] == 0) {
                    // Add this vertex into subgraph
                    tmpNode.selfIdxInVecNode = _globalTime;
                    tmpNode.mappedIdx = childNode;
                    tmpNode.value = vecNode[childNode].value;
                    tmpSubgraph.vecNode.push_back(tmpNode);
                    //
                    oldMappedToNew[childNode] = _globalTime++;
                    insideColor[childNode] = 1;

                    nodeStack.push(childNode);

                    nodeFinished = false;

                    break;
                }
            }
        }
        for (int u = 0; u < vecNode[curNode].upValence; u++) {
            int childNode = vecNode[curNode].upNode[u];
            if (vecNode[childNode].value >= LowValue &&
                vecNode[childNode].value <= HighValue) {

                if (insideColor[childNode] == 0) {
                    // Add this vertex into subgraph
                    tmpNode.selfIdxInVecNode = _globalTime;
                    tmpNode.mappedIdx = childNode;
                    tmpNode.value = vecNode[childNode].value;
                    tmpSubgraph.vecNode.push_back(tmpNode);
                    //

                    oldMappedToNew[childNode] = _globalTime++;
                    insideColor[childNode] = 1;

                    nodeStack.push(childNode);

                    nodeFinished = false;

                    break;
                }
            }

        }
        if (nodeFinished) {// all incident nodes are discovered
            // this node is finised
            insideColor[curNode] = 2;
            //
            nodeStack.pop();
        }

    }
    // assemble the edges
    for (unsigned int i = 0; i < tmpSubgraph.vecNode.size(); i++) {
        curNode = tmpSubgraph.vecNode[i].mappedIdx;
        int validNodeIdx = 0;
        for (int d = 0; d < vecNode[curNode].downValence; d++) {
            if (insideColor[vecNode[curNode].downNode[d]]) {// reachable from source
                tmpSubgraph.vecNode[i].downNode[validNodeIdx] = oldMappedToNew[vecNode[curNode].downNode[d]];
                tmpSubgraph.vecNode[i].downEdgeLabel[validNodeIdx] = vecNode[curNode].downEdgeLabel[d];
                validNodeIdx++;
            }
        }
        tmpSubgraph.vecNode[i].downValence = validNodeIdx;
        //
        validNodeIdx = 0;
        for (int u = 0; u < vecNode[curNode].upValence; u++) {
            if (insideColor[vecNode[curNode].upNode[u]]) {// reachable from source
                tmpSubgraph.vecNode[i].upNode[validNodeIdx] = oldMappedToNew[vecNode[curNode].upNode[u]];
                tmpSubgraph.vecNode[i].upEdgeLabel[validNodeIdx] = vecNode[curNode].upEdgeLabel[u];
                validNodeIdx++;
            }
        }
        tmpSubgraph.vecNode[i].upValence = validNodeIdx;
        //
        //record the index in reeb graph node array
        tmpSubgraph.vecNode[i].mappedIdx = vecNode[curNode].mappedIdx;

    }
    // simplified the subgraph by cutting non-loop branches
    // selfIdxInVecNode = -1 means this vertex is deleted
    while (1) {// THERE IS AT LEAST ONE LOOP THERE
        bool exitLoop = true;
        for (unsigned int i = 0; i < tmpSubgraph.vecNode.size(); i++) {
            if (tmpSubgraph.vecNode[i].selfIdxInVecNode >= 0) {
                // sum the degree of this vertex
                int degree = tmpSubgraph.vecNode[i].downValence + tmpSubgraph.vecNode[i].upValence;
                //as each vertex is reachable from source
                // degree has at least 1
                if (degree == 1) {// delete this vertex and trace out
                    int leafNodeIdx = i;
                    int OtherEndPointIdx = 0;
                    int OtherDegree = 0;
                    while (1) {
                        //compute the vertex connected to this vertex
                        if (tmpSubgraph.vecNode[leafNodeIdx].downValence == 1) {
                            OtherEndPointIdx = tmpSubgraph.vecNode[leafNodeIdx].downNode[0];
                        } else {
                            OtherEndPointIdx = tmpSubgraph.vecNode[leafNodeIdx].upNode[0];
                        }

                        // delete leaf node
                        tmpSubgraph.vecNode[leafNodeIdx].selfIdxInVecNode = -1;
                        //
                        OtherDegree = tmpSubgraph.vecNode[OtherEndPointIdx].downValence +
                                      tmpSubgraph.vecNode[OtherEndPointIdx].upValence;
                        //compute the degree of the other vertex  and delete the edge
                        int vPos = 0;
                        for (vPos = 0; vPos < tmpSubgraph.vecNode[OtherEndPointIdx].downValence; vPos++) {
                            if (tmpSubgraph.vecNode[OtherEndPointIdx].downNode[vPos] == leafNodeIdx)
                                break;
                        }
                        if (vPos < tmpSubgraph.vecNode[OtherEndPointIdx].downValence) {
                            if (vPos + 1 == tmpSubgraph.vecNode[OtherEndPointIdx].downValence) {
                                tmpSubgraph.vecNode[OtherEndPointIdx].downNode[vPos] = -1;
                                tmpSubgraph.vecNode[OtherEndPointIdx].downEdgeLabel[vPos] = -1;
                            } else {
                                tmpSubgraph.vecNode[OtherEndPointIdx].downNode[0] = tmpSubgraph.vecNode[OtherEndPointIdx].downNode[1];
                                tmpSubgraph.vecNode[OtherEndPointIdx].downEdgeLabel[0] = tmpSubgraph.vecNode[OtherEndPointIdx].downEdgeLabel[1];
                            }
                            tmpSubgraph.vecNode[OtherEndPointIdx].downValence--;
                        } else {
                            for (vPos = 0; vPos < tmpSubgraph.vecNode[OtherEndPointIdx].upValence; vPos++)
                                if (tmpSubgraph.vecNode[OtherEndPointIdx].upNode[vPos] == leafNodeIdx)
                                    break;
                            if (vPos + 1 == tmpSubgraph.vecNode[OtherEndPointIdx].upValence) {
                                tmpSubgraph.vecNode[OtherEndPointIdx].upNode[vPos] = -1;
                                tmpSubgraph.vecNode[OtherEndPointIdx].upEdgeLabel[vPos] = -1;
                            } else {
                                tmpSubgraph.vecNode[OtherEndPointIdx].upNode[0] = tmpSubgraph.vecNode[OtherEndPointIdx].upNode[1];
                                tmpSubgraph.vecNode[OtherEndPointIdx].upEdgeLabel[0] = tmpSubgraph.vecNode[OtherEndPointIdx].upEdgeLabel[1];
                            }
                            tmpSubgraph.vecNode[OtherEndPointIdx].upValence--;
                        }
                        OtherDegree--;
                        //
                        if (OtherDegree != 1)
                            break;
                        else
                            leafNodeIdx = OtherEndPointIdx;
                    }
                    exitLoop = false;
                }// if == 1
            }// IF >= 0

        }// FOR I
        if (exitLoop)
            break;
    }// while(1)
//
    tmpNode.Clear();
    oldMappedToNew.clear();
    //
    int nNodeCounter = 0;
    for (unsigned int i = 0; i < tmpSubgraph.vecNode.size(); i++) {
        if (tmpSubgraph.vecNode[i].selfIdxInVecNode >= 0) {
            tmpNode.selfIdxInVecNode = nNodeCounter;
            //record the index in non-simplified temp graph
            tmpNode.mappedIdx = tmpSubgraph.vecNode[i].selfIdxInVecNode;
            tmpNode.value = tmpSubgraph.vecNode[i].value;
            tmpNode.downValence = tmpSubgraph.vecNode[i].downValence;
            tmpNode.upValence = tmpSubgraph.vecNode[i].upValence;
            //
            retGraph.vecNode.push_back(tmpNode);
            //
            oldMappedToNew[tmpSubgraph.vecNode[i].selfIdxInVecNode] = nNodeCounter++;
        }
    }
    // assemble edges
    for (unsigned int i = 0; i < retGraph.vecNode.size(); i++) {
        int tmpGraphNodeIdx = retGraph.vecNode[i].mappedIdx;
        for (int v = 0; v < retGraph.vecNode[i].downValence; v++) {
            // reachable from target
            retGraph.vecNode[i].downNode[v] = oldMappedToNew[tmpSubgraph.vecNode[tmpGraphNodeIdx].downNode[v]];
            retGraph.vecNode[i].downEdgeLabel[v] = tmpSubgraph.vecNode[tmpGraphNodeIdx].downEdgeLabel[v];
        }
        for (int v = 0; v < retGraph.vecNode[i].upValence; v++) {
            // reachable from target
            retGraph.vecNode[i].upNode[v] = oldMappedToNew[tmpSubgraph.vecNode[tmpGraphNodeIdx].upNode[v]];
            retGraph.vecNode[i].upEdgeLabel[v] = tmpSubgraph.vecNode[tmpGraphNodeIdx].upEdgeLabel[v];
        }
        retGraph.vecNode[i].mappedIdx = tmpSubgraph.vecNode[tmpGraphNodeIdx].mappedIdx;
    }
    //delete temp graph
    tmpSubgraph.Clear();
    return;
}

bool __psbmInterimGraph::DFS_Reachable_STACK(
        psbmReebGraph &rbGrpah,
        const int SourceNodeIdx,
        const int TargetNodeIdx,
        char *Color,
        const float HighValue,
        const float LowValue) {
    std::stack<int> nodeStack;

    Color[SourceNodeIdx] = 1;

    nodeStack.push(SourceNodeIdx);

    bool nodeFinished = true; // to add the possibility that this node is isolated
    bool ret = false;

    std::list<int>::iterator sIter;

    int curNode = -1;

    while (!nodeStack.empty()) {
        curNode = nodeStack.top();
        nodeFinished = true;

        for (int d = 0; d < 2; d++) {
            if (vecNode[curNode].downNode[d] >= 0) {
                int childNode = vecNode[curNode].downNode[d];
                if (vecNode[childNode].value >= LowValue &&
                    vecNode[childNode].value <= HighValue) {
                    if (Color[childNode] == 0) {
                        Color[childNode] = 1;
                        if (childNode == TargetNodeIdx)
                            ret = true;
                        nodeStack.push(childNode);

                        nodeFinished = false;

                        break;
                    }
                }
            }
        }
        for (int u = 0; u < 2; u++) {
            if (vecNode[curNode].upNode[u] >= 0) {
                int childNode = vecNode[curNode].upNode[u];
                if (vecNode[childNode].value >= LowValue &&
                    vecNode[childNode].value <= HighValue) {
                    if (Color[childNode] == 0) {
                        Color[childNode] = 1;
                        if (childNode == TargetNodeIdx)
                            ret = true;
                        nodeStack.push(childNode);

                        nodeFinished = false;

                        break;
                    }
                }
            }
        }
        if (nodeFinished) {// all incident nodes are discovered
            // this node is finised
            Color[curNode] = 2;

            //
            nodeStack.pop();
        }
        if (ret) {
            while (!nodeStack.empty())
                nodeStack.pop();
        }
    }
    return ret;
}

bool __psbmInterimGraph::DFS_UP_Reachable_STACK(
        // psbmReebGraph &rbGraph,
        const int SourceNodeIdx,
        const int TargetNodeIdx,
        const float HighValue,
        const float LowValue) {
    std::stack<int> nodeStack;
    //
    char *Color = new char[vecNode.size()];
    memset(Color, 0, vecNode.size() * sizeof(char));
    //
    Color[SourceNodeIdx] = 1;

    nodeStack.push(SourceNodeIdx);

    bool nodeFinished = true; // to add the possibility that this node is isolated
    bool ret = false;

    std::list<int>::iterator sIter;

    int curNode = -1;

    float curLowValue = 0.f;

    while (!nodeStack.empty()) {
        curNode = nodeStack.top();
        nodeFinished = true;

        curLowValue = vecNode[curNode].value;

        for (int d = 0; d < 2; d++) {
            if (vecNode[curNode].downNode[d] >= 0) {
                int childNode = vecNode[curNode].downNode[d];
                if (vecNode[childNode].value >= curLowValue &&
                    vecNode[childNode].value <= HighValue) {
                    if (Color[childNode] == 0) {
                        Color[childNode] = 1;
                        if (childNode == TargetNodeIdx)
                            ret = true;
                        nodeStack.push(childNode);

                        nodeFinished = false;

                        break;
                    }
                }
            }
        }
        for (int u = 0; u < 2; u++) {
            if (vecNode[curNode].upNode[u] >= 0) {
                int childNode = vecNode[curNode].upNode[u];
                if (vecNode[childNode].value >= curLowValue &&
                    vecNode[childNode].value <= HighValue) {
                    if (Color[childNode] == 0) {
                        Color[childNode] = 1;
                        if (childNode == TargetNodeIdx)
                            ret = true;
                        nodeStack.push(childNode);

                        nodeFinished = false;

                        break;
                    }
                }
            }
        }
        if (nodeFinished) {// all incident nodes are discovered
            // this node is finised
            Color[curNode] = 2;

            //
            nodeStack.pop();
        }
        if (ret) {
            while (!nodeStack.empty())
                nodeStack.pop();
        }
    }
    return ret;
}
