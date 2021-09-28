/*
(c) 2012 Fengtao Fan
*/
#include "MapLoopsBackToMeshLevelSetAndArc.h"
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

void MapLoopsBackToMeshLevelSetAndArc::FindEdgeConnectingTwoVertices(const int src_x, const int dist_y, const int eid_y,
                                                                     int &outEdge) {
    int edge_u = 0;
    int edge_v = 0;
    int triangle_id = TriangleAdjacentToOneEdgesOneVertex(eid_y, src_x);
    //
    edge_u = (*inMeshPtr).vecTriangle[triangle_id].e01;
    if ((*inMeshPtr).vecTriangle[triangle_id].e01 == eid_y)
        edge_u = (*inMeshPtr).vecTriangle[triangle_id].e02;
    //
    edge_v = (*inMeshPtr).vecTriangle[triangle_id].e01 +
             (*inMeshPtr).vecTriangle[triangle_id].e02 +
             (*inMeshPtr).vecTriangle[triangle_id].e12 -
             (eid_y + edge_u);
    //
    if ((*inMeshPtr).vecEdge[edge_u].v0 == dist_y ||
        (*inMeshPtr).vecEdge[edge_u].v1 == dist_y) {
        outEdge = edge_u;
    } else {
        outEdge = edge_v;
    }
    return;
}

void MapLoopsBackToMeshLevelSetAndArc::ComputeLevelCycle() {
    // preq: level cycle pointType is initialized
//	//////////////////////////
    /*
    The level set can be treated as a closed triangle strip
    So one can trace the triangle strip
    */
    /*
       boost::unordered_set<int> triangleDegree;
       std::vector<int> triangleOrder ;
       for (unsigned int i = 0; i < _levelCyclePointTypePtr->size() - 1; i++)
       {
           int triangle_id = 0;
           int edge_a = (*_levelCyclePointTypePtr)[i].first;
           int edge_b = (*_levelCyclePointTypePtr)[i + 1].first;
           //if (i == 0)
           //{
           //	triangle_id = TriangleAdjacentToOneEdgesOneVertex(edge_b, edge_a);
           //}
           //else
           //{
           //	if (i == _arcCyclePointTypePtr->size() - 2)
           //	{
           //		triangle_id = TriangleAdjacentToOneEdgesOneVertex(edge_a, edge_b);
           //	}
           //	else
           //	{
           //		triangle_id = TriangleAdjacentToTwoEdges(edge_a, edge_b);
           //	}
           //}
           //
           triangle_id = TriangleAdjacentToTwoEdges(edge_a, edge_b);
           //
           boost::unordered_set<int>::iterator sIter = triangleDegree.find(triangle_id);
           if (sIter == triangleDegree.end())
           {
               triangleDegree.insert(triangle_id);
               triangleOrder.push_back(triangle_id);
           }
           else
           {//
               std::cout << "NOT a triangle trip" << std::endl;
               exit(9);
           }
           //
       }
       //
       _simpGraph triangleStripGraph;
       boost::unordered_map<int, int> verToGraphNode;
       boost::unordered_map<int, int>::iterator mIter;
       int nGrapNodeCounter = 0;
       for (unsigned int i = 0; i < triangleOrder.size(); i++)
       {
           int node[3] = {	(*inMeshPtr).vecTriangle[triangleOrder[i]].v0,
                           (*inMeshPtr).vecTriangle[triangleOrder[i]].v1,
                           (*inMeshPtr).vecTriangle[triangleOrder[i]].v2};
           int gNode[3] = {0};
           for (int j = 0; j < 3; j++)
           {
               mIter = verToGraphNode.find(node[j]);
               if (mIter == verToGraphNode.end())
               {
                   triangleStripGraph.AddNode(nGrapNodeCounter, node[j]);
                   gNode[j] = nGrapNodeCounter;
                   verToGraphNode[node[j]] = nGrapNodeCounter;
                   nGrapNodeCounter++;
               }
               else
               {
                   gNode[j] = mIter->second;
               }
           }
           triangleStripGraph.AddEdge(gNode[0], gNode[1]); // since graph node use set to store adjacent list
           triangleStripGraph.AddEdge(gNode[0], gNode[2]); // no need to check the edge exists or not
           triangleStripGraph.AddEdge(gNode[1], gNode[2]);
       }
       //
       for (unsigned int i = 0; i < _levelCyclePointTypePtr->size() - 1; i++)
       {
           int edge_idx = (*_levelCyclePointTypePtr)[i].first;
           triangleStripGraph.RemoveEdge( verToGraphNode[(*inMeshPtr).vecEdge[edge_idx].v0],  verToGraphNode[(*inMeshPtr).vecEdge[edge_idx].v1]);
       }
       //there exist two disjont two components
       std::vector<bool> colors(triangleStripGraph.vecNode.size(), false);
       std::vector<int> parents(triangleStripGraph.vecNode.size(), -1);
       int components = 0;
       int currentNode = -1;
       int oppositeNode = -1;
       for (unsigned int i = 0; i < colors.size(); i++)
       {
           if (!colors[i])
           {// visit the components
               std::queue<int> Q;
               Q.push(i);
               colors[i] = true;
               parents[i] = -1;
               while (!Q.empty())
               {
                   int node = Q.front();
                   Q.pop();
                   for (std::set<int>::iterator sIter = triangleStripGraph.vecNode[node].adjList.begin();
                       sIter != triangleStripGraph.vecNode[node].adjList.end(); sIter++)
                   {
                       if (!colors[*sIter])
                       {
                           colors[*sIter] = true;
                           Q.push(*sIter);
                           parents[*sIter] = node;
                       }
                       else
                       {
                           if (parents[node] != *sIter)
                           {// non-tree edge
                               currentNode = node;
                               oppositeNode = *sIter;
                           }
                       }
                   }
               }
               components++;
           }
       }
       if (components != 2)
       {
           std::cout << "wrong method again" << std::endl;
           exit(4);
       }
       outLevelCycle_VertexOnMesh.reserve(triangleStripGraph.vecNode.size());
       outLevelCycle_EdgeOnMesh.reserve(triangleStripGraph.vecNode.size());
       //

       //for (int i = 0; i < triangleStripGraph.vecNode.size(); i++)
       //{
       //	if (parents[i] < 0 && triangleStripGraph.vecNode[i].adjList.size() > 1)
       //	{
       //		currentNode = i;
       //		break;
       //	}
       //}
       if (currentNode < 0 || oppositeNode < 0)
       {
           std::cout << "CANNOT find intial node for level cycle " << std::endl;
           exit(4);
       }

       //for (std::set<int>::iterator sIter = triangleStripGraph.vecNode[currentNode].adjList.begin();
       //	sIter != triangleStripGraph.vecNode[currentNode].adjList.end(); sIter++)
       //{
       //	if (currentNode != parents[*sIter])
       //	{
       //		oppositeNode = *sIter;
       //		break;
       //	}
       //}
       //
       std::vector<int> leftPath;
       int iterNode = currentNode;
       while (iterNode >= 0)
       {
           leftPath.push_back(triangleStripGraph.vecNode[iterNode].color);
           iterNode = parents[iterNode];
       }
       for (std::vector<int>::reverse_iterator rvIter = leftPath.rbegin();
           rvIter != leftPath.rend(); rvIter++)
       {
           outLevelCycle_VertexOnMesh.push_back(*rvIter);
       }
       //
       iterNode = oppositeNode;
       while (iterNode >= 0)
       {
           outLevelCycle_VertexOnMesh.push_back(triangleStripGraph.vecNode[iterNode].color);
           iterNode = parents[iterNode];
       }
       if (outLevelCycle_VertexOnMesh.back() != leftPath.back())
       {
           std::cout << "in different compoents" << std::endl;
           exit(4);
       }
       //
       for (unsigned int i = 0; i < outLevelCycle_VertexOnMesh.size() - 1; i++)
       {
           int curNode = outLevelCycle_VertexOnMesh[i];
           int oppNode = outLevelCycle_VertexOnMesh[i + 1];
           for (std::vector<int>::iterator vIter = (*inMeshPtr).vecVertex[curNode].adjEdges.begin();
               vIter != (*inMeshPtr).vecVertex[curNode].adjEdges.end(); vIter++)
           {
               if ((*inMeshPtr).vecEdge[*vIter].v0 == oppNode ||
                   (*inMeshPtr).vecEdge[*vIter].v1 == oppNode)
               {
                   outLevelCycle_EdgeOnMesh.push_back(*vIter);
                   break;
               }
           }
       }
       //
       //
       return;
   */

//	// critical point
    //std::set<int> processed_vertex;
    //std::set<int>::iterator findIter;
    //
    int pot_vertex = 0;
    int pot_edge_id = 0;
    int current_edge_id = 0;
    int current_active_vertex = 0;
    int shared_vertex = 0;
    /*
    Initializing current active vertex
    */
    if ((*_levelCyclePointTypePtr)[0].second) {// it is an edge
        current_edge_id = (*_levelCyclePointTypePtr)[0].first;
        current_active_vertex = (*inMeshPtr).vecEdge[current_edge_id].v0;
        //
        //processed_vertex.insert(current_active_vertex);
        //
        outLevelCycle_VertexOnMesh.push_back(current_active_vertex);
    } else {// it is only a vertex
        current_active_vertex = (*_levelCyclePointTypePtr)[0].first;
        //
        //processed_vertex.insert(current_active_vertex);
        //
        outLevelCycle_VertexOnMesh.push_back(current_active_vertex);
    }
    // boundary check is done YET !!!!!!!!!!!!!!!
    for (unsigned int i = 1; i < _levelCyclePointTypePtr->size() - 1; i++) {
        /*
        Invariant : current_active_vertex is always on current edge if current array element is an edge
                    or it is the current vertex if current array element is a vertex
        */
        if ((*_levelCyclePointTypePtr)[i].second) {// this is an edge type vertex
            current_edge_id = (*_levelCyclePointTypePtr)[i].first;
            if (inMeshPtr->vecEdge[current_edge_id].v0 != current_active_vertex &&
                inMeshPtr->vecEdge[current_edge_id].v1 !=
                current_active_vertex) {// if current active vertex is on this edge, then do nothing
                // else add one more edge and update current_active_vertex
                // since the levele set cross the edge only once, so there must exist at least one vertex on the edge not visited yet
                if ((*_levelCyclePointTypePtr)[i - 1].second) {// previous one is an edge
                    shared_vertex = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i - 1].first].v0 +
                                    inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i - 1].first].v1 -
                                    current_active_vertex;
                    pot_vertex = inMeshPtr->vecEdge[current_edge_id].v0;
                    if (shared_vertex == pot_vertex)
                        pot_vertex = inMeshPtr->vecEdge[current_edge_id].v1;
                    //
                    FindEdgeConnectingTwoVertices(current_active_vertex, pot_vertex, current_edge_id, pot_edge_id);
                    //
                    current_active_vertex = pot_vertex;
                    //
                    outLevelCycle_VertexOnMesh.push_back(current_active_vertex);
                    outLevelCycle_EdgeOnMesh.push_back(pot_edge_id);
                    //
                    //processed_vertex.insert(current_active_vertex);
                } else {
                    std::cout << "error in levelset mapping, exist an vertex?" << std::endl;
                    exit(2);
                }
            }
        } else {// it is a vertex
            std::cout << "error in levelset mapping, exist an vertex?" << std::endl;
            exit(2);
            pot_vertex = (*_levelCyclePointTypePtr)[i].first;
            // check if this vertex is there or not
            // if it is there, ignore it
            // otherwise, push it
            //
            pot_edge_id = -1;
            //if (processed_vertex.find(pot_vertex) == processed_vertex.end())
            {
                for (std::vector<int>::iterator vIter = inMeshPtr->vecVertex[current_active_vertex].adjEdges.begin();
                     vIter != inMeshPtr->vecVertex[current_active_vertex].adjEdges.end(); vIter++) {
                    if (inMeshPtr->vecEdge[*vIter].v0 == pot_vertex ||
                        inMeshPtr->vecEdge[*vIter].v1 == pot_vertex) {
                        pot_edge_id = *vIter;
                        break;
                    }
                }
                if (pot_edge_id > -1) {
                    outLevelCycle_VertexOnMesh.push_back(pot_vertex);
                    outLevelCycle_EdgeOnMesh.push_back(pot_edge_id);
                } else {
                    std::cout << "No edge connecting two vertices " << std::endl;
                    exit(0);
                }
                //
                current_active_vertex = pot_vertex;
                //
                //processed_vertex.insert(current_active_vertex);
            }
        }// else it is a vertex
    }
    // Handle the boundary case
    pot_vertex = outLevelCycle_VertexOnMesh.front();
    if ((*_levelCyclePointTypePtr).back().second) {
        if (pot_vertex != current_active_vertex) {
            current_edge_id = (*_levelCyclePointTypePtr).back().first;
            if (inMeshPtr->vecEdge[current_edge_id].v0 == current_active_vertex ||
                inMeshPtr->vecEdge[current_edge_id].v1 == current_active_vertex) {
                outLevelCycle_EdgeOnMesh.push_back(current_edge_id);
                outLevelCycle_VertexOnMesh.push_back(pot_vertex);
            } else {
                FindEdgeConnectingTwoVertices(current_active_vertex, pot_vertex, current_edge_id, pot_edge_id);
                //
                outLevelCycle_EdgeOnMesh.push_back(pot_edge_id);
                outLevelCycle_VertexOnMesh.push_back(pot_vertex);
            }
        }
    } else {
        std::cout << "error in levelset mapping, exist an vertex?" << std::endl;
        exit(2);
        pot_edge_id = -1;
        for (std::vector<int>::iterator vIter = inMeshPtr->vecVertex[current_active_vertex].adjEdges.begin();
             vIter != inMeshPtr->vecVertex[current_active_vertex].adjEdges.end(); vIter++) {
            if (inMeshPtr->vecEdge[*vIter].v0 == pot_vertex ||
                inMeshPtr->vecEdge[*vIter].v1 == pot_vertex) {
                pot_edge_id = *vIter;
                break;
            }
        }
        if (pot_edge_id > -1) {
            outLevelCycle_VertexOnMesh.push_back(pot_vertex);
            outLevelCycle_EdgeOnMesh.push_back(pot_edge_id);
        } else {
            std::cout << "No edge connecting two vertices " << std::endl;
            exit(0);
        }
    }
    return;
    //
    //// Trace a vertex-edge path on the mesh which is homotopic to the input cycle on mesh
    //// the second part is the pair to its parent and the edge connecting it to its parent
    //std::map<int, std::pair<int,int>> processed_vertex;
    //int current_edge_a = 0;
    //int current_edge_b = 0;
    //int current_edge_c = 0;
    //int opposite_vertex = 0;
    //int tracing_termination_pt = 0;
    //int current_active_node = 0;
    //if ((*_levelCyclePointTypePtr)[0].second)
    //{// it is an edge
    //	current_edge_a = (*_levelCyclePointTypePtr)[0].first;
    //	tracing_termination_pt = (*inMeshPtr).vecEdge[current_edge_a].v0;
    //	processed_vertex[(*inMeshPtr).vecEdge[current_edge_a].v0] = std::pair<int,int>(-1, -1);
    //	processed_vertex[(*inMeshPtr).vecEdge[current_edge_a].v1] = std::pair<int,int>((*inMeshPtr).vecEdge[current_edge_a].v0, current_edge_a);
    //	current_active_node = (*inMeshPtr).vecEdge[current_edge_a].v1;
    //}
    //else
    //{// it is only a vertex
    //	current_edge_a = (*_levelCyclePointTypePtr)[0].first;
    //	tracing_termination_pt = current_edge_a;
    //	processed_vertex[current_edge_a] = std::pair<int,int>(-1, -1);
    //	current_active_node = current_edge_a;
    //}
    //////////////////////////////
    //for (unsigned int i = 0; i < _levelCyclePointTypePtr->size() - 2; i++)
    //{
    //	// find triangle conntaining the edge [i, i+1]
    //	int triangle_or_edge_id = 0;

    //	if ((*_levelCyclePointTypePtr)[i].second)
    //	{
    //		if ((*_levelCyclePointTypePtr)[i+1].second)
    //		{// both points are on edges
    //		 // there exists an triangle containing them
    //			current_edge_a = (*_levelCyclePointTypePtr)[i].first;
    //			current_edge_b = (*_levelCyclePointTypePtr)[i+1].first;
    //			triangle_or_edge_id = TriangleAdjacentToTwoEdges(current_edge_a, current_edge_b);
    //			opposite_vertex = (*inMeshPtr).vecTriangle[triangle_or_edge_id].v0 +
    //							(*inMeshPtr).vecTriangle[triangle_or_edge_id].v1 +
    //							(*inMeshPtr).vecTriangle[triangle_or_edge_id].v2 -
    //							((*inMeshPtr).vecEdge[current_edge_a].v0 + (*inMeshPtr).vecEdge[current_edge_a].v1);
    //			//
    //			current_edge_c = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e01 +
    //							 (*inMeshPtr).vecTriangle[triangle_or_edge_id].e02 +
    //							 (*inMeshPtr).vecTriangle[triangle_or_edge_id].e12 -
    //							 (current_edge_a + current_edge_b);
    //			//
    //			if (processed_vertex.find(opposite_vertex) == processed_vertex.end())
    //			{
    //				if ((*inMeshPtr).vecEdge[current_edge_b].v1 == current_active_node ||
    //					(*inMeshPtr).vecEdge[current_edge_b].v0 == current_active_node)
    //				{
    //					processed_vertex[opposite_vertex] = std::pair<int,int>(((*inMeshPtr).vecEdge[current_edge_b].v0 + (*inMeshPtr).vecEdge[current_edge_b].v1) - opposite_vertex, current_edge_b);
    //				}
    //				else
    //				{
    //					processed_vertex[opposite_vertex] = std::pair<int,int>(((*inMeshPtr).vecEdge[current_edge_c].v0 + (*inMeshPtr).vecEdge[current_edge_c].v1) - opposite_vertex, current_edge_c);
    //				}
    //			}
    //			if (opposite_vertex == tracing_termination_pt)
    //			{
    //				current_active_node = (*inMeshPtr).vecEdge[current_edge_b].v0;
    //				if ((*inMeshPtr).vecEdge[current_edge_b].v0 == tracing_termination_pt)
    //					current_active_node = (*inMeshPtr).vecEdge[current_edge_b].v1;
    //			}
    //			else
    //			{
    //				current_active_node = opposite_vertex;
    //			}
    //

    //		}
    //		else
    //		{// first is on edge , the second is a vertex
    //			opposite_vertex = (*_levelCyclePointTypePtr)[i+1].first;
    //			if (processed_vertex.find(opposite_vertex) == processed_vertex.end())
    //			{
    //				current_edge_a = (*_levelCyclePointTypePtr)[i].first;
    //				FindEdgeConnectingTwoVertices(opposite_vertex, current_active_node, current_edge_a, current_edge_c);
    //				processed_vertex[opposite_vertex] = std::pair<int, int>( current_active_node, current_edge_c);
    //
    //			}
    //			current_active_node = opposite_vertex;
    //		}
    //	}
    //	else
    //	{
    //		if((*_levelCyclePointTypePtr)[i+1].second)
    //		{//
    //			opposite_vertex = (*_levelCyclePointTypePtr)[i].first;
    //			current_edge_c = (*_levelCyclePointTypePtr)[i+1].first;
    //			//
    //			triangle_or_edge_id = TriangleAdjacentToOneEdgesOneVertex(current_edge_c, opposite_vertex);
    //			//
    //			current_edge_a = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e01;
    //			if ((*inMeshPtr).vecTriangle[triangle_or_edge_id].e01 == current_edge_c)
    //				current_edge_a = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e02;
    //			//
    //			current_edge_b = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e01 +
    //							 (*inMeshPtr).vecTriangle[triangle_or_edge_id].e02 +
    //							 (*inMeshPtr).vecTriangle[triangle_or_edge_id].e12 -
    //							 (current_edge_a  + current_edge_c);
    //			//
    //			int pot_vertex = (*inMeshPtr).vecEdge[current_edge_a].v0 + (*inMeshPtr).vecEdge[current_edge_a].v1 - opposite_vertex;
    //			if (processed_vertex.find(pot_vertex) == processed_vertex.end())
    //			{
    //				processed_vertex[pot_vertex] = std::pair<int,int>(opposite_vertex, current_edge_a);
    //			}
    //			//
    //			pot_vertex = (*inMeshPtr).vecEdge[current_edge_b].v0 + (*inMeshPtr).vecEdge[current_edge_b].v1 - opposite_vertex;
    //			if (processed_vertex.find(pot_vertex) == processed_vertex.end())
    //			{
    //				processed_vertex[pot_vertex] = std::pair<int,int>(opposite_vertex, current_edge_b);
    //
    //			}
    //			current_active_node = pot_vertex;
    //		}
    //		else
    //		{
    //			int src_node = (*_levelCyclePointTypePtr)[i].first;
    //			opposite_vertex = (*_levelCyclePointTypePtr)[i+1].first;
    //			if (processed_vertex.find(opposite_vertex) == processed_vertex.end())
    //			{
    //				for (std::vector<int>::iterator vIter = (*inMeshPtr).vecVertex[src_node].adjEdges.begin();
    //					vIter != (*inMeshPtr).vecVertex[src_node].adjEdges.begin();
    //					vIter++)
    //				{
    //					if ((*inMeshPtr).vecEdge[*vIter].v0 == opposite_vertex ||
    //						(*inMeshPtr).vecEdge[*vIter].v1 == opposite_vertex)
    //					{
    //						processed_vertex[opposite_vertex] = std::pair<int,int>(src_node,  *vIter);
    //						break;
    //					}
    //				}
    //			}
    //			current_active_node = opposite_vertex;
    //		}
    //	}

    //}
    //// find the starting point for tracing the closed loop
    //int starting_point = 0;
    //int starting_edge = 0;
    //if (_levelCyclePointTypePtr->back().second)
    //{
    //	int triangle_id = 0;
    //	current_edge_a = (*_levelCyclePointTypePtr).back().first;
    //	if ((*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].second)
    //	{
    //		current_edge_b = (*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].first;
    //		triangle_id = TriangleAdjacentToTwoEdges(current_edge_a, current_edge_b);
    //		starting_point = (*inMeshPtr).vecTriangle[triangle_id].v0 +
    //						 (*inMeshPtr).vecTriangle[triangle_id].v1 +
    //						 (*inMeshPtr).vecTriangle[triangle_id].v2 -
    //						 ((*inMeshPtr).vecEdge[current_edge_a].v0 + (*inMeshPtr).vecEdge[current_edge_a].v1);
    //		//
    //		FindEdgeConnectingTwoVertices(starting_point, tracing_termination_pt, current_edge_a, starting_edge);
    //	}
    //	else
    //	{
    //		starting_point = (*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].first;
    //		FindEdgeConnectingTwoVertices(starting_point, tracing_termination_pt, current_edge_a, starting_edge);
    //	}
    //}
    //else
    //{//
    //	if ((*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].second)
    //	{
    //		current_edge_a = (*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].first;
    //		starting_point = (*inMeshPtr).vecEdge[current_edge_a].v0;
    //		FindEdgeConnectingTwoVertices(tracing_termination_pt, starting_point, current_edge_a, starting_edge);
    //	}
    //	else
    //	{
    //		starting_point = (*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].first;
    //		for (std::vector<int>::iterator vIter = (*inMeshPtr).vecVertex[starting_point].adjEdges.begin();
    //			vIter != (*inMeshPtr).vecVertex[starting_point].adjEdges.begin();
    //			vIter++)
    //		{
    //			if ((*inMeshPtr).vecEdge[*vIter].v0 == tracing_termination_pt ||
    //				(*inMeshPtr).vecEdge[*vIter].v1 == tracing_termination_pt)
    //			{
    //				starting_edge = *vIter;
    //				break;
    //			}
    //		}
    //	}
    //}
    ////
    //outLevelCycle_VertexOnMesh.push_back(tracing_termination_pt);
    //outLevelCycle_VertexOnMesh.push_back(starting_point);
    //outLevelCycle_EdgeOnMesh.push_back(starting_edge);
    //int current_vertex = starting_point;
    //std::pair<int, int> current_vertex_parents = processed_vertex[current_vertex];
    //while (current_vertex_parents.first != -1)
    //{
    //	outLevelCycle_VertexOnMesh.push_back(current_vertex_parents.first);
    //	outLevelCycle_EdgeOnMesh.push_back(current_vertex_parents.second);
    //	current_vertex_parents = processed_vertex[current_vertex_parents.first];
    //}
    ////////////////////////////
    //return;
}
//void MapLoopsBackToMeshLevelSetAndArc::ComputeLevelCycle()
//{
//	// preq: both level cycle and arc cycle , point and pointType are initialized
//
//	//
//	// Trace a vertex-edge path on the mesh which is homotopic to the input cycle on mesh
//	// the second part is the pair to its parent and the edge connecting it to its parent
//	std::map<int, std::pair<int,int>> processed_vertex;
//	int current_edge_a = 0;
//	int current_edge_b = 0;
//	int current_edge_c = 0;
//	int opposite_vertex = 0;
//	int tracing_termination_pt = 0;
//	int current_active_node = 0;
//	if ((*_levelCyclePointTypePtr)[0].second)
//	{// it is an edge
//		current_edge_a = (*_levelCyclePointTypePtr)[0].first;
//		tracing_termination_pt = (*inMeshPtr).vecEdge[current_edge_a].v0;
//		processed_vertex[(*inMeshPtr).vecEdge[current_edge_a].v0] = std::pair<int,int>(-1, -1);
//		processed_vertex[(*inMeshPtr).vecEdge[current_edge_a].v1] = std::pair<int,int>((*inMeshPtr).vecEdge[current_edge_a].v0, current_edge_a);
//		current_active_node = (*inMeshPtr).vecEdge[current_edge_a].v1;
//	}
//	else
//	{// it is only a vertex
//		current_edge_a = (*_levelCyclePointTypePtr)[0].first;
//		tracing_termination_pt = current_edge_a;
//		processed_vertex[current_edge_a] = std::pair<int,int>(-1, -1);
//		current_active_node = current_edge_a;
//	}
//	////////////////////////////
//	for (unsigned int i = 0; i < _levelCyclePointTypePtr->size() - 2; i++)
//	{
//		// find triangle conntaining the edge [i, i+1]
//		int triangle_or_edge_id = 0;
//
//		if ((*_levelCyclePointTypePtr)[i].second)
//		{
//			if ((*_levelCyclePointTypePtr)[i+1].second)
//			{// both points are on edges
//			 // there exists an triangle containing them
//				current_edge_a = (*_levelCyclePointTypePtr)[i].first;
//				current_edge_b = (*_levelCyclePointTypePtr)[i+1].first;
//				triangle_or_edge_id = TriangleAdjacentToTwoEdges(current_edge_a, current_edge_b);
//				opposite_vertex = (*inMeshPtr).vecTriangle[triangle_or_edge_id].v0 +
//								(*inMeshPtr).vecTriangle[triangle_or_edge_id].v1 +
//								(*inMeshPtr).vecTriangle[triangle_or_edge_id].v2 -
//								((*inMeshPtr).vecEdge[current_edge_a].v0 + (*inMeshPtr).vecEdge[current_edge_a].v1);
//				//
//				current_edge_c = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e01 +
//								 (*inMeshPtr).vecTriangle[triangle_or_edge_id].e02 +
//								 (*inMeshPtr).vecTriangle[triangle_or_edge_id].e12 -
//								 (current_edge_a + current_edge_b);
//				//
//				if (processed_vertex.find(opposite_vertex) == processed_vertex.end())
//				{
//					if ((*inMeshPtr).vecEdge[current_edge_b].v1 == current_active_node ||
//						(*inMeshPtr).vecEdge[current_edge_b].v0 == current_active_node)
//					{
//						processed_vertex[opposite_vertex] = std::pair<int,int>(((*inMeshPtr).vecEdge[current_edge_b].v0 + (*inMeshPtr).vecEdge[current_edge_b].v1) - opposite_vertex, current_edge_b);
//					}
//					else
//					{
//						processed_vertex[opposite_vertex] = std::pair<int,int>(((*inMeshPtr).vecEdge[current_edge_c].v0 + (*inMeshPtr).vecEdge[current_edge_c].v1) - opposite_vertex, current_edge_c);
//					}
//				}
//				if (opposite_vertex == tracing_termination_pt)
//				{
//					current_active_node = (*inMeshPtr).vecEdge[current_edge_b].v0;
//					if ((*inMeshPtr).vecEdge[current_edge_b].v0 == tracing_termination_pt)
//						current_active_node = (*inMeshPtr).vecEdge[current_edge_b].v1;
//				}
//				else
//				{
//					current_active_node = opposite_vertex;
//				}
//
//
//			}
//			else
//			{// first is on edge , the second is a vertex
//				opposite_vertex = (*_levelCyclePointTypePtr)[i+1].first;
//				if (processed_vertex.find(opposite_vertex) == processed_vertex.end())
//				{
//					current_edge_a = (*_levelCyclePointTypePtr)[i].first;
//					FindEdgeConnectingTwoVertices(opposite_vertex, current_active_node, current_edge_a, current_edge_c);
//					processed_vertex[opposite_vertex] = std::pair<int, int>( current_active_node, current_edge_c);
//
//				}
//				current_active_node = opposite_vertex;
//			}
//		}
//		else
//		{
//			if((*_levelCyclePointTypePtr)[i+1].second)
//			{//
//				opposite_vertex = (*_levelCyclePointTypePtr)[i].first;
//				current_edge_c = (*_levelCyclePointTypePtr)[i+1].first;
//				//
//				triangle_or_edge_id = TriangleAdjacentToOneEdgesOneVertex(current_edge_c, opposite_vertex);
//				//
//				current_edge_a = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e01;
//				if ((*inMeshPtr).vecTriangle[triangle_or_edge_id].e01 == current_edge_c)
//					current_edge_a = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e02;
//				//
//				current_edge_b = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e01 +
//								 (*inMeshPtr).vecTriangle[triangle_or_edge_id].e02 +
//								 (*inMeshPtr).vecTriangle[triangle_or_edge_id].e12 -
//								 (current_edge_a  + current_edge_c);
//				//
//				int pot_vertex = (*inMeshPtr).vecEdge[current_edge_a].v0 + (*inMeshPtr).vecEdge[current_edge_a].v1 - opposite_vertex;
//				if (processed_vertex.find(pot_vertex) == processed_vertex.end())
//				{
//					processed_vertex[pot_vertex] = std::pair<int,int>(opposite_vertex, current_edge_a);
//				}
//				//
//				pot_vertex = (*inMeshPtr).vecEdge[current_edge_b].v0 + (*inMeshPtr).vecEdge[current_edge_b].v1 - opposite_vertex;
//				if (processed_vertex.find(pot_vertex) == processed_vertex.end())
//				{
//					processed_vertex[pot_vertex] = std::pair<int,int>(opposite_vertex, current_edge_b);
//
//				}
//				current_active_node = pot_vertex;
//			}
//			else
//			{
//				int src_node = (*_levelCyclePointTypePtr)[i].first;
//				opposite_vertex = (*_levelCyclePointTypePtr)[i+1].first;
//				if (processed_vertex.find(opposite_vertex) == processed_vertex.end())
//				{
//					for (std::vector<int>::iterator vIter = (*inMeshPtr).vecVertex[src_node].adjEdges.begin();
//						vIter != (*inMeshPtr).vecVertex[src_node].adjEdges.begin();
//						vIter++)
//					{
//						if ((*inMeshPtr).vecEdge[*vIter].v0 == opposite_vertex ||
//							(*inMeshPtr).vecEdge[*vIter].v1 == opposite_vertex)
//						{
//							processed_vertex[opposite_vertex] = std::pair<int,int>(src_node,  *vIter);
//							break;
//						}
//					}
//				}
//				current_active_node = opposite_vertex;
//			}
//		}
//
//	}
//	// find the starting point for tracing the closed loop
//	int starting_point = 0;
//	int starting_edge = 0;
//	if (_levelCyclePointTypePtr->back().second)
//	{
//		int triangle_id = 0;
//		current_edge_a = (*_levelCyclePointTypePtr).back().first;
//		if ((*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].second)
//		{
//			current_edge_b = (*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].first;
//			triangle_id = TriangleAdjacentToTwoEdges(current_edge_a, current_edge_b);
//			starting_point = (*inMeshPtr).vecTriangle[triangle_id].v0 +
//							 (*inMeshPtr).vecTriangle[triangle_id].v1 +
//							 (*inMeshPtr).vecTriangle[triangle_id].v2 -
//							 ((*inMeshPtr).vecEdge[current_edge_a].v0 + (*inMeshPtr).vecEdge[current_edge_a].v1);
//			//
//			FindEdgeConnectingTwoVertices(starting_point, tracing_termination_pt, current_edge_a, starting_edge);
//		}
//		else
//		{
//			starting_point = (*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].first;
//			FindEdgeConnectingTwoVertices(starting_point, tracing_termination_pt, current_edge_a, starting_edge);
//		}
//	}
//	else
//	{//
//		if ((*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].second)
//		{
//			current_edge_a = (*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].first;
//			starting_point = (*inMeshPtr).vecEdge[current_edge_a].v0;
//			FindEdgeConnectingTwoVertices(tracing_termination_pt, starting_point, current_edge_a, starting_edge);
//		}
//		else
//		{
//			starting_point = (*_levelCyclePointTypePtr)[(*_levelCyclePointTypePtr).size() - 2].first;
//			for (std::vector<int>::iterator vIter = (*inMeshPtr).vecVertex[starting_point].adjEdges.begin();
//				vIter != (*inMeshPtr).vecVertex[starting_point].adjEdges.begin();
//				vIter++)
//			{
//				if ((*inMeshPtr).vecEdge[*vIter].v0 == tracing_termination_pt ||
//					(*inMeshPtr).vecEdge[*vIter].v1 == tracing_termination_pt)
//				{
//					starting_edge = *vIter;
//					break;
//				}
//			}
//		}
//	}
//	//
//	outLevelCycle_VertexOnMesh.push_back(tracing_termination_pt);
//	outLevelCycle_VertexOnMesh.push_back(starting_point);
//	outLevelCycle_EdgeOnMesh.push_back(starting_edge);
//	int current_vertex = starting_point;
//	std::pair<int, int> current_vertex_parents = processed_vertex[current_vertex];
//	while (current_vertex_parents.first != -1)
//	{
//		outLevelCycle_VertexOnMesh.push_back(current_vertex_parents.first);
//		outLevelCycle_EdgeOnMesh.push_back(current_vertex_parents.second);
//		current_vertex_parents = processed_vertex[current_vertex_parents.first];
//	}
//	//////////////////////////
//	// critical point
//	//int crit_vertex_a = _arcCyclePointTypePtr->front().first;
//	//int crit_vertex_b = _arcCyclePointTypePtr->back().first;
//	//std::set<int> processed_vertex;
//	//std::set<int>::iterator findIter;
//	////
//	//int pot_vertex = 0;
//	//// boundary check is done YET !!!!!!!!!!!!!!!
//	//for (unsigned int i = 0; i < _levelCyclePointTypePtr->size() - 1; i++)
//	//{
// //
//	//	if ((*_levelCyclePointTypePtr)[i].second)
//	//	{// this is an edge type vertex
//	//		pot_vertex = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v0;
//	//		//
//	//		if (i != 0 && pot_vertex == outLevelCycle_VertexOnMesh.back())
//	//			continue;
//	//		// check if this vertex is there or not
//	//		// if it is there, ignore it
//	//		// otherwise, push it
//	//		//
//	//			//
//	//		outLevelCycle_VertexOnMesh.push_back(pot_vertex);
//	//			//
//	//			if (i != 0)
//	//			{// find the edge connecting (i-1)--(i)
//	//				// First, if the (i-1)-vertex is on the edge, as the other vertex of this edge is taken
//	//				// so the edge connecting (i-1)--(i) is just the processing edge
//	//				if (inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v0 == outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2] ||
//	//					inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v1 == outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2] )
//	//				{
//	//					outLevelCycle_EdgeOnMesh.push_back((*_levelCyclePointTypePtr)[i].first);
//	//				}
//	//				else
//	//				{
//	//				// find the triangle first containing the edge and prev vertex
//	//					int triangle_id = 0;
//	//					for (unsigned int nTri = 0; nTri < 1; nTri++)
//	//					{
//	//						int tot_sum = inMeshPtr->vecTriangle[inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].AdjTri[nTri]].v0 +
//	//								inMeshPtr->vecTriangle[inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].AdjTri[nTri]].v1 +
//	//								inMeshPtr->vecTriangle[inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].AdjTri[nTri]].v2;
//	//						int tot_edge_sum = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v0 +
//	//								inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v1;
//
//	//						if (tot_sum - tot_edge_sum == outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2])
//	//						{
//	//							triangle_id = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].AdjTri[nTri];
//	//						}
//	//						else
//	//						{
//	//							triangle_id = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].AdjTri[nTri + 1];
//	//						}
//	//					}
//	//					//
//	//					int edge_vec[3] = {inMeshPtr->vecTriangle[triangle_id].e01, inMeshPtr->vecTriangle[triangle_id].e02,
//	//										inMeshPtr->vecTriangle[triangle_id].e12};
//	//					for (int eid = 0; eid < 3; eid++)
//	//					{
//	//						if ((	inMeshPtr->vecEdge[edge_vec[eid]].v0 == pot_vertex ||
//	//								inMeshPtr->vecEdge[edge_vec[eid]].v1 == pot_vertex) &&
//	//							(inMeshPtr->vecEdge[edge_vec[eid]].v0 == outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2] ||
//	//							inMeshPtr->vecEdge[edge_vec[eid]].v1 == outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2]))
//	//						{
//	//							outLevelCycle_EdgeOnMesh.push_back(edge_vec[eid]);
//	//							break;
//	//						}
//	//					}
//	//					//
//	//				}
//	//			}
//	//			//
//
//	//
//	//	}
//	//	else
//	//	{// it is a vertex
//	//		pot_vertex = (*_levelCyclePointTypePtr)[i].first;
//	//		// check if this vertex is there or not
//	//		// if it is there, ignore it
//	//		// otherwise, push it
//	//		//
//	//		if (i != 0 && pot_vertex == outLevelCycle_VertexOnMesh.back())
//	//			continue;
//	//			//
//	//			outLevelCycle_VertexOnMesh.push_back(pot_vertex);
//	//			//
//	//			if (i != 0)
//	//			{
//	//				for (unsigned int k = 0; k < inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2]].adjEdges.size();
//	//					k++)
//	//				{
//	//					if (inMeshPtr->vecEdge[inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2]].adjEdges[k]].v0 == pot_vertex ||
//	//						inMeshPtr->vecEdge[inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2]].adjEdges[k]].v1 == pot_vertex)
//	//					{
//	//						outLevelCycle_EdgeOnMesh.push_back(inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2]].adjEdges[k]);
//	//						break;
//	//					}
//	//				}
//	//			}// if i != 0
//	//
//	//	}// else it is a vertex
//	//}
//	//// Handle the boundary case
//	//if (outLevelCycle_VertexOnMesh.back() != outLevelCycle_VertexOnMesh.front())
//	//{
//	//	for (unsigned int i = 0; i < inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh.front()].adjEdges.size();
//	//		i++)
//	//	{
//	//		if (inMeshPtr->vecEdge[inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh.front()].adjEdges[i]].v0 == outLevelCycle_VertexOnMesh.back() ||
//	//			inMeshPtr->vecEdge[inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh.front()].adjEdges[i]].v1 == outLevelCycle_VertexOnMesh.back())
//	//		{
//	//			outLevelCycle_EdgeOnMesh.push_back(inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh.front()].adjEdges[i]);
//	//			break;
//	//		}
//	//	}
//	//	// add the first vertex into the back to close the loop
//	//	outLevelCycle_VertexOnMesh.push_back(outLevelCycle_VertexOnMesh.front());
//	//}
//	//
//	return;
//}

//void MapLoopsBackToMeshLevelSetAndArc::ComputeLevelCycle()
//{
//	// preq: both level cycle and arc cycle , point and pointType are initialized
//	// critical point
//	int crit_vertex_a = _arcCyclePointTypePtr->front().first;
//	int crit_vertex_b = _arcCyclePointTypePtr->back().first;
//	//
//	std::set<int> processed_vertex;
//	std::set<int>::iterator findIter;
//	//
//	int pot_vertex = 0;
//	// boundary check is done YET !!!!!!!!!!!!!!!
//	for (unsigned int i = 0; i < _levelCyclePointTypePtr->size() - 1; i++)
//	{
//		if (outLevelCycle_VertexOnMesh.size() == 53)
//			i = i;
//		if ((*_levelCyclePointTypePtr)[i].second)
//		{// this is an edge type vertex
//			pot_vertex = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v0;
//			//
//			if (i != 0 && pot_vertex == outLevelCycle_VertexOnMesh.back())
//				continue;
//			findIter = processed_vertex.find(pot_vertex);
//
//			if (findIter != processed_vertex.end())
//			{
//				pot_vertex = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v1;
//				findIter = processed_vertex.find(pot_vertex);
//			}
//			else
//			{
//				//
//				if (inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v1 == crit_vertex_a ||
//					inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v1 == crit_vertex_b)
//				{
//
//					findIter = processed_vertex.find(inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v1);
//					if (findIter == processed_vertex.end())
//						pot_vertex = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v1;
//				}
//			}
//			//
//			if (i != 0 && pot_vertex == outLevelCycle_VertexOnMesh.back())
//				continue;
//			// check if this vertex is there or not
//			// if it is there, ignore it
//			// otherwise, push it
//			//
//			//findIter = processed_vertex.find(pot_vertex);
//			if (i == 0 || pot_vertex != outLevelCycle_VertexOnMesh.back())//findIter == processed_vertex.end())
//			{
//				processed_vertex.insert(pot_vertex);
//				//
//				outLevelCycle_VertexOnMesh.push_back(pot_vertex);
//				//
//				if (i != 0)
//				{// find the edge connecting (i-1)--(i)
//					// First, if the (i-1)-vertex is on the edge, as the other vertex of this edge is taken
//					// so the edge connecting (i-1)--(i) is just the processing edge
//					if (inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v0 == outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2] ||
//						inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v1 == outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2] )
//					{
//						outLevelCycle_EdgeOnMesh.push_back((*_levelCyclePointTypePtr)[i].first);
//					}
//					else
//					{
//					// find the triangle first containing the edge and prev vertex
//						int triangle_id = 0;
//						for (unsigned int nTri = 0; nTri < 1; nTri++)
//						{
//							int tot_sum = inMeshPtr->vecTriangle[inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].AdjTri[nTri]].v0 +
//									inMeshPtr->vecTriangle[inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].AdjTri[nTri]].v1 +
//									inMeshPtr->vecTriangle[inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].AdjTri[nTri]].v2;
//							int tot_edge_sum = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v0 +
//									inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].v1;
//
//							if (tot_sum - tot_edge_sum == outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2])
//							{
//								triangle_id = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].AdjTri[nTri];
//							}
//							else
//							{
//								triangle_id = inMeshPtr->vecEdge[(*_levelCyclePointTypePtr)[i].first].AdjTri[nTri + 1];
//							}
//						}
//						//
//						int edge_vec[3] = {inMeshPtr->vecTriangle[triangle_id].e01, inMeshPtr->vecTriangle[triangle_id].e02,
//											inMeshPtr->vecTriangle[triangle_id].e12};
//						for (int eid = 0; eid < 3; eid++)
//						{
//							if ((	inMeshPtr->vecEdge[edge_vec[eid]].v0 == pot_vertex ||
//									inMeshPtr->vecEdge[edge_vec[eid]].v1 == pot_vertex) &&
//								(inMeshPtr->vecEdge[edge_vec[eid]].v0 == outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2] ||
//								inMeshPtr->vecEdge[edge_vec[eid]].v1 == outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2]))
//							{
//								outLevelCycle_EdgeOnMesh.push_back(edge_vec[eid]);
//								break;
//							}
//						}
//						//
//					}
//				}// if findIter
//				//
//
//			}
//			// else
//			// do nothing
//		}
//		else
//		{// it is a vertex
//			pot_vertex = (*_levelCyclePointTypePtr)[i].first;
//			// check if this vertex is there or not
//			// if it is there, ignore it
//			// otherwise, push it
//			//
//			if (i != 0 && pot_vertex == outLevelCycle_VertexOnMesh.back())
//				continue;
//			findIter = processed_vertex.find(pot_vertex);
//			if (findIter == processed_vertex.end())
//			{
//				processed_vertex.insert(pot_vertex);
//				//
//				outLevelCycle_VertexOnMesh.push_back(pot_vertex);
//				//
//				if (i != 0)
//				{
//					for (unsigned int k = 0; k < inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2]].adjEdges.size();
//						k++)
//					{
//						if (inMeshPtr->vecEdge[inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2]].adjEdges[k]].v0 == pot_vertex ||
//							inMeshPtr->vecEdge[inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2]].adjEdges[k]].v1 == pot_vertex)
//						{
//							outLevelCycle_EdgeOnMesh.push_back(inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh[outLevelCycle_VertexOnMesh.size() - 2]].adjEdges[k]);
//							break;
//						}
//					}
//				}// if i != 0
//			}// if findIter
//		}// else it is a vertex
//	}
//	// Handle the boundary case
//	for (unsigned int i = 0; i < inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh.front()].adjEdges.size();
//		i++)
//	{
//		if (inMeshPtr->vecEdge[inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh.front()].adjEdges[i]].v0 == outLevelCycle_VertexOnMesh.back() ||
//			inMeshPtr->vecEdge[inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh.front()].adjEdges[i]].v1 == outLevelCycle_VertexOnMesh.back())
//		{
//			outLevelCycle_EdgeOnMesh.push_back(inMeshPtr->vecVertex[outLevelCycle_VertexOnMesh.front()].adjEdges[i]);
//			break;
//		}
//	}
//	// add the first vertex into the back to close the loop
//	outLevelCycle_VertexOnMesh.push_back(outLevelCycle_VertexOnMesh.front());
//	//
//	return;
//}
//
void MapLoopsBackToMeshLevelSetAndArc::WalkAlongArcWithInitialSide_InverseArcOrder(const int init_vertex_for_side,
                                                                                   int &current_active_vertex,
                                                                                   int &refEdge, int &outTriagleId,
                                                                                   const std::vector<std::pair<int, int> > &arcPointType,
                                                                                   std::vector<int> &vecArc_vertex,
                                                                                   std::vector<int> &vecArc_edge) {//
    // preq: the arc is not an edge,
    // current_active_vertex is the lowest critical point on the mesh
    int pot_vertex = 0;
    int pot_edge_id = 0;
    int current_edge_id = 0;
    //init_vertex_for_side is on the first edge and it is NOT in the vertex_array
    int shared_vertex = 0;
    if (arcPointType.size() == 2) {// init_vertex_for_side is the opposite vertex to current_ative_vertex
        vecArc_vertex.push_back(init_vertex_for_side);
        //
        int three_edges[3] = {inMeshPtr->vecTriangle[outTriagleId].e01,
                              inMeshPtr->vecTriangle[outTriagleId].e02,
                              inMeshPtr->vecTriangle[outTriagleId].e12};
        for (int teid = 0; teid < 3; teid++) {
            if ((inMeshPtr->vecEdge[three_edges[teid]].v0 == current_active_vertex ||
                 inMeshPtr->vecEdge[three_edges[teid]].v1 == current_active_vertex) &&
                (inMeshPtr->vecEdge[three_edges[teid]].v0 == init_vertex_for_side ||
                 inMeshPtr->vecEdge[three_edges[teid]].v1 == init_vertex_for_side)) {
                pot_edge_id = three_edges[teid];
                break;
            }
        }
        //
        vecArc_edge.push_back(pot_edge_id);
        //
        // update current_active_vertex and refEdge
        int third_vertex = inMeshPtr->vecTriangle[outTriagleId].v0 +
                           inMeshPtr->vecTriangle[outTriagleId].v1 +
                           inMeshPtr->vecTriangle[outTriagleId].v2 -
                           (current_active_vertex + init_vertex_for_side);
        //
        current_active_vertex = init_vertex_for_side;
        for (int teid = 0; teid < 3; teid++) {
            if ((inMeshPtr->vecEdge[three_edges[teid]].v0 == current_active_vertex ||
                 inMeshPtr->vecEdge[three_edges[teid]].v1 == current_active_vertex) &&
                (inMeshPtr->vecEdge[three_edges[teid]].v0 == third_vertex ||
                 inMeshPtr->vecEdge[three_edges[teid]].v1 == third_vertex)) {
                refEdge = three_edges[teid];
                break;
            }
        }
    } else {
        /*
        Initializing current active vertex
        */
        if (!arcPointType[arcPointType.size() - 2].second) {//  it is only a vertex
            std::cout << "vertex except critical points on the arc ?" << std::endl;
            exit(0);
        }
        // Handle the intial case
        current_edge_id = arcPointType[arcPointType.size() - 2].first;
        pot_vertex = init_vertex_for_side;
        //
        FindEdgeConnectingTwoVertices(current_active_vertex, pot_vertex, current_edge_id, pot_edge_id);
        //
        current_active_vertex = pot_vertex;
        //
        vecArc_vertex.push_back(current_active_vertex);
        vecArc_edge.push_back(pot_edge_id);
        //
        for (unsigned int i = arcPointType.size() - 3; i > 0; i--) {
            /*
            Invariant : current_active_vertex is always on current edge if current array element is an edge
                        or it is the current vertex if current array element is a vertex
            */
            if (arcPointType[i].second) {// this is an edge type vertex
                current_edge_id = arcPointType[i].first;
                if (arcPointType[i + 1].second) {// previous has one edge
                    if (inMeshPtr->vecEdge[current_edge_id].v0 != current_active_vertex &&
                        inMeshPtr->vecEdge[current_edge_id].v1 !=
                        current_active_vertex) {// if current active vertex is on this edge, then do nothing
                        // else add one more edge and update current_active_vertex
                        // since the levele set cross the edge only once, so there must exist at least one vertex on the edge not visited yet
                        shared_vertex = inMeshPtr->vecEdge[arcPointType[i + 1].first].v0 +
                                        inMeshPtr->vecEdge[arcPointType[i + 1].first].v1 - current_active_vertex;
                        pot_vertex = inMeshPtr->vecEdge[current_edge_id].v0;
                        if (shared_vertex == pot_vertex)
                            pot_vertex = inMeshPtr->vecEdge[current_edge_id].v1;
                        //
                        FindEdgeConnectingTwoVertices(current_active_vertex, pot_vertex, current_edge_id, pot_edge_id);
                        //
                        current_active_vertex = pot_vertex;
                        //
                        vecArc_vertex.push_back(current_active_vertex);
                        vecArc_edge.push_back(pot_edge_id);
                        //
                    }
                } else {// previous one is an vertex
                    std::cout << "vertex except critical points on the arc ?" << std::endl;
                    exit(0);
                }
            } else {// it is a vertex
                std::cout << "vertex except critical points on the arc ?" << std::endl;
                exit(0);
            }// else it is a vertex
        }
        // Handle the boundary case
        pot_vertex = arcPointType.front().first; // the last vertex
        if (current_active_vertex != pot_vertex) {
            current_edge_id = arcPointType[1].first;
            FindEdgeConnectingTwoVertices(pot_vertex, current_active_vertex, current_edge_id, pot_edge_id);
            //
            vecArc_vertex.push_back(pot_vertex);
            vecArc_edge.push_back(pot_edge_id);
            //
            outTriagleId = TriangleAdjacentToOneEdgesOneVertex(current_edge_id, pot_vertex);
            refEdge = pot_edge_id;
            current_active_vertex = pot_vertex;
        } else {//
            std::cout << "the one before downforking pt is a vertex ?" << std::endl;
            exit(1);
        }
    }
    // invariant : outTrinagleId is the one triangle which is adjacent to the last edge in vecArc_edge
    return;
}

void MapLoopsBackToMeshLevelSetAndArc::WalkAlongArcWithInitialSide_ArcOrder(const int init_vertex_for_side,
                                                                            int &current_active_vertex, int &refEdge,
                                                                            int &outTriagleId,
                                                                            const std::vector<std::pair<int, int> > &arcPointType,
                                                                            std::vector<int> &vecArc_vertex,
                                                                            std::vector<int> &vecArc_edge) {//
    // preq: the arc is not an edge,
    // current_active_vertex is the lowest critical point on the mesh
    int pot_vertex = 0;
    int pot_edge_id = 0;
    int current_edge_id = 0;
    //init_vertex_for_side is on the first edge and it is NOT in the vertex_array
    int shared_vertex = 0;
    if (arcPointType.size() == 2) {// init_vertex_for_side is the opposite vertex to current_ative_vertex
        vecArc_vertex.push_back(init_vertex_for_side);
        //
        int three_edges[3] = {inMeshPtr->vecTriangle[outTriagleId].e01,
                              inMeshPtr->vecTriangle[outTriagleId].e02,
                              inMeshPtr->vecTriangle[outTriagleId].e12};
        for (int teid = 0; teid < 3; teid++) {
            if ((inMeshPtr->vecEdge[three_edges[teid]].v0 == current_active_vertex ||
                 inMeshPtr->vecEdge[three_edges[teid]].v1 == current_active_vertex) &&
                (inMeshPtr->vecEdge[three_edges[teid]].v0 == init_vertex_for_side ||
                 inMeshPtr->vecEdge[three_edges[teid]].v1 == init_vertex_for_side)) {
                pot_edge_id = three_edges[teid];
                break;
            }
        }
        //
        vecArc_edge.push_back(pot_edge_id);
        //
        // update current_active_vertex and refEdge
        int third_vertex = inMeshPtr->vecTriangle[outTriagleId].v0 +
                           inMeshPtr->vecTriangle[outTriagleId].v1 +
                           inMeshPtr->vecTriangle[outTriagleId].v2 -
                           (current_active_vertex + init_vertex_for_side);
        //
        current_active_vertex = init_vertex_for_side;
        for (int teid = 0; teid < 3; teid++) {
            if ((inMeshPtr->vecEdge[three_edges[teid]].v0 == current_active_vertex ||
                 inMeshPtr->vecEdge[three_edges[teid]].v1 == current_active_vertex) &&
                (inMeshPtr->vecEdge[three_edges[teid]].v0 == third_vertex ||
                 inMeshPtr->vecEdge[three_edges[teid]].v1 == third_vertex)) {
                refEdge = three_edges[teid];
                break;
            }
        }
    } else {
        /*
        Initializing current active vertex
        */
        if (!arcPointType[1].second) {//  it is only a vertex
            std::cout << "vertex except critical points on the arc ?" << std::endl;
            exit(0);
        }
        // Handle the intial case
        current_edge_id = arcPointType[1].first;
        pot_vertex = init_vertex_for_side;
        //
        FindEdgeConnectingTwoVertices(current_active_vertex, pot_vertex, current_edge_id, pot_edge_id);
        //
        current_active_vertex = pot_vertex;
        //
        vecArc_vertex.push_back(current_active_vertex);
        vecArc_edge.push_back(pot_edge_id);
        //
        for (unsigned int i = 2; i < arcPointType.size() - 1; i++) {
            /*
            Invariant : current_active_vertex is always on current edge if current array element is an edge
                        or it is the current vertex if current array element is a vertex
            */
            if (arcPointType[i].second) {// this is an edge type vertex
                current_edge_id = arcPointType[i].first;
                if (arcPointType[i - 1].second) {// previous has one edge
                    if (inMeshPtr->vecEdge[current_edge_id].v0 != current_active_vertex &&
                        inMeshPtr->vecEdge[current_edge_id].v1 !=
                        current_active_vertex) {// if current active vertex is on this edge, then do nothing
                        // else add one more edge and update current_active_vertex
                        // since the levele set cross the edge only once, so there must exist at least one vertex on the edge not visited yet
                        shared_vertex = inMeshPtr->vecEdge[arcPointType[i - 1].first].v0 +
                                        inMeshPtr->vecEdge[arcPointType[i - 1].first].v1 - current_active_vertex;
                        pot_vertex = inMeshPtr->vecEdge[current_edge_id].v0;
                        if (shared_vertex == pot_vertex)
                            pot_vertex = inMeshPtr->vecEdge[current_edge_id].v1;
                        //
                        FindEdgeConnectingTwoVertices(current_active_vertex, pot_vertex, current_edge_id, pot_edge_id);
                        //
                        current_active_vertex = pot_vertex;
                        //
                        vecArc_vertex.push_back(current_active_vertex);
                        vecArc_edge.push_back(pot_edge_id);
                        //
                    }
                } else {// previous one is an vertex
                    std::cout << "vertex except critical points on the arc ?" << std::endl;
                    exit(0);
                }
            } else {// it is a vertex
                std::cout << "vertex except critical points on the arc ?" << std::endl;
                exit(0);
            }// else it is a vertex
        }
        // Handle the boundary case
        pot_vertex = arcPointType.back().first; // the last vertex
        if (current_active_vertex != pot_vertex) {
            current_edge_id = arcPointType[arcPointType.size() - 2].first;
            FindEdgeConnectingTwoVertices(pot_vertex, current_active_vertex, current_edge_id, pot_edge_id);
            //
            vecArc_vertex.push_back(pot_vertex);
            vecArc_edge.push_back(pot_edge_id);
            //
            outTriagleId = TriangleAdjacentToOneEdgesOneVertex(current_edge_id, pot_vertex);
            refEdge = pot_edge_id;
            current_active_vertex = pot_vertex;
        } else {//
            std::cout << "the one before downforking pt is a vertex ?" << std::endl;
            exit(1);
        }
        // invariant : outTrinagleId is the one triangle which is adjacent to the last edge in vecArc_edge
    }
    return;
}

void MapLoopsBackToMeshLevelSetAndArc::ComputeArcCycle() {
//	// critical point
    if ((*_arcCyclePointTypePtr).size() == 2) {// this is arc corresponding to an edge on the mesh
        int dst_node = (*_arcCyclePointTypePtr).back().first;
        int src_node = (*_arcCyclePointTypePtr).front().first;
        for (std::vector<int>::iterator vIter = (*inMeshPtr).vecVertex[src_node].adjEdges.begin();
             vIter != (*inMeshPtr).vecVertex[src_node].adjEdges.end();
             vIter++) {
            if ((*inMeshPtr).vecEdge[*vIter].v0 == dst_node ||
                (*inMeshPtr).vecEdge[*vIter].v1 == dst_node) {
                outArcCycle_EdgeOnMesh.push_back(*vIter);
                break;
            }
        }
        //
        outArcCycle_VertexOnMesh.push_back(src_node);
        outArcCycle_VertexOnMesh.push_back(dst_node);
    } else {
        std::set<int> processed_vertex;
        std::set<int>::iterator findIter;
        //
        int pot_vertex = 0;
        int pot_edge_id = 0;
        int current_edge_id = 0;
        int current_active_vertex = 0;
        int shared_vertex = 0;
        /*
        Initializing current active vertex
        */
        if ((*_arcCyclePointTypePtr)[0].second) {// it is an edge
            current_edge_id = (*_arcCyclePointTypePtr)[0].first;
            current_active_vertex = (*inMeshPtr).vecEdge[current_edge_id].v0;
            //
            processed_vertex.insert(current_active_vertex);
            //
            outArcCycle_VertexOnMesh.push_back(current_active_vertex);
        } else {// it is only a vertex
            current_active_vertex = (*_arcCyclePointTypePtr)[0].first;
            //
            processed_vertex.insert(current_active_vertex);
            //
            outArcCycle_VertexOnMesh.push_back(current_active_vertex);
        }
        // boundary check is done YET !!!!!!!!!!!!!!!
        for (unsigned int i = 1; i < _arcCyclePointTypePtr->size() - 1; i++) {
            /*
            Invariant : current_active_vertex is always on current edge if current array element is an edge
                        or it is the current vertex if current array element is a vertex
            */
            if ((*_arcCyclePointTypePtr)[i].second) {// this is an edge type vertex
                current_edge_id = (*_arcCyclePointTypePtr)[i].first;
                if ((*_arcCyclePointTypePtr)[i - 1].second) {// previous has one edge
                    if (inMeshPtr->vecEdge[current_edge_id].v0 != current_active_vertex &&
                        inMeshPtr->vecEdge[current_edge_id].v1 !=
                        current_active_vertex) {// if current active vertex is on this edge, then do nothing
                        // else add one more edge and update current_active_vertex
                        // since the levele set cross the edge only once, so there must exist at least one vertex on the edge not visited yet
                        shared_vertex = inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i - 1].first].v0 +
                                        inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i - 1].first].v1 -
                                        current_active_vertex;
                        pot_vertex = inMeshPtr->vecEdge[current_edge_id].v0;
                        if (shared_vertex == pot_vertex)
                            pot_vertex = inMeshPtr->vecEdge[current_edge_id].v1;
                        //
                        FindEdgeConnectingTwoVertices(current_active_vertex, pot_vertex, current_edge_id, pot_edge_id);
                        //
                        current_active_vertex = pot_vertex;
                        //
                        outArcCycle_VertexOnMesh.push_back(current_active_vertex);
                        outArcCycle_EdgeOnMesh.push_back(pot_edge_id);
                        //
                        processed_vertex.insert(current_active_vertex);
                    }
                } else {// previous one is an vertex
                    if (inMeshPtr->vecEdge[current_edge_id].v0 != current_active_vertex &&
                        inMeshPtr->vecEdge[current_edge_id].v1 != current_active_vertex) {// take arbitrary decision
                        pot_vertex = inMeshPtr->vecEdge[current_edge_id].v0;
                        //
                        FindEdgeConnectingTwoVertices(current_active_vertex, pot_vertex, current_edge_id, pot_edge_id);
                        //
                        current_active_vertex = pot_vertex;
                        //
                        outArcCycle_VertexOnMesh.push_back(current_active_vertex);
                        outArcCycle_EdgeOnMesh.push_back(pot_edge_id);
                        //
                        processed_vertex.insert(current_active_vertex);
                    } else {//
                        std::cout << "potential error in arc computation" << std::endl;
                        exit(0);
                    }
                }
            } else {// it is a vertex
                pot_vertex = (*_arcCyclePointTypePtr)[i].first;
                // check if this vertex is there or not
                // if it is there, ignore it
                // otherwise, push it
                //
                pot_edge_id = -1;
                if (processed_vertex.find(pot_vertex) == processed_vertex.end()) {
                    for (std::vector<int>::iterator vIter = inMeshPtr->vecVertex[current_active_vertex].adjEdges.begin();
                         vIter != inMeshPtr->vecVertex[current_active_vertex].adjEdges.end(); vIter++) {
                        if (inMeshPtr->vecEdge[*vIter].v0 == pot_vertex ||
                            inMeshPtr->vecEdge[*vIter].v1 == pot_vertex) {
                            pot_edge_id = *vIter;
                            break;
                        }
                    }
                    if (pot_edge_id > -1) {
                        outArcCycle_VertexOnMesh.push_back(pot_vertex);
                        outArcCycle_EdgeOnMesh.push_back(pot_edge_id);
                    } else {
                        std::cout << "No edge connecting two vertices " << std::endl;
                        exit(0);
                    }
                    //
                    current_active_vertex = pot_vertex;
                    //
                    processed_vertex.insert(current_active_vertex);
                }
            }// else it is a vertex
        }
        // Handle the boundary case
        pot_vertex = (*_arcCyclePointTypePtr).back().first; // the last vertex
        if (current_active_vertex != pot_vertex) {
            pot_edge_id = -1;
            for (std::vector<int>::iterator vIter = inMeshPtr->vecVertex[current_active_vertex].adjEdges.begin();
                 vIter != inMeshPtr->vecVertex[current_active_vertex].adjEdges.end(); vIter++) {
                if (inMeshPtr->vecEdge[*vIter].v0 == pot_vertex ||
                    inMeshPtr->vecEdge[*vIter].v1 == pot_vertex) {
                    pot_edge_id = *vIter;
                    break;
                }
            }
            if (pot_edge_id > -1) {
                outArcCycle_VertexOnMesh.push_back(pot_vertex);
                outArcCycle_EdgeOnMesh.push_back(pot_edge_id);
            } else {
                std::cout << "No edge connecting two vertices " << std::endl;
                exit(0);
            }
        }

    }
    return;
/*
	if ((*_arcCyclePointTypePtr).size() == 2)
	{// this is arc corresponding to an edge on the mesh
		int dst_node = (*_arcCyclePointTypePtr).back().first;
		int src_node = (*_arcCyclePointTypePtr).front().first;
		for (std::vector<int>::iterator vIter = (*inMeshPtr).vecVertex[src_node].adjEdges.begin();
			vIter != (*inMeshPtr).vecVertex[src_node].adjEdges.end();
			vIter++)
		{
			if ((*inMeshPtr).vecEdge[*vIter].v0 == dst_node ||
				(*inMeshPtr).vecEdge[*vIter].v1 == dst_node)
			{
				outArcCycle_EdgeOnMesh.push_back(*vIter);
				break;
			}
		}
		//
		outArcCycle_VertexOnMesh.push_back(dst_node);
		outArcCycle_VertexOnMesh.push_back(src_node);
	}
	else
	{
		boost::unordered_set<int> triangleDegree;
		std::vector<int> triangleOrder ;
		for (unsigned int i = 0; i < _arcCyclePointTypePtr->size() - 1; i++)
		{
			int triangle_id = 0;
			int edge_a = (*_arcCyclePointTypePtr)[i].first;
			int edge_b = (*_arcCyclePointTypePtr)[i + 1].first;
			if (i == 0)
			{
				triangle_id = TriangleAdjacentToOneEdgesOneVertex(edge_b, edge_a);
			}
			else
			{
				if (i == _arcCyclePointTypePtr->size() - 2)
				{
					triangle_id = TriangleAdjacentToOneEdgesOneVertex(edge_a, edge_b);
				}
				else
				{
					triangle_id = TriangleAdjacentToTwoEdges(edge_a, edge_b);
				}
			}
			//
			//
			boost::unordered_set<int>::iterator sIter = triangleDegree.find(triangle_id);
			if (sIter == triangleDegree.end())
			{
				triangleDegree.insert(triangle_id);
				triangleOrder.push_back(triangle_id);
			}
			else
			{//
				std::cout << "NOT a triangle trip" << std::endl;
				exit(9);
			}
			//
		}
		//
		_simpGraph triangleStripGraph;
		boost::unordered_map<int, int> verToGraphNode;
		boost::unordered_map<int, int>::iterator mIter;
		int nGrapNodeCounter = 0;
		for (unsigned int i = 0; i < triangleOrder.size(); i++)
		{
			int node[3] = {	(*inMeshPtr).vecTriangle[triangleOrder[i]].v0,
							(*inMeshPtr).vecTriangle[triangleOrder[i]].v1,
							(*inMeshPtr).vecTriangle[triangleOrder[i]].v2};
			int gNode[3] = {0};
			for (int j = 0; j < 3; j++)
			{
				mIter = verToGraphNode.find(node[j]);
				if (mIter == verToGraphNode.end())
				{
					triangleStripGraph.AddNode(nGrapNodeCounter, node[j]);
					gNode[j] = nGrapNodeCounter;
					verToGraphNode[node[j]] = nGrapNodeCounter;
					nGrapNodeCounter++;
				}
				else
				{
					gNode[j] = mIter->second;
				}
			}
			triangleStripGraph.AddEdge(gNode[0], gNode[1]); // since graph node use set to store adjacent list
			triangleStripGraph.AddEdge(gNode[0], gNode[2]); // no need to check the edge exists or not
			triangleStripGraph.AddEdge(gNode[1], gNode[2]);
		}
		//
		for (unsigned int i = 1; i < _arcCyclePointTypePtr->size() - 1; i++)
		{
			int edge_idx = (*_arcCyclePointTypePtr)[i].first;
			triangleStripGraph.RemoveEdge( verToGraphNode[(*inMeshPtr).vecEdge[edge_idx].v0],  verToGraphNode[(*inMeshPtr).vecEdge[edge_idx].v1]);
		}
		//there exist two disjont two components
		std::vector<bool> colors(triangleStripGraph.vecNode.size(), false);
		std::vector<int> parents(triangleStripGraph.vecNode.size(), -1);
		int components = 1;
		int currentNode = verToGraphNode[_arcCyclePointTypePtr->front().first];
		int oppositeNode = verToGraphNode[_arcCyclePointTypePtr->back().first];
	// visit the components
		std::queue<int> Q;
		Q.push(currentNode);
		colors[currentNode] = true;
		parents[currentNode] = -1;
		while (!Q.empty())
		{
			int node = Q.front();
			Q.pop();
			for (std::set<int>::iterator sIter = triangleStripGraph.vecNode[node].adjList.begin();
				sIter != triangleStripGraph.vecNode[node].adjList.end(); sIter++)
			{
				if (!colors[*sIter])
				{
					colors[*sIter] = true;
					Q.push(*sIter);
					parents[*sIter] = node;
				}
			}
		}
		for (unsigned int i = 0; i < colors.size(); i++)
		{
			if (!colors[i])
			{
				components++;
			}
		}
		if (components != 1)
		{
			std::cout << "wrong method again for arc" << std::endl;
			exit(4);
		}
		outLevelCycle_VertexOnMesh.reserve(triangleStripGraph.vecNode.size());
		outLevelCycle_EdgeOnMesh.reserve(triangleStripGraph.vecNode.size());
		//

		//
		//
		int iterNode = oppositeNode;
		while (iterNode >= 0)
		{
			outArcCycle_VertexOnMesh.push_back(triangleStripGraph.vecNode[iterNode].color);
			iterNode = parents[iterNode];
		}
		if (outArcCycle_VertexOnMesh.back() != triangleStripGraph.vecNode[currentNode].color)
		{
			std::cout << "in different compoents for arc" << std::endl;
			exit(4);
		}
		//
		for (unsigned int i = 0; i < outArcCycle_VertexOnMesh.size() - 1; i++)
		{
			int curNode = outArcCycle_VertexOnMesh[i];
			int oppNode = outArcCycle_VertexOnMesh[i + 1];
			for (std::vector<int>::iterator vIter = (*inMeshPtr).vecVertex[curNode].adjEdges.begin();
				vIter != (*inMeshPtr).vecVertex[curNode].adjEdges.end(); vIter++)
			{
				if ((*inMeshPtr).vecEdge[*vIter].v0 == oppNode ||
					(*inMeshPtr).vecEdge[*vIter].v1 == oppNode)
				{
					outArcCycle_EdgeOnMesh.push_back(*vIter);
					break;
				}
			}
		}
	}
	return;
*/
/////////////////////////////////////////////

/*	std::map<int, std::pair<int,int>> processed_vertex;
	int current_edge_a = 0;
	int current_edge_b = 0;
	int current_edge_c = 0;
	int opposite_vertex = 0;
	int tracing_termination_pt = 0;
	//
	if ((*_arcCyclePointTypePtr).size() == 2)
	{// this is arc corresponding to an edge on the mesh
		int dst_node = (*_arcCyclePointTypePtr).back().first;
		int src_node = (*_arcCyclePointTypePtr).front().first;
		for (std::vector<int>::iterator vIter = (*inMeshPtr).vecVertex[src_node].adjEdges.begin();
			vIter != (*inMeshPtr).vecVertex[src_node].adjEdges.end();
			vIter++)
		{
			if ((*inMeshPtr).vecEdge[*vIter].v0 == dst_node ||
				(*inMeshPtr).vecEdge[*vIter].v1 == dst_node)
			{
				outArcCycle_EdgeOnMesh.push_back(*vIter);
				break;
			}
		}
		//
		outArcCycle_VertexOnMesh.push_back(dst_node);
		outArcCycle_VertexOnMesh.push_back(src_node);
	}
	else
	{ // it always starts with a vertex
		////////////////////////////
		current_edge_a = (*_arcCyclePointTypePtr)[0].first;
		tracing_termination_pt = current_edge_a;
		processed_vertex[current_edge_a] = std::pair<int,int>(-1, -1);

		////////////////////////////
		for (unsigned int i = 0; i < _arcCyclePointTypePtr->size() - 1; i++)
		{
			// find triangle conntaining the edge [i, i+1]
			int triangle_or_edge_id = 0;
			if ((*_arcCyclePointTypePtr)[i].second)
			{
				if ((*_arcCyclePointTypePtr)[i+1].second)
				{// both points are on edges
				 // there exists an triangle containing them
					current_edge_a = (*_arcCyclePointTypePtr)[i].first;
					current_edge_b = (*_arcCyclePointTypePtr)[i+1].first;
					triangle_or_edge_id = TriangleAdjacentToTwoEdges(current_edge_a, current_edge_b);
					opposite_vertex = (*inMeshPtr).vecTriangle[triangle_or_edge_id].v0 +
									  (*inMeshPtr).vecTriangle[triangle_or_edge_id].v1 +
									  (*inMeshPtr).vecTriangle[triangle_or_edge_id].v2 -
									  ((*inMeshPtr).vecEdge[current_edge_a].v0 + (*inMeshPtr).vecEdge[current_edge_a].v1);
					//
					if (processed_vertex.find(opposite_vertex) == processed_vertex.end())
					{
						processed_vertex[opposite_vertex] = std::pair<int,int>(((*inMeshPtr).vecEdge[current_edge_b].v0 + (*inMeshPtr).vecEdge[current_edge_b].v1) - opposite_vertex, current_edge_b);
					}
				}
				else
				{// first is on edge , the second is a vertex
					opposite_vertex = (*_arcCyclePointTypePtr)[i+1].first;
					if (processed_vertex.find(opposite_vertex) == processed_vertex.end())
					{
						current_edge_a = (*_arcCyclePointTypePtr)[i].first;
						FindEdgeConnectingTwoVertices(opposite_vertex, (*inMeshPtr).vecEdge[current_edge_a].v0, current_edge_a, current_edge_c);
						processed_vertex[opposite_vertex] = std::pair<int, int>( (*inMeshPtr).vecEdge[current_edge_a].v0, current_edge_c);
					}
				}
			}
			else
			{
				if ((*_arcCyclePointTypePtr)[i+1].second)
				{//
					opposite_vertex = (*_arcCyclePointTypePtr)[i].first;
					current_edge_c = (*_arcCyclePointTypePtr)[i+1].first;
					//
					triangle_or_edge_id = TriangleAdjacentToOneEdgesOneVertex(current_edge_c, opposite_vertex);
					//
					current_edge_a = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e01;
					if ((*inMeshPtr).vecTriangle[triangle_or_edge_id].e01 == current_edge_c)
						current_edge_a = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e02;
					//
					current_edge_b = (*inMeshPtr).vecTriangle[triangle_or_edge_id].e01 +
										(*inMeshPtr).vecTriangle[triangle_or_edge_id].e02 +
										(*inMeshPtr).vecTriangle[triangle_or_edge_id].e12 -
									(current_edge_a  + current_edge_c);
					//
					int pot_vertex = (*inMeshPtr).vecEdge[current_edge_a].v0 + (*inMeshPtr).vecEdge[current_edge_a].v1 - opposite_vertex;
					if (processed_vertex.find(pot_vertex) == processed_vertex.end())
					{
						processed_vertex[pot_vertex] = std::pair<int,int>(opposite_vertex, current_edge_a);
					}
					//
					pot_vertex = (*inMeshPtr).vecEdge[current_edge_b].v0 + (*inMeshPtr).vecEdge[current_edge_b].v1 - opposite_vertex;
					if (processed_vertex.find(pot_vertex) == processed_vertex.end())
					{
						processed_vertex[pot_vertex] = std::pair<int,int>(opposite_vertex, current_edge_b);
					}
				}
				else
				{
					int src_node = (*_arcCyclePointTypePtr)[i].first;
					opposite_vertex = (*_arcCyclePointTypePtr)[i+1].first;
					if (processed_vertex.find(opposite_vertex) == processed_vertex.end())
					{
						for (std::vector<int>::iterator vIter = (*inMeshPtr).vecVertex[src_node].adjEdges.begin();
							vIter != (*inMeshPtr).vecVertex[src_node].adjEdges.begin();
							vIter++)
						{
							if ((*inMeshPtr).vecEdge[*vIter].v0 == opposite_vertex ||
								(*inMeshPtr).vecEdge[*vIter].v1 == opposite_vertex)
							{
								processed_vertex[opposite_vertex] = std::pair<int,int>(src_node,  *vIter);
								break;
							}
						}
					}
				}
			}// else i
		}
		//
		int current_vertex = (*_arcCyclePointTypePtr).back().first;
		outArcCycle_VertexOnMesh.push_back(current_vertex);
		std::pair<int, int> current_vertex_parents = processed_vertex[current_vertex];
		while (current_vertex_parents.first != -1)
		{
			outArcCycle_VertexOnMesh.push_back(current_vertex_parents.first);
			outArcCycle_EdgeOnMesh.push_back(current_vertex_parents.second);
			current_vertex_parents = processed_vertex[current_vertex_parents.first];
		}
	}// else == 2

 */
    //////////////////////////
    ////////////////////////////////////////////////////
    // preq: arc cycle point and pointType are initialized
    ////std::cout << "IN ComputeArcCycle" << std::endl;
    //std::set<int> processed_vertex;
    //std::set<int>::iterator findIter;
    ////
    //int pot_vertex = 0;
    ////
    //if (_arcCyclePointTypePtr->size() == 3)
    //{
    //	if (	(*_arcCyclePointTypePtr)[1].second &&
    //			!(*_arcCyclePointTypePtr)[0].second &&
    //			!(*_arcCyclePointTypePtr)[2].second &&
    //			(
    //			(inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[1].first].v0 == (*_arcCyclePointTypePtr)[0].first ||
    //			 inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[1].first].v0 == (*_arcCyclePointTypePtr)[2].first) &&
    //			(inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[1].first].v1 == (*_arcCyclePointTypePtr)[0].first ||
    //			 inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[1].first].v1 == (*_arcCyclePointTypePtr)[2].first)
    //			)
    //		)
    //	{// this is the 2 points case
    //		outArcCycle_VertexOnMesh.push_back((*_arcCyclePointTypePtr)[0].first);
    //		outArcCycle_VertexOnMesh.push_back((*_arcCyclePointTypePtr)[2].first);
    //		outArcCycle_EdgeOnMesh.push_back((*_arcCyclePointTypePtr)[1].first);
    //		return;
    //	}
    //}
    //// boundary check is done YET !!!!!!!!!!!!!!!
    //for (unsigned int i = 0; i < _arcCyclePointTypePtr->size(); i++)
    //{
    //	if ((*_arcCyclePointTypePtr)[i].second)
    //	{// this is an edge type vertex
    //		pot_vertex = inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].v0;
    //		if (i != 0 && pot_vertex == outArcCycle_VertexOnMesh.back())
    //			continue;
    //		findIter = processed_vertex.find(pot_vertex);
    //		if (findIter != processed_vertex.end())
    //		{
    //			pot_vertex = inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].v1;
    //			//if (i != 0 && pot_vertex == outArcCycle_VertexOnMesh.back())
    //			//continue;
    //		}
    //	}
    //	else
    //	{// this is a vertex
    //		pot_vertex = (*_arcCyclePointTypePtr)[i].first;
    //	}
    //	if (i != 0 && pot_vertex == outArcCycle_VertexOnMesh.back())
    //		continue;
    //	// check if this vertex is there or not
    //	// if it is there, ignore it
    //	// otherwise, push it
    //	//
    //	findIter = processed_vertex.find(pot_vertex);
    //	if (i == 0 || pot_vertex != outArcCycle_VertexOnMesh.back())//if (findIter == processed_vertex.end())
    //	{
    //		processed_vertex.insert(pot_vertex);
    //		//
    //		outArcCycle_VertexOnMesh.push_back(pot_vertex);
    //		//
    //		if (i != 0)
    //		{// find the edge connecting (i-1)--(i)
    //			// First, if the (i-1)-vertex is on the edge, as the other vertex of this edge is taken
    //			// so the edge connecting (i-1)--(i) is just the processing edge
    //			if ((*_arcCyclePointTypePtr)[i].second)
    //			{
    //				if (inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].v0 == outArcCycle_VertexOnMesh[outArcCycle_VertexOnMesh.size() - 2] ||
    //					inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].v1 == outArcCycle_VertexOnMesh[outArcCycle_VertexOnMesh.size() - 2] )
    //				{
    //					outArcCycle_EdgeOnMesh.push_back((*_arcCyclePointTypePtr)[i].first);
    //				}
    //				else
    //				{
    //				// find the triangle first containing the edge and prev vertex
    //					int triangle_id = 0;
    //					for (unsigned int nTri = 0; nTri < 1; nTri++)
    //					{
    //						int tot_sum = inMeshPtr->vecTriangle[inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].AdjTri[nTri]].v0 +
    //								inMeshPtr->vecTriangle[inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].AdjTri[nTri]].v1 +
    //								inMeshPtr->vecTriangle[inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].AdjTri[nTri]].v2;
    //						int tot_edge_sum = inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].v0 +
    //								inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].v1;

    //						if (tot_sum - tot_edge_sum == outArcCycle_VertexOnMesh[outArcCycle_VertexOnMesh.size() - 2])
    //						{
    //							triangle_id = inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].AdjTri[nTri];
    //						}
    //						else
    //						{
    //							triangle_id = inMeshPtr->vecEdge[(*_arcCyclePointTypePtr)[i].first].AdjTri[nTri + 1];
    //						}
    //					}
    //					//
    //					int edge_vec[3] = {inMeshPtr->vecTriangle[triangle_id].e01, inMeshPtr->vecTriangle[triangle_id].e02,
    //										inMeshPtr->vecTriangle[triangle_id].e12};
    //					for (int eid = 0; eid < 3; eid++)
    //					{
    //						if ((	inMeshPtr->vecEdge[edge_vec[eid]].v0 == pot_vertex ||
    //								inMeshPtr->vecEdge[edge_vec[eid]].v1 == pot_vertex) &&
    //							(inMeshPtr->vecEdge[edge_vec[eid]].v0 == outArcCycle_VertexOnMesh[outArcCycle_VertexOnMesh.size() - 2] ||
    //							inMeshPtr->vecEdge[edge_vec[eid]].v1 == outArcCycle_VertexOnMesh[outArcCycle_VertexOnMesh.size() - 2]))
    //						{
    //							outArcCycle_EdgeOnMesh.push_back(edge_vec[eid]);
    //							break;
    //						}
    //					}
    //					//
    //				}
    //			}//.second
    //			else
    //			{
    //				// Handle the vertex case
    //				for (unsigned int k = 0; k < inMeshPtr->vecVertex[pot_vertex].adjEdges.size();
    //					k++)
    //				{
    //					if (inMeshPtr->vecEdge[inMeshPtr->vecVertex[pot_vertex].adjEdges[k]].v0 == outArcCycle_VertexOnMesh[outArcCycle_VertexOnMesh.size() - 2] ||
    //						inMeshPtr->vecEdge[inMeshPtr->vecVertex[pot_vertex].adjEdges[k]].v1 == outArcCycle_VertexOnMesh[outArcCycle_VertexOnMesh.size() - 2])
    //					{
    //						outArcCycle_EdgeOnMesh.push_back(inMeshPtr->vecVertex[pot_vertex].adjEdges[k]);
    //						break;
    //					}
    //				}	// for
    //			}
    //		} // i != 0
    //		// else
    //		// do nothing for i == 0
    //	}// if findIter
    //}
    //std::cout << "OUT ComputeArcCycle" << std::endl;
    return;
}

////
//void MapLoopsBackToMeshLevelSetAndArc::ComputeLoopOnMeshPolygon()
//{
//	// preq: all arcs are processed before
//// path is encoded in this way
//// v0--e0--v1--e1--....vn--en--v(n+1)[==v0]
//// each cycle are passing from high to low
//// each edge is taken in this way [a...b)-[b...c)-[c...a)
//	int curArcId = 0;
//	int startNodeId = 0;
//	int endNodeId = 0;
//	Vector3 translatedPt;
//	Vector3 prePt, nxtPt;
//	// walk around the path in terms of simplified arcs
//	for (unsigned int i = 0; i < _meshLoopPtr->simpArcIndexSet.size(); i++)
//	{
//		curArcId = _meshLoopPtr->simpArcIndexSet[i];
//		startNodeId = _meshLoopPtr->nodeIndexSet[i];
//		endNodeId = _meshLoopPtr->nodeIndexSet[i+1];
//		//
//		if ((*_pVecSimplifiedArc)[curArcId].nCriticalNode0 == startNodeId)
//		{// since the arc are traversed from high to low,
//		// need to travel in the opposite direction
//			//
//			for ( int iarc = (*_arcVerOnMeshPtr)[curArcId].size() - 1; iarc > 0 ; iarc--)
//			{// be aware of which arc to take like _offsetPathArcOnMeshPtr or _pathArcOnMeshPtr
//				outArcCycle_VertexOnMesh.push_back((*_arcVerOnMeshPtr)[curArcId][iarc]);
//				outArcCycle_EdgeOnMesh.push_back((*_arcEdgOnMeshPtr)[curArcId][iarc - 1]);
//			}
//		}
//		else
//		{// ignore the last point
//			for (unsigned int iarc = 0; iarc < (*_arcVerOnMeshPtr)[curArcId].size() - 1 ; iarc++)
//			{
//				outArcCycle_VertexOnMesh.push_back((*_arcVerOnMeshPtr)[curArcId][iarc]);
//				outArcCycle_EdgeOnMesh.push_back((*_arcEdgOnMeshPtr)[curArcId][iarc]);
//			}
//		}
//	}
//	// push the first point into the loop to close it
//	outArcCycle_VertexOnMesh.push_back(outArcCycle_VertexOnMesh[0]);
//	//
//	return;
//}
bool MapLoopsBackToMeshLevelSetAndArc::VerifyingEmbeddedPath(std::vector<int> &vertices, std::vector<int> &edges) {
    for (unsigned int i = 0; i < vertices.size() - 1; i++) {
        if (!((inMeshPtr->vecEdge[edges[i]].v0 == vertices[i] &&
               inMeshPtr->vecEdge[edges[i]].v1 == vertices[i + 1]) ||
              (inMeshPtr->vecEdge[edges[i]].v1 == vertices[i] &&
               inMeshPtr->vecEdge[edges[i]].v0 == vertices[i + 1]))
                ) {
            return false;
        }
    }
    return true;
}
