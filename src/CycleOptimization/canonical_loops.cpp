/*
(c) 2012 Fengtao Fan
*/
#include "canonical_loops.h"

#include <algorithm>
#include <queue>
#include <map>
#include <limits>

void canonical_loops_computing::PreprocessingForEdges(std::vector<std::vector<float> > &adjListWithWeight) {
    // edge is represented by an array
    //pair_to_edge_index_mapping.resize(_graphPtr->vecNode.size());
    //
    //int edge_counter = 0;
    //
    //for (unsigned int i = 0; i < _graphPtr->vecNode.size(); i++)
    //{
    //	int eid = 0;
    //	for (std::set<int>::iterator sIter = _graphPtr->vecNode[i].adjList.begin();
    //		sIter != _graphPtr->vecNode[i].adjList.end();
    //		sIter++)
    //	{
    //		if (*sIter > int(i))
    //		{
    //			pair_to_edge_index_mapping[i].push_back(std::pair<int, int>(*sIter, edge_counter++));
    //			edge_weight_vec.push_back(adjListWithWeight[i][eid]);
    //		}
    //		//
    //		eid++;
    //	}
    //	std::vector<bool> tmpVec(pair_to_edge_index_mapping[i].size(), false);
    //	edge_flag_list.push_back(tmpVec);
    //}
    //
    edge_weight_vec.resize(_graphPtr->edge_size);
    edge_flag_list.resize(_graphPtr->edge_size);
    //for (unsigned int i = 0; i < _graphPtr->vecNode.size(); i++)
    //{
    //	int eid = 0;
    //	for (unsigned int eid = 0; eid < _graphPtr->vecNode[i].adjVertexIdVec.size(); eid++)
    //	{
    //		int v_ano_id = _graphPtr->vecNode[i].adjVertexIdVec[eid];
    //		if (v_ano_id > int(i))
    //		{
    //			edge_weight_vec[_graphPtr->vecNode[i].adjEdgeIdVec[eid]] = adjListWithWeight[i][eid];
    //		}
    //		//
    //	}
    //}
    //
    edge_flag_list.assign(_graphPtr->edge_size, false);
    //
    edge_weight_vec[0] = 7;
    edge_weight_vec[1] = 9;
    edge_weight_vec[2] = 14;
    edge_weight_vec[3] = 10;
    edge_weight_vec[4] = 15;
    edge_weight_vec[5] = 11;
    edge_weight_vec[6] = 2;
    edge_weight_vec[7] = 6;
    edge_weight_vec[8] = 9;
    return;
}

void canonical_loops_computing::SetEdgeWeights(std::vector<float> &in_weight) {
    edge_weight_vec.resize(_graphPtr->vecNode.size());
    edge_weight_vec.assign(in_weight.begin(), in_weight.end());
    //
    edge_flag_list.resize(_graphPtr->edge_size);
    edge_flag_list.assign(_graphPtr->edge_size, false);
    return;
}

void canonical_loops_computing::Dijkstra_shortest_path_tree(const int source,
                                                            _simpGraph_vec &graph,
                                                            std::vector<float> &weight,
                                                            std::set<std::pair<float, int>, myFloatIntPairLessThan> &remaining_vertex,
                                                            std::vector<int> &parents_vertex,
                                                            std::vector<int> &parents_edge,
                                                            std::vector<float> &distance) {

    /*initialize*/
    std::set<std::pair<float, int>, myFloatIntPairLessThan>::iterator distIter;
    //
    int vec_index = 0;
    for (std::vector<_simpGraphNode_vec>::iterator graphNodeIter = graph.vecNode.begin();
         graphNodeIter != graph.vecNode.end();
         graphNodeIter++) {
        /*
        initializing the distance
        */
        if (graphNodeIter->selfIndex == source) {
            remaining_vertex.insert(std::pair<float, int>(0.0f, graphNodeIter->selfIndex));
            //
            parents_vertex[vec_index] = -1;
            parents_edge[vec_index] = -1;
            distance[vec_index++] = 0.0;
        } else {//
            remaining_vertex.insert(std::pair<float, int>(std::numeric_limits<float>::max(), graphNodeIter->selfIndex));
            //
            parents_vertex[vec_index] = -1;
            parents_edge[vec_index] = -1;
            distance[vec_index++] = std::numeric_limits<float>::max();
        }
    }
    //
    while (!remaining_vertex.empty()) {
        std::pair<float, int> current_node_pair = *remaining_vertex.begin();
        // check it is stil connected or not
        if (current_node_pair.first == std::numeric_limits<float>::max())
            break;
        // remove it from the remaining vertex set
        remaining_vertex.erase(remaining_vertex.begin());
        //
        int current_node = current_node_pair.second;
        //
        for (unsigned int eid = 0; eid < graph.vecNode[current_node].adjVertexIdVec.size(); eid++) {
            int adjNode = graph.vecNode[current_node].adjVertexIdVec[eid];
            int edge_id = graph.vecNode[current_node].adjEdgeIdVec[eid];
            float alt = current_node_pair.first + weight[edge_id];
            //
            if (alt < distance[adjNode]) {
                //
                std::pair<float, int> tempPair(distance[adjNode], adjNode);
                distIter = remaining_vertex.find(tempPair);
                remaining_vertex.erase(distIter);
                //
                tempPair.first = alt;
                remaining_vertex.insert(tempPair);
                //
                distance[adjNode] = alt;
                //
                parents_vertex[adjNode] = current_node;
                parents_edge[adjNode] = edge_id;
            }
        }
    }

}

//void canonical_loops_computing::Dijkstra_shortest_path_tree(const int source,
//							_simpGraph_vec &graph,
//							std::vector<std::pair<std::set<int>, double>> & path_to_source_in_tree,
//							std::vector<double> & weight,
//							std::set<std::pair<double, int>, myDoubleIntPairLessThan> & remaining_vertex,
//							std::vector<std::pair<int, double>> &parents)
//{
//
//	/*initialize*/
//	std::set<std::pair<double, int>, myDoubleIntPairLessThan>::iterator distIter;
//	//std::vector<double> distance(graph.vecNode.size(), std::numeric_limits<double>::max());
//	//
//	distIter = remaining_vertex.begin();
//	//
//	int vec_index = 0;
//	for (std::vector<_simpGraphNode_vec>::iterator graphNodeIter = graph.vecNode.begin();
//		graphNodeIter != graph.vecNode.end();
//		graphNodeIter++)
//	{
//		/*
//		initializing the distance
//		*/
//		if (graphNodeIter->selfIndex == source)
//		{
//			path_to_source_in_tree[graphNodeIter->selfIndex].second = 0.0;
//			remaining_vertex.insert(std::pair<double, int>(0.0, graphNodeIter->selfIndex));
//		}
//		else
//		{
//			//distance[graphNodeIter->selfIndex] = std::numeric_limits<double>::max();
//			path_to_source_in_tree[graphNodeIter->selfIndex].second = std::numeric_limits<double>::max();
//			remaining_vertex.insert(std::pair<double, int>(std::numeric_limits<double>::max(), graphNodeIter->selfIndex));
//		}
//		// parents is intialized before
//		parents[vec_index++] = std::pair<int, double>(-1, 0.0);
//	}
//	//
//	while (!remaining_vertex.empty())
//	{
//		int cur_node = remaining_vertex.begin()->second;
//		// check it is stil connected or not
//		if (path_to_source_in_tree[cur_node].second == std::numeric_limits<double>::max())
//			break;
//		// remove it from the remaining vertex set
//		remaining_vertex.erase(remaining_vertex.begin());
//		//
//		for (unsigned int eid = 0; eid < graph.vecNode[cur_node].adjVertexIdVec.size(); eid++)
//		{
//			int v_ano_id = graph.vecNode[cur_node].adjVertexIdVec[eid];
//			int edge_id = graph.vecNode[cur_node].adjEdgeIdVec[eid];
//			double alt = path_to_source_in_tree[cur_node].second + weight[edge_id];
//			//
//			if (alt < path_to_source_in_tree[v_ano_id].second)
//			{
//				//
//				std::pair<double, int> tempPair(path_to_source_in_tree[v_ano_id].second, v_ano_id);
//				distIter = remaining_vertex.find(tempPair);
//				remaining_vertex.erase(distIter);
//				//
//				tempPair.first = alt;
//				remaining_vertex.insert(tempPair);
//				//
//				path_to_source_in_tree[v_ano_id].second = alt;
//				//
//				parents[v_ano_id].first = cur_node;
//				parents[v_ano_id].second = weight[edge_id];
//			}
//		}
//	}
//
//}
void canonical_loops_computing::compute_canonical_loops_fibo(std::vector<std::set<int> > &canonical_loops,
                                                             std::vector<int> &non_tree_edges,
                                                             std::vector<Annotation_Type> &canonical_loops_annotation,
                                                             std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                                             edge_annotation_computing &edge_annotations_proxy) {
    std::vector<int> parents_vertex(_graphPtr->vecNode.size());
    std::vector<int> parents_edge(_graphPtr->vecNode.size());
    std::vector<float> distance(_graphPtr->vecNode.size());
    std::vector<unsigned char> colors(_graphPtr->vecNode.size());
    std::vector<Annotation_Type> accumulated_parent_annotations(_graphPtr->vecNode.size(),
                                                                Annotation_Type(edge_annotations_proxy.vec_size, 0));
    //
    std::vector<FIBO::Node *> vecNode;
    std::vector<FIBO::Edge *> vecEdge;
    // create vecNode
    for (int i = 0; i < (*_graphPtr).vecNode.size(); i++) {
        FIBO::Node *tempVer = new FIBO::Node(i, std::numeric_limits<float>::max());
        vecNode.push_back(tempVer);
    }
    // create vecEdge
    for (unsigned int i = 0; i < (*_meshPtr).vecEdge.size(); i++) {
        const int v0 = (*_meshPtr).vecEdge[i].v0;
        const int v1 = (*_meshPtr).vecEdge[i].v1;
        //
        FIBO::Edge *tempEdge = new FIBO::Edge(vecNode[v0], vecNode[v1], edge_weight_vec[i]);
        tempEdge->head->addIncomingEdge(tempEdge);
        tempEdge->tail->addOutgoingEdge(tempEdge);
        tempEdge->edge_idx = i;
        vecEdge.push_back(tempEdge);

        tempEdge = new FIBO::Edge(vecNode[v1], vecNode[v0], edge_weight_vec[i]);
        tempEdge->head->addIncomingEdge(tempEdge);
        tempEdge->tail->addOutgoingEdge(tempEdge);
        tempEdge->edge_idx = i;
        vecEdge.push_back(tempEdge);
    }
    //
    enum FIBO::State unTouchedState, scannedState; // used for Dijkstra
    // set up the appropriate states
    unTouchedState = FIBO::UNLABELED;
    scannedState = FIBO::SCANNED;
    //
    for (unsigned int i = 0; i < base_pt_index.size(); i++)//base_pt_index.size(); i++)
    {
        //
        int source_node_id = base_pt_index[i];
        // parents,  path_to_source_in_tree[i].second and remaining_vertex will be intitalized in Dijkstra
        for (int vid = 0; vid < parents_edge.size(); vid++) {
            //parents_edge[vid] = -1;
            //parents_vertex[vid] = -1;
            //vecNode[vid]->key = std::numeric_limits<float>::max();
            vecNode[vid]->state = unTouchedState;
        }
        DijkstraAlgorithm::DijkstraShortestPath(
                source_node_id,
                vecNode,
                unTouchedState,
                scannedState,
                parents_vertex,
                parents_edge,
                distance);
        //

        //
        for (unsigned int vid = 0; vid < parents_edge.size(); vid++) {
            if (parents_edge[vid] > -1) {
                edge_flag_list[parents_edge[vid]] = true;
            }
            //
            colors[vid] = 0;
        }
        //
        collect_non_tree_edge_cycles(source_node_id,
                                     parents_vertex,
                                     parents_edge,
                                     distance,
                                     colors,
                                     accumulated_parent_annotations,
                                     canonical_loops,
                                     non_tree_edges,
                                     canonical_loops_annotation,
                                     sorted_loops,
                                     edge_annotations_proxy);
        //
    }
}

void canonical_loops_computing::compute_canonical_loops(std::vector<std::set<int> > &canonical_loops,
                                                        std::vector<Annotation_Type> &canonical_loops_annotation,
                                                        std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                                        edge_annotation_computing &edge_annotations_proxy) {
    std::set<std::pair<float, int>, myFloatIntPairLessThan> remaining_vertex;
    std::vector<int> parents_vertex(_graphPtr->vecNode.size());
    std::vector<int> parents_edge(_graphPtr->vecNode.size());
    std::vector<float> distance(_graphPtr->vecNode.size());
    std::vector<unsigned char> colors(_graphPtr->vecNode.size());
    std::vector<Annotation_Type> accumulated_parent_annotations(_graphPtr->vecNode.size(),
                                                                Annotation_Type(edge_annotations_proxy.vec_size, 0));
    //
    for (unsigned int i = 0; i < base_pt_index.size(); i++)//base_pt_index.size(); i++)
    {
        //
        int source_node_id = base_pt_index[i];
        // parents,  path_to_source_in_tree[i].second and remaining_vertex will be intitalized in Dijkstra
        Dijkstra_shortest_path_tree(source_node_id,
                                    *_graphPtr,
                                    edge_weight_vec,
                                    remaining_vertex,
                                    parents_vertex,
                                    parents_edge,
                                    distance);
        //

        //
        for (unsigned int i = 0; i < parents_edge.size(); i++) {
            if (parents_edge[i] > -1) {
                edge_flag_list[parents_edge[i]] = true;
            }
            //
            colors[i] = 0;
            //

        }
        //
        /*
        collect_non_tree_edge_cycles(source_node_id,
                                    parents_vertex,
                                    parents_edge,
                                    distance,
                                    colors,
                                    accumulated_parent_annotations,
                                    canonical_loops,
                                    canonical_loops_annotation,
                                    sorted_loops,
                                    edge_annotations_proxy);
        *///
    }
}

//void canonical_loops_computing::compute_canonical_loops(std::vector<std::set<int>> &canonical_loops,
//										std::vector<Annotation_Type> &canonical_loops_annotation,
//										std::set<std::pair<double, int>, myDoubleIntPairLessThan> &sorted_loops,
//										edge_annotation_computing &edge_annotations_proxy)
//{
//	std::set<std::pair<double, int>, myDoubleIntPairLessThan> remaining_vertex;
//	std::vector<std::pair<int, double>> parents(_graphPtr->vecNode.size());
//	//
//	std::vector<std::pair<std::set<int>, double>> path_to_source_in_tree;
//	path_to_source_in_tree.resize(_graphPtr->vecNode.size());
//	//
//	std::vector<std::vector<Annotation_Type>> accumulated_parent_annotations;
//	//
//	for (unsigned int i = 0; i < base_pt_index.size(); i++)//base_pt_index.size(); i++)
//	{
//		//
//		int source_node_id = base_pt_index[i];
//		// parents,  path_to_source_in_tree[i].second and remaining_vertex will be intitalized in Dijkstra
//		Dijkstra_shortest_path_tree(source_node_id,
//									*_graphPtr,
//									path_to_source_in_tree,
//									edge_weight_vec,
//									remaining_vertex,
//									parents);
//		//
//
//		//
//		for (unsigned int i = 0; i < _graphPtr->vecNode.size(); i++)
//		{
//			path_to_source_in_tree[i].first.clear();
//			int current  = i;
//			while (current != -1)
//			{
//				int v_id = current < parents[current].first ? current : parents[current].first;
//				int the_other_v_id = parents[current].first + current - v_id;
//				if (v_id == -1)
//					break;
//				for (unsigned int eid = 0; eid < _graphPtr->vecNode[v_id].adjVertexIdVec.size(); eid++)
//				{
//					int current_ano_v_id = _graphPtr->vecNode[v_id].adjVertexIdVec[eid];
//					if (current_ano_v_id == the_other_v_id)
//					{
//						int edge_id = _graphPtr->vecNode[v_id].adjEdgeIdVec[eid];
//						edge_flag_list[edge_id] = true;
//						//
//						//path_to_source_in_tree[i].second += parents[current].second;
//						path_to_source_in_tree[i].first.insert(edge_id);
//					}
//				}
//				//
//				// path_to_source_in_tree[i].push_back(current);
//				current = parents[current].first;
//			}
//		}
//		//
//		collect_non_tree_edge_cycles(path_to_source_in_tree,
//								canonical_loops,
//								canonical_loops_annotation,
//								sorted_loops,
//								edge_annotations_proxy);
//	}
//}
//bool canonical_loops_computing::is_cycle_nontrivial(std::set<int> &cycle, Annotation_Type &cur_edge_annotation)
//{
//	return true;
//}
void canonical_loops_computing::collect_non_tree_edge_cycles(
        std::vector<std::pair<std::set<int>, float> > &path_to_source_in_tree,
        std::vector<std::set<int> > &canonical_loops,
        std::vector<Annotation_Type> &canonical_loops_annotation,
        std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
        edge_annotation_computing &edge_annotations_proxy) {

    // for each non-tree edge output the cycle
    //
    for (unsigned int vid = 0; vid < _graphPtr->vecNode.size(); vid++) {
        //
        for (unsigned int eid = 0; eid < _graphPtr->vecNode[vid].adjVertexIdVec.size(); eid++) {//
            int ano_node_id = _graphPtr->vecNode[vid].adjVertexIdVec[eid];
            if (ano_node_id > int(vid)) {
                int edge_id = _graphPtr->vecNode[vid].adjEdgeIdVec[eid];
                if (!edge_flag_list[edge_id]) {
                    // non-tree edge
                    //
                    std::vector<int> tempCycle(path_to_source_in_tree[ano_node_id].first.begin(),
                                               path_to_source_in_tree[ano_node_id].first.end());
                    tempCycle.push_back(edge_id);
                    //
                    std::copy(path_to_source_in_tree[vid].first.begin(),
                              path_to_source_in_tree[vid].first.end(), std::inserter(tempCycle, tempCycle.end()));
                    //
                    Annotation_Type cur_edge_anno;
                    if (!edge_annotations_proxy.CheckZeroAnnotationsAndReturnAnnotation(tempCycle,
                                                                                        cur_edge_anno)) {// only collect nontrivial cycle
                        // need to unique the cycle
                        float cycle_weight = path_to_source_in_tree[ano_node_id].second +
                                             path_to_source_in_tree[vid].second + edge_weight_vec[edge_id];
                        //
                        std::vector<int> intersection(path_to_source_in_tree[vid].first.size());
                        std::vector<int>::iterator vecIter;
                        vecIter = std::set_intersection(path_to_source_in_tree[vid].first.begin(),
                                                        path_to_source_in_tree[vid].first.end(),
                                                        path_to_source_in_tree[ano_node_id].first.begin(),
                                                        path_to_source_in_tree[ano_node_id].first.end(),
                                                        intersection.begin());
                        std::set<int> clean_set(tempCycle.begin(), tempCycle.end());
                        if (vecIter - intersection.begin() > 0) {// remove extra copies
                            for (std::vector<int>::iterator it = intersection.begin();
                                 it != vecIter;
                                 it++) {
                                clean_set.erase(clean_set.find(*it));
                                cycle_weight -= 2 * edge_weight_vec[*it];
                            }
                        }
                        //
                        sorted_loops.insert(std::pair<float, int>(cycle_weight, (int) canonical_loops.size()));
                        canonical_loops.push_back(clean_set);
                        canonical_loops_annotation.push_back(cur_edge_anno);
                    }
                } else {// reset the edge flag bits
                    edge_flag_list[edge_id] = false;
                }
            }
        }
    }
}

void canonical_loops_computing::ResetBaseVertices(std::vector<std::vector<int> > &basis_loops,
                                                  edge_annotation_computing &edge_annotations_proxy) {
    std::vector<int> parents_vertex(_graphPtr->vecNode.size());
    std::vector<int> parents_edge(_graphPtr->vecNode.size());
    std::vector<int> distance(_graphPtr->vecNode.size());
    std::vector<unsigned char> colors(_graphPtr->vecNode.size());
    std::vector<Annotation_Type> accumulated_parent_annotations(_graphPtr->vecNode.size(),
                                                                Annotation_Type(edge_annotations_proxy.vec_size, 0));
    //
    std::vector<FIBO::Node *> vecNode;
    std::vector<FIBO::Edge *> vecEdge;
    // create vecNode
    for (int i = 0; i < (*_graphPtr).vecNode.size(); i++) {
        FIBO::Node *tempVer = new FIBO::Node(i, std::numeric_limits<int>::max());
        vecNode.push_back(tempVer);
    }
    // create vecEdge
    for (unsigned int i = 0; i < (*_meshPtr).vecEdge.size(); i++) {
        const int v0 = (*_meshPtr).vecEdge[i].v0;
        const int v1 = (*_meshPtr).vecEdge[i].v1;
        //
        FIBO::Edge *tempEdge = new FIBO::Edge(vecNode[v0], vecNode[v1], 1);
        tempEdge->head->addIncomingEdge(tempEdge);
        tempEdge->tail->addOutgoingEdge(tempEdge);
        tempEdge->edge_idx = i;
        vecEdge.push_back(tempEdge);

        tempEdge = new FIBO::Edge(vecNode[v1], vecNode[v0], 1);
        tempEdge->head->addIncomingEdge(tempEdge);
        tempEdge->tail->addOutgoingEdge(tempEdge);
        tempEdge->edge_idx = i;
        vecEdge.push_back(tempEdge);
    }
    //
    enum FIBO::State unTouchedState, scannedState; // used for Dijkstra
    // set up the appropriate states
    unTouchedState = FIBO::UNLABELED;
    scannedState = FIBO::SCANNED;
    //
    Annotation_Type org_loop_annotation;
    int org_loop_size = 0;
    std::set<int> base_pts_unique;
    for (unsigned int i = 0; i < basis_loops.size(); i++)//base_pt_index.size(); i++)
    {
        //
        org_loop_size = basis_loops[i].size();
        int source_node_id = (*_meshPtr).vecEdge[basis_loops[i][basis_loops[i].size() / 2]].v0;
        // parents,  path_to_source_in_tree[i].second and remaining_vertex will be intitalized in Dijkstra
        for (int vid = 0; vid < parents_edge.size(); vid++) {
            //parents_edge[vid] = -1;
            //parents_vertex[vid] = -1;
            //vecNode[vid]->key = std::numeric_limits<float>::max();
            vecNode[vid]->state = unTouchedState;
        }
        DijkstraAlgorithm::DijkstraShortestPath(
                source_node_id,
                vecNode,
                unTouchedState,
                scannedState,
                parents_vertex,
                parents_edge,
                distance);
        //

        //
        for (unsigned int vid = 0; vid < parents_edge.size(); vid++) {
            if (parents_edge[vid] > -1) {
                edge_flag_list[parents_edge[vid]] = true;
            }
            //
            colors[vid] = 0;
        }
        //
        edge_annotations_proxy.CheckZeroAnnotationsAndReturnAnnotation(basis_loops[i], org_loop_annotation);
        if (!Find_shortes_cycle_with_same_annotation(source_node_id,
                                                     parents_vertex,
                                                     parents_edge,
                                                     distance,
                                                     colors,
                                                     accumulated_parent_annotations,
                                                     org_loop_annotation,
                                                     org_loop_size,
                                                     edge_annotations_proxy,
                                                     base_pts_unique)) {
            for (std::vector<int>::iterator sIter = basis_loops[i].begin();
                 sIter != basis_loops[i].end();
                 sIter++) {
                base_pts_unique.insert((*_meshPtr).vecEdge[*sIter].v0);
                base_pts_unique.insert((*_meshPtr).vecEdge[*sIter].v1);
            }
        }
        //
    }
    //
    if (base_pt_index.size() > base_pts_unique.size()) {
        base_pt_index.clear();
        for (std::set<int>::iterator sIter = base_pts_unique.begin();
             sIter != base_pts_unique.end();
             sIter++) {
            base_pt_index.push_back(*sIter);
        }
    }
    return;
}

void canonical_loops_computing::ExtractShortestBasisFromCycleSets(const int genus,
                                                                  std::vector<int> &shortest_basis_indices,
                                                                  std::vector<Annotation_Type> &canonical_loops_annotation,
                                                                  std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
                                                                  edge_annotation_computing &edge_annotations_proxy) {
    std::vector<Annotation_Type> basis_annotations;
    for (std::set<std::pair<float, int>, myFloatIntPairLessThan>::iterator sIter = sorted_loops.begin();
         sIter != sorted_loops.end();
         sIter++) {
        basis_annotations.push_back(canonical_loops_annotation[sIter->second]);
        //
        if (!edge_annotations_proxy.GaussianElimination(
                basis_annotations)) {// this vector is independent from previous the basis
            shortest_basis_indices.push_back(sIter->second);
            if (shortest_basis_indices.size() == 2 * genus)
                break;
        } else {
            basis_annotations.pop_back();
        }
    }
    return;
}

bool canonical_loops_computing::Find_shortes_cycle_with_same_annotation(
        const int src_node,
        std::vector<int> &parents_vertex,
        std::vector<int> &parents_edge,
        std::vector<int> &distance,
        std::vector<unsigned char> &colors,
        std::vector<Annotation_Type> &accumulated_parent_annotations,
        const Annotation_Type &org_loop_annotation,
        const int org_loop_weight,
        edge_annotation_computing &edge_annotations_proxy,
        std::set<int> &base_pt_set_unique) {
    bool findNewCycle = false;
    //
    // accumulate the edge annotations along the path from each node to the source node
    //std::vector<Annotation_Type> accumulated_parent_annotations(_graphPtr->vecNode.size(), Annotation_Type(edge_annotations_proxy.vec_size, 0));
    std::queue<int> Q;
    //std::vector<bool> colors(parents_vertex.size(), false);
    Q.push(src_node);
    //
    colors[src_node] = 1;
    //
    for (unsigned int i = 0; i < accumulated_parent_annotations[src_node].size(); i++)
        accumulated_parent_annotations[src_node][i] = 0;
    //
    while (!Q.empty()) {
        int current_node = Q.front();
        Q.pop();
        //
        for (unsigned int i = 0; i < (*_graphPtr).vecNode[current_node].adjVertexIdVec.size(); i++) {
            int adj_vertex_idx = (*_graphPtr).vecNode[current_node].adjVertexIdVec[i];
            int adj_edge_idx = (*_graphPtr).vecNode[current_node].adjEdgeIdVec[i];
            if (edge_flag_list[adj_edge_idx]) {// only visit tree edges
                if (!colors[adj_vertex_idx]) {// it is not visited
                    edge_annotations_proxy.Add(edge_annotations_proxy.edge_annotations[adj_edge_idx],
                                               accumulated_parent_annotations[current_node],
                                               accumulated_parent_annotations[adj_vertex_idx]); //
                    colors[adj_vertex_idx] = 1;
                    Q.push(adj_vertex_idx);
                }
            }
        }
        //
    }
    //
    // for each non-tree edge output the cycle
    //
    std::pair<std::pair<int, int>, int> shortest_loop_mark;
    int shortest_loop_weight = org_loop_weight;
    Annotation_Type cycleAnnotation(edge_annotations_proxy.vec_size, 0);
    for (unsigned int eid = 0; eid < edge_flag_list.size(); eid++) {
        if (!edge_flag_list[eid]) {
            //
            int left_node = (*_meshPtr).vecEdge[eid].v0;
            int right_node = (*_meshPtr).vecEdge[eid].v1;
            //
            edge_annotations_proxy.Add(accumulated_parent_annotations[left_node],
                                       accumulated_parent_annotations[right_node],
                                       cycleAnnotation);
            edge_annotations_proxy.Add(edge_annotations_proxy.edge_annotations[eid],
                                       cycleAnnotation,
                                       cycleAnnotation);
            //
            if (edge_annotations_proxy.IsEqual(cycleAnnotation, org_loop_annotation)) {
                int weight = 1 + distance[left_node] + distance[right_node];
                if (shortest_loop_weight > weight) {// update
                    shortest_loop_weight = weight;
                    //
                    shortest_loop_mark = std::pair<std::pair<int, int>, int>(std::pair<int, int>(left_node, right_node),
                                                                             eid);
                    findNewCycle = true;
                }
            }
        } else {
            edge_flag_list[eid] = false;
        }
    }
    //
    //std::set<int> tempCycle;
    int left_node = shortest_loop_mark.first.first;
    int right_node = shortest_loop_mark.first.second;
    int eid = shortest_loop_mark.second;
    //
    base_pt_set_unique.insert(src_node);
    while (parents_vertex[left_node] != -1) {
        base_pt_set_unique.insert(left_node);
        left_node = parents_vertex[left_node];
    }
    while (parents_vertex[right_node] != -1) {
        base_pt_set_unique.insert(right_node);
        right_node = parents_vertex[right_node];
    }
    //
    return findNewCycle;
}

void canonical_loops_computing::collect_non_tree_edge_cycles(
        const int src_node,
        std::vector<int> &parents_vertex,
        std::vector<int> &parents_edge,
        std::vector<float> &distance,
        std::vector<unsigned char> &colors,
        std::vector<Annotation_Type> &accumulated_parent_annotations,
        std::vector<std::set<int> > &canonical_loops,
        std::vector<int> &non_tree_edges,
        std::vector<Annotation_Type> &canonical_loops_annotation,
        std::set<std::pair<float, int>, myFloatIntPairLessThan> &sorted_loops,
        edge_annotation_computing &edge_annotations_proxy) {
    //
    // accumulate the edge annotations along the path from each node to the source node
    //std::vector<Annotation_Type> accumulated_parent_annotations(_graphPtr->vecNode.size(), Annotation_Type(edge_annotations_proxy.vec_size, 0));
    std::queue<int> Q;
    //std::vector<bool> colors(parents_vertex.size(), false);
    Q.push(src_node);
    //
    colors[src_node] = 1;
    //
    for (unsigned int i = 0; i < accumulated_parent_annotations[src_node].size(); i++)
        accumulated_parent_annotations[src_node][i] = 0;
    //
    while (!Q.empty()) {
        int current_node = Q.front();
        Q.pop();
        //
        for (unsigned int i = 0; i < (*_graphPtr).vecNode[current_node].adjVertexIdVec.size(); i++) {
            int adj_vertex_idx = (*_graphPtr).vecNode[current_node].adjVertexIdVec[i];
            int adj_edge_idx = (*_graphPtr).vecNode[current_node].adjEdgeIdVec[i];
            if (edge_flag_list[adj_edge_idx]) {// only visit tree edges
                if (!colors[adj_vertex_idx]) {// it is not visited
                    edge_annotations_proxy.Add(edge_annotations_proxy.edge_annotations[adj_edge_idx],
                                               accumulated_parent_annotations[current_node],
                                               accumulated_parent_annotations[adj_vertex_idx]); //
                    colors[adj_vertex_idx] = 1;
                    Q.push(adj_vertex_idx);
                }
            }
        }
        //
    }
    //
    // for each non-tree edge output the cycle
    //
    std::vector<std::pair<std::pair<int, int>, int> > ver_canonical_loops;
    std::vector<Annotation_Type> ver_canonical_loops_annotation;
    //std::set<std::pair<float, int>, myFloatIntPairLessThan> ver_sorted_loops;
    std::vector<std::pair<float, int> > ver_sorted_loops;
    std::map<Annotation_Type, std::pair<float, int>, Annotation_Type_LessThan> uniqueBitMap;
    std::map<Annotation_Type, std::pair<float, int>, Annotation_Type_LessThan>::iterator mIter;
    Annotation_Type cycleAnnotation(edge_annotations_proxy.vec_size, 0);
    for (unsigned int eid = 0; eid < edge_flag_list.size(); eid++) {
        if (!edge_flag_list[eid]) {
            //
            int left_node = (*_meshPtr).vecEdge[eid].v0;
            int right_node = (*_meshPtr).vecEdge[eid].v1;
            //
            edge_annotations_proxy.Add(accumulated_parent_annotations[left_node],
                                       accumulated_parent_annotations[right_node],
                                       cycleAnnotation);
            edge_annotations_proxy.Add(edge_annotations_proxy.edge_annotations[eid],
                                       cycleAnnotation,
                                       cycleAnnotation);
            //
            if (!edge_annotations_proxy.IsZeroAnnotation(cycleAnnotation)) {
                float weight = edge_weight_vec[eid] + distance[left_node] + distance[right_node];
                mIter = uniqueBitMap.find(cycleAnnotation);
                if (mIter == uniqueBitMap.end()) {
                    ver_canonical_loops.push_back(
                            std::pair<std::pair<int, int>, int>(std::pair<int, int>(left_node, right_node), eid));
                    ver_canonical_loops_annotation.push_back(cycleAnnotation);
                    ver_sorted_loops.push_back(std::pair<float, int>(weight, (int) ver_canonical_loops.size() - 1));
                    //ver_sorted_loops.insert(std::pair<float, int>(weight, (int)ver_canonical_loops.size() - 1));
                    //
                    uniqueBitMap[cycleAnnotation] = std::pair<float, int>(weight, (int) ver_canonical_loops.size() - 1);
                } else {
                    if (mIter->second.first > weight) {// update
                        //std::set<std::pair<float, int>, myFloatIntPairLessThan>::iterator sIter = ver_sorted_loops.find(mIter->second);
                        //ver_sorted_loops.erase(sIter);
                        //ver_sorted_loops.insert(std::pair<float, int>(weight, mIter->second.second));
                        ver_sorted_loops[mIter->second.second] = std::pair<float, int>(weight, mIter->second.second);
                        //
                        ver_canonical_loops[mIter->second.second] = std::pair<std::pair<int, int>, int>(
                                std::pair<int, int>(left_node, right_node), eid);
                        //
                        mIter->second = ver_sorted_loops[mIter->second.second];//std::pair<float, int>(weight, mIter->second.second);
                    }
                }
            }
        } else {
            edge_flag_list[eid] = false;
        }
    }
    //
    // extract the shortest basis
    std::vector<int> ver_shortest_basis;
    std::set<std::pair<float, int>, myFloatIntPairLessThan> ver_sorted_loops_sorted;
    for (unsigned int i = 0; i < ver_sorted_loops.size(); i++) {
        ver_sorted_loops_sorted.insert(ver_sorted_loops[i]);
        //ver_shortest_basis.push_back(i);
    }
    ExtractShortestBasisFromCycleSets(_genus, ver_shortest_basis, ver_canonical_loops_annotation,
                                      ver_sorted_loops_sorted, edge_annotations_proxy);

    // store these cycles in the shortest basis
    //
    for (std::vector<int>::iterator vIter = ver_shortest_basis.begin();
         vIter != ver_shortest_basis.end(); vIter++) {
        //
        std::set<int> tempCycle;
        int left_node = ver_canonical_loops[*vIter].first.first;
        int right_node = ver_canonical_loops[*vIter].first.second;
        int eid = ver_canonical_loops[*vIter].second;
        //
        float cycle_weight = edge_weight_vec[eid] + distance[left_node] + distance[right_node];
        //
        tempCycle.insert(eid);
        while (parents_vertex[left_node] != -1) {
            tempCycle.insert(parents_edge[left_node]);
            left_node = parents_vertex[left_node];
        }
        while (parents_vertex[right_node] != -1) {
            tempCycle.insert(parents_edge[right_node]);
            right_node = parents_vertex[right_node];
        }
        //
        sorted_loops.insert(std::pair<float, int>(cycle_weight, (int) canonical_loops.size()));
        canonical_loops_annotation.push_back(ver_canonical_loops_annotation[*vIter]);
        canonical_loops.push_back(tempCycle);
        non_tree_edges.push_back(eid);
        //
    }
    return;
}
