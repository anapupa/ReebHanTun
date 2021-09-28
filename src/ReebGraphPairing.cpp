/*
(c) 2012 Fengtao Fan
*/
#include "ReebGraphPairing.h"
#include <vector>
#include <algorithm>
#include <iostream>

namespace RG_PAIR_ALG = ReebGraphPairingAlgorithm;

void RG_PAIR_ALG::ReebGraphPairing::ReportPairs() {
    for (unsigned int i = 0; i < resPairing.size(); i++) {
        std::cout << "< " << resPairing[i].first << " , " << resPairing[i].second << " > " << std::endl;
    }
    return;
}

void RG_PAIR_ALG::ReebGraphPairing::AddMinimumReebNode(const int critNodeId, const double height) {
    BinaryTreeNode *root = new BinaryTreeNode;
    //
    root->critNodeId = critNodeId;
    root->height = height;
    root->paired = false;
    // set link structures
    root->parent = NULL; // as the parent
    root->rightChild = NULL;
    root->leftChild = NULL;
    //
    // add to forest
    forest.push_back(root);
    // add to mapping
    critNodeToTreeNodeMapping[critNodeId] = root;
    return;
}

void RG_PAIR_ALG::ReebGraphPairing::AddUpForkingReebNode(const int critNodeId, const double height,
                                                         const int critNodeId_down) {
    // find the leaf representing the upforking node
    BinaryTreeNode *upforkingParentNode = critNodeToTreeNodeMapping[critNodeId_down];
    // create two new leaves add them this node
    BinaryTreeNode *tmp = new BinaryTreeNode;
    // set data
    tmp->critNodeId = critNodeId;
    tmp->height = height;
    tmp->paired = false;
    // set the links
    tmp->parent = upforkingParentNode;
    if (!upforkingParentNode->leftChild)
        upforkingParentNode->leftChild = tmp;
    else
        upforkingParentNode->rightChild = tmp;
    //
    tmp->leftChild = tmp->rightChild = NULL;

    // add to the mapping
    critNodeToTreeNodeMapping[critNodeId] = tmp;
    //
    return;
}

void RG_PAIR_ALG::ReebGraphPairing::AddDownForkingReebNode(const int critNodeId, const double height,
                                                           const int critNodeId_down_0,
                                                           const int critNodeId_down_1) {
    // there are two leaf nodes corresponding to the downforking node
    // thus, the mapping captures only one of it
    // find the tree nodes corresponding to the two downforked reeb nodes
    bool ancestor_type = true; // true for common ancestor
    // false for different roots on two paths
    BinaryTreeNode *downForkingPtr_0 = critNodeToTreeNodeMapping[critNodeId_down_0];
    BinaryTreeNode *downForkingPtr_1 = critNodeToTreeNodeMapping[critNodeId_down_1];
    //
    /**************************************************/
    // create the node and added it into the tree
    BinaryTreeNode *tmpNewNode = new BinaryTreeNode;
    // set data
    tmpNewNode->critNodeId = critNodeId;
    tmpNewNode->height = height;
    tmpNewNode->paired = true;
    // set the links
    tmpNewNode->parent = NULL;
    if (!downForkingPtr_1->leftChild)
        downForkingPtr_1->leftChild = tmpNewNode;
    else
        downForkingPtr_1->rightChild = tmpNewNode;
    if (!downForkingPtr_0->leftChild)
        downForkingPtr_0->leftChild = tmpNewNode;
    else
        downForkingPtr_0->rightChild = tmpNewNode;
    //
    tmpNewNode->leftChild = NULL;
    tmpNewNode->rightChild = NULL;
    // note that the downforking has two parents
    // so the parent bit is not set here
    // only children pointer are set
    // add to the mapping
    critNodeToTreeNodeMapping[critNodeId] = tmpNewNode;
    //
    /**************************************************/
    std::vector<BinaryTreeNode *> path_0;
    std::vector<BinaryTreeNode *> path_1;
    // push the first element of each path
    path_0.push_back(downForkingPtr_0);
    path_1.push_back(downForkingPtr_1);
    // store the pairing information
    std::pair<int, int> downForkingPair;
    downForkingPair.first = critNodeId;

    //
    do {
        if (path_0.back()->parent == NULL && path_1.back()->parent == NULL) {
            ancestor_type = false;
            //
            if (path_0.back()->height < path_1.back()->height)
                downForkingPair.second = path_1.back()->critNodeId;
            else
                downForkingPair.second = path_0.back()->critNodeId;
            // store to the vector of pairing
            resPairing.push_back(downForkingPair);
            break;
        }
        if (path_0.back()->parent == NULL && path_1.back()->height <= path_0.back()->height ||
            path_1.back()->parent == NULL && path_0.back()->height <= path_1.back()->height) {
            ancestor_type = false;
            if (path_0.back()->parent == NULL)
                downForkingPair.second = path_0.back()->critNodeId;
            else
                downForkingPair.second = path_1.back()->critNodeId;
            // store to the vector of pairing
            resPairing.push_back(downForkingPair);
            break;
        }

        // check the different cases to go down
        if (path_0.back()->parent == NULL) {// no need to proceed for path_0
            // by if condition, path_1 has height > path_0's height
            while (path_1.back()->height > path_0.back()->height && path_1.back()->parent) {
                path_1.push_back(path_1.back()->parent);
            }
        } else {
            if (path_1.back()->parent == NULL) {
                while (path_0.back()->height > path_1.back()->height && path_0.back()->parent) {
                    path_0.push_back(path_0.back()->parent);
                }
            } else {// both of them are not roots
                if (path_0.back()->critNodeId == path_1.back()->critNodeId) {
                    ancestor_type = true;
                    downForkingPair.second = path_0.back()->critNodeId;
                    // store to the vector of pairing
                    resPairing.push_back(downForkingPair);
                    break;
                }
                if (path_0.back()->height > path_1.back()->height) {
                    while (path_0.back()->height > path_1.back()->height && path_0.back()->parent) {
                        path_0.push_back(path_0.back()->parent);
                    }
                } else {
                    while (path_1.back()->height > path_0.back()->height && path_1.back()->parent) {
                        path_1.push_back(path_1.back()->parent);
                    }
                }
            }
        }

    } while (1);
    //
    critNodeToTreeNodeMapping[downForkingPair.second]->paired = true;
    // check which case it is
    if (ancestor_type) {
        std::pair<BinaryTreeNode *, BinaryTreeNode *> leftEdge;
        std::pair<BinaryTreeNode *, BinaryTreeNode *> rightEdge;
        // leftEdge and rightEdge have the common piont in higher height
        // leftEdge.first == rightEdge.first
        // glue the left edge and right edge
        leftEdge.first = rightEdge.first = tmpNewNode;
        leftEdge.second = downForkingPtr_0;
        rightEdge.second = downForkingPtr_1;
        //
        while (leftEdge.second->critNodeId != rightEdge.second->critNodeId) {
            if (leftEdge.second->height < rightEdge.second->height) {// swap left and right
                std::pair<BinaryTreeNode *, BinaryTreeNode *> tmpEdge = leftEdge;
                leftEdge = rightEdge;
                rightEdge = tmpEdge;
            }
            //left edge height > right edge height
            // cut right edge
            if (rightEdge.second->leftChild == rightEdge.first) {
                rightEdge.second->leftChild = leftEdge.second;
            } else {
                rightEdge.second->rightChild = leftEdge.second;
            }
            leftEdge.first->parent = leftEdge.second;
            //
            leftEdge.first = rightEdge.first = leftEdge.second;
            //
            leftEdge.second = leftEdge.second->parent;
        }
        // ancestor is true
        leftEdge.first->parent = leftEdge.second;
        leftEdge.second->leftChild = NULL;
    } else {// common ancestor
        std::pair<BinaryTreeNode *, BinaryTreeNode *> leftEdge;
        std::pair<BinaryTreeNode *, BinaryTreeNode *> rightEdge;
        // leftEdge and rightEdge have the common piont in higher height
        // leftEdge.first == rightEdge.first
        // glue the left edge and right edge
        leftEdge.first = rightEdge.first = tmpNewNode;
        leftEdge.second = downForkingPtr_0;
        rightEdge.second = downForkingPtr_1;
        while (leftEdge.second->critNodeId != downForkingPair.second &&
               rightEdge.second->critNodeId != downForkingPair.second) {
            if (leftEdge.second->height < rightEdge.second->height) {// swap left and right
                std::pair<BinaryTreeNode *, BinaryTreeNode *> tmpEdge = leftEdge;
                leftEdge = rightEdge;
                rightEdge = tmpEdge;
            }
            //left edge height > right edge height
            // cut right edge
            if (rightEdge.second->leftChild == rightEdge.first) {
                rightEdge.second->leftChild = leftEdge.second;
            } else {
                rightEdge.second->rightChild = leftEdge.second;
            }
            leftEdge.first->parent = leftEdge.second;
            //
            leftEdge.first = rightEdge.first = leftEdge.second;
            //
            leftEdge.second = leftEdge.second->parent;
        }

        if (leftEdge.second->critNodeId == downForkingPair.second) {//
            std::pair<BinaryTreeNode *, BinaryTreeNode *> tmpEdge = leftEdge;
            leftEdge = rightEdge;
            rightEdge = tmpEdge;
        }
        // rightEdge has end point
        //cut the root out
        if (rightEdge.second->leftChild = rightEdge.first)
            rightEdge.second->leftChild = NULL;
        else
            rightEdge.second->rightChild = NULL;
        //
        rightEdge.first->parent = leftEdge.second;
        //
        while (leftEdge.second->height > rightEdge.second->height) {//
            //
            leftEdge.first = leftEdge.second;
            leftEdge.second = leftEdge.second->parent;
        }
        //
        leftEdge.first->parent = NULL;
        if (leftEdge.second->leftChild == leftEdge.first)
            leftEdge.second->leftChild = rightEdge.second;
        else
            leftEdge.second->rightChild = rightEdge.second;
        //
        leftEdge.first->parent = rightEdge.second;
        rightEdge.second->parent = leftEdge.second;
        rightEdge.second->leftChild = leftEdge.first;
        // remove it from the forest
        std::list<BinaryTreeNode *>::iterator listIter = forest.begin();
        for (; listIter != forest.end(); listIter++) {
            if (*listIter == rightEdge.second)
                break;
        }
        forest.erase(listIter);

    }

    return;
}

void RG_PAIR_ALG::ReebGraphPairing::AddMaximumReebNode(const int critNodeId, const double height,
                                                       const int critNodeId_down) {
    // find the leaf corresponding to this critical node first
    std::pair<int, int> maxPair;
    maxPair.first = critNodeId; // the maximun Reeb node
    /**************************************************/
    // create the node and added it into the tree
    BinaryTreeNode *tmp = new BinaryTreeNode;
    // set data
    tmp->critNodeId = critNodeId;
    tmp->height = height;
    tmp->paired = true;
    // set the links
    BinaryTreeNode *currentParent = critNodeToTreeNodeMapping[critNodeId_down];
    tmp->parent = currentParent;
    if (!currentParent->leftChild)
        currentParent->leftChild = tmp;
    else
        currentParent->rightChild = tmp;
    //
    tmp->leftChild = NULL;
    tmp->rightChild = NULL;
    // add to the mapping
    critNodeToTreeNodeMapping[critNodeId] = tmp;
    //
    /**************************************************/
    BinaryTreeNode *travNodePtr = tmp;//critNodeToTreeNodeMapping[critNodeId];
    BinaryTreeNode *delPtr = NULL;
    while (travNodePtr->paired)
        travNodePtr = travNodePtr->parent;
    travNodePtr->paired = true;
    //
    maxPair.second = travNodePtr->critNodeId; // the paired node
    // store to the vector of pairing
    resPairing.push_back(maxPair);
    return;
}