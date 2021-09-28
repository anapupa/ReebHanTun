/*
(c) 2012 Fengtao Fan
*/
#ifndef MY_REEB_GRAPH_PAIRING_ALG_H
#define MY_REEB_GRAPH_PAIRING_ALG_H

#include <vector>
#include <map>
#include <list>
#include <queue>
#include <iostream>

namespace ReebGraphPairingAlgorithm {
    // a binary tree is represented as a pointer to the root node
    // T-->[first element];
    // That is, T == NULL is empty;
    // tree root has null parent;
    // leaf has the null leftChild and null rightChild;
    // its parent always has lower height
    struct BinaryTreeNode {
        // link structure
        BinaryTreeNode *parent;
        BinaryTreeNode *leftChild;
        BinaryTreeNode *rightChild;
        // identify the input index
        int critNodeId;
        double height;
        bool paired;
    };

    class ReebGraphPairing {
    public:
        // constructor
        ReebGraphPairing() {
        }

        ReebGraphPairing(const ReebGraphPairing &rhs) {
            // never use it
            std::cout << "NEVER USE copy constructor" << std::endl;
        }

        // assign operator
        ReebGraphPairing &operator=(const ReebGraphPairing &rhs) {
            // never use it
            std::cout << "NEVER USE assignment operator" << std::endl;
            return *this;
        }

        // deconstructor
        ~ReebGraphPairing() {
            Clear();
        }

        void clean_sequence_tree(BinaryTreeNode *root) {//
            if (!root)
                return;
            BinaryTreeNode *leaf = root;
            while (leaf->leftChild || leaf->rightChild) {
                if (leaf->leftChild)
                    leaf = leaf->leftChild;
                else
                    leaf = leaf->rightChild;
            }
            //
            BinaryTreeNode *delPtr = NULL;
            while (leaf) {
                delPtr = leaf;
                leaf = leaf->parent;
                delete delPtr;
                delPtr = NULL;
            }
        }

        void clean_binary_tree(BinaryTreeNode *root) {
            // empty tree
            if (!root)
                return;
            //
            std::queue<BinaryTreeNode *> Q;
            Q.push(root);
            root = NULL;
            while (!Q.empty()) {
                BinaryTreeNode *tmp = Q.front();
                Q.pop();
                //
                //tmp->paired = false;
                //
                if (tmp->leftChild) {
                    tmp->leftChild->parent = NULL;
                    Q.push(tmp->leftChild);/*
					if (!tmp->leftChild->paired)
					{
						std::cout << "this is a loop" << std::endl;
					}*/
                }
                if (tmp->rightChild) {
                    tmp->rightChild->parent = NULL;
                    Q.push(tmp->rightChild);/*
					if (!tmp->rightChild->paired)
					{
						std::cout << "this is a loop" << std::endl;
					}*/
                }
                //
                tmp->parent = NULL;
                tmp->leftChild = NULL;
                tmp->rightChild = NULL;

                delete tmp;
                tmp = NULL;
            }
            //}
            //// clean left child tree
            //clean_binary_tree(root->leftChild);
            //// clean right child tree
            //clean_binary_tree(root->rightChild);
            ////clean the node itself
            //BinaryTreeNode * tmp = root;
            //// set its parent as null pointer
            //if (root->parent)
            //{
            //	//if (root->parent->leftChild == root)
            //	//	root->parent->leftChild = NULL;
            //	//else
            //	//{
            //	//	root->parent->rightChild = NULL;
            //	//}
            //	root->parent = NULL;
            //}
            //// delete data
            //delete tmp;
        }

        void Clear() {// clean all dynamically allocated memory
            //std::cout << "Clearing tree" << std::endl;
            // clean binary trees
            {
                std::list<BinaryTreeNode *> tmpTreeList;
                std::list<BinaryTreeNode *>::iterator listIter;
                for (listIter = forest.begin(); listIter != forest.end(); listIter++) {
                    clean_binary_tree(*listIter);
                    //clean_sequence_tree(*listIter);
                    *listIter = NULL;
                }
                forest.swap(tmpTreeList);
            }
            // clean the mapping
            {
                std::map<int, BinaryTreeNode *> tmpMap;
                std::map<int, BinaryTreeNode *>::iterator mIter;
                for (mIter = critNodeToTreeNodeMapping.begin(); mIter != critNodeToTreeNodeMapping.end(); mIter++)
                    mIter->second = NULL;
                critNodeToTreeNodeMapping.swap(tmpMap);
            }
            // clean the vector data
            {
                std::vector<std::pair<int, int> > tmpPairing;
                resPairing.swap(tmpPairing);
            }
        }

        //
        void AddMinimumReebNode(const int critNodeId, const double height);

        void AddUpForkingReebNode(const int critNodeId, const double height,
                                  const int critNodeId_down);

        void AddDownForkingReebNode(const int critNodeId, const double height,
                                    const int critNodeId_down_0,
                                    const int critNodeId_down_1);

        void AddMaximumReebNode(const int critNodeId, const double height,
                                const int critNodeId_down);

        void ReportPairs();

    public:
        std::list<BinaryTreeNode *> forest; // a list of tree
        std::map<int, BinaryTreeNode *> critNodeToTreeNodeMapping;
        std::vector<std::pair<int, int> > resPairing; // pairing results in a pair of crit index
    };
};

#endif //MY_REEB_GRAPH_PAIRING_ALG_H
