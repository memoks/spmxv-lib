
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

#include "include/data_structure/tree.h"


// Utility functions
// --------------------------------------------------------------------------------------------------------

extern void __tree_node_getLeafs(tree_node_t* curr, int* currIndex, tree_node_t** leafs_out);

// --------------------------------------------------------------------------------------------------------

void tree_node_init(tree_node_t* node)
{
	node->left = NULL;
	node->right = NULL;
	node->parent = NULL;
}

void tree_node_terminateSingle(tree_node_t* node)
{
	if(node->parent != NULL)
	{
		if(node == node->parent->left)
			node->parent->left = NULL;
		else
			node->parent->right = NULL;
	}

	node->parent = NULL;
}

// ----------------------------------------------------------------------------------------------------------

tree_node_t* tree_node_new(tree_node_t* parent,
		tree_node_t* left, tree_node_t* right)
{
	tree_node_t* newNode = (tree_node_t*) malloc(sizeof(tree_node_t));
	newNode->left = left;
	newNode->right = right;
	newNode->parent = parent;

	return newNode;
}

inline int tree_node_isLeaf(tree_node_t* node)
{
	return node->left == NULL && node->right == NULL;
}

inline int tree_node_isFull(tree_node_t* node)
{
	return node->left != NULL && node->right != NULL;
}

void tree_node_delete(tree_node_t* head)
{
	if(head == NULL)
		return;

	tree_node_delete(head->left);
	tree_node_delete(head->right);

	free(head);
}

void tree_node_deleteSingle(tree_node_t* node)
{
	free(node);
}

int tree_node_getLeafCount(tree_node_t* node)
{
	if(node == NULL)
		return 0;

	if(tree_node_isLeaf(node))
		return 1;

	return tree_node_getLeafCount(node->left) + tree_node_getLeafCount(node->right);
}

int tree_node_getNodeCount(tree_node_t* node)
{
	if(node == NULL)
		return 0;

	return 1 + tree_node_getNodeCount(node->left) + tree_node_getNodeCount(node->right);
}

int tree_node_getParentCount(tree_node_t* node)
{
	int parentCount = 0;
	tree_node_t* parent = node->parent;

	while(parent != NULL)
	{
		++parentCount;
		parent = parent->parent;
	}

	return parentCount;
}

int tree_node_isLeftChild(tree_node_t* node)
{
	if(node->parent == NULL)
		return FALSE;

	if(node->parent->left == node)
		return TRUE;

	return FALSE;
}

int tree_node_isRightChild(tree_node_t* node)
{
	if(node->parent == NULL)
		return FALSE;

	if(node->parent->right == node)
		return TRUE;

	return FALSE;
}

tree_node_t* tree_node_getSibling(tree_node_t* node)
{
	if(node->parent == NULL)
		return NULL;

	if(tree_node_isLeftChild(node))
		return node->parent->right;

	return node->parent->left;
}

tree_node_t** tree_node_getLeafs(tree_node_t* head, int* leafCount_out)
{
	*leafCount_out = tree_node_getLeafCount(head);

	tree_node_t** leafs = malloc(sizeof(tree_node_t*) * (*leafCount_out));
	int currIndex = 0;
	__tree_node_getLeafs(head, &currIndex, leafs);

	return leafs;
}

inline void tree_node_addLeft(tree_node_t* parent, tree_node_t* left)
{
	parent->left = left;
	left->parent = parent;
}

inline void tree_node_addRight(tree_node_t* parent, tree_node_t* right)
{
	parent->right = right;
	right->parent = parent;
}

inline tree_node_t* tree_node_getHead(tree_node_t* node)
{
	if(node == NULL)
		return NULL;

	while(node->parent != NULL)
		node = node->parent;

	return node;
}

// Utility functions
// --------------------------------------------------------------------------------------------------------

void __tree_node_getLeafs(tree_node_t* curr, int* currIndex, tree_node_t** leafs_out)
{
	if(curr == NULL)
		return;

	if(tree_node_isLeaf(curr))
	{
		leafs_out[*currIndex] = curr;
		++(*currIndex);
	}
	else
	{
		__tree_node_getLeafs(curr->left, currIndex, leafs_out);
		__tree_node_getLeafs(curr->right, currIndex, leafs_out);
	}
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
