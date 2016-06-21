
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

#include "include/config.h"
#include "include/data_structure/comm_tree.h"

// Helper Functions
// --------------------------------------------------------------------------------------------------

extern void __comm_tree_delete(tree_node_t* head);
extern comm_tree_t* __comm_tree_createCommTree(int numBlocks);

// --------------------------------------------------------------------------------------------------

comm_tree_t* comm_tree_new(int procId)
{
	comm_tree_t* newComm = (comm_tree_t*) malloc(sizeof(comm_tree_t));
	newComm->executionContextId = procId;
	newComm->batchCount = 0;
	newComm->hasWork = TRUE;

	tree_node_init(&newComm->node);

	return newComm;
}

inline void comm_tree_delete(comm_tree_t* head)
{
	__comm_tree_delete(&head->node);
}

extern void comm_tree_deleteSingle(comm_tree_t* commNode)
{
	free(commNode);
}

extern void comm_tree_addLeft(comm_tree_t* parent, comm_tree_t* left)
{
	parent->node.left = &left->node;
	left->node.parent = &parent->node;
}

extern void comm_tree_addRight(comm_tree_t* parent, comm_tree_t* right)
{
	parent->node.right = &right->node;
	right->node.parent = &parent->node;
}

inline comm_tree_t* comm_tree_getLeft(comm_tree_t* parent)
{
	return comm_tree_getComm(parent->node.left);
}

inline comm_tree_t* comm_tree_getRight(comm_tree_t* parent)
{
	return comm_tree_getComm(parent->node.right);
}

inline comm_tree_t* comm_tree_getComm(tree_node_t* node)
{
	if(node == NULL)
		return NULL;

	return tree_entry(node, comm_tree_t, node);
}

inline int comm_tree_isFull(comm_tree_t* commNode)
{
	return tree_node_isFull(&commNode->node);
}

inline int comm_tree_isLeaf(comm_tree_t* commNode)
{
	return tree_node_isLeaf(&commNode->node);
}

int comm_tree_getChildCount(tree_node_t* node)
{
	if(node == NULL)
		return 0;

	if(tree_node_isLeaf(node))
		return 1;

	return comm_tree_getChildCount(node->left) + comm_tree_getChildCount(node->right);
}

inline void comm_tree_print(comm_tree_t* commNode)
{
	PRINTF("ProcId: %d, hasWork: %d, batch count: %d\n", commNode->executionContextId, commNode->hasWork, commNode->batchCount);
}

void comm_tree_printInOrder(tree_node_t* node)
{
	if(node == NULL)
		return;

	comm_tree_printInOrder(node->left);
	comm_tree_print(comm_tree_getComm(node));
	comm_tree_printInOrder(node->right);
}

void comm_tree_printPreOrder(tree_node_t* node)
{
	if(node == NULL)
		return;

	comm_tree_print(comm_tree_getComm(node));
	comm_tree_printPreOrder(node->left);
	comm_tree_printPreOrder(node->right);
}

void comm_tree_printLeafs(tree_node_t* node)
{
	if(node == NULL)
		return;

	if(tree_node_isLeaf(node))
	{
		comm_tree_print(comm_tree_getComm(node));
	}
	else
	{
		comm_tree_printLeafs(node->left);
		comm_tree_printLeafs(node->right);
	}
}

comm_tree_t* comm_tree_createCommTree(int numBlocks)
{
	if(numBlocks <= 0)
		return NULL;

	return __comm_tree_createCommTree(numBlocks);
}

// Helper Functions
// --------------------------------------------------------------------------------------------------

comm_tree_t* __comm_tree_createCommTree(int numBlocks)
{
	if(numBlocks <= 0)
		return NULL;

	comm_tree_t* commTree = comm_tree_new(-1);

	if(numBlocks >= 1)
	{
		comm_tree_addLeft(commTree, __comm_tree_createCommTree(numBlocks / 2 + numBlocks % 2));
		comm_tree_addRight(commTree, __comm_tree_createCommTree(numBlocks / 2));
	}

	return commTree;
}

void __comm_tree_delete(tree_node_t* head)
{
	if(head == NULL)
		return;

	__comm_tree_delete(head->left);
	__comm_tree_delete(head->right);

	comm_tree_t* commNode = comm_tree_getComm(head);
	comm_tree_deleteSingle(commNode);
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
