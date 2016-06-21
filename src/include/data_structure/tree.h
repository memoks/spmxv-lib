
#ifndef TREE_H_
#define TREE_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

#include "include/config.h"

struct tree_node
{
	struct tree_node* parent;
	struct tree_node* left;
	struct tree_node* right;
};

typedef struct tree_node tree_node_t;

// Different methods
// (intended to be used with thread_info data structure)
// ------------------------------------------------------------
extern void tree_node_init(tree_node_t* node);
extern void tree_node_terminateSingle(tree_node_t* node);
// ------------------------------------------------------------

extern tree_node_t* tree_node_new(tree_node_t* parent,
		tree_node_t* left, tree_node_t* right);
extern int tree_node_isLeaf(tree_node_t* node);
extern int tree_node_isFull(tree_node_t* node);
extern void tree_node_delete(tree_node_t* head);
extern void tree_node_deleteSingle(tree_node_t* node);
extern int tree_node_getLeafCount(tree_node_t* node);
extern int tree_node_getNodeCount(tree_node_t* node);
extern int tree_node_getParentCount(tree_node_t* node);
extern void tree_node_addLeft(tree_node_t* parent, tree_node_t* left);
extern void tree_node_addRight(tree_node_t* parent, tree_node_t* right);
extern int tree_node_isLeftChild(tree_node_t* node);
extern int tree_node_isRightChild(tree_node_t* node);
extern tree_node_t* tree_node_getSibling(tree_node_t* node);
extern tree_node_t* tree_node_getHead(tree_node_t* node);
extern tree_node_t** tree_node_getLeafs(
		tree_node_t* node, int* leafCount_out);

/**
 * tree_entry - get the struct for this entry
 * @ptr:	the &struct tree pointer.
 * @type:	the type of the struct this is embedded in.
 * @member:	the name of the tree_struct within the struct.
 */
#define tree_entry(ptr, type, member) \
	((type *)((char *)(ptr)-(unsigned long)(&((type *)0)->member)))

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* TREE_H_ */
