
#ifndef COMM_TREE_H_
#define COMM_TREE_H_

#include "include/data_structure/tree.h"
#include "include/data_structure/sub_mtx.h"
#include "include/scheduler/job_queue.h"

#define COMM_TREE_DEFAULT_EXECUTION_CONTEXT_ID -1

struct comm_tree
{
	tree_node_t node;

	int executionContextId;
	int hasWork;
	int batchCount;
};

typedef struct comm_tree comm_tree_t;


extern comm_tree_t* comm_tree_new(int procId);
extern void comm_tree_delete(comm_tree_t* head);
extern void comm_tree_deleteSingle(comm_tree_t* commNode);
extern void comm_tree_addLeft(comm_tree_t* parent, comm_tree_t* left);
extern void comm_tree_addRight(comm_tree_t* parent, comm_tree_t* right);
extern comm_tree_t* comm_tree_getLeft(comm_tree_t* parent);
extern comm_tree_t* comm_tree_getRight(comm_tree_t* parent);
extern comm_tree_t* comm_tree_getComm(tree_node_t* node);
extern int comm_tree_isFull(comm_tree_t* commNode);
extern int comm_tree_isLeaf(comm_tree_t* commNode);
extern int comm_tree_getChildCount(tree_node_t* node);
extern void comm_tree_print(comm_tree_t* commNode);
extern void comm_tree_printPreOrder(tree_node_t* node);
extern void comm_tree_printInOrder(tree_node_t* node);
extern void comm_tree_printLeafs(tree_node_t* node);
extern comm_tree_t* comm_tree_createCommTree(int numBlocks);

#endif // COMM_TREE_H_
