
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <omp.h>

#include "include/config.h"
#include "include/timer/custom_timer.h"
#include "include/data_structure/tree.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/comm_tree.h"
#include "include/scheduler/tree_sched.h"
#include "include/scheduler/block_info.h"

#include "include/scheduler/tree_scheduler.h"

// TODO comment all

// Some helper functions
// -------------------------------------------------------------------------------------------------------------------
static void __tree_scheduler_traverse(block_info_t* myInfo, comm_tree_t* affNode, int numBlocks);
static void __tree_scheduler_expand(tree_node_t* node, int* currIndex, int* victims);
static tree_node_t* __tree_scheduler_genereateCommTree(
		int leftBlockFirstId, int numBlocks, tree_node_t* parent, int jobBatchCount);
static void __tree_scheduler_printCommTreePreOrder(tree_node_t* head, int* currNodeIndex);
static void __tree_scheduler_deleteCommTree(tree_node_t* head);
// -------------------------------------------------------------------------------------------------------------------

// TODO test
void tree_scheduler_init(block_info_t** blockInfos, int numBlocks, int jobBatchCount)
{
#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("Tree-scheduler working\n");
#endif

	// int blockId;
	// for(blockId = 0; blockId < numBlocks; ++blockId)
	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();
		block_info_t* myBlock = blockInfos[blockId];
		tree_sched_init(myBlock->id, numBlocks, &myBlock->sched.treeSched);
	}

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("Generating comm-tree... ");
#endif
	tree_node_t* commHead = __tree_scheduler_genereateCommTree(0, numBlocks, NULL, jobBatchCount);


#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("done\n");
	PRINTF("Finding dummy leafs... ");
#endif

	int leafCount = 0;
	tree_node_t** commTreeLeafs = tree_node_getLeafs(commHead, &leafCount);


#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("done\n");
#endif

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("done\n");
#endif

	/*
	#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Printing communication tree in pre-order...\n");
		int currNodeIndex = 0;
		tree_node_t* debugHead = tree_node_getHead(&blockInfos[0]->treeSched.affNode);
		__tree_scheduler_printCommTreePreOrder(debugHead, &currNodeIndex);
	#endif
	*/

	// TODO parallel implementation
	// for(blockId = 0; blockId < numBlocks; ++blockId)
	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();
		block_info_t* myBlock = blockInfos[blockId];
		comm_tree_t* affNode = comm_tree_getComm(commTreeLeafs[blockId]);
		__tree_scheduler_traverse(myBlock, affNode, numBlocks);
	}

	free(commTreeLeafs);

	#if PROGRAM_MODE >= TRACE_MODE
		int i;
		int j;
		for(i = 0; i < numBlocks; ++i)
		{
			PRINTF("Block-%d\n", i);
			tree_sched_print(&blockInfos[i]->sched.treeSched);
			PRINTF("\n");
		}
	#endif
}

// TODO test
void tree_scheduler_terminate(block_info_t** blockInfos, int numBlocks)
{
	int i;
	for(i = 0; i < numBlocks; ++i)
	{
		block_info_t* myBlock = blockInfos[i];
		tree_sched_terminate(&myBlock->sched.treeSched);
	}
}

int tree_scheduler_findVictim(block_info_t* myBlock, block_info_t** blockInfos, int numBlocks)
{
	if(myBlock->sched.treeSched.victimIndex >= myBlock->sched.treeSched.victimCount)
		return FALSE;

	myBlock->sched.treeSched.lastVictim = myBlock->sched.treeSched.victims[myBlock->sched.treeSched.victimIndex];
	++myBlock->sched.treeSched.victimIndex;

	return TRUE;
}

// Helper Functions
// ------------------------------------------------------------------------------------------------------

static void __tree_scheduler_traverse(block_info_t* myInfo, comm_tree_t* affNode, int numBlocks)
{
	int currIndex = 0;
	int* victims = myInfo->sched.treeSched.victims;

	tree_node_t* curr = &affNode->node;
	tree_node_t* parent = curr->parent;

	while(parent != NULL)
	{
		if(tree_node_isLeftChild(curr))
			__tree_scheduler_expand(parent->right, &currIndex, victims);
		else
			__tree_scheduler_expand(parent->left, &currIndex, victims);

		curr = parent;
		parent = curr->parent;
	}

	myInfo->sched.treeSched.victims = victims;
}

// TODO test
// TODO might be better if we steal from left most child
static void __tree_scheduler_expand(tree_node_t* node, int* currIndex, int* victims)
{
	if(tree_node_isLeaf(node))
	{
		comm_tree_t* commTree = comm_tree_getComm(node);

		victims[*currIndex] = commTree->executionContextId;
		*currIndex = *currIndex + 1;
	}
	else
	{
		// code for stealing from leftmost child first
		// __tree_scheduler_expand(node->right, currIndex, victims);
		// __tree_scheduler_expand(node->left, currIndex, victims);
		__tree_scheduler_expand(node->right, currIndex, victims);
		__tree_scheduler_expand(node->left, currIndex, victims);
	}
}

// TODO test
static tree_node_t* __tree_scheduler_genereateCommTree(
		int leftBlockFirstId, int numBlocks, tree_node_t* parent, int jobBatchCount)
{
	if(numBlocks < 0)
		return NULL;

	comm_tree_t* commNode = comm_tree_new(COMM_TREE_DEFAULT_EXECUTION_CONTEXT_ID);
	tree_node_t* node = &commNode->node;
	node->parent = parent;

	if(numBlocks > 1)
	{
		// As long as # of blocks is greater than 1, create non-leaf nodes
		int leftBlockCount = numBlocks / 2 + numBlocks % 2;
		int rightBlockCount = numBlocks / 2;

		int leftJobBatchCount = jobBatchCount / 2;
		int rightJobBatchCount = jobBatchCount / 2 + jobBatchCount % 2;

		int rightBlockFirstId = leftBlockFirstId + leftBlockCount;

		node->left = __tree_scheduler_genereateCommTree(
				leftBlockFirstId, leftBlockCount, node, leftJobBatchCount);
		node->right = __tree_scheduler_genereateCommTree(
				rightBlockFirstId, rightBlockCount, node, rightJobBatchCount);
	}
	else // if(numBlocks == 1)
	{
		// Can't go further so just set the id of execution-context.
		commNode->executionContextId = leftBlockFirstId;
	}

	return &commNode->node;
}

// TODO test
static void __tree_scheduler_printCommTreePreOrder(tree_node_t* head, int* currNodeIndex)
{
	if(head == NULL)
		return;

	if(tree_node_isLeaf(head))
	{
		comm_tree_t* commNode = comm_tree_getComm(head);
		comm_tree_print(commNode);
	}
	else
	{
		PRINTF("Node-%d\n", *currNodeIndex);
		++(*currNodeIndex);

		__tree_scheduler_printCommTreePreOrder(head->left, currNodeIndex);
		__tree_scheduler_printCommTreePreOrder(head->right, currNodeIndex);
	}
}

static void __tree_scheduler_deleteCommTree(tree_node_t* head)
{
	if(head == NULL)
		return;

	if(tree_node_isLeaf(head))
		return;

	__tree_scheduler_deleteCommTree(head->left);
	__tree_scheduler_deleteCommTree(head->right);
	tree_node_deleteSingle(head);
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
