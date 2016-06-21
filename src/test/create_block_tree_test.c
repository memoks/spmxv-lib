

#include <stdlib.h>
#include <stdio.h>

#include "include/input_parser.h"
#include "include/data_structure/block_tree.h"
#include "include/data_structure/tree.h"
#include "include/scheduler/job_queue.h"


void postOrderTraverse(tree_node_t* head)
{
	if(head == NULL)
		return;

	postOrderTraverse(head->left);
	postOrderTraverse(head->right);

	block_tree_t* block = block_getBlock(head);
	block_print(block);
}
/*
void distributePEs(tree_node_t* head, int* peList, int startInd, int peCount)
{
	if(peCount <= 1)
	{
		printf("pe-%d is assigned an internal node. Child count: %d\n", startInd, block_getChildCount(head));
		block_tree_t* block = block_getBlock(head);
		return;
	}

	int leftPeCount = peCount / 2 + peCount % 2;
	int rightPeCount = peCount / 2;

	if(head->left != NULL)
		distributePEs(head->left, peList, startInd, leftPeCount);
	else
	{
		block_tree_t* block = block_getBlock(head);
		printf("block%d assigned to p%d\n", startInd, peList[startInd]);
	}

	if(head->right != NULL)
		distributePEs(head->right, peList, startInd + leftPeCount, rightPeCount);
	else
	{
		block_tree_t* block = block_getBlock(head);
		printf("block%d assigned to p%d\n", startInd, peList[startInd + leftPeCount]);
	}
}
*/
int main(void)
{
	vector_int_t* blockData = input_readBlockFile("../../input/wheel_mtx/dim_wheel_256K");

	printf("Creating block tree... ");
	block_tree_t* head = block_createBlockTree(blockData, 902103, 723605);
	printf("done\n");

	job_queue_t* jq = job_queue_new(0);
	job_queue_fillRowParallel(jq, &head->node);

	job_queue_print(jq);

	return EXIT_SUCCESS;
}
