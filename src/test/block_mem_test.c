/*

#include <stdlib.h>
#include <stdio.h>

#include "include/data_structure/tree.h"
#include "include/data_structure/block_tree.h"


void inOrderTraverse(tree_node_t* head)
{
	if(head == NULL)
		return;

	inOrderTraverse(head->left);
	inOrderTraverse(head->right);

	block_t* block = block_getBlock(head);
	block_print(block);
}

int main(void)
{
	printf("Creating blocks...\n");
	block_t* head = block_new(0, 0, 0, 0);
	block_print(head);
	printf("testing address calculation...\n");
	block_t* sameHead = block_getBlock(&head->node);
	block_print(sameHead);

	block_t* l1 = block_new(1, 1, 1, 1);
	block_t* l1l = block_new(2, 2, 2, 2);
	block_t* l1r = block_new(3, 3, 3, 3);
	block_t* r1 = block_new(4, 4, 4, 4);



	block_print(l1);
	block_print(l1l);
	block_print(l1r);
	block_print(r1);

	printf("Creating tree...\n");

	block_addRight(head, r1);
	block_addLeft(head, l1);
	block_addLeft(l1, l1l);
	block_addRight(l1, l1r);

	printf("traversing tree...\n");

	inOrderTraverse(&head->node);

	//block_deleteSingle(head);
	block_delete(&head->node);

	return EXIT_SUCCESS;
}
*/
