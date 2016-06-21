
#include <stdlib.h>
#include <stdio.h>

#include <include/data_structure/tree.h>
/*
void preOrderTraverse(tree_node_t* head)
{
	printf("p\n");

	if(head->left != NULL)
	{
		printf("l");
		preOrderTraverse(head->left);
	}

	if(head->right != NULL)
	{
		printf("r");
		preOrderTraverse(head->right);
	}
}

int main(void)
{
	tree_node_t* head = tree_node_new(NULL, NULL, NULL);
	tree_node_t* l1 = tree_node_new(head, NULL, NULL);
	head->left = l1;

	tree_node_t* l1l = tree_node_new(l1, NULL, NULL);
	l1->left = l1l;
	tree_node_t* l1r = tree_node_new(l1, NULL, NULL);
	l1->right = l1r;

	tree_node_t* r1 = tree_node_new(head, NULL, NULL);
	head->right = r1;


	preOrderTraverse(head);

	tree_node_delete(head);

	return EXIT_SUCCESS;
}
*/
