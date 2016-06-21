/*

#include <stdio.h>
#include <stdlib.h>

#include "include/data_structure/tree.h"
#include "include/data_structure/block_tree.h"


struct demo
{
	index_t start;
	index_t size;
	tree_node_t node;
};

typedef struct demo demo_t;

void fillDemo(demo_t* demo, int startRow, int startCol, int rowCount, int colCount)
{
	demo->start.i = startRow;
	demo->start.j = startCol;
	demo->size.i = rowCount;
	demo->size.j = colCount;
	demo->node.left = NULL;
	demo->node.right = NULL;
	demo->node.parent = NULL;
}

void printDemo(demo_t* demo)
{
	printf("start(%d, %d), size(%d, %d)\n",
			demo->start.i, demo->start.j, demo->size.i, demo->size.j);
}

demo_t* createDemo(void)
{
	demo_t* demo = (demo_t*) malloc(sizeof(demo_t));
	demo->start.i = 12;
	demo->start.j = 12;
	demo->size.i = 123;
	demo->size.j = 123;
	demo->node.left = NULL;
	demo->node.right = NULL;
	demo->node.parent = NULL;

	return demo;
}


int main(void)
{
	demo_t demo;
	fillDemo(&demo, 1, 1, 2, 2);
	printDemo(&demo);

	demo_t* mofo = createDemo();
	printDemo(mofo);

	printf("address of demo= %p\n", &demo);
	printf("address of demo->start= %p\t size: %ld\t offset: %lu\n", &demo.start, sizeof(index_t*), (unsigned long) &((struct block *)0)->start);
	printf("address of demo->size= %p\t size: %ld\t offset: %lu\n", &demo.size, sizeof(index_t*), (unsigned long) &((struct block *)0)->size);
	printf("address of demo->node= %p\t size: %ld\t offset: %lu\n", &demo.node, sizeof(tree_node_t*), (unsigned long) &((struct block *)0)->node);
	printf("address calculated using node: %p\n", (demo_t*) (((char *) &demo.node) - ((unsigned long) &((demo_t*)0)->node)));
	printf("\n");


	printf("address of demo= %p\n", mofo);
	printf("address of demo->start= %p\t size: %ld\t offset: %lu\n", &mofo->start, sizeof(index_t*), (unsigned long) &((struct block *)0)->start);
	printf("address of demo->size= %p\t size: %ld\t offset: %lu\n", &mofo->size, sizeof(index_t*), (unsigned long) &((struct block *)0)->size);
	printf("address of demo->node= %p\t size: %ld\t offset: %lu\n", &mofo->node, sizeof(tree_node_t*), (unsigned long) &((struct block *)0)->node);
	printf("address calculated using node: %p\n", (demo_t*) (((char *) &mofo->node) - ((unsigned long) &((demo_t*)0)->node)));
	printf("\n");

	free(mofo);

	block_t* blockPtr = block_new(0, 0, 1, 1);
	block_print(blockPtr);
	printf("address of blockPtr= %p\n", blockPtr);
	printf("address of blockPtr->start= %p\t size: %ld\t offset: %lu\n", &blockPtr->start, sizeof(index_t*), (unsigned long) &((struct block *)0)->start);
	printf("address of blockPtr->size= %p\t size: %ld\t offset: %lu\n", &blockPtr->size, sizeof(index_t*), (unsigned long) &((struct block *)0)->size);
	printf("address of blockPtr->node= %p\t size: %ld\t offset: %lu\n", &blockPtr->node, sizeof(tree_node_t*), (unsigned long) &((struct block *)0)->node);
	printf("\n");

	printf("address calculated from node: %p\n", block_getBlock(&blockPtr->node));

	block_t* sameBlock = block_getBlock(&blockPtr->node);
	sameBlock->size.i = 10;
	block_print(sameBlock);
	block_print(blockPtr);

	block_deleteSingle(blockPtr);

	return EXIT_SUCCESS;
}
*/
