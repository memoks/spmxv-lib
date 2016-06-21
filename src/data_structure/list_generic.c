
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

#include "include/data_structure/sub_mtx.h"

#include "include/data_structure/list_generic.h"

// List-Node-Generic Functions
// -------------------------------------------------------------------------------------------------------

lng_t* lng_new(void *dataPtr)
{
	lng_t* node = (lng_t*) malloc(sizeof(lng_t));
	INIT_LIST_HEAD(&node->listHead);
	node->dataPtr = dataPtr;

	return node;
}

lng_t* lng_getListGeneric(list_head_t* listHead)
{
	return list_entry(listHead, lng_t, listHead);
}

lng_t* lng_next(lng_t* node)
{
	return list_entry(&node->listHead.next, lng_t, listHead);
}

lng_t* lng_prev(lng_t* node)
{
	return list_entry(&node->listHead.prev, lng_t, listHead);
}

void lng_addTailData(void* dataPtr, lng_t* list)
{
	lng_t* node = lng_new(dataPtr);
	lng_addTail(node, list);
}

void lng_addTail(lng_t* new, lng_t* headNode)
{
	list_add_tail(&new->listHead, &headNode->listHead);
}

void lng_addFrontData(void* dataPtr, lng_t* headNode)
{
	lng_t* new = lng_new(dataPtr);
	lng_addFront(new, headNode);
}

void lng_addFront(lng_t* new, lng_t* headNode)
{
	list_add(&new->listHead, &headNode->listHead);
}

void* lng_remove(lng_t* node)
{
	void* dataPtr = node->dataPtr;
	lng_deleteShallow(node);
	return dataPtr;
}

void lng_deleteShallow(lng_t* node)
{
	list_del(&node->listHead);
	free(node);
}

void lng_deleteAllShallow(lng_t* headNode)
{
	if(headNode == NULL)
		return;

	list_head_t* curr;
	list_head_t* temp;

	list_for_each_safe(curr, temp, &headNode->listHead)
	{
		lng_t* temp = lng_getListGeneric(curr);
		lng_deleteShallow(temp);
	}
}

int lng_isEmpty(lng_t* headNode)
{
	if(headNode == NULL)
		return TRUE;

	return list_empty(&headNode->listHead);
}

int lng_size(lng_t* headNode)
{
	if(headNode == NULL)
		return 0;

	int size = 0;
	list_head_t* curr = NULL;
	list_head_t* temp = NULL;

	list_for_each_safe(curr, temp, &headNode->listHead)
		++size;

	return size;
}

// List-Generic Structure
// -------------------------------------------------------------------------------------------------------

static void __lg_split(lg_t* source, int totalItemCount, int totalSplitCount, lg_t** lists);

// -------------------------------------------------------------------------------------------------------

lg_t* lg_new(void)
{
	lg_t* new = (lg_t*) malloc(sizeof(lg_t));
	new->headNode = lng_new(NULL);
	new->size = 0;

	return new;
}

lng_t* lg_tail(lg_t* list)
{
	return lng_getListGeneric(list->headNode->listHead.prev);
}

lng_t* lg_front(lg_t* list)
{
	return lng_getListGeneric(list->headNode->listHead.next);
}

void lg_addTailData(void* dataPtr, lg_t* list)
{
	lng_addTailData(dataPtr, list->headNode);
	++list->size;
}

void lg_addTail(lng_t* newNode, lg_t* list)
{
	lng_addTail(newNode, list->headNode);
	++list->size;
}

void lg_addFrontData(void* dataPtr, lg_t* list)
{
	lng_addFrontData(dataPtr, list->headNode);
	++list->size;
}

void lg_addFront(lng_t* newNode, lg_t* list)
{
	lng_addFront(newNode, list->headNode);
	++list->size;
}

int lg_size(lg_t* list)
{
	return list->size;
}

int lg_isEmpty(lg_t* list)
{
	return lng_isEmpty(list->headNode);
}

void* lg_remove(lg_t* list, lng_t* node)
{
	--list->size;
	return lng_remove(node);
}

// TODO test remove method changed
void* lg_removeTail(lg_t* list)
{
	--list->size;
	lng_t* tailNode = lng_getListGeneric(list->headNode->listHead.prev);
	return lg_remove(list, tailNode);
}

// TODO test remove method changed
void* lg_removeFront(lg_t* list)
{
	--list->size;
	lng_t* frontNode = lng_getListGeneric(list->headNode->listHead.next);
	return lg_remove(list, frontNode);
}

void lg_removeAll(lg_t* list)
{
	list->size = 0;
	lng_deleteAllShallow(list->headNode);
}

void lg_deleteShallow(lg_t* list)
{
	lng_deleteAllShallow(list->headNode);
	free(list->headNode);
	free(list);
}

void lg_deleteShallowMultiple(lg_t** lists, int numLists)
{
	if(lists == NULL)
		return;

	int i;
	for(i = 0; i < numLists; ++i)
	{
		lg_deleteShallow(lists[i]);
	}

	free(lists);
}

void lg_decompose(lg_t* source, int partCount, lg_t*** lists_out)
{
	if(source == NULL)
		return;

	int totalItemCount = lg_size(source);
	int minItemCountPerPart = totalItemCount / partCount;
	int extraItemCount = totalItemCount % partCount;

	lg_t** lists = (lg_t**) malloc(sizeof(lg_t*) * partCount);
	int i;
	for(i = 0; i < partCount; ++i)
	{
		int currListItemCount = minItemCountPerPart;
		if(i < extraItemCount)
			++currListItemCount;

		lists[i] = lg_new();

		int j;
		for(j = 0; j < currListItemCount; ++j)
		{
			void* data = lg_removeFront(source);
			lg_addTailData(data, lists[i]);
		}
	}

	*lists_out = lists;
}

void lg_scatter(lg_t* source, int partCount, lg_t*** lists_out)
{
	if(source == NULL)
		return;

	int totalItemCount = lg_size(source);

	int i;
	lg_t** lists = (lg_t**) malloc(sizeof(lg_t*) * partCount);
	for(i = 0; i < partCount; ++i)
		lists[i] = lg_new();

	int currListIndex = 0;
	for(i = 0; i < totalItemCount; ++i)
	{
		void* data = lg_removeFront(source);
		lg_addTailData(data, lists[currListIndex]);

		++currListIndex;
		currListIndex = currListIndex % partCount;
	}

	*lists_out = lists;
}

void lg_split(lg_t* source, int partCount, lg_t*** lists_out)
{
	if(source == NULL)
		return;

	lg_t** lists = (lg_t**) malloc(sizeof(lg_t*) * partCount);
	int i;
	for(i = 0; i < partCount; ++i)
		lists[i] = lg_new();

	int itemCount = lg_size(source);
	__lg_split(source, itemCount, partCount, lists);


	*lists_out = lists;
}

// Helper Functions
// --------------------------------------------------------------------------------

static void __lg_split(lg_t* source, int totalItemCount, int totalSplitCount, lg_t** lists)
{
	if(totalSplitCount < 1 || totalItemCount < 1)
		return;

	if(totalSplitCount > 1)
	{
		int itemCountLeft = totalItemCount / 2 + totalItemCount % 2;
		int itemCountRight = totalItemCount - itemCountLeft;
		int splitCountLeft = totalSplitCount / 2 + totalSplitCount % 2;
		int splitCountRight = totalSplitCount - splitCountLeft;

		// go left first because we will use remove front function
		__lg_split(source, itemCountLeft, splitCountLeft, lists);
		__lg_split(source, itemCountRight, splitCountRight, &lists[splitCountLeft]);
	}
	else // if(totalSplitCount == 1)
	{
		int i;
		for(i = 0; i < totalItemCount; ++i)
		{
			void* data = lg_removeFront(source);
			lg_addTailData(data, lists[0]);
		}
	}
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
