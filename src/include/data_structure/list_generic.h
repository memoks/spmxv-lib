/*
 * list_generic.h
 *
 *  Created on: Jan 20, 2015
 *      Author: endoplasmic
 */

#ifndef LIST_GENERIC_H_
#define LIST_GENERIC_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/data_structure/list.h"

// List-Node-Generic Structure
// ---------------------------------------------------------------------------------

struct list_node_generic
{
	list_head_t listHead;
	void* dataPtr;
};

typedef struct list_node_generic lng_t;

extern lng_t* lng_new(void *dataPtr);
extern lng_t* lng_getListGeneric(list_head_t* head);
extern lng_t* lng_next(lng_t* node);
extern lng_t* lng_prev(lng_t* node);
extern void lng_addTailData(void* dataPtr, lng_t* node);
extern void lng_addTail(lng_t* new, lng_t* list);
extern void lng_addFrontData(void* dataPtr, lng_t* node);
extern void lng_addFront(lng_t* new, lng_t* node);
extern void* lng_remove(lng_t* node);
extern void lng_deleteShallow(lng_t* node);
extern void lng_deleteAllShallow(lng_t* headNode);
extern int lng_isEmpty(lng_t* list);
extern int lng_size(lng_t* list);


// List-Generic Structure
// ---------------------------------------------------------------------------------

struct list_generic
{
	lng_t* headNode;
	size_t size;
};

typedef struct list_generic lg_t;

extern lg_t* lg_new(void);
extern lng_t* lg_tail(lg_t* list);
extern lng_t* lg_front(lg_t* list);
extern void lg_addTailData(void* dataPtr, lg_t* list);
extern void lg_addTail(lng_t* newNode, lg_t* list);
extern void lg_addFrontData(void* dataPtr, lg_t* list);
extern void lg_addFront(lng_t* newNode, lg_t* list);
extern int lg_size(lg_t* list);
extern int lg_isEmpty(lg_t* list);
extern void* lg_remove(lg_t* list, lng_t* node);
extern void* lg_removeTail(lg_t* list);
extern void* lg_removeFront(lg_t* list);
extern void lg_removeAll(lg_t* list);
extern void lg_deleteShallow(lg_t* list);
extern void lg_deleteShallowMultiple(lg_t** lists, int numLists);

/**
 * Partitions a given list into multiple lists where each output list
 * will have chunks of original data as they were ordered in source list.
 *
 * @source (lg_t*):
 * @partCount (int):
 * @lists_out (lg_t***):
 *
 * Suppose we have;
 * source = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
 * partCount = 4
 *
 * Then lists_out will be;
 * *lists_out[0]: 0 1 2
 * *lists_out[1]: 3 4 9
 * *lists_out[2]: 5 6
 * *lists_out[3]: 7 8
 */
extern void lg_decompose(lg_t* source, int partCount, lg_t*** lists_out);

/**
 * Scatters all elements in a given list into multiple lists one by one.
 *
 *
 * @source (lg_t*):
 * @partCount (int):
 * @lists_out (lg_t***):
 *
 * Suppose we have;
 * source = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
 * partCount = 4
 *
 * Then lists_out will be;
 * *lists_out[0]: 0 4 8
 * *lists_out[1]: 1 5 9
 * *lists_out[2]: 2 6
 * *lists_out[3]: 3 7
 */
extern void lg_scatter(lg_t* source, int partCount, lg_t*** lists_out);

/**
 * Recursively splits a list into 2 lists, eventually list count reaches
 * "partCount" parameter.
 *
 * @source (lg_t*):
 * @partCount (int):
 * @lists_out (lg_t***):
 *
 * Suppose we have;
 * source = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
 * partCount = 4
 *
 * Then lists_out will be;
 * *lists_out[0]: 0 1 2
 * *lists_out[1]: 3 4
 * *lists_out[2]: 5 6 7
 * *lists_out[3]: 8 9
 */
extern void lg_split(lg_t* source, int partCount, lg_t*** lists_out);

// Prototypes & Definitions of data type dependent functions.
// ---------------------------------------------------------------------------------

/**
 * lg_??_toArray function.
 * Converts a given generic list into an array of "structs" (not pointers).
 *
 * @list: generic list data structure.
 * @arr_out: (return value) array of structures (not pointers) for a given list.
 * @length_out: (return value) length of the returned array.
 */
#define FUNC_PROTOTYPE_LG_TOARRAY(												\
		/* function name prefix */ PREFIX, 										\
		/* type to operate on */ type_t)										\
	void lg_##PREFIX##_toArray(												\
			lg_t* list, type_t** arr_out, DECIMAL* length_out)

#define FUNC_DEFINITION_LG_TOARRAY(											\
		/* function name prefix */ PREFIX, 										\
		/* type to operate on */ type_t)										\
	FUNC_PROTOTYPE_LG_TOARRAY(PREFIX, type_t)									\
	{																			\
		size_t length = lg_size(list);											\
		if(length <= 0)															\
		{																		\
			*arr_out = NULL;													\
			return;																\
		}																		\
																				\
		type_t* arr = (type_t*) malloc(sizeof(type_t) * length);				\
		int index = 0;															\
																				\
		list_head_t* currListHead = NULL;										\
		list_for_each(currListHead, &list->headNode->listHead)					\
		{																		\
			/* Copy from list node to an array index. */ 						\
			lng_t* currListNode = lng_getListGeneric(currListHead);				\
																				\
			/* ASSUMPTION: copy function with the following name must be */		\
			/* defined. */														\
			PREFIX##_copy((type_t*) currListNode->dataPtr, &arr[index]);		\
			++index;															\
		}																		\
																				\
		/* return values */														\
		*arr_out = arr;															\
		*length_out = length;													\
	}

/**
 * lg_??_toArrayMultiple function.
 * Converts multiple instance of generic list into an arrays of
 * "structs" (not pointers).
 *
 * @lists: array of generic list data structures.
 * @numLists: length of given list array.
 * @arrs_out: (return value) arrays of structures (not pointers) for a given lists.
 * @lengths_out: (return value) lengths of each returned array of structures.
 */
#define FUNC_PROTOTYPE_LG_TOARRAY_MULTIPLE(										\
		/* function name prefix */ PREFIX, 										\
		/* type to operate on */ type_t)										\
		void lg_##PREFIX##_toArrayMultiple										\
			(lg_t** lists, int numLists, type_t*** arrs_out, int** lengths_out)

#define FUNC_DEFINITION_LG_TOARRAY_MULTIPLE(									\
		/* function name prefix */ PREFIX, 										\
		/* type to operate on */ type_t)										\
	FUNC_PROTOTYPE_LG_TOARRAY_MULTIPLE(PREFIX, type_t)							\
	{																			\
		type_t** arrs = (type_t**) malloc(sizeof(type_t*) * numLists);			\
		DECIMAL* lengths = (DECIMAL*) malloc(sizeof(DECIMAL) * numLists);		\
																				\
		int i;																	\
		for(i = 0; i < numLists; ++i)											\
		{																		\
			lengths[i] = 0;														\
			arrs[i] = NULL;														\
			lg_##PREFIX##_toArray(lists[i], &arrs[i], &lengths[i]);				\
		}																		\
																				\
		/* return values */														\
		*arrs_out = arrs;														\
		*lengths_out = lengths;													\
	}

/**
 * lg_??_copyAdd function.
 * Copies given "source" list elements and inserts them into "destination"
 * list using given add function. There are 2 add functions.
 * <ul>
 * <li>lg_addTail</li>
 * <li>lg_addFront</li>
 * </ul>
 *
 * @source: generic list whose elements will be copied.
 * @destination: generic list in which copied elements will be put.
 */
#define FUNC_PROTOTYPE_LG_COPYADD(												\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t,										\
		/* "Tail" or "Front" */ LG_ADD_TYPE)									\
	void lg_##PREFIX##_copyAdd##LG_ADD_TYPE##Data(lg_t* source, lg_t* destination)

#define FUNC_COMMON_BODY_LG_COPYADD(											\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t,										\
		/* "Tail" or "Front" */ LG_ADD_TYPE)									\
	if(source == NULL || destination == NULL)									\
		return;																	\
																				\
	list_head_t* currListHead = NULL;											\
	list_for_each(currListHead, &source->headNode->listHead)					\
	{																			\
		lng_t* currLng = lng_getListGeneric(currListHead);						\
																				\
		/* ASSUMPTION: this is a deep copy since function copies data also, */	\
		/* but it doesn't copy the actual links or whatsoever connection */		\
		/* list-node-generic has. This is important... only the data. Copy (*/	\
		/* function with the following name must be defined somewhere else. */	\
		type_t* copyData = PREFIX##_copyToPtr(currLng->dataPtr);				\
																				\
		lg_add##LG_ADD_TYPE##Data(copyData, destination);						\
	}

#define FUNC_DEFINITION_LG_COPYADD(												\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t,										\
		/* "Tail" or "Front" */ LG_ADD_TYPE)									\
	FUNC_PROTOTYPE_LG_COPYADD(PREFIX, type_t, LG_ADD_TYPE)						\
	{																			\
		FUNC_COMMON_BODY_LG_COPYADD(PREFIX, type_t, LG_ADD_TYPE)				\
	}

/**
 * lg_??_print function.
 * Prints all elements in the parameter list.
 *
 * @message: message to be printed before list contents.
 * @list: generic list whose contents will be printed.
 */
#define FUNC_PROTOTYPE_LG_PRINT(												\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t)										\
	void lg_##PREFIX##_print(char* message, lg_t* list)

#define FUNC_DEFINITION_LG_PRINT(												\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t)										\
	FUNC_PROTOTYPE_LG_PRINT(PREFIX, type_t)										\
	{																			\
		int itemCount = lg_size(list);											\
		PRINTF("%s: listing %d items...\n", message, itemCount);				\
																				\
		list_head_t* currListHead = NULL;										\
		list_for_each(currListHead, &list->headNode->listHead)					\
		{																		\
			lng_t* currListNode = lng_getListGeneric(currListHead);				\
																				\
			/* ASSUMPTION: print function with the following name */			\
			/* must be defined. */												\
			PREFIX##_print((type_t*) currListNode->dataPtr);					\
		}																		\
	}

/**
 * lg_??_print function.
 * Prints all elements in given set of lists.
 *
 * @lists: generic list array whose contents will be printed.
 * @numLists: number of generic lists.
 */
#define FUNC_PROTOTYPE_LG_PRINT_MULTIPLE(										\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t)										\
	void lg_##PREFIX##_printMultiple(lg_t** lists, int numLists)

#define FUNC_DEFINITION_LG_PRINT_MULTIPLE(										\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t)										\
	FUNC_PROTOTYPE_LG_PRINT_MULTIPLE(PREFIX, type_t)							\
	{																			\
		if(lists == NULL)														\
			return;																\
																				\
		int totalNoOfElements = 0;												\
		int i;																	\
		for(i = 0; i < numLists; ++i)											\
		{																		\
			if(lists[i] != NULL)												\
			{																	\
				char message[DEFAULT_STR_BUFF_SIZE];							\
				sprintf(message, "list %d", i);									\
				lg_##PREFIX##_print(message, lists[i]);							\
																				\
				totalNoOfElements += lg_size(lists[i]);							\
			}																	\
		}																		\
		PRINTF("\n");															\
		PRINTF("Total Number of Elements: %d\n", totalNoOfElements);			\
	}

/**
 * lg_??_deleteDeep function.
 * Deletes not only the list but also all of elements in the it.
 *
 * @list: generic list structure which is going to be "completely" deleted.
 */
#define FUNC_PROTOTYPE_LG_DELETEDEEP(											\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t)										\
	void lg_##PREFIX##_deleteDeep(lg_t* list)

#define FUNC_DEFINITION_LG_DELETEDEEP(											\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t)										\
	FUNC_PROTOTYPE_LG_DELETEDEEP(PREFIX, type_t)								\
	{																			\
		if(list == NULL)														\
			return;																\
																				\
		list_head_t* curr;														\
		list_head_t* temp;														\
																				\
		list_for_each_safe(curr, temp, &list->headNode->listHead)				\
		{																		\
			lng_t* currListNode = lng_getListGeneric(curr);						\
																				\
			PREFIX##_delete((type_t*) currListNode->dataPtr);					\
			lng_deleteShallow(currListNode);									\
		}																		\
																				\
		free(list->headNode);													\
		free(list);																\
	}

/**
 * lg_??_deleteDeepMultiple function.
 * Deletes not only multiple lists (given as parameter) but also all
 * of elements in them.
 *
 * @lists: an array of generic lists whose contents will be "completely" deleted.
 * @numLists: length of generic list array.
 */
#define FUNC_PROTOTYPE_LG_DELETEDEEP_MULTIPLE(									\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t)										\
	void lg_##PREFIX##_deleteDeepMultiple(lg_t** lists, int numLists)

#define FUNC_DEFINITION_LG_DELETEDEEP_MULTIPLE(									\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t)										\
	FUNC_PROTOTYPE_LG_DELETEDEEP_MULTIPLE(PREFIX, type_t)						\
	{																			\
		if(lists == NULL)														\
			return;																\
																				\
		for(int i = 0; i < numLists; ++i)										\
			lg_##PREFIX##_deleteDeep(lists[i]);									\
																				\
		free(lists);															\
	}

/**
 * lg_??_partitionDeep function.
 * Divides a given list into specified number of partition using one of
 * the following algorithms.
 * <ul>
 *   <li>decompose</li>
 *   <li>scatter</li>
 *   <li>split</li>
 * </ul>
 * However, after function call structure of the list remains unchanged.
 *
 * @source: list structure to divide.
 * @partCount: data type that is stored in list structure.
 * @lists_out: (return value) "partCount" number of generic lists.
 */
#define FUNC_PROTOTYPE_LG_PARTITIONDEEP(										\
		/* function name prefix */ PREFIX,										\
		/* type to operate on */ type_t,										\
		/* partition function */ SUFFIX)										\
	void lg_##PREFIX##_##SUFFIX##Deep(											\
			lg_t* source, int partCount, lg_t*** lists_out)

#define FUNC_DEFINITION_LG_PARTITIONDEEP(										\
		/* function name prefix */ PREFIX,										\
 		/* type to operate on */ type_t,										\
 		/* partition function */ SUFFIX)										\
	FUNC_PROTOTYPE_LG_PARTITIONDEEP(PREFIX, type_t, SUFFIX)						\
	{																			\
 		if(source == NULL)														\
 			return;																\
 																				\
 		/* Create a deep copy of source list */									\
		lg_t* copy = lg_new();													\
		lg_##PREFIX##_copyAdd##Tail##Data(source, copy)	;						\
																				\
		/* Partition copy list */												\
		lg_##SUFFIX(copy, partCount, lists_out);								\
																				\
		/* Delete (by now) empty list */										\
		lg_deleteShallow(copy);													\
	}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* LIST_GENERIC_H_ */

