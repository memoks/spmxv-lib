
#include <stdlib.h>
#include <stdio.h>

#include "include/data_structure/stack_generic.h"

// Generic Functions
// -----------------------------------------------------------------------------------------

stack_generic_t* stack_generic_new(void)
{
	stack_generic_t* stack = (stack_generic_t*) malloc(sizeof(stack_generic_t));
	stack->head = lng_new(NULL);
	stack->curr = stack->head;
	stack->size = 0;

	return stack;
}

void* stack_generic_pop(stack_generic_t* stack)
{
	if(stack->curr == stack->head)
		return NULL;

	void* dataPtr = stack->curr->dataPtr;
	lng_t* newCurr = lng_prev(stack->curr);
	lng_deleteShallow(stack->curr);
	stack->curr = newCurr;

	return dataPtr;
}

void stack_generic_push(void* dataPtr, stack_generic_t* stack)
{
	lng_t* new = lng_new(dataPtr);
	lng_addFront(new, stack->curr);
	stack->curr = new;
}

int stack_generic_isEmpty(stack_generic_t* stack)
{
	if(stack == NULL)
		return TRUE;

	if(stack->curr == stack->head)
		return TRUE;

	return FALSE;
}

void stack_generic_deleteAllShallow(stack_generic_t* stack)
{
	lng_deleteAllShallow(stack->head);
}
