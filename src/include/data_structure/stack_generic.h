
#ifndef STACK_GENERIC_H_
#define STACK_GENERIC_H_

#include <stdlib.h>

#include "include/config.h"
#include "include/data_structure/list_generic.h"

// TODO test all these

struct stack_generic
{
	lng_t* head;
	lng_t* curr;
	size_t size;
};

typedef struct stack_generic stack_generic_t;

// Generic functions
// ---------------------------------------------------------------------

extern stack_generic_t* stack_generic_new(void);
extern void* stack_generic_pop(stack_generic_t* stack);
extern void stack_generic_push(void* dataPtr, stack_generic_t* stack);
extern int stack_generic_isEmpty(stack_generic_t* stack);
extern void stack_generic_deleteAllShallow(stack_generic_t* stack);


#endif /* STACK_GENERIC_H_ */
