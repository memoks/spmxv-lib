
#ifndef TREE_SCHED_H
#define TREE_SCHED_H

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/data_structure/tree.h"

struct tree_sched
{
	// current index of victim in victims array
	int victimIndex;

	// victims[victimIndex]
	int lastVictim;

	// How many jobs are stolen until now?
	int stealCount;

	// length of the victims array
	int victimCount;

	// victims array.
	// Stores the order in which victims are visited.
	int* victims;
};

typedef struct tree_sched tree_sched_t;

extern void tree_sched_init(int blockId,
		int numBlocks, void* treeSched);
extern void tree_sched_terminate(void* ringSched);
extern tree_sched_t* tree_sched_new(void);
extern void tree_sched_delete(tree_sched_t* ringSched);
extern void tree_sched_reset(int blockId, void* ringSched);
extern void tree_sched_print(tree_sched_t* ringSched);

/**
 * tree_sched - get the struct for this entry
 * @ptr:	the &struct tree pointer.
 * @type:	the type of the struct this is embedded in.
 * @member:	the name of the tree_struct within the struct.
 */
#define tree_sched_entry(ptr, type, member) \
	((type *)((char *)(ptr)-(unsigned long)(&((type *)0)->member)))

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* TREE_SCHED_H */
