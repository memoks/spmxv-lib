
#ifndef BLOCK_INFO_H_
#define BLOCK_INFO_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"

#include "include/scheduler/job_queue.h"
#include "include/scheduler/job_batch.h"
#include "include/scheduler/ring_sched.h"
#include "include/scheduler/tree_sched.h"
#include "include/scheduler/ring_scheduler.h"
#include "include/scheduler/tree_scheduler.h"

// TODO comment all

// Scheduler data structure
// ------------------------------------------
// ASSUMPTION: there should not be any other structure data in
// sched (such as int type...) , other than scheduler type.
// Because of block_info's sched_terminate function pointer.
union sched
{
	ring_sched_t ringSched;
	tree_sched_t treeSched;
};

typedef union sched sched_t;

// ------------------------------------------

struct block_info
{
	int id;
	job_queue_t queue;

	int schedType;
	sched_t sched;

	int (*scheduler_findVictim) \
			(struct block_info* myBlock, struct block_info** blockInfos, int numBlocks);

	void (*sched_init) (int blockId, int numBlocks, void* sched);
	void (*sched_reset) (int blockId, void* schedAddr);
	void (*sched_terminate) (void* sched);
};

typedef struct block_info block_info_t;

extern void block_init(
		block_info_t* blockInfo, int blockId, int numBlocks,
		int stealTreshold, int stealingScheme);
extern void block_deleteNonPtr(block_info_t* blockInfo);
extern block_info_t* block_new(int blockId, int stealTreshold, int stealingScheme, int numBlocks);
extern block_info_t** block_newMultiple(int numBlocks, int stealTreshold, int stealingScheme);
extern void block_delete(block_info_t* blockInfo);
extern void block_deleteMultiple(block_info_t** blockInfos, int numBlocks);
extern block_info_t* block_getBlockFromRingSched(ring_sched_t* ringSched);
extern block_info_t* block_getBlockFromTreeSched(tree_sched_t* treeSched);
extern void block_reset(block_info_t* blockInfo);
extern void block_resetMultiple(block_info_t** blockInfos, int numBlocks);
extern void block_reform(
		block_info_t* blockInfo, lg_t* executedJobBatchList);
extern void block_reformMultiple(
		block_info_t** blockInfos,
		lg_t** executedJobBatchListPerBlock, int numBlocks);
extern void block_print(block_info_t* blockInfo);
extern void block_printDetailed(block_info_t* blockInfo, void* spm);
extern void block_printMultiple(
		block_info_t** blockInfos, int numBlocks);
extern void block_printDetailedMultiple(
		block_info_t** blockInfos, int numBlocks, void* spm);

extern void block_convertToJDSPartialMultiple(
		spm_cmp_t* spmCsrCounterpart, DECIMAL* permutation,
		block_info_t** blockInfos_inout, int numBlocks);

extern void block_convertFromPartialJDSToHybridJDSCSRMultiple(
		block_info_t** blockInfos_inout, int numBlocks);

// TODO test
extern 	int block_localityAwareStealHalf(
		block_info_t* currBlock, block_info_t** allBlocks, int blockCount);

/**
 * block_info_t - get the struct for this entry
 * @ptr:	the &treeSched / ringSched pointer.
 * @type:	the type of the struct this is embedded in.
 * @member:	the name of the struct within the struct.
 */
#define block_entry(ptr, type, member) \
	((type *) ((char *)(ptr) - (unsigned long)(&((type *)0)->member)))

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* BLOCK_INFO_H_ */
