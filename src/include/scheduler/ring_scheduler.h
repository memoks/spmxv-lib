
#ifndef RING_SCHEDULER_H_
#define RING_SCHEDULER_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

struct block_info;
typedef struct block_info block_info_t;

extern void ring_scheduler_init(
		block_info_t** blockInfos, int numBlocks,
		int jobBatchCount);
extern void ring_scheduler_terminate(
		block_info_t** blockInfos, int numBlocks);
extern int ring_scheduler_findVictim(
		block_info_t* myBlock,
		block_info_t** blockInfos, int numBlocks);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* RING_SCHEDULER_H_ */
