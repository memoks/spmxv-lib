
#ifndef SCHEDULER_H_
#define SCHEDULER_H_

struct block_info;
typedef struct block_info block_info_t;

extern void tree_scheduler_init(
		block_info_t** blockInfos, int numBlocks,
		int jobBatchCount);
extern void tree_scheduler_terminate(
		block_info_t** blockInfos, int numBlocks);
extern int tree_scheduler_findVictim(
		block_info_t* myBlock,
		block_info_t** blockInfos, int numBlocks);

#endif /* SCHEDULER_H_ */
