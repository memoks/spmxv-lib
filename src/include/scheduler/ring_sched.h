
#ifndef RING_SCHED_H_
#define RING_SCHED_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

// TODO comment
struct ring_sched
{
	int victimId;
	int stealRange;
};

typedef struct ring_sched ring_sched_t;

extern void ring_sched_init(
		int blockId, int numBlocks, void* ringSched);
extern void ring_sched_terminate(void* ringSched);
extern ring_sched_t* ring_sched_new(int blockId, int numBlocks);
extern void ring_sched_delete(ring_sched_t* ringSched);
extern void ring_sched_reset(int blockId, void* ringSched);
extern void ring_sched_print(ring_sched_t* ringSched);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* RING_SCHED_H_ */
