
#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include "include/io/input_parser.h"

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"

extern void algorithm_parallelMergeQuickSort(
		quintet_t* quintets_inout,
		DECIMAL length, int threadCount,
		int (*quintet_cmpFunc) (const void*, const void*));

extern void algorithm_parallelMerge(
		quintet_t* quintets_inout,
		DECIMAL* borders, DECIMAL borderCount,
		int (*quintet_cmpFunc) (const void*, const void*));

extern void algorithm_merge(
		quintet_t* quintets1, DECIMAL length1,
		quintet_t* quintets2, DECIMAL length2,
		quintet_t* quintetsMerged_inout,
		int (*quintet_cmpFunc) (const void*, const void*));

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* ALGORITHM_H_ */
