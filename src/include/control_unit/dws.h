
#ifndef SRC_INCLUDE_CONTROL_UNIT_DWS_H_
#define SRC_INCLUDE_CONTROL_UNIT_DWS_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/data_structure/quintet.h"
#include "include/data_structure/sub_mtx.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/list_generic.h"
#include "include/data_structure/vector.h"
#include "include/scheduler/block_info.h"
#include "include/control_unit/static.h"
#include "include/control_unit/fast_run.h"

struct dws_args
{
	static_args_t staticArgs;

	// These will be passed as parameter
	// -----------------------------------------------------------------------------------------
	int algorithmType;
	int stealingScheme;
	int stealTreshold;

	spm_cmp_t* spmCsr;
	DECIMAL* permutation;
	vector_int_t* rowPartitioning;
	vector_int_t* columnPartitioning;
	int subMtxCount;

	// These will be calculated from parameters above
	// -----------------------------------------------------------------------------------------
	block_info_t** warps;
	int numWarps;
	int threadsPerWarp;
	lg_t** execHistoryPerWarp;
	lg_t** execHistoryPerThread;

	// function pointers
	// -----------------------------------------------------------------------------------------

	// for spmxv multiply, warm-up, and data-preparation routines
	void (*prep) (struct dws_args* args);
	void (*spmxv_dws_warmup) (struct dws_args* args, vector_real_t* x, vector_real_t* y_out);
	void (*spmxv_dws) (static_args_t* args, vector_real_t* x, vector_real_t* y_out);

	// for scheduling algorithm
	void (*scheduler_init) (block_info_t** localExecutionGroup, int numBlocks, int subMtxCount);
	void (*scheduler_terminate) (block_info_t** localExecutionGroup, int numBlocks);
	int (*scheduler_findVictim) (block_info_t* myBlock, block_info_t** localExecutionGroup, int numBlocks);
};

typedef struct dws_args dws_args_t;

extern dws_args_t* dws_args_new();
extern void dws_args_initDefault(dws_args_t* args);
extern void dws_args_init(
		dws_args_t* args, fast_run_t* fr,
		int algorithmType, int stealingScheme, int stealTreshold);
extern void dws_args_print(dws_args_t* args);
extern void dws_args_warmup(dws_args_t* args, vector_real_t* x, vector_real_t* y_inout);
extern void dws_args_setup(dws_args_t* args);
extern void dws_args_reset(dws_args_t* args);
extern void dws_args_terminate(dws_args_t* args);
extern void dws_args_delete(dws_args_t* args);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SRC_INCLUDE_CONTROL_UNIT_DWS_H_ */
