
#ifndef SRC_INCLUDE_CONTROL_UNIT_OMP_TASK_H_
#define SRC_INCLUDE_CONTROL_UNIT_OMP_TASK_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/control_unit/static.h"
#include "include/control_unit/fast_run.h"

struct omp_task_args
{
	static_args_t staticArgs;

	lg_t** execHistPerThread;

	// function pointers
	// -----------------------------------------------------------------------------------------

	// for spmxv multiply and data-preparation routines
	// (omp_loop routine doesn't have a different warm-up routine)
	void (*prep) (struct omp_task_args* args);
	void (*spmxv_omp_task_warmup) (struct omp_task_args* args, vector_real_t* x, vector_real_t* y_out);
	void (*spmxv_omp_task) (static_args_t* args, vector_real_t* x, vector_real_t* y_out);
};

typedef struct omp_task_args omp_task_args_t;


extern omp_task_args_t* omp_task_args_new();
extern void omp_task_args_initDefault(omp_task_args_t* args);
extern void omp_task_args_init(omp_task_args_t* args, fast_run_t* fr);
// TODO omp_task_args_terminate()

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SRC_INCLUDE_CONTROL_UNIT_OMP_TASK_H_ */
