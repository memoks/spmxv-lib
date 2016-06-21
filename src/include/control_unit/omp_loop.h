
#ifndef SRC_INCLUDE_CONTROL_UNIT_OMP_LOOP_H_
#define SRC_INCLUDE_CONTROL_UNIT_OMP_LOOP_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/control_unit/static.h"
#include "include/control_unit/fast_run.h"

struct omp_loop_args
{
	// ASSUMPTION this must be the first element of this struct
	static_args_t staticArgs;

	int chunkSize;

	// function pointers
	// -----------------------------------------------------------------------------------------

	// for spmxv multiply and data-preparation routines
	// (omp_loop routine doesn't have a different warm-up routine)
	void (*prep) (struct omp_loop_args* args);
	void (*spmxv_omp_loop_dynamic) (struct omp_loop_args* args, vector_real_t* x, vector_real_t* y_out);
	void (*spmxv_omp_loop_guided) (struct omp_loop_args* args, vector_real_t* x, vector_real_t* y_out);
};

typedef struct omp_loop_args omp_loop_args_t;

extern omp_loop_args_t* omp_loop_args_new();
extern void omp_loop_args_initDefault(omp_loop_args_t* args);
extern void omp_loop_args_init(
		omp_loop_args_t* args, fast_run_t* fr, int chunkSize);
extern void omp_loop_args_print(omp_loop_args_t* args);
extern void omp_loop_args_terminate(omp_loop_args_t* args);
extern void omp_loop_args_delete(omp_loop_args_t* args);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SRC_INCLUDE_CONTROL_UNIT_OMP_LOOP_H_ */
