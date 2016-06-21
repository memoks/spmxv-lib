
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/control_unit/omp_task.h"

omp_task_args_t* omp_task_args_new()
{
	omp_task_args_t* args = (omp_task_args_t*) malloc(sizeof(omp_task_args_t));
	omp_task_args_initDefault(args);

	return args;
}

void omp_task_args_initDefault(omp_task_args_t* args)
{
	static_args_initDefault(&args->staticArgs);

	args->execHistPerThread = NULL;

	// function pointers
	// -----------------------------------------------------------------------------------------

	args->prep = NULL;
	args->spmxv_omp_task_warmup = NULL;
	args->spmxv_omp_task = NULL;
}

void omp_task_args_init(	omp_task_args_t* args, fast_run_t* fr)
{
	static_args_init(&args->staticArgs, fr);
	args->spmxv_omp_task = args->staticArgs.spmxv_static;
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
