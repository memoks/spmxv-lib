#ifndef SRC_INCLUDE_CONTROL_UNIT_GWS_H_
#define SRC_INCLUDE_CONTROL_UNIT_GWS_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/scheduler/job_queue.h"
#include "include/control_unit/static.h"
#include "include/control_unit/fast_run.h"

struct gws_args
{
	static_args_t staticArgs;

	job_queue_t* globalQueue;
	lg_t** execHistoryPerThread;

	// function pointers
	// -----------------------------------------------------------------------------------------

	// for spmxv multiply, warm-up, and data-preparation routines
	void (*prep) (struct gws_args* args);
	void (*spmxv_gws_warmup) (struct gws_args* args, vector_real_t* x, vector_real_t* y_out);
	void (*spmxv_gws) (static_args_t* args, vector_real_t* x, vector_real_t* y_out);
};

typedef struct gws_args gws_args_t;

extern gws_args_t* gws_args_new();
extern void gws_args_initDefault(gws_args_t* args);
extern void gws_args_init(gws_args_t* args, fast_run_t* fr);
extern void gws_args_print(gws_args_t* args);
extern void gws_args_warmup(gws_args_t* args, vector_real_t* x, vector_real_t* y_inout);
extern void gws_args_setup(gws_args_t* args);
extern void gws_args_terminate(gws_args_t* args);
extern void gws_args_delete(gws_args_t* args);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SRC_INCLUDE_CONTROL_UNIT_GWS_H_ */
