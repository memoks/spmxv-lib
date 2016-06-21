
#ifndef DWS_H_
#define DWS_H_

#include "include/config.h"
#include "include/io/cli.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/list_generic.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/sub_mtx.h"
#include "include/scheduler/job_queue.h"
#include "include/scheduler/job_batch.h"
#include "include/scheduler/block_info.h"

#include "include/interface/common.h"

struct dws
{
	common_routine* cr;

	spm_cmp_t* spm;
	block_info_t** blockInfos;
	sub_mtx_tree_t* subMtxHead;
	lg_t** execHistoryPerBlock;
	sub_mtx_dim_t** subMtxArrayPerBlock;
	DECIMAL* subMtxCountPerBlock;

	void (*initialize)(struct dws* self);
	void (*warmUp)(struct dws* self, vector_real_t* x);
	void (*mult)(struct dws* self,
			vector_real_t* x, vector_real_t* y_out);
	void (*scalarMult)(struct dws* self,
			vector_real_t* x, REAL scalar, vector_real_t* y_out);
	REAL (*innerProd)(struct dws* self,
			vector_real_t* r1, vector_real_t* r2);
	void (*terminate)(struct dws* self);
};

typedef struct dws dws_routine;

extern dws_routine* dws_routine_new(cli_options_t* options);
extern void dws_routine_cleanUp(dws_routine* dr);


#endif /* DWS_H_ */
