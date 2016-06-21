
#ifndef COMMON_H_
#define COMMON_H_

#include "include/config.h"
#include "include/io/cli.h"
#include "include/io/input_parser.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/spm_storage.h"


struct common
{
	cli_options_t* options;

	DECIMAL nnz;
	DECIMAL rowCount;
	DECIMAL colCount;
	quintet_t* triplets;
	vector_int_t* rowOrderLookup;
	vector_int_t* colOrderLookup;
	vector_int_t* subMtxData;
	spm_cmp_t* spm;

	void (*initialize)(struct common* self);
	void (*reorderVector)(struct common* self, vector_real_t* v);
	void (*terminate)(struct common* self);
};

typedef struct common common_routine;

extern common_routine* common_routine_new(cli_options_t* options);
extern void common_routine_cleanUp(common_routine* cr);


#endif /* COMMON_H_ */
