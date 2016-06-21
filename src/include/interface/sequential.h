
#ifndef SEQUENTIAL_H_
#define SEQUENTIAL_H_

#include "include/config.h"
#include "include/io/cli.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/vector.h"

#include "include/interface/common.h"

struct sequential
{
	common_routine* cr;

	spm_cmp_t* spm;

	void (*initialize)(struct sequential* self);
	void (*warmUp)(struct sequential* self, vector_real_t* x);
	void (*terminate)(struct sequential* self);

	void (*copy)(struct sequential* self,
			vector_real_t* dest, vector_real_t* source);
	void (*add)(struct sequential* self,
			vector_real_t* r, vector_real_t* lhs, vector_real_t* rhs);
	void (*subtract)(struct sequential* self,
			vector_real_t* r, vector_real_t* lhs, vector_real_t* rhs);
	void (*accumulate)(struct sequential* self, vector_real_t* res,
			vector_real_t* lhs, REAL scalar, vector_real_t* rhs);
	void (*matrixMult)(struct sequential* self,
			vector_real_t* x, vector_real_t* y_out);
	void (*scalarMult)(struct sequential* self,
			vector_real_t* x, REAL scalar, vector_real_t* y_out);
	REAL (*innerProd)(struct sequential* self,
			vector_real_t* r1, vector_real_t* r2);

	void (*asolve)(struct sequential* self,
			vector_real_t* res, vector_real_t* b);
};

typedef struct sequential sequential_routine;

extern sequential_routine* sequential_routine_new(cli_options_t* options);
extern void sequential_routine_cleanUp(sequential_routine* sr);

#endif /* SEQUENTIAL_H_ */
