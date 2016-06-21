
#ifndef SPMXV_SEQUENTIAL_H_
#define SPMXV_SEQUENTIAL_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <omp.h>

#include "include/config.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/spm_storage.h"


#define SPMXV_SEQUENTIAL_CSR(		\
		/* spm_cmp_t* */ spm,		\
		/* vector_real_t* */ x,		\
		/* vector_real_t* */ y)		\
	{								\
		SPMXV_CSR_PARTIAL(			\
				0, 					\
				spm->rowCount, 		\
				spm->ptr, 			\
				spm->ind, 			\
				spm->values, 		\
				x->data, 			\
				y->data)			\
	}

#define SPMXV_SEQUENTIAL_JDS(		\
		/* spm_jds_t* */ spm,		\
		/* vector_real_t* */ x,		\
		/* vector_real_t* */ y)		\
	{								\
		SPMXV_JDS_PARTIAL(			\
				spm->rowCount, 		\
				spm->colCount, 		\
				spm->nnz, 			\
				0,					\
				spm->idiagLength, 	\
				spm->idiag, 		\
				spm->jdiag, 		\
				spm->dj, 			\
				spm->permutation, 	\
				x->data, 			\
				y_inout->data)		\
	}

extern void spmxv_sequential_CSR(
		spm_cmp_t* spmCsr, vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_sequential_CSC(
		spm_cmp_t* spmCsc, vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_sequential_ICSR(
		spm_inc_t* spmIcsr, vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_sequential_JDS(
		spm_jds_t* spmJds, vector_real_t* x, vector_real_t* y_inout);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SPMXV_SEQUENTIAL_H_ */
