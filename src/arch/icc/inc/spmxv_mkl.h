
#ifndef SPMXV_MKL_H_
#define SPMXV_MKL_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/data_structure/spm_storage.h"
#include "include/data_structure/vector.h"

extern void spmxv_mkl(
		spm_cmp_t* spmCsr,
		vector_real_t* x, vector_real_t* y_out);

extern void spmxv_mklPartial(
		DECIMAL startRowInd,
		DECIMAL endRowInd,
		DECIMAL* RESTRICT rowPtr,
		DECIMAL* RESTRICT colInd,
		REAL* RESTRICT values,
		REAL* RESTRICT x,
		REAL* RESTRICT y_out);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SPMXV_MKL_H_ */
