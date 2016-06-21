
#include "inc/spmxv_sequential.h"
#include "inc/kernel.h"

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

inline void spmxv_sequential_CSR(spm_cmp_t* spmCsr, vector_real_t* x, vector_real_t* y_inout)
{
	SPMXV_CSR_PARTIAL(
			0, spmCsr->rowCount, spmCsr->ptr, spmCsr->ind, spmCsr->values, x->data, y_inout->data)
}

inline void spmxv_sequential_CSC(spm_cmp_t* spmCsc, vector_real_t* x, vector_real_t* y_inout)
{
	SPMXV_CSC_PARTIAL(
			0, spmCsc->rowCount, spmCsc->ptr, spmCsc->ind, spmCsc->values, x->data, y_inout->data);
}

inline void spmxv_sequential_ICSR(spm_inc_t* spmIcsr, vector_real_t* x, vector_real_t* y_inout)
{
	SPMXV_ICSR(
			spmIcsr->rowCount, spmIcsr->colCount, spmIcsr->nnz,
			spmIcsr->ptr, spmIcsr->ind, spmIcsr->values, x->data, y_inout->data)
}

inline void spmxv_sequential_JDS(spm_jds_t* spmJds, vector_real_t* x, vector_real_t* y_inout)
{
	SPMXV_JDS_PARTIAL(
			0, spmJds->idiagLength, spmJds->idiag, spmJds->jdiag,
			spmJds->dj, spmJds->permutation, x->data, y_inout->data)
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
