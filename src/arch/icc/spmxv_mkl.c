
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

#include "mkl.h"
#include "mkl_types.h"
#include "mkl_cblas.h"

#include "include/config.h"
#include "include/timer/custom_timer.h"
#include "arch/icc/inc/spmxv_mkl.h"

void spmxv_mkl(spm_cmp_t* spm, vector_real_t* x, vector_real_t* y_out)
{

#ifdef SINGLE_PRECISION
	mkl_cspblas_scsrgemv("N", &spm->rowCount, spm->values, spm->ptr, spm->ind, x->data, y_out->data);
#else
	mkl_cspblas_dcsrgemv("N", &spm->rowCount, spm->values, spm->ptr, spm->ind, x->data, y_out->data);
#endif

}

void spmxv_mklPartial(
		DECIMAL startRowInd,
		DECIMAL endRowInd,
		DECIMAL* RESTRICT rowPtr,
		DECIMAL* RESTRICT colInd,
		REAL* RESTRICT values,
		REAL* RESTRICT x,
		REAL* RESTRICT y_out)
{
	DECIMAL rowCount = endRowInd - startRowInd;

#ifdef SINGLE_PRECISION
	mkl_cspblas_scsrgemv("N", &rowCount, &values[rowPtr[startRowInd]], &rowPtr[startRowInd],
				&colInd[rowPtr[startRowInd]], x, &y_out[startRowInd]);
#else
	mkl_cspblas_dcsrgemv("N", &rowCount, &values[rowPtr[startRowInd]], &rowPtr[startRowInd],
				&colInd[rowPtr[startRowInd]], x, &y_out[startRowInd]);
#endif

}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
