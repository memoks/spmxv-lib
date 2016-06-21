
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "inc/kernel.h"

/*
void SPMXV_ICSR_ATOMIC(
		const INTEGER rowCount, const INTEGER colCount, const INTEGER nnz,
		const INTEGER* RESTRICT ptr, const INTEGER* RESTRICT inc, const REAL* RESTRICT values,
		const REAL* RESTRICT x, REAL* RESTRICT y)
{
	if(nnz <= 0)
		return;

	INTEGER r = 0;
	INTEGER t = 0;
	INTEGER i = 0;
	INTEGER j = inc[t];
	while(t < nnz)
	{
		REGISTER REAL temp = 0.0;

		i = i + ptr[r];
		while(j < colCount)
		{
			temp = temp + values[t] * x[j];

			t = t + 1;
			j = j + inc[t];
		}

		#pragma omp atomic
		y[i] += temp;

		j = j - colCount;
		r = r + 1;
	}
}


void SPMXV_ICSR(INTEGER rowCount, INTEGER colCount, INTEGER nnz,
		INTEGER* RESTRICT ptr, INTEGER* RESTRICT inc, REAL* RESTRICT values,
		REAL* RESTRICT x, REAL* RESTRICT y)
{
	INTEGER r = 0;	// rowInc index
	INTEGER t = 0;	// non-zero index
	INTEGER i = 0;
	INTEGER j = inc[t];
	while(t < nnz)
	{
		REGISTER REAL temp = 0.0;

		i = i + ptr[r];
		while(j < colCount)
		{
			// y[i] = y[i] + values[t] * x[j];
			temp = temp + values[t] * x[j];

			t = t + 1;
			j = j + inc[t];
		}

		y[i] = temp;

		j = j - colCount;
		r = r + 1;
	}
}

*/

void spmxv_csr_partial(
		const DECIMAL startRowInd, const DECIMAL endRowInd,
		const DECIMAL* RESTRICT rowPtr, const DECIMAL* RESTRICT colInd,
		const REAL* RESTRICT values, const REAL* RESTRICT x, REAL* RESTRICT y)
{

#ifdef __ICC
	__assume_aligned(rowPtr, ALIGNMENT);
	__assume_aligned(colInd, ALIGNMENT);
	__assume_aligned(values, ALIGNMENT);
	__assume_aligned(x, ALIGNMENT);
	__assume_aligned(y, ALIGNMENT);
#endif

	DECIMAL i;
	DECIMAL j;
	for(i = startRowInd; i < endRowInd; ++i)
	{
		REGISTER REAL temp = 0.0;
		#pragma simd reduction(+:temp)
		for(j = rowPtr[i]; j < rowPtr[i + 1]; ++j)
			temp = temp + values[j] * x[colInd[j]];

		y[i] = temp;
	}
}

void spmxv_csr_add(
		const DECIMAL startRowInd, const DECIMAL endRowInd,
		const DECIMAL* RESTRICT rowPtr, const DECIMAL* RESTRICT colInd,
		const REAL* RESTRICT values, const REAL* RESTRICT x, REAL* RESTRICT y)
{

#ifdef __ICC
	__assume_aligned(rowPtr, ALIGNMENT);
	__assume_aligned(colInd, ALIGNMENT);
	__assume_aligned(values, ALIGNMENT);
	__assume_aligned(x, ALIGNMENT);
	__assume_aligned(y, ALIGNMENT);
#endif

	DECIMAL i;
	DECIMAL j;
	for(i = startRowInd; i < endRowInd; ++i)
	{
		REGISTER REAL temp = 0.0;
		#pragma simd reduction(+:temp)
		for(j = rowPtr[i]; j < rowPtr[i + 1]; ++j)
			temp = temp + values[j] * x[colInd[j]];

		y[i] += temp;
	}
}

void spmxv_jds_partial(
		const DECIMAL idiagStart, const DECIMAL idiagEnd,
		const DECIMAL* RESTRICT idiag, const DECIMAL* RESTRICT jdiag,
		const REAL* RESTRICT dj, const DECIMAL* RESTRICT perm,
		const REAL* RESTRICT x, REAL* RESTRICT y)
{

#ifdef __ICC
	__assume_aligned(idiag, ALIGNMENT);
	__assume_aligned(jdiag, ALIGNMENT);
	__assume_aligned(dj, ALIGNMENT);
	__assume_aligned(perm, ALIGNMENT);
	__assume_aligned(x, ALIGNMENT);
	__assume_aligned(y, ALIGNMENT);
#endif

	DECIMAL i;
	DECIMAL j;
	for(i = idiagStart; i < idiagEnd; ++i)
	{
		#pragma simd
		for(j = idiag[i]; j < idiag[i + 1]; ++j)
		{
			int rowIndex = j - idiag[i];
			y[perm[rowIndex]] += dj[j] * x[jdiag[j]];
		}
	}
}

void spmxv_jds_partial_optimized(
		const DECIMAL rowStart,	const DECIMAL idiagStart, const DECIMAL idiagEnd,
		const DECIMAL* RESTRICT idiag, const DECIMAL* RESTRICT jdiag,
		const REAL* RESTRICT dj, const REAL* RESTRICT x, REAL* RESTRICT y)
{

#ifdef __ICC
	__assume_aligned(idiag, ALIGNMENT);
	__assume_aligned(jdiag, ALIGNMENT);
	__assume_aligned(dj, ALIGNMENT);
	__assume_aligned(x, ALIGNMENT);
	__assume_aligned(y, ALIGNMENT);
#endif

	DECIMAL i;
	DECIMAL j;
	for(i = idiagStart; i < idiagEnd; ++i)
	{
		#pragma simd
		for(j = idiag[i]; j < idiag[i + 1]; ++j)
		{
			int rowIndex = j - idiag[i];
			y[rowStart + rowIndex] += dj[j] * x[jdiag[j]];
		}
	}
}

/*
for(i = startRowInd; i < endRowInd; ++i)				\
{														\
	REGISTER REAL temp = 0.0;							\
	for(j = rowPtr[i]; j < rowPtr[i + 1]; ++j)			\
		temp += values[j] * x[colInd[j]];				\
														\
	y_out[i] = temp;									\
}
*/

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
