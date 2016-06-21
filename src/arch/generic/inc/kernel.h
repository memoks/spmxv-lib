
#ifndef KERNEL_H_
#define KERNEL_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"

// CSR based partial macro
// ----------------------------------------------------------------

#define SPMXV_CSR_PARTIAL(									\
	/* INTEGER */ startRowInd, 								\
	/* INTEGER */ endRowInd, 								\
	/* const INTEGER* RESTRICT */ rowPtr, 					\
	/* const INTEGER* RESTRICT */ colInd, 					\
	/* const REAL* RESTRICT */ values, 						\
	/* const REAL* RESTRICT */ x, 							\
	/* REAL* RESTRICT y_out */ y_out)						\
	{														\
		DECIMAL i;											\
		DECIMAL j;											\
		SPMXV_CSR_PARTIAL_FOR_LOOP(							\
				startRowInd,								\
				endRowInd,									\
				rowPtr,										\
				colInd,										\
				values,										\
				x,											\
				y_out)										\
	}


// CSC based partial macro
// ----------------------------------------------------------------

#define SPMXV_CSC_PARTIAL(									\
	/* INTEGER */ startColInd, 								\
	/* INTEGER */ endColInd, 								\
	/* const INTEGER* RESTRICT */ colPtr, 					\
	/* const INTEGER* RESTRICT */ rowInd, 					\
	/* const REAL* RESTRICT */ values, 						\
	/* const REAL* RESTRICT */ x, 							\
	/* REAL* RESTRICT */ y_out)								\
	{														\
		DECIMAL i;											\
		DECIMAL j;											\
		for(i = startColInd; i < endColInd; ++i)			\
		{													\
			REGISTER REAL xi = x[i];						\
			for(j = colPtr[i]; j < colPtr[i + 1]; ++j)		\
				y_out[rowInd[j]] += values[j] * xi;			\
				/* y_out[rowInd[j]] += values[j] * x[i]; */	\
		}													\
	}



// Column parallel functions
// Designed to be used on borders which requires locking
// ----------------------------------------------------------------

#define SPMXV_CSR_PARTIAL_ATOMIC(							\
	/* INTEGER */ startRowInd, 								\
	/* INTEGER */ endRowInd, 								\
	/* const INTEGER* RESTRICT */ rowPtr, 					\
	/* const INTEGER* RESTRICT */ colInd, 					\
	/* const REAL* RESTRICT */ values, 						\
	/* const REAL* RESTRICT */ x, 							\
	/* REAL* RESTRICT */ y_out)								\
	{														\
		DECIMAL i;											\
		DECIMAL j;											\
		for(i = startRowInd; i < endRowInd; ++i)			\
		{													\
			REGISTER REAL temp = 0.0;						\
			for(j = rowPtr[i]; j < rowPtr[i + 1]; ++j)		\
				temp += values[j] * x[colInd[j]];			\
															\
			PRAGMA_OMP_ATOMIC								\
			y_out[i] += temp;								\
		}													\
	}

// Designed to be used inside omp_parallel_for pragma
// indexes i and j must be declared beforehand and
// marked as private inside the loop
// Sample;
// INTEGER i;
// INTEGER j;
// pragma omp parallel for private(i, j)
// SPMXV_CSR_PARTIAL_FOR_LOOP(...)
// ----------------------------------------------------------------

#define SPMXV_CSR_PARTIAL_FOR_LOOP(							\
	/* INTEGER */ startRowInd, 								\
	/* INTEGER */ endRowInd, 								\
	/* const INTEGER* RESTRICT */ rowPtr, 					\
	/* const INTEGER* RESTRICT */ colInd, 					\
	/* const REAL* RESTRICT */ values, 						\
	/* const REAL* RESTRICT */ x, 							\
	/* REAL* RESTRICT y_out */ y_out)						\
	for(i = startRowInd; i < endRowInd; ++i)				\
	{														\
		REGISTER REAL temp = 0.0;							\
		for(j = rowPtr[i]; j < rowPtr[i + 1]; ++j)			\
			temp += values[j] * x[colInd[j]];				\
															\
		y_out[i] = temp;									\
	}


// Designed to be used inside omp_parallel_for pragma
// indexes i and j must be declared beforehand and
// marked as private inside the loop
// Sample;
// INTEGER i;
// INTEGER j;
// pragma omp parallel for private(i, j)
// SPMXV_CSC_PARTIAL_FOR_LOOP_ATOMIC(...)
// ----------------------------------------------------------------

#define SPMXV_CSC_PARTIAL_FOR_LOOP_ATOMIC(					\
	/* INTEGER */ startColInd, 								\
	/* INTEGER */ endColInd, 								\
	/* const INTEGER* RESTRICT */ colPtr, 					\
	/* const INTEGER* RESTRICT */ rowInd, 					\
	/* const REAL* RESTRICT */ values, 						\
	/* const REAL* RESTRICT */ x, 							\
	/* REAL* RESTRICT */ y_out)								\
	for(i = startColInd; i < endColInd; ++i)				\
	{														\
		REGISTER REAL xi = x[i];							\
		for(j = colPtr[i]; j < colPtr[i + 1]; ++j)			\
		{													\
			PRAGMA_OMP_ATOMIC								\
			y_out[rowInd[j]] += values[j] * xi;				\
															\
			/* y_out[rowInd[j]] += values[j] * x[i]; */		\
		}													\
	}

#define SPMXV_ICSR(											\
	/* const INTEGER */ rowCount,							\
	/* const INTEGER */ colCount,							\
	/* const INTEGER */ nnz,								\
	/* const INTEGER* RESTRICT */ ptr,	 					\
	/* const INTEGER* RESTRICT */ inc,	 					\
	/* const REAL* RESTRICT */ values, 						\
	/* const REAL* RESTRICT */ x, 							\
	/* REAL* RESTRICT y_out */ y)							\
	{														\
		DECIMAL r = 0;	/* rowInc index */					\
		DECIMAL t = 0;	/* non-zero index */				\
		DECIMAL i = 0;										\
		DECIMAL j = inc[t];									\
		while(t < nnz)										\
		{													\
			REGISTER REAL temp = 0.0;						\
															\
			i = i + ptr[r];									\
			while(j < colCount)								\
			{												\
				/* y[i] = y[i] + values[t] * x[j]; */		\
				temp = temp + values[t] * x[j];				\
															\
				t = t + 1;									\
				j = j + inc[t];								\
			}												\
															\
			y[i] = temp;									\
															\
			j = j - colCount;								\
			r = r + 1;										\
		}													\
	}

#define SPMXV_ICSR_ATOMIC(									\
	/* const INTEGER */ rowCount,							\
	/* const INTEGER */ colCount,							\
	/* const INTEGER */ nnz,								\
	/* const INTEGER* RESTRICT */ ptr,	 					\
	/* const INTEGER* RESTRICT */ inc,	 					\
	/* const REAL* RESTRICT */ values, 						\
	/* const REAL* RESTRICT */ x, 							\
	/* REAL* RESTRICT y */ y)								\
	{														\
		DECIMAL r = 0;	/* rowInc index */					\
		DECIMAL t = 0;	/* non-zero index */				\
		DECIMAL i = 0;										\
		DECIMAL j = inc[t];									\
		while(t < nnz)										\
		{													\
			REGISTER REAL temp = 0.0;						\
															\
			i = i + ptr[r];									\
			while(j < colCount)								\
			{												\
				/* y[i] = y[i] + values[t] * x[j]; */		\
				temp = temp + values[t] * x[j];				\
															\
				t = t + 1;									\
				j = j + inc[t];								\
			}												\
															\
			PRAGMA_OMP_ATOMIC								\
			y[i] += temp;									\
															\
			j = j - colCount;								\
			r = r + 1;										\
		}													\
	}

// JDS (jagged diagonal storage) format SPMXV macro
// ----------------------------------------------------------------

#define SPMXV_JDS_PARTIAL(											\
	/* const INTEGER */ idiagStart,									\
	/* const INTEGER */ idiagEnd,									\
	/* const INTEGER* RESTRICT */ idiag,							\
	/* const INTEGER* RESTRICT */ jdiag,							\
	/* REAL* RESTRICT */ dj,										\
	/* const INTEGER* RESTRICT */ perm,								\
	/* const REAL* RESTRICT */ x,									\
	/* REAL* RESTRICT */ y_inout)									\
																	\
	DECIMAL i;														\
	DECIMAL j;														\
	for(i = idiagStart; i < idiagEnd; ++i)							\
	{																\
		PRAGMA_IVDEP												\
		PRAGMA_SIMD													\
		for(j = idiag[i]; j < idiag[i + 1];	++j)					\
		{															\
			int rowIndex = j - idiag[i];							\
			y_inout[perm[rowIndex]] += dj[j] * x[jdiag[j]];			\
		}															\
	}

// To enable vectorization
// -----------------------------------------------------------------------

extern void spmxv_csr_partial(
		const DECIMAL startRowInd, const DECIMAL endRowInd,
		const DECIMAL* RESTRICT rowPtr, const DECIMAL* RESTRICT colInd,
		const REAL* RESTRICT values, const REAL* RESTRICT x,
		REAL* RESTRICT y);

extern void spmxv_csr_add(
		const DECIMAL startRowInd, const DECIMAL endRowInd,
		const DECIMAL* RESTRICT rowPtr, const DECIMAL* RESTRICT colInd,
		const REAL* RESTRICT values, const REAL* RESTRICT x,
		REAL* RESTRICT y);

extern void spmxv_jds_partial(
		const DECIMAL idiagStart, const DECIMAL idiagEnd,
		const DECIMAL* RESTRICT idiag, const DECIMAL* RESTRICT jdiag,
		const REAL* RESTRICT dj, const DECIMAL* RESTRICT perm,
		const REAL* RESTRICT x, REAL* RESTRICT y);

extern void spmxv_jds_partial_optimized(
		const DECIMAL rowStart,
		const DECIMAL idiagStart, const DECIMAL idiagEnd,
		const DECIMAL* RESTRICT idiag, const DECIMAL* RESTRICT jdiag,
		const REAL* RESTRICT dj,
		const REAL* RESTRICT x, REAL* RESTRICT y);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* KERNEL_H_ */
