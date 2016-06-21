
#ifndef SPM_STORAGE_H_
#define SPM_STORAGE_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"

// traditional CSR / CSC storage format
// -------------------------------------------------------------------------

struct spm_cmp
{
	REAL* values;
	DECIMAL* ind;
	DECIMAL* ptr;

	DECIMAL nnz;
	DECIMAL rowCount;
	DECIMAL colCount;
};

typedef struct spm_cmp spm_cmp_t;

// functions for traditional CSR / CSC storage format
// -------------------------------------------------------------------------

extern spm_cmp_t* spm_cmp_new(DECIMAL rowCount, DECIMAL colCount, DECIMAL nnz);
extern void spm_cmp_initDefault(spm_cmp_t* spm);

/**
 * ASSUMPTION: function copies the values of values, ptr, and ind
 * arrays. So they should have matching length in source and
 * destination.
 */
extern void spm_cmp_copy(spm_cmp_t* source, spm_cmp_t* destination);

/**
 * Returns the elements at i th row and j th column of give spm. Returns
 * zero if there is no element at that location.
 */
extern REAL spm_cmp_get(spm_cmp_t* spm, DECIMAL ptrIndex, DECIMAL indIndex);

extern void spm_cmp_print(spm_cmp_t* spm);
extern void spm_cmp_printMultiple(char* m, spm_cmp_t** spms, int numSpms);
extern void spm_cmp_print2DFormat(char* m, spm_cmp_t* spm);
extern void spm_cmp_print2DFormatMultiple(char* m, spm_cmp_t** spm, int numSpms);
extern void spm_cmp_printToFile(char* filePath, char* m, spm_cmp_t* spm);

extern void spm_cmp_delete(spm_cmp_t* spm);
extern void spm_cmp_deleteNonPtr(spm_cmp_t* spm);
extern void spm_cmp_extractNonZeroStatistics_CSR(
		spm_cmp_t* spmCsr, REAL* rowMinNNZ_out,
		REAL* rowAvgNNZ_out, REAL* rowMaxNNZ_out,
		REAL* columnMinNNZ_out, REAL* columnAvgNNZ_out,
		REAL* columnMaxNNZ_out);

extern void spm_cmp_plotShovedLeft(spm_cmp_t* spm);
extern void spm_cmp_plotToBufferShovedLeft(
		spm_cmp_t* spm, char* buff_inout);

extern void spm_cmp_printAttributes(spm_cmp_t* spm);
extern void spm_cmp_toStringAttributes(
		spm_cmp_t* spm, char* buff_inout);
extern DECIMAL spm_cmp_findMaxPtrGap(
		spm_cmp_t* spm, int startPtrIndex, int endPtrIndex);

// ICSR / ICSC storage format
// ---------------------------------------------------------------------------

struct spm_inc
{
	REAL* values;
	DECIMAL* ind;
	DECIMAL* ptr;
	DECIMAL ptrLength;

	DECIMAL nnz;
	DECIMAL rowCount;
	DECIMAL colCount;
};

typedef struct spm_inc spm_inc_t;

// functions for ICSR / ICSC storage format
// ---------------------------------------------------------------------------

extern spm_inc_t* spm_inc_new(
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount, DECIMAL ptrLength);
extern void spm_inc_init(spm_inc_t* spmInc,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount, DECIMAL ptrLength);
extern void spm_inc_initDefault(spm_inc_t* spmInc);
extern void spm_inc_copy(spm_inc_t* source, spm_inc_t* destination);
extern void spm_inc_delete(spm_inc_t* spm);
extern void spm_inc_deleteNonPtr(spm_inc_t* spm);
extern void spm_inc_deleteMultiple(spm_inc_t** spms, int numSpms);
extern void spm_inc_print(spm_inc_t* spmInc);
extern void spm_inc_printMultiple(char* m, spm_inc_t** spmIncs, int numSpms);
extern void spm_inc_printAttributes(spm_inc_t* spm);
extern void spm_inc_toStringAttributes(spm_inc_t* spm, char* buff_inout);

// JDS (jagged diagonal storage) format
// ---------------------------------------------------------------------------

struct spm_jds
{
	REAL* dj;
	DECIMAL* jdiag;
	DECIMAL* idiag;
	DECIMAL* permutation;

	DECIMAL idiagLength;
	DECIMAL nnz;
	DECIMAL rowCount;
	DECIMAL colCount;

	spm_cmp_t* csrCounterpart;
};

// functions for JDS (jagged diagonal storage)
// ---------------------------------------------------------------------------

typedef struct spm_jds spm_jds_t;

extern spm_jds_t* spm_jds_new(
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount, DECIMAL idiagLength);
extern void spm_jds_initDefault(spm_jds_t* spmJds);
extern void spm_jds_copy(spm_jds_t* source, spm_jds_t* destination);
extern void spm_jds_printAttributes(spm_jds_t* spmJds);
extern void spm_jds_toStringAttributes(spm_jds_t* spmJds, char* buff_inout);
extern void spm_jds_print(spm_jds_t* spmJds);
extern void spm_jds_deleteNonPtr(spm_jds_t* spmJds);
extern void spm_jds_delete(spm_jds_t* spmJds);

/**
 * ASSUMPTION: To add a permutation, use this function. It
 * handles aligning. DO NOT use any other means.
 */
extern void spm_jds_addPermutation(spm_jds_t* spmJds, DECIMAL* permutation);
extern void spm_jds_extractNonZeroStatistics_JDS(
		spm_jds_t* spmJds,
		REAL* rowMinNNZ_out, REAL* rowAvgNNZ_out, REAL* rowMaxNNZ_out,
		REAL* columnMinNNZ_out, REAL* columnAvgNNZ_out, REAL* columnMaxNNZ_out);
extern void spm_jds_plot_JDS(spm_jds_t* spmJds);
extern void spm_jds_plotToBuffer_JDS(
		spm_jds_t* spmJds, char* buff_inout);

/**
 * <p>
 * Finds the best suitable cut to partition a given matrix
 * in JDS format into 2 matrices in;
 * <ul>
 * <li>JDS format.</li>
 * <li>CSR format.</li>
 * </ul>
 * in an effort to create better environment for
 * vectorization.
 * </p>
 */
// TODO this is naive implementation, optimization is possible.
extern void spm_jds_findOptimumCut_JDS(
		spm_jds_t* spmJds,
		DECIMAL* optColInd_out, DECIMAL* optRowInd_out);

extern void spm_jds_findOptimumCut_CSR(
		spm_cmp_t* spmCsr,
		DECIMAL* optColInd_out, DECIMAL* optRowInd_out);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SPM_STORAGE_H_ */
