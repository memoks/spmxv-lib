
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "include/config.h"
#include "include/io/input_parser.h"
#include "include/io/converter.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/sub_mtx.h"
#include "include/scheduler/job_batch.h"

// all
int main(int argsc, char* argsv[])
{
	{
		sub_mtx_dim_t* subMtx = sub_mtx_new();
		sub_mtx_init(subMtx, 0, 0, 10, 10);

		spm_cmp_t* spmCsr = spm_cmp_new(8, 118, 30);
		spmCsr->values[0] = 0;		spmCsr->ind[0] = 0;		spmCsr->ptr[0] = 0;
		spmCsr->values[1] = 7;		spmCsr->ind[1] = 1;		spmCsr->ptr[1] = 7;
		spmCsr->values[2] = 14;		spmCsr->ind[2] = 2;		spmCsr->ptr[2] = 13;
		spmCsr->values[3] = 20;		spmCsr->ind[3] = 3;		spmCsr->ptr[3] = 17;
		spmCsr->values[4] = 25;		spmCsr->ind[4] = 4;		spmCsr->ptr[4] = 21;
		spmCsr->values[5] = 27;		spmCsr->ind[5] = 5;		spmCsr->ptr[5] = 25;
		spmCsr->values[6] = 29;		spmCsr->ind[6] = 6;		spmCsr->ptr[6] = 28;
		spmCsr->values[7] = 1;		spmCsr->ind[7] = 0;		spmCsr->ptr[7] = 30;
		spmCsr->values[8] = 8;		spmCsr->ind[8] = 1;		spmCsr->ptr[8] = 30;
		spmCsr->values[9] = 15;		spmCsr->ind[9] = 2;
		spmCsr->values[10] = 21;	spmCsr->ind[10] = 3;
		spmCsr->values[11] = 26;	spmCsr->ind[11] = 4;
		spmCsr->values[12] = 28;	spmCsr->ind[12] = 5;
		spmCsr->values[13] = 2;		spmCsr->ind[13] = 0;
		spmCsr->values[14] = 9;		spmCsr->ind[14] = 1;
		spmCsr->values[15] = 16;	spmCsr->ind[15] = 2;
		spmCsr->values[16] = 22;	spmCsr->ind[16] = 3;
		spmCsr->values[17] = 3;		spmCsr->ind[17] = 0;
		spmCsr->values[18] = 10;	spmCsr->ind[18] = 1;
		spmCsr->values[19] = 17;	spmCsr->ind[19] = 2;
		spmCsr->values[20] = 23;	spmCsr->ind[20] = 3;
		spmCsr->values[21] = 4;		spmCsr->ind[21] = 0;
		spmCsr->values[22] = 11;	spmCsr->ind[22] = 1;
		spmCsr->values[23] = 18;	spmCsr->ind[23] = 2;
		spmCsr->values[24] = 24;	spmCsr->ind[24] = 3;
		spmCsr->values[25] = 5;		spmCsr->ind[25] = 0;
		spmCsr->values[26] = 12;	spmCsr->ind[26] = 1;
		spmCsr->values[27] = 19;	spmCsr->ind[27] = 2;
		spmCsr->values[28] = 6;		spmCsr->ind[28] = 0;
		spmCsr->values[29] = 13;	spmCsr->ind[29] = 1;

		spm_jds_t* spmJds = spm_jds_new(30, 8, 118, 7);
		spmJds->dj[0] = 0;		spmJds->jdiag[0] = 0;		spmJds->idiag[0] = 0;
		spmJds->dj[1] = 1;		spmJds->jdiag[1] = 0;		spmJds->idiag[1] = 7;
		spmJds->dj[2] = 2;		spmJds->jdiag[2] = 0;		spmJds->idiag[2] = 14;
		spmJds->dj[3] = 3;		spmJds->jdiag[3] = 0;		spmJds->idiag[3] = 20;
		spmJds->dj[4] = 4;		spmJds->jdiag[4] = 0;		spmJds->idiag[4] = 25;
		spmJds->dj[5] = 5;		spmJds->jdiag[5] = 0;		spmJds->idiag[5] = 27;
		spmJds->dj[6] = 6;		spmJds->jdiag[6] = 0;		spmJds->idiag[6] = 29;
		spmJds->dj[7] = 7;		spmJds->jdiag[7] = 1;		spmJds->idiag[7] = 30;
		spmJds->dj[8] = 8;		spmJds->jdiag[8] = 1;
		spmJds->dj[9] = 9;		spmJds->jdiag[9] = 1;
		spmJds->dj[10] = 10;	spmJds->jdiag[10] = 1;
		spmJds->dj[11] = 11;	spmJds->jdiag[11] = 1;
		spmJds->dj[12] = 12;	spmJds->jdiag[12] = 1;
		spmJds->dj[13] = 13;	spmJds->jdiag[13] = 1;
		spmJds->dj[14] = 14;	spmJds->jdiag[14] = 2;
		spmJds->dj[15] = 15;	spmJds->jdiag[15] = 2;
		spmJds->dj[16] = 16;	spmJds->jdiag[16] = 2;
		spmJds->dj[17] = 17;	spmJds->jdiag[17] = 2;
		spmJds->dj[18] = 18;	spmJds->jdiag[18] = 2;
		spmJds->dj[19] = 19;	spmJds->jdiag[19] = 2;
		spmJds->dj[20] = 20;	spmJds->jdiag[20] = 3;
		spmJds->dj[21] = 21;	spmJds->jdiag[21] = 3;
		spmJds->dj[22] = 22;	spmJds->jdiag[22] = 3;
		spmJds->dj[23] = 23;	spmJds->jdiag[23] = 3;
		spmJds->dj[24] = 24;	spmJds->jdiag[24] = 3;
		spmJds->dj[25] = 25;	spmJds->jdiag[25] = 4;
		spmJds->dj[26] = 26;	spmJds->jdiag[26] = 4;
		spmJds->dj[27] = 27;	spmJds->jdiag[27] = 5;
		spmJds->dj[28] = 28;	spmJds->jdiag[28] = 5;
		spmJds->dj[29] = 29;	spmJds->jdiag[29] = 6;

		DECIMAL* permutation = (DECIMAL*) malloc(sizeof(DECIMAL) * spmJds->rowCount);
		permutation[0] = 6;
		permutation[1] = 5;
		permutation[2] = 8;
		permutation[3] = 0;
		permutation[4] = 1;
		permutation[5] = 3;
		permutation[6] = 2;

		spm_jds_addPermutation(spmJds, permutation);

		spmJds->csrCounterpart = spmCsr;

		PRINTF("Whole JDS...\n");
		spm_jds_print(spmJds);
		spm_jds_plot_JDS(spmJds);

		PRINTF("CSR counterpart...\n");
		spm_cmp_print(spmCsr);

		// conversion
		// ----------------------------------------------------------------------------------

		job_batch_t* partialJdsBatch = job_batch_newPartialJds(subMtx, spmJds);
		job_batch_t* hybridJdsCsrBatch = job_batch_new();

		converter_extractHybridJDSCSR(partialJdsBatch, hybridJdsCsrBatch);

		PRINTF("Plotting partitioned JDS...\n");
		spm_jds_print(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);
		spm_jds_plot_JDS(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);

		PRINTF("Plotting partitioned CSR...\n");
		spm_cmp_print(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);
		spm_cmp_plotShovedLeft(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);

		sub_mtx_delete(subMtx);
		free(permutation);
		job_batch_delete(partialJdsBatch);
		job_batch_delete(hybridJdsCsrBatch);
	}
	PRINTF("=============================================================\n");
	{
		sub_mtx_dim_t* subMtx = sub_mtx_new();
		sub_mtx_init(subMtx, 0, 0, 10, 10);

		spm_cmp_t* spmCsr = spm_cmp_new(7, 118, 29);
		spmCsr->values[0] = 0;		spmCsr->ind[0] = 0;		spmCsr->ptr[0] = 0;
		spmCsr->values[1] = 7;		spmCsr->ind[1] = 1;		spmCsr->ptr[1] = 5;
		spmCsr->values[2] = 14;		spmCsr->ind[2] = 2;		spmCsr->ptr[2] = 10;
		spmCsr->values[3] = 21;		spmCsr->ind[3] = 3;		spmCsr->ptr[3] = 15;
		spmCsr->values[4] = 26;		spmCsr->ind[4] = 4;		spmCsr->ptr[4] = 19;
		spmCsr->values[5] = 1;		spmCsr->ind[5] = 0;		spmCsr->ptr[5] = 23;
		spmCsr->values[6] = 8;		spmCsr->ind[6] = 1;		spmCsr->ptr[6] = 26;
		spmCsr->values[7] = 15;		spmCsr->ind[7] = 2;		spmCsr->ptr[7] = 29;
		spmCsr->values[8] = 22;		spmCsr->ind[8] = 3;
		spmCsr->values[9] = 27;		spmCsr->ind[9] = 4;
		spmCsr->values[10] = 2;		spmCsr->ind[10] = 0;
		spmCsr->values[11] = 9;		spmCsr->ind[11] = 1;
		spmCsr->values[12] = 16;	spmCsr->ind[12] = 2;
		spmCsr->values[13] = 23;	spmCsr->ind[13] = 3;
		spmCsr->values[14] = 28;	spmCsr->ind[14] = 4;
		spmCsr->values[15] = 3;		spmCsr->ind[15] = 0;
		spmCsr->values[16] = 10;	spmCsr->ind[16] = 1;
		spmCsr->values[17] = 17;	spmCsr->ind[17] = 2;
		spmCsr->values[18] = 24;	spmCsr->ind[18] = 3;
		spmCsr->values[19] = 4;		spmCsr->ind[19] = 0;
		spmCsr->values[20] = 11;	spmCsr->ind[20] = 1;
		spmCsr->values[21] = 18;	spmCsr->ind[21] = 2;
		spmCsr->values[22] = 25;	spmCsr->ind[22] = 3;
		spmCsr->values[23] = 5;		spmCsr->ind[23] = 0;
		spmCsr->values[24] = 12;	spmCsr->ind[24] = 1;
		spmCsr->values[25] = 19;	spmCsr->ind[25] = 2;
		spmCsr->values[26] = 6;		spmCsr->ind[26] = 0;
		spmCsr->values[27] = 13;	spmCsr->ind[27] = 1;
		spmCsr->values[28] = 20;	spmCsr->ind[28] = 2;


		spm_cmp_plotShovedLeft(spmCsr);

		spm_jds_t* spmJds = spm_jds_new(29, 7, 118, 6);
		spmJds->dj[0] = 0;		spmJds->jdiag[0] = 0;		spmJds->idiag[0] = 0;
		spmJds->dj[1] = 1;		spmJds->jdiag[1] = 0;		spmJds->idiag[1] = 7;
		spmJds->dj[2] = 2;		spmJds->jdiag[2] = 0;		spmJds->idiag[2] = 14;
		spmJds->dj[3] = 3;		spmJds->jdiag[3] = 0;		spmJds->idiag[3] = 21;
		spmJds->dj[4] = 4;		spmJds->jdiag[4] = 0;		spmJds->idiag[4] = 26;
		spmJds->dj[5] = 5;		spmJds->jdiag[5] = 0;		spmJds->idiag[5] = 29;
		spmJds->dj[6] = 6;		spmJds->jdiag[6] = 0;
		spmJds->dj[7] = 7;		spmJds->jdiag[7] = 1;
		spmJds->dj[8] = 8;		spmJds->jdiag[8] = 1;
		spmJds->dj[9] = 9;		spmJds->jdiag[9] = 1;
		spmJds->dj[10] = 10;	spmJds->jdiag[10] = 1;
		spmJds->dj[11] = 11;	spmJds->jdiag[11] = 1;
		spmJds->dj[12] = 12;	spmJds->jdiag[12] = 1;
		spmJds->dj[13] = 13;	spmJds->jdiag[13] = 1;
		spmJds->dj[14] = 14;	spmJds->jdiag[14] = 2;
		spmJds->dj[15] = 15;	spmJds->jdiag[15] = 2;
		spmJds->dj[16] = 16;	spmJds->jdiag[16] = 2;
		spmJds->dj[17] = 17;	spmJds->jdiag[17] = 2;
		spmJds->dj[18] = 18;	spmJds->jdiag[18] = 2;
		spmJds->dj[19] = 19;	spmJds->jdiag[19] = 2;
		spmJds->dj[20] = 20;	spmJds->jdiag[20] = 2;
		spmJds->dj[21] = 21;	spmJds->jdiag[21] = 3;
		spmJds->dj[22] = 22;	spmJds->jdiag[22] = 3;
		spmJds->dj[23] = 23;	spmJds->jdiag[23] = 3;
		spmJds->dj[24] = 24;	spmJds->jdiag[24] = 3;
		spmJds->dj[25] = 25;	spmJds->jdiag[25] = 3;
		spmJds->dj[26] = 26;	spmJds->jdiag[26] = 4;
		spmJds->dj[27] = 27;	spmJds->jdiag[27] = 4;
		spmJds->dj[28] = 28;	spmJds->jdiag[28] = 4;


		DECIMAL* permutation = (DECIMAL*) malloc(sizeof(DECIMAL) * spmJds->rowCount);
		permutation[0] = 6;
		permutation[1] = 5;
		permutation[2] = 8;
		permutation[3] = 0;
		permutation[4] = 1;
		permutation[5] = 3;
		permutation[6] = 2;

		spm_jds_addPermutation(spmJds, permutation);

		spmJds->csrCounterpart = spmCsr;

		PRINTF("Whole JDS...\n");
		spm_jds_print(spmJds);
		spm_jds_plot_JDS(spmJds);

		PRINTF("CSR counterpart...\n");
		spm_cmp_print(spmCsr);

		// conversion
		// ----------------------------------------------------------------------------------

		job_batch_t* partialJdsBatch = job_batch_newPartialJds(subMtx, spmJds);
		job_batch_t* hybridJdsCsrBatch = job_batch_new();

		converter_extractHybridJDSCSR(partialJdsBatch, hybridJdsCsrBatch);

		PRINTF("Plotting partitioned JDS...\n");
		spm_jds_print(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);
		spm_jds_plot_JDS(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);

		PRINTF("Plotting partitioned CSR...\n");
		spm_cmp_print(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);
		spm_cmp_plotShovedLeft(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);

		sub_mtx_delete(subMtx);
		free(permutation);
		job_batch_delete(partialJdsBatch);
		job_batch_delete(hybridJdsCsrBatch);
	}
	PRINTF("=============================================================\n");
	{
		sub_mtx_dim_t* subMtx = sub_mtx_new();
		sub_mtx_init(subMtx, 0, 0, 10, 10);

		spm_cmp_t* spmCsr = spm_cmp_new(10, 10, 55);
		spmCsr->values[0] = 0;		spmCsr->ind[0] = 0;		spmCsr->ptr[0] = 0;
		spmCsr->values[1] = 10;		spmCsr->ind[1] = 1;		spmCsr->ptr[1] = 10;
		spmCsr->values[2] = 19;		spmCsr->ind[2] = 2;		spmCsr->ptr[2] = 19;
		spmCsr->values[3] = 27;		spmCsr->ind[3] = 3;		spmCsr->ptr[3] = 27;
		spmCsr->values[4] = 34;		spmCsr->ind[4] = 4;		spmCsr->ptr[4] = 34;
		spmCsr->values[5] = 40;		spmCsr->ind[5] = 5;		spmCsr->ptr[5] = 40;
		spmCsr->values[6] = 45;		spmCsr->ind[6] = 6;		spmCsr->ptr[6] = 45;
		spmCsr->values[7] = 49;		spmCsr->ind[7] = 7;		spmCsr->ptr[7] = 49;
		spmCsr->values[8] = 52;		spmCsr->ind[8] = 8;		spmCsr->ptr[8] = 52;
		spmCsr->values[9] = 54;		spmCsr->ind[9] = 9;		spmCsr->ptr[9] = 54;
		spmCsr->values[10] = 1;		spmCsr->ind[10] = 10;	spmCsr->ptr[10] = 55;
		spmCsr->values[11] = 11;	spmCsr->ind[11] = 11;
		spmCsr->values[12] = 20;	spmCsr->ind[12] = 12;
		spmCsr->values[13] = 28;	spmCsr->ind[13] = 13;
		spmCsr->values[14] = 35;	spmCsr->ind[14] = 14;
		spmCsr->values[15] = 41;	spmCsr->ind[15] = 15;
		spmCsr->values[16] = 46;	spmCsr->ind[16] = 16;
		spmCsr->values[17] = 50;	spmCsr->ind[17] = 17;
		spmCsr->values[18] = 53;	spmCsr->ind[18] = 18;
		spmCsr->values[19] = 2;		spmCsr->ind[19] = 19;
		spmCsr->values[20] = 12;	spmCsr->ind[20] = 20;
		spmCsr->values[21] = 21;	spmCsr->ind[21] = 21;
		spmCsr->values[22] = 29;	spmCsr->ind[22] = 22;
		spmCsr->values[23] = 36;	spmCsr->ind[23] = 23;
		spmCsr->values[24] = 42;	spmCsr->ind[24] = 24;
		spmCsr->values[25] = 47;	spmCsr->ind[25] = 25;
		spmCsr->values[26] = 51;	spmCsr->ind[26] = 26;
		spmCsr->values[27] = 3;		spmCsr->ind[27] = 27;
		spmCsr->values[28] = 13;	spmCsr->ind[28] = 28;
		spmCsr->values[29] = 22;	spmCsr->ind[29] = 29;
		spmCsr->values[30] = 30;	spmCsr->ind[30] = 30;
		spmCsr->values[31] = 37;	spmCsr->ind[31] = 31;
		spmCsr->values[32] = 43;	spmCsr->ind[32] = 32;
		spmCsr->values[33] = 48;	spmCsr->ind[33] = 33;
		spmCsr->values[34] = 4;		spmCsr->ind[34] = 34;
		spmCsr->values[35] = 14;	spmCsr->ind[35] = 35;
		spmCsr->values[36] = 23;	spmCsr->ind[36] = 36;
		spmCsr->values[37] = 31;	spmCsr->ind[37] = 37;
		spmCsr->values[38] = 38;	spmCsr->ind[38] = 38;
		spmCsr->values[39] = 44;	spmCsr->ind[39] = 39;
		spmCsr->values[40] = 5;		spmCsr->ind[40] = 40;
		spmCsr->values[41] = 15;	spmCsr->ind[41] = 41;
		spmCsr->values[42] = 24;	spmCsr->ind[42] = 42;
		spmCsr->values[43] = 32;	spmCsr->ind[43] = 43;
		spmCsr->values[44] = 39;	spmCsr->ind[44] = 44;
		spmCsr->values[45] = 6;		spmCsr->ind[45] = 45;
		spmCsr->values[46] = 16;	spmCsr->ind[46] = 46;
		spmCsr->values[47] = 25;	spmCsr->ind[47] = 47;
		spmCsr->values[48] = 33;	spmCsr->ind[48] = 48;
		spmCsr->values[49] = 7;		spmCsr->ind[49] = 49;
		spmCsr->values[50] = 17;	spmCsr->ind[50] = 50;
		spmCsr->values[51] = 26;	spmCsr->ind[51] = 51;
		spmCsr->values[52] = 8;		spmCsr->ind[52] = 52;
		spmCsr->values[53] = 18;	spmCsr->ind[53] = 53;
		spmCsr->values[54] = 9;		spmCsr->ind[54] = 54;


		spm_jds_t* spmJds = spm_jds_new(55, 10, 10, 11);
		spmJds->dj[0] = 0;		spmJds->jdiag[0] = 0;		spmJds->idiag[0] = 0;
		spmJds->dj[1] = 1;		spmJds->jdiag[1] = 0;		spmJds->idiag[1] = 10;
		spmJds->dj[2] = 2;		spmJds->jdiag[2] = 0;		spmJds->idiag[2] = 19;
		spmJds->dj[3] = 3;		spmJds->jdiag[3] = 0;		spmJds->idiag[3] = 27;
		spmJds->dj[4] = 4;		spmJds->jdiag[4] = 0;		spmJds->idiag[4] = 34;
		spmJds->dj[5] = 5;		spmJds->jdiag[5] = 0;		spmJds->idiag[5] = 40;
		spmJds->dj[6] = 6;		spmJds->jdiag[6] = 0;		spmJds->idiag[6] = 45;
		spmJds->dj[7] = 7;		spmJds->jdiag[7] = 0;		spmJds->idiag[7] = 49;
		spmJds->dj[8] = 8;		spmJds->jdiag[8] = 0;		spmJds->idiag[8] = 52;
		spmJds->dj[9] = 9;		spmJds->jdiag[9] = 0;		spmJds->idiag[9] = 54;
		spmJds->dj[10] = 10;	spmJds->jdiag[10] = 1;		spmJds->idiag[10] = 55;
		spmJds->dj[11] = 11;	spmJds->jdiag[11] = 1;
		spmJds->dj[12] = 12;	spmJds->jdiag[12] = 1;
		spmJds->dj[13] = 13;	spmJds->jdiag[13] = 1;
		spmJds->dj[14] = 14;	spmJds->jdiag[14] = 1;
		spmJds->dj[15] = 15;	spmJds->jdiag[15] = 1;
		spmJds->dj[16] = 16;	spmJds->jdiag[16] = 1;
		spmJds->dj[17] = 17;	spmJds->jdiag[17] = 1;
		spmJds->dj[18] = 18;	spmJds->jdiag[18] = 1;
		spmJds->dj[19] = 19;	spmJds->jdiag[19] = 2;
		spmJds->dj[20] = 20;	spmJds->jdiag[20] = 2;
		spmJds->dj[21] = 21;	spmJds->jdiag[21] = 2;
		spmJds->dj[22] = 22;	spmJds->jdiag[22] = 2;
		spmJds->dj[23] = 23;	spmJds->jdiag[23] = 2;
		spmJds->dj[24] = 24;	spmJds->jdiag[24] = 2;
		spmJds->dj[25] = 25;	spmJds->jdiag[25] = 2;
		spmJds->dj[26] = 26;	spmJds->jdiag[26] = 2;
		spmJds->dj[27] = 27;	spmJds->jdiag[27] = 3;
		spmJds->dj[28] = 28;	spmJds->jdiag[28] = 3;
		spmJds->dj[29] = 29;	spmJds->jdiag[29] = 3;
		spmJds->dj[30] = 30;	spmJds->jdiag[30] = 3;
		spmJds->dj[31] = 31;	spmJds->jdiag[31] = 3;
		spmJds->dj[32] = 32;	spmJds->jdiag[32] = 3;
		spmJds->dj[33] = 33;	spmJds->jdiag[33] = 3;
		spmJds->dj[34] = 34;	spmJds->jdiag[34] = 4;
		spmJds->dj[35] = 35;	spmJds->jdiag[35] = 4;
		spmJds->dj[36] = 36;	spmJds->jdiag[36] = 4;
		spmJds->dj[37] = 37;	spmJds->jdiag[37] = 4;
		spmJds->dj[38] = 38;	spmJds->jdiag[38] = 4;
		spmJds->dj[39] = 39;	spmJds->jdiag[39] = 4;
		spmJds->dj[40] = 40;	spmJds->jdiag[40] = 5;
		spmJds->dj[41] = 41;	spmJds->jdiag[41] = 5;
		spmJds->dj[42] = 42;	spmJds->jdiag[42] = 5;
		spmJds->dj[43] = 43;	spmJds->jdiag[43] = 5;
		spmJds->dj[44] = 44;	spmJds->jdiag[44] = 5;
		spmJds->dj[45] = 45;	spmJds->jdiag[45] = 6;
		spmJds->dj[46] = 46;	spmJds->jdiag[46] = 6;
		spmJds->dj[47] = 47;	spmJds->jdiag[47] = 6;
		spmJds->dj[48] = 48;	spmJds->jdiag[48] = 6;
		spmJds->dj[49] = 49;	spmJds->jdiag[49] = 7;
		spmJds->dj[50] = 50;	spmJds->jdiag[50] = 7;
		spmJds->dj[51] = 51;	spmJds->jdiag[51] = 7;
		spmJds->dj[52] = 52;	spmJds->jdiag[52] = 8;
		spmJds->dj[53] = 53;	spmJds->jdiag[53] = 8;
		spmJds->dj[54] = 54;	spmJds->jdiag[54] = 9;

		DECIMAL* permutation = (DECIMAL*) malloc(sizeof(DECIMAL) * spmJds->rowCount);
		permutation[0] = 6;
		permutation[1] = 5;
		permutation[2] = 8;
		permutation[3] = 0;
		permutation[4] = 1;
		permutation[5] = 3;
		permutation[6] = 2;
		permutation[7] = 4;
		permutation[8] = 9;
		permutation[9] = 7;

		spm_jds_addPermutation(spmJds, permutation);

		spmJds->csrCounterpart = spmCsr;

		PRINTF("Whole JDS...\n");
		spm_jds_print(spmJds);
		spm_jds_plot_JDS(spmJds);

		PRINTF("CSR counterpart...\n");
		spm_cmp_print(spmCsr);

		// conversion
		// ----------------------------------------------------------------------------------

		job_batch_t* partialJdsBatch = job_batch_newPartialJds(subMtx, spmJds);
		job_batch_t* hybridJdsCsrBatch = job_batch_new();

		converter_extractHybridJDSCSR(partialJdsBatch, hybridJdsCsrBatch);

		PRINTF("Plotting partitioned JDS...\n");
		spm_jds_print(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);
		spm_jds_plot_JDS(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);

		PRINTF("Plotting partitioned CSR...\n");
		spm_cmp_print(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);
		spm_cmp_plotShovedLeft(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);

		sub_mtx_delete(subMtx);
		free(permutation);
		job_batch_delete(partialJdsBatch);
		job_batch_delete(hybridJdsCsrBatch);
	}
	PRINTF("=============================================================\n");
	{
		sub_mtx_dim_t* subMtx = sub_mtx_new();
		sub_mtx_init(subMtx, 0, 0, 10, 10);

		spm_cmp_t* spmCsr = spm_cmp_new(10, 10, 15);
		spmCsr->values[0] = 0;		spmCsr->ind[0] = 0;		spmCsr->ptr[0] = 0;
		spmCsr->values[1] = 9;		spmCsr->ind[1] = 1;		spmCsr->ptr[1] = 3;
		spmCsr->values[2] = 14;		spmCsr->ind[2] = 2;		spmCsr->ptr[2] = 5;
		spmCsr->values[3] = 1;		spmCsr->ind[3] = 0;		spmCsr->ptr[3] = 7;
		spmCsr->values[4] = 10;		spmCsr->ind[4] = 1;		spmCsr->ptr[4] = 9;
		spmCsr->values[5] = 2;		spmCsr->ind[5] = 0;		spmCsr->ptr[5] = 11;
		spmCsr->values[6] = 11;		spmCsr->ind[6] = 1;		spmCsr->ptr[6] = 12;
		spmCsr->values[7] = 3;		spmCsr->ind[7] = 0;		spmCsr->ptr[7] = 13;
		spmCsr->values[8] = 12;		spmCsr->ind[8] = 1;		spmCsr->ptr[8] = 14;
		spmCsr->values[9] = 4;		spmCsr->ind[9] = 0;		spmCsr->ptr[9] = 15;
		spmCsr->values[10] = 13;	spmCsr->ind[10] = 1;
		spmCsr->values[11] = 5;		spmCsr->ind[11] = 0;
		spmCsr->values[12] = 6;		spmCsr->ind[12] = 0;
		spmCsr->values[13] = 7;		spmCsr->ind[13] = 0;
		spmCsr->values[14] = 8;		spmCsr->ind[14] = 0;


		spm_jds_t* spmJds = spm_jds_new(15, 10, 10, 4);
		spmJds->dj[0] = 0;		spmJds->jdiag[0] = 0;		spmJds->idiag[0] = 0;
		spmJds->dj[1] = 1;		spmJds->jdiag[1] = 0;		spmJds->idiag[1] = 9;
		spmJds->dj[2] = 2;		spmJds->jdiag[2] = 0;		spmJds->idiag[2] = 14;
		spmJds->dj[3] = 3;		spmJds->jdiag[3] = 0;		spmJds->idiag[3] = 15;
		spmJds->dj[4] = 4;		spmJds->jdiag[4] = 0;
		spmJds->dj[5] = 5;		spmJds->jdiag[5] = 0;
		spmJds->dj[6] = 6;		spmJds->jdiag[6] = 0;
		spmJds->dj[7] = 7;		spmJds->jdiag[7] = 0;
		spmJds->dj[8] = 8;		spmJds->jdiag[8] = 0;
		spmJds->dj[9] = 9;		spmJds->jdiag[9] = 1;
		spmJds->dj[10] = 10;	spmJds->jdiag[10] = 1;
		spmJds->dj[11] = 11;	spmJds->jdiag[11] = 1;
		spmJds->dj[12] = 12;	spmJds->jdiag[12] = 1;
		spmJds->dj[13] = 13;	spmJds->jdiag[13] = 1;
		spmJds->dj[14] = 14;	spmJds->jdiag[14] = 2;

		DECIMAL* permutation = (DECIMAL*) malloc(sizeof(DECIMAL) * spmJds->rowCount);
		permutation[0] = 6;
		permutation[1] = 5;
		permutation[2] = 8;
		permutation[3] = 0;
		permutation[4] = 1;
		permutation[5] = 3;
		permutation[6] = 2;
		permutation[7] = 4;
		permutation[8] = 9;
		permutation[9] = 7;

		spm_jds_addPermutation(spmJds, permutation);

		spmJds->csrCounterpart = spmCsr;

		PRINTF("Whole JDS...\n");
		spm_jds_print(spmJds);
		spm_jds_plot_JDS(spmJds);

		PRINTF("CSR counterpart...\n");
		spm_cmp_print(spmCsr);

		// conversion
		// ----------------------------------------------------------------------------------

		job_batch_t* partialJdsBatch = job_batch_newPartialJds(subMtx, spmJds);
		job_batch_t* hybridJdsCsrBatch = job_batch_new();

		converter_extractHybridJDSCSR(partialJdsBatch, hybridJdsCsrBatch);

		PRINTF("Plotting partitioned JDS...\n");
		spm_jds_print(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);
		spm_jds_plot_JDS(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);

		PRINTF("Plotting partitioned CSR...\n");
		spm_cmp_print(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);
		spm_cmp_plotShovedLeft(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);

		sub_mtx_delete(subMtx);
		free(permutation);
		job_batch_delete(partialJdsBatch);
		job_batch_delete(hybridJdsCsrBatch);
	}
	PRINTF("=============================================================\n");
	{
		sub_mtx_dim_t* subMtx = sub_mtx_new();
		sub_mtx_init(subMtx, 0, 0, 10, 10);

		spm_cmp_t* spmCsr = spm_cmp_new(10, 10, 17);
		spmCsr->values[0] = 0;		spmCsr->ind[0] = 0;		spmCsr->ptr[0] = 0;
		spmCsr->values[1] = 3;		spmCsr->ind[1] = 1;		spmCsr->ptr[1] = 9;
		spmCsr->values[2] = 5;		spmCsr->ind[2] = 2;		spmCsr->ptr[2] = 16;
		spmCsr->values[3] = 7;		spmCsr->ind[3] = 3;		spmCsr->ptr[3] = 17;
		spmCsr->values[4] = 9;		spmCsr->ind[4] = 4;
		spmCsr->values[5] = 11;		spmCsr->ind[5] = 5;
		spmCsr->values[6] = 13;		spmCsr->ind[6] = 6;
		spmCsr->values[7] = 15;		spmCsr->ind[7] = 7;
		spmCsr->values[8] = 16;		spmCsr->ind[8] = 8;
		spmCsr->values[9] = 1;		spmCsr->ind[9] = 0;
		spmCsr->values[10] = 4;		spmCsr->ind[10] = 1;
		spmCsr->values[11] = 6;		spmCsr->ind[11] = 2;
		spmCsr->values[12] = 8;		spmCsr->ind[12] = 3;
		spmCsr->values[13] = 10;	spmCsr->ind[13] = 4;
		spmCsr->values[14] = 12;	spmCsr->ind[14] = 5;
		spmCsr->values[15] = 14;	spmCsr->ind[15] = 6;
		spmCsr->values[16] = 2;		spmCsr->ind[16] = 0;


		spm_jds_t* spmJds = spm_jds_new(17, 10, 10, 9);
		spmJds->dj[0] = 0;		spmJds->jdiag[0] = 0;		spmJds->idiag[0] = 0;
		spmJds->dj[1] = 1;		spmJds->jdiag[1] = 0;		spmJds->idiag[1] = 3;
		spmJds->dj[2] = 2;		spmJds->jdiag[2] = 0;		spmJds->idiag[2] = 5;
		spmJds->dj[3] = 3;		spmJds->jdiag[3] = 1;		spmJds->idiag[3] = 7;
		spmJds->dj[4] = 4;		spmJds->jdiag[4] = 1;		spmJds->idiag[4] = 9;
		spmJds->dj[5] = 5;		spmJds->jdiag[5] = 2;		spmJds->idiag[5] = 11;
		spmJds->dj[6] = 6;		spmJds->jdiag[6] = 2;		spmJds->idiag[6] = 13;
		spmJds->dj[7] = 7;		spmJds->jdiag[7] = 3;		spmJds->idiag[7] = 15;
		spmJds->dj[8] = 8;		spmJds->jdiag[8] = 3;		spmJds->idiag[8] = 16;
		spmJds->dj[9] = 9;		spmJds->jdiag[9] = 4;		spmJds->idiag[9] = 17;
		spmJds->dj[10] = 10;	spmJds->jdiag[10] = 4;
		spmJds->dj[11] = 11;	spmJds->jdiag[11] = 5;
		spmJds->dj[12] = 12;	spmJds->jdiag[12] = 5;
		spmJds->dj[13] = 13;	spmJds->jdiag[13] = 6;
		spmJds->dj[14] = 14;	spmJds->jdiag[14] = 6;
		spmJds->dj[15] = 15;	spmJds->jdiag[15] = 7;
		spmJds->dj[16] = 16;	spmJds->jdiag[16] = 8;

		DECIMAL* permutation = (DECIMAL*) malloc(sizeof(DECIMAL) * spmJds->rowCount);
		permutation[0] = 6;
		permutation[1] = 5;
		permutation[2] = 8;
		permutation[3] = 0;
		permutation[4] = 1;
		permutation[5] = 3;
		permutation[6] = 2;
		permutation[7] = 4;
		permutation[8] = 9;
		permutation[9] = 7;

		spm_jds_addPermutation(spmJds, permutation);

		spmJds->csrCounterpart = spmCsr;

		PRINTF("Whole JDS...\n");
		spm_jds_print(spmJds);
		spm_jds_plot_JDS(spmJds);

		PRINTF("CSR counterpart...\n");
		spm_cmp_print(spmCsr);

		// conversion
		// ----------------------------------------------------------------------------------

		job_batch_t* partialJdsBatch = job_batch_newPartialJds(subMtx, spmJds);
		job_batch_t* hybridJdsCsrBatch = job_batch_new();

		converter_extractHybridJDSCSR(partialJdsBatch, hybridJdsCsrBatch);

		PRINTF("Plotting partitioned JDS...\n");
		spm_jds_print(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);
		spm_jds_plot_JDS(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);

		PRINTF("Plotting partitioned CSR...\n");
		spm_cmp_print(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);
		spm_cmp_plotShovedLeft(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);

		sub_mtx_delete(subMtx);
		free(permutation);
		job_batch_delete(partialJdsBatch);
		job_batch_delete(hybridJdsCsrBatch);
	}
	PRINTF("=============================================================\n");
	{
		sub_mtx_dim_t* subMtx = sub_mtx_new();
		sub_mtx_init(subMtx, 0, 0, 10, 9);

		spm_cmp_t* spmCsr = spm_cmp_new(10, 10, 18);
		spmCsr->values[0] = 0;		spmCsr->ind[0] = 0;		spmCsr->ptr[0] = 0;
		spmCsr->values[1] = 3;		spmCsr->ind[1] = 1;		spmCsr->ptr[1] = 9;
		spmCsr->values[2] = 6;		spmCsr->ind[2] = 2;		spmCsr->ptr[2] = 15;
		spmCsr->values[3] = 9;		spmCsr->ind[3] = 3;		spmCsr->ptr[3] = 18;
		spmCsr->values[4] = 11;		spmCsr->ind[4] = 4;
		spmCsr->values[5] = 13;		spmCsr->ind[5] = 5;
		spmCsr->values[6] = 15;		spmCsr->ind[6] = 6;
		spmCsr->values[7] = 16;		spmCsr->ind[7] = 7;
		spmCsr->values[8] = 17;		spmCsr->ind[8] = 8;
		spmCsr->values[9] = 1;		spmCsr->ind[9] = 0;
		spmCsr->values[10] = 4;		spmCsr->ind[10] = 1;
		spmCsr->values[11] = 7;		spmCsr->ind[11] = 2;
		spmCsr->values[12] = 10;	spmCsr->ind[12] = 3;
		spmCsr->values[13] = 12;	spmCsr->ind[13] = 4;
		spmCsr->values[14] = 14;	spmCsr->ind[14] = 5;
		spmCsr->values[15] = 2;		spmCsr->ind[15] = 0;
		spmCsr->values[16] = 5;		spmCsr->ind[16] = 1;
		spmCsr->values[17] = 8;		spmCsr->ind[17] = 2;


		spm_jds_t* spmJds = spm_jds_new(18, 10, 10, 9);
		spmJds->dj[0] = 0;		spmJds->jdiag[0] = 0;		spmJds->idiag[0] = 0;
		spmJds->dj[1] = 1;		spmJds->jdiag[1] = 0;		spmJds->idiag[1] = 3;
		spmJds->dj[2] = 2;		spmJds->jdiag[2] = 0;		spmJds->idiag[2] = 6;
		spmJds->dj[3] = 3;		spmJds->jdiag[3] = 1;		spmJds->idiag[3] = 9;
		spmJds->dj[4] = 4;		spmJds->jdiag[4] = 1;		spmJds->idiag[4] = 11;
		spmJds->dj[5] = 5;		spmJds->jdiag[5] = 1;		spmJds->idiag[5] = 13;
		spmJds->dj[6] = 6;		spmJds->jdiag[6] = 2;		spmJds->idiag[6] = 15;
		spmJds->dj[7] = 7;		spmJds->jdiag[7] = 2;		spmJds->idiag[7] = 16;
		spmJds->dj[8] = 8;		spmJds->jdiag[8] = 2;		spmJds->idiag[8] = 17;
		spmJds->dj[9] = 9;		spmJds->jdiag[9] = 3;		spmJds->idiag[9] = 18;
		spmJds->dj[10] = 10;	spmJds->jdiag[10] = 3;
		spmJds->dj[11] = 11;	spmJds->jdiag[11] = 4;
		spmJds->dj[12] = 12;	spmJds->jdiag[12] = 4;
		spmJds->dj[13] = 13;	spmJds->jdiag[13] = 5;
		spmJds->dj[14] = 14;	spmJds->jdiag[14] = 5;
		spmJds->dj[15] = 15;	spmJds->jdiag[15] = 6;
		spmJds->dj[16] = 16;	spmJds->jdiag[16] = 7;
		spmJds->dj[17] = 17;	spmJds->jdiag[17] = 8;

		DECIMAL* permutation = (DECIMAL*) malloc(sizeof(DECIMAL) * spmJds->rowCount);
		permutation[0] = 6;
		permutation[1] = 5;
		permutation[2] = 8;
		permutation[3] = 0;
		permutation[4] = 1;
		permutation[5] = 3;
		permutation[6] = 2;
		permutation[7] = 4;
		permutation[8] = 9;
		permutation[9] = 7;

		spm_jds_addPermutation(spmJds, permutation);

		spmJds->csrCounterpart = spmCsr;

		PRINTF("Plotting whole JDS...\n");
		spm_jds_print(spmJds);
		spm_jds_plot_JDS(spmJds);

		PRINTF("CSR counterpart...\n");
		spm_cmp_print(spmCsr);

		// conversion
		// ----------------------------------------------------------------------------------

		job_batch_t* partialJdsBatch = job_batch_newPartialJds(subMtx, spmJds);
		job_batch_t* hybridJdsCsrBatch = job_batch_new();

		converter_extractHybridJDSCSR(partialJdsBatch, hybridJdsCsrBatch);

		PRINTF("Plotting partitioned JDS...\n");
		spm_jds_print(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);
		spm_jds_plot_JDS(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);

		PRINTF("Plotting partitioned CSR...\n");
		spm_cmp_print(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);
		spm_cmp_plotShovedLeft(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);

		sub_mtx_delete(subMtx);
		free(permutation);
		job_batch_delete(partialJdsBatch);
		job_batch_delete(hybridJdsCsrBatch);
	}
	PRINTF("=============================================================\n");
	{
		sub_mtx_dim_t* subMtx = sub_mtx_new();
		sub_mtx_init(subMtx, 0, 0, 10, 10);

		spm_cmp_t* spmCsr = spm_cmp_new(10, 10, 13);
		spmCsr->values[0] = 0;	spmCsr->ind[0] = 0;		spmCsr->ptr[0] = 0;
		spmCsr->values[1] = 2;	spmCsr->ind[1] = 1;		spmCsr->ptr[1] = 8;
		spmCsr->values[2] = 4;	spmCsr->ind[2] = 2;		spmCsr->ptr[2] = 13;
		spmCsr->values[3] = 6;	spmCsr->ind[3] = 3;
		spmCsr->values[4] = 8;	spmCsr->ind[4] = 4;
		spmCsr->values[5] = 10;	spmCsr->ind[5] = 5;
		spmCsr->values[6] = 11;	spmCsr->ind[6] = 6;
		spmCsr->values[7] = 12;	spmCsr->ind[7] = 7;
		spmCsr->values[8] = 1;	spmCsr->ind[8] = 0;
		spmCsr->values[9] = 3;	spmCsr->ind[9] = 1;
		spmCsr->values[10] = 5;	spmCsr->ind[10] = 2;
		spmCsr->values[11] = 7;	spmCsr->ind[11] = 3;
		spmCsr->values[12] = 9;	spmCsr->ind[12] = 4;

		spm_jds_t* spmJds = spm_jds_new(13, 10, 10, 8);
		spmJds->dj[0] = 0;		spmJds->jdiag[0] = 0;		spmJds->idiag[0] = 0;
		spmJds->dj[1] = 1;		spmJds->jdiag[1] = 0;		spmJds->idiag[1] = 2;
		spmJds->dj[2] = 2;		spmJds->jdiag[2] = 1;		spmJds->idiag[2] = 4;
		spmJds->dj[3] = 3;		spmJds->jdiag[3] = 1;		spmJds->idiag[3] = 6;
		spmJds->dj[4] = 4;		spmJds->jdiag[4] = 2;		spmJds->idiag[4] = 8;
		spmJds->dj[5] = 5;		spmJds->jdiag[5] = 2;		spmJds->idiag[5] = 10;
		spmJds->dj[6] = 6;		spmJds->jdiag[6] = 3;		spmJds->idiag[6] = 11;
		spmJds->dj[7] = 7;		spmJds->jdiag[7] = 3;		spmJds->idiag[7] = 12;
		spmJds->dj[8] = 8;		spmJds->jdiag[8] = 4;		spmJds->idiag[8] = 13;
		spmJds->dj[9] = 9;		spmJds->jdiag[9] = 4;
		spmJds->dj[10] = 10;	spmJds->jdiag[10] = 5;
		spmJds->dj[11] = 11;	spmJds->jdiag[11] = 6;
		spmJds->dj[12] = 12;	spmJds->jdiag[12] = 7;

		DECIMAL* permutation = (DECIMAL*) malloc(sizeof(DECIMAL) * spmJds->rowCount);
		permutation[0] = 6;
		permutation[1] = 5;
		permutation[2] = 8;
		permutation[3] = 0;
		permutation[4] = 1;
		permutation[5] = 3;
		permutation[6] = 2;
		permutation[7] = 4;
		permutation[8] = 9;
		permutation[9] = 7;

		spm_jds_addPermutation(spmJds, permutation);

		spmJds->csrCounterpart = spmCsr;

		PRINTF("Whole JDS...\n");
		spm_jds_print(spmJds);
		spm_jds_plot_JDS(spmJds);

		PRINTF("CSR counterpart...\n");
		spm_cmp_print(spmCsr);

		// conversion
		// ----------------------------------------------------------------------------------

		job_batch_t* partialJdsBatch = job_batch_newPartialJds(subMtx, spmJds);
		job_batch_t* hybridJdsCsrBatch = job_batch_new();

		converter_extractHybridJDSCSR(partialJdsBatch, hybridJdsCsrBatch);

		PRINTF("Plotting partitioned JDS...\n");
		spm_jds_print(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);
		spm_jds_plot_JDS(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);

		PRINTF("Plotting partitioned CSR...\n");
		spm_cmp_print(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);
		spm_cmp_plotShovedLeft(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);

		sub_mtx_delete(subMtx);
		free(permutation);
		job_batch_delete(partialJdsBatch);
		job_batch_delete(hybridJdsCsrBatch);
	}
	PRINTF("=============================================================\n");
	{
		sub_mtx_dim_t* subMtx = sub_mtx_new();
		sub_mtx_init(subMtx, 0, 0, 10, 10);

		spm_cmp_t* spmCsr = spm_cmp_new(10, 10, 14);
		spmCsr->values[0] = 0.0;		spmCsr->ind[0] = 0;		spmCsr->ptr[0] = 0;
		spmCsr->values[1] = 9.0;		spmCsr->ind[1] = 1;		spmCsr->ptr[1] = 2;
		spmCsr->values[2] = 1.0;		spmCsr->ind[2] = 0;		spmCsr->ptr[2] = 4;
		spmCsr->values[3] = 10.0;		spmCsr->ind[3] = 1;		spmCsr->ptr[3] = 6;
		spmCsr->values[4] = 2.0;		spmCsr->ind[4] = 0;		spmCsr->ptr[4] = 8;
		spmCsr->values[5] = 11.0;		spmCsr->ind[5] = 1;		spmCsr->ptr[5] = 9;
		spmCsr->values[6] = 3.0;		spmCsr->ind[6] = 0;		spmCsr->ptr[6] = 10;
		spmCsr->values[7] = 12.0;		spmCsr->ind[7] = 1;		spmCsr->ptr[7] = 11;
		spmCsr->values[8] = 4.0;		spmCsr->ind[8] = 0;		spmCsr->ptr[8] = 12;
		spmCsr->values[9] = 13.0;		spmCsr->ind[9] = 1;		spmCsr->ptr[9] = 13;
		spmCsr->values[10] = 5.0;		spmCsr->ind[10] = 0;	spmCsr->ptr[10] = 14;
		spmCsr->values[11] = 6.0;		spmCsr->ind[11] = 1;
		spmCsr->values[12] = 7.0;		spmCsr->ind[12] = 0;
		spmCsr->values[13] = 8.0;		spmCsr->ind[13] = 1;

		spm_jds_t* spmJds = spm_jds_new(14, 10, 10, 2);
		spmJds->dj[0] = 0;		spmJds->jdiag[0] = 0;	spmJds->idiag[0] = 0;
		spmJds->dj[1] = 1;		spmJds->jdiag[1] = 0;	spmJds->idiag[1] = 9;
		spmJds->dj[2] = 2;		spmJds->jdiag[2] = 0;	spmJds->idiag[2] = 14;
		spmJds->dj[3] = 3;		spmJds->jdiag[3] = 0;
		spmJds->dj[4] = 4;		spmJds->jdiag[4] = 0;
		spmJds->dj[5] = 5;		spmJds->jdiag[5] = 0;
		spmJds->dj[6] = 6;		spmJds->jdiag[6] = 0;
		spmJds->dj[7] = 7;		spmJds->jdiag[7] = 0;
		spmJds->dj[8] = 8;		spmJds->jdiag[8] = 0;
		spmJds->dj[9] = 9;		spmJds->jdiag[9] = 1;
		spmJds->dj[10] = 10;	spmJds->jdiag[10] = 1;
		spmJds->dj[11] = 11;	spmJds->jdiag[11] = 1;
		spmJds->dj[12] = 12;	spmJds->jdiag[12] = 1;
		spmJds->dj[13] = 13;	spmJds->jdiag[13] = 1;

		DECIMAL* permutation = (DECIMAL*) malloc(sizeof(DECIMAL) * spmJds->rowCount);
		permutation[0] = 6;
		permutation[1] = 5;
		permutation[2] = 8;
		permutation[3] = 0;
		permutation[4] = 1;
		permutation[5] = 3;
		permutation[6] = 2;
		permutation[7] = 4;
		permutation[8] = 9;
		permutation[9] = 7;

		spm_jds_addPermutation(spmJds, permutation);

		spmJds->csrCounterpart = spmCsr;

		PRINTF("Plotting whole JDS...\n");
		spm_jds_print(spmJds);
		spm_jds_plot_JDS(spmJds);

		PRINTF("CSR counterpart...\n");
		spm_cmp_print(spmCsr);

		// conversion
		// ----------------------------------------------------------------------------------

		job_batch_t* partialJdsBatch = job_batch_newPartialJds(subMtx, spmJds);
		job_batch_t* hybridJdsCsrBatch = job_batch_new();

		converter_extractHybridJDSCSR(partialJdsBatch, hybridJdsCsrBatch);

		PRINTF("Plotting partitioned JDS...\n");
		spm_jds_print(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);
		spm_jds_plot_JDS(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);

		PRINTF("Plotting partitioned CSR...\n");
		spm_cmp_print(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);
		spm_cmp_plotShovedLeft(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);


		sub_mtx_delete(subMtx);
		free(permutation);
		job_batch_delete(partialJdsBatch);
		job_batch_delete(hybridJdsCsrBatch);
	}
	PRINTF("=============================================================\n");
	{
		sub_mtx_dim_t* subMtx = sub_mtx_new();
		sub_mtx_init(subMtx, 0, 0, 10, 10);

		spm_cmp_t* spmCsr = spm_cmp_new(10, 10, 25);
		spmCsr->values[0] = 1.0;	spmCsr->ind[0] = 0;		spmCsr->ptr[0] = 0;
		spmCsr->values[1] = 10.0;	spmCsr->ind[1] = 1;		spmCsr->ptr[1] = 8;
		spmCsr->values[2] = 16.0;	spmCsr->ind[2] = 4;		spmCsr->ptr[2] = 13;
		spmCsr->values[3] = 19.0;	spmCsr->ind[3] = 5;		spmCsr->ptr[3] = 16;
		spmCsr->values[4] = 21.0;	spmCsr->ind[4] = 6;		spmCsr->ptr[4] = 18;
		spmCsr->values[5] = 23.0;	spmCsr->ind[5] = 7;		spmCsr->ptr[5] = 20;
		spmCsr->values[6] = 24.0;	spmCsr->ind[6] = 8;		spmCsr->ptr[6] = 22;
		spmCsr->values[7] = 25.0;	spmCsr->ind[7] = 9;		spmCsr->ptr[7] = 23;
		spmCsr->values[8] = 2.0;	spmCsr->ind[8] = 1;		spmCsr->ptr[8] = 24;
		spmCsr->values[9] = 11.0;	spmCsr->ind[9] = 2;		spmCsr->ptr[9] = 25;
		spmCsr->values[10] = 17.0;	spmCsr->ind[10] = 5;	spmCsr->ptr[10] = 25;
		spmCsr->values[11] = 20.0;	spmCsr->ind[11] = 8;
		spmCsr->values[12] = 22.0;	spmCsr->ind[12] = 9;
		spmCsr->values[13] = 3.0;	spmCsr->ind[13] = 1;
		spmCsr->values[14] = 12.0;	spmCsr->ind[14] = 3;
		spmCsr->values[15] = 18.0;	spmCsr->ind[15] = 4;
		spmCsr->values[16] = 4.0;	spmCsr->ind[16] = 1;
		spmCsr->values[17] = 13.0;	spmCsr->ind[17] = 7;
		spmCsr->values[18] = 5.0;	spmCsr->ind[18] = 4;
		spmCsr->values[19] = 14.0;	spmCsr->ind[19] = 9;
		spmCsr->values[20] = 6.0;	spmCsr->ind[20] = 1;
		spmCsr->values[21] = 15.0;	spmCsr->ind[21] = 4;
		spmCsr->values[22] = 7.0;	spmCsr->ind[22] = 7;
		spmCsr->values[23] = 8.0;	spmCsr->ind[23] = 7;
		spmCsr->values[24] = 9.0;	spmCsr->ind[24] = 6;

		spm_jds_t* spmJds = spm_jds_new(25, 10, 10, 8);
		spmJds->dj[0] = 1.0;	spmJds->jdiag[0] = 0;	spmJds->idiag[0] = 0;
		spmJds->dj[1] = 2.0;	spmJds->jdiag[1] = 1;	spmJds->idiag[1] = 9;
		spmJds->dj[2] = 3.0;	spmJds->jdiag[2] = 4;	spmJds->idiag[2] = 15;
		spmJds->dj[3] = 4.0;	spmJds->jdiag[3] = 5;	spmJds->idiag[3] = 18;
		spmJds->dj[4] = 5.0;	spmJds->jdiag[4] = 6;	spmJds->idiag[4] = 20;
		spmJds->dj[5] = 6.0;	spmJds->jdiag[5] = 7;	spmJds->idiag[5] = 22;
		spmJds->dj[6] = 7.0;	spmJds->jdiag[6] = 8;	spmJds->idiag[6] = 23;
		spmJds->dj[7] = 8.0;	spmJds->jdiag[7] = 9;	spmJds->idiag[7] = 24;
		spmJds->dj[8] = 9.0;	spmJds->jdiag[8] = 1;	spmJds->idiag[8] = 25;
		spmJds->dj[9] = 10.0;	spmJds->jdiag[9] = 2;
		spmJds->dj[10] = 11.0;	spmJds->jdiag[10] = 5;
		spmJds->dj[11] = 12.0;	spmJds->jdiag[11] = 8;
		spmJds->dj[12] = 13.0;	spmJds->jdiag[12] = 9;
		spmJds->dj[13] = 14.0;	spmJds->jdiag[13] = 1;
		spmJds->dj[14] = 15.0;	spmJds->jdiag[14] = 3;
		spmJds->dj[15] = 16.0;	spmJds->jdiag[15] = 4;
		spmJds->dj[16] = 17.0;	spmJds->jdiag[16] = 1;
		spmJds->dj[17] = 18.0;	spmJds->jdiag[17] = 7;
		spmJds->dj[18] = 19.0;	spmJds->jdiag[18] = 4;
		spmJds->dj[19] = 20.0;	spmJds->jdiag[19] = 9;
		spmJds->dj[20] = 21.0;	spmJds->jdiag[20] = 1;
		spmJds->dj[21] = 22.0;	spmJds->jdiag[21] = 4;
		spmJds->dj[22] = 23.0;	spmJds->jdiag[22] = 7;
		spmJds->dj[23] = 24.0;	spmJds->jdiag[23] = 7;
		spmJds->dj[24] = 25.0;	spmJds->jdiag[24] = 6;

		DECIMAL* permutation = (DECIMAL*) malloc(sizeof(DECIMAL) * spmJds->rowCount);
		permutation[0] = 6;
		permutation[1] = 5;
		permutation[2] = 8;
		permutation[3] = 0;
		permutation[4] = 1;
		permutation[5] = 3;
		permutation[6] = 2;
		permutation[7] = 4;
		permutation[8] = 9;
		permutation[9] = 7;

		spm_jds_addPermutation(spmJds, permutation);

		spmJds->csrCounterpart = spmCsr;

		PRINTF("Whole JDS...\n");
		spm_jds_print(spmJds);
		spm_jds_plot_JDS(spmJds);

		PRINTF("CSR counterpart...\n");
		spm_cmp_print(spmCsr);

		// conversion
		// ----------------------------------------------------------------------------------

		job_batch_t* partialJdsBatch = job_batch_newPartialJds(subMtx, spmJds);
		job_batch_t* hybridJdsCsrBatch = job_batch_new();

		converter_extractHybridJDSCSR(partialJdsBatch, hybridJdsCsrBatch);

		PRINTF("Plotting partitioned JDS...\n");
		spm_jds_print(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);
		spm_jds_plot_JDS(hybridJdsCsrBatch->data.hybrid_jds_csr.jds);

		PRINTF("Plotting partitioned CSR...\n");
		spm_cmp_print(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);
		spm_cmp_plotShovedLeft(hybridJdsCsrBatch->data.hybrid_jds_csr.csr);

		sub_mtx_delete(subMtx);
		free(permutation);
		job_batch_delete(partialJdsBatch);
		job_batch_delete(hybridJdsCsrBatch);
	}

	return EXIT_SUCCESS;
}

/*
int main(int argsc, char* argsv[])
{

	return EXIT_SUCCESS;
}
*/
