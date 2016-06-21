#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>


#include "patoh.h"
struct mm_data {

	int N, M, NNZ;
	int *x;
	int *y;
	double *v;

	int is_symmetric;
	int is_binary;
	
	int ndiagonal;
	int realnnz;
};

extern void convertPartVecToOrdering(int* partVec, int** ordering_out, int** outputDim_out, int length, int kway);
extern void printToFile(int* arr, int length, char* filename);
extern void printSubMatrices(int* partVec, int partVecLength, int kway, struct mm_data* mtx, char* mtxPath, int* outdimVec);
extern void printToConsole(int* arr, int length);
extern void applyOrdering(struct mm_data* mtx, int* ordering);

#define STRINGSIZE 1024
void mtxread(struct mm_data *mm, char* filename){
	FILE* fp = fopen(filename, "r");
	if(fp == NULL) {
		fprintf(stderr, "%s %d cannot open matrix file %s\n", __FILE__, __LINE__, filename);
		exit(2);
	};
	char strtmp[100];
	char* strread;
	do {
		strread = fgets(strtmp, 100, fp); //printf("%s", strtmp);
		if(strread == NULL){
			fprintf(stderr, "%s %d cannot read matrix file %s\n", __FILE__, __LINE__, filename);
			exit(3);
		}
	} while(strtmp[0] == '%');
	int nr, nc, nnz;
	sscanf(strtmp, "%d %d %d", &nr, &nc, &nnz); //printf("%d %d %d\n", nr, nc, nnz);
	mm->N=nr;mm->M=nc;mm->NNZ=nnz;
	mm->x = (int*)calloc(nnz, sizeof(int));	
	mm->y = (int*)calloc(nnz, sizeof(int));	
	mm->v = (double*)calloc(nnz, sizeof(double));	
	int i;
	for(i=0; i<nnz; i++){
		int r, c;
		float v = 1;

		if(mm->is_binary)
			fscanf(fp, "%d %d", &r, &c, &v); //printf("%d %d %g\n", r, c, v);
		else
			fscanf(fp, "%d %d %g", &r, &c, &v); //printf("%d %d %g\n", r, c, v);

		mm->x[i]=r-1;
		mm->y[i]=c-1;
		mm->v[i]=v;
	}
	fclose(fp);

	// dublicate elements if mtx is symmetric
	if(mm->is_symmetric)
	{
		// find symmetric NNZ
		int newNNZ = mm->NNZ;
		for(i = 0; i < nnz; ++i)
		{
			if(mm->x[i] != mm->y[i])
				++newNNZ;
		}

		int* newX = (int*)calloc(newNNZ, sizeof(int));
		int* newY = (int*)calloc(newNNZ, sizeof(int));
		double* newV = (double*)calloc(newNNZ, sizeof(double));
		for(i = 0; i < mm->NNZ; ++i)
		{
			newX[i] = mm->x[i];
			newY[i] = mm->y[i];
			newV[i] = mm->v[i];
		}

		int index = mm->NNZ;
		for(i = 0; i < mm->NNZ; ++i)
		{
			if(newX[i] != newY[i])
			{
				newX[index] = newY[i];
				newY[index] = newX[i];
				newV[index] = newV[i];
				++index;
			}
		}

		free(mm->x);
		free(mm->y);
		free(mm->v);

		mm->x = newX;
		mm->y = newY;
		mm->v = newV;
		mm->NNZ = newNNZ;
	}
}

int parseinput(char *parameter, const char *format, ...) {

	va_list ap;

	char *saveptr, *token;
	int v, flag = 0;

	va_start(ap, format);

	char *fmt = (char *) malloc(sizeof(char) * STRINGSIZE);
	sprintf(fmt, "%s", format);
	token = strtok_r(fmt, "|", &saveptr);

	while (token != NULL) {

		v = va_arg(ap, int);
		if (!strcmp(parameter, token)) {
			flag = 1;
			break;
		}

		token = strtok_r(NULL, "|", &saveptr);
	}

	va_end(ap);
	free(fmt);

	if (!flag) {}
//		printusage();

	return v;
}

void mm2csr(int nr, int nc, int nnz, int* x, int* y, double* v, int** ptr, int** ind, double** values)
{
	int n = nr;

	int *xpins;  xpins = calloc(nr + 2, sizeof(int));
	int *pins;   pins = calloc(nnz,  sizeof(int));
			double *vals = calloc(nnz,  sizeof(double));

	*ptr = xpins; // IROW
	*ind = pins;  // ICOL
	*values = vals;

	int i;
	for(i=0; i < nnz; i++) {
			xpins[x[i]+2] ++;
	}

	for(i=2; i<=n; i++)
			xpins[i] += xpins[i-1];

	for(i=0; i<nnz; i++) {
			pins[ xpins[ x[i]+1 ] ] = y[i];
			vals[ xpins[ x[i]+1 ] ] = v[i];
							xpins[ x[i]+1 ] ++;
	}
}

int main( int argc, char** argv )
{
	if(argc < 5)
	{
		printf("Usage: ./bin mtx_path isSymmetric isBinary k\n");
		return 1;
	}

//	printf("this is spmxv\n");exit(1);
	register int i, j;
	char order = 0;
//	if(argc != 7) {
//		printusage();
//		return 1;
//	}
	char *inputfilename = argv[1];
	// char* strmethod = argv[2];
#define ROWBASED 1
	int method = parseinput("ROW", "ROW|COL|asd", ROWBASED,2,3);
	int isSymmetric = atoi(argv[2]);
	int isBinary = atoi(argv[3]);
	int kway = atoi(argv[4]);

	char strtmp[1000];
	char strtmp2[1000];

	int *Lptr, *Lind;	
	double* Lval;	
	struct mm_data *mmL = calloc(1, sizeof(struct mm_data)); sprintf(strtmp, "%s", inputfilename);
	mmL->is_symmetric = isSymmetric;
	mmL->is_binary = isBinary;

	mtxread(mmL, strtmp);

	int nrow = mmL->N;
	int ncol = mmL->M;
//printf("%d\n", __LINE__);
	mm2csr(ncol, nrow, mmL->NNZ, mmL->y, mmL->x, mmL->v, &Lptr, &Lind, &Lval);
//printf("%d\n", __LINE__);

/*		for (i=0; i<ncol; i++){
			printf("%d: ", i);
			int k;
			for(k=Lptr[i]; k<Lptr[i+1]; k++) {
				printf("%d:%g ", Lind[k], Lval[k]);
			}
			printf("\n");
		}
*/
	PaToH_Parameters args;

/***** DO NOT CHANGE *****/
//	PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);
	PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_SPEED);
	args.final_imbal = 0.10; // %1
/***** ENDOF DO NOT CHANGE *****/

	args._k = kway;
	int _c = nrow; // vertex
	int _n = ncol; // net
	int nconst = 1;
	int* cwghts = calloc(_c, sizeof(int));
	int* partvec = calloc(_c, sizeof(int));
	int* partweights = calloc(_c, sizeof(int));
	int* nwghts = calloc(_n, sizeof(int));
	int* xpins = Lptr;
	int* pins = Lind;
	int cutsize = 0;
	for(i=0; i<_c;i++){
		cwghts[i]=1;
	}
	for(i=0; i<_n;i++){
		nwghts[i]=1;
	}

	// printf("number of cells:\t%d\n", _c);
	// printf("number of nets:\t%d\n", _n);
	PaToH_Alloc(&args, _c, _n, nconst, cwghts, nwghts, xpins, pins);
	PaToH_Partition(&args, _c, _n, cwghts, nwghts, xpins, pins, partvec, partweights, &cutsize);
	PaToH_Free();

//	printf("cutsize:\t%d\n", cutsize);

//	printf("partvec...\n");
//	printToConsole(partvec, _c);

	int* ordering = NULL;
	int* outputDim = NULL;
	convertPartVecToOrdering(partvec, &ordering, &outputDim, _c, kway);
//	printf("Ordering...\n");
//	printToConsole(ordering, _c);

	char buff[1024];

	sprintf(buff, "%s_ordering", inputfilename);
	printToFile(ordering, _c, buff);

	sprintf(buff, "%s_outputdim", inputfilename);
	printToFile(outputDim, kway, buff);

	applyOrdering(mmL, ordering);
	printSubMatrices(partvec, _c, kway, mmL, inputfilename, outputDim);

	free(ordering);
	free(outputDim);

	free(Lptr);
	free(Lind);
	free(Lval);

	free(mmL->x);
	free(mmL->y);
	free(mmL->v);

	printf("done!!\n");
}

void convertPartVecToOrdering(int* partVec, int** ordering_out, int** outputDim_out, int length, int kway)
{
	int* ordering = (int*) malloc(sizeof(int) * length);
	int* outputDim = (int*) malloc(sizeof(int) * kway);
	int i;
	int j;

	int ind = 0;
	for(i = 0; i < kway; ++i)
	{
		for(j = 0; j < length; ++j)
		{
			if(partVec[j] == i)
				ordering[j] = ind++;
		}

		outputDim[i] = ind;
	}

	for(i = kway - 1 ; 0 < i; --i)
	{
		outputDim[i] = outputDim[i] - outputDim[i - 1];
	}

	*ordering_out = ordering;
	*outputDim_out = outputDim;
}

void printToFile(int* arr, int length, char* filename)
{
	FILE *f = fopen(filename, "w");

	// fprintf(f, "%d\n", length);

	int i;
	for(i = 0; i < length; ++i)
		fprintf(f, "%d\n", arr[i]);

	fclose(f);
}

void printToConsole(int* arr, int length)
{
	int i;
	for(i = 0; i < length; ++i)
		printf("%d\n", arr[i]);
	printf("\n");
}

void applyOrdering(struct mm_data* mtx, int* ordering)
{
	int i;
	for(i = 0; i < mtx->NNZ; ++i)
	{
		// printf("%d => %d to %d\n", i, mtx->x[i], ordering[mtx->x[i]]);
		mtx->x[i] = ordering[mtx->x[i]];
	}
}

void printSubMatrices(int* partVec, int partVecLength, int kway, struct mm_data* mtx, char* mtxPath, int* outdimVec)
{
	int i;
	int j;

	int outdimMinus[kway + 1];
	outdimMinus[0] = 0;
	for(i = 1; i < kway; ++i)
	{
		outdimMinus[i] = outdimMinus[i - 1] + outdimVec[i - 1];
		// printf("%d\n", outdimMinus[i]);
	}
	outdimMinus[i] = outdimMinus[i - 1] + outdimVec[i - 1];

	// open a file stream for each submtx
	FILE** files = (FILE**) malloc(sizeof(FILE*) * kway);
	for(i = 0; i < kway; ++i)
	{
		char subMtxPath[2048];
		subMtxPath[0] = '\0';
		sprintf(subMtxPath, "%s_core%d", mtxPath, i);

		files[i] = fopen(subMtxPath, "w");
	}

	int nnzCountPerSubMtx[kway];
	for(i = 0; i < kway; ++i)
		nnzCountPerSubMtx[i] = 0;

	for(i = 0; i < kway; ++i)
	{
		int lower = outdimMinus[i];
		int upper = outdimMinus[i + 1];

		for(j = 0; j < mtx->NNZ; ++j)
		{
			if(mtx->x[j] >= lower && mtx->x[j] < upper)
			{
				++nnzCountPerSubMtx[i];
			}
		}
	}

	for(i = 0; i < kway; ++i)
	{
		int lower = outdimMinus[i];
		int upper = outdimMinus[i + 1];

		fprintf(files[i], "%s %s_core_%d\n", "%", mtxPath, i);
		fprintf(files[i], "%d %d %d\n", upper - lower, mtx->M, nnzCountPerSubMtx[i]);

		printf("Core_%d => rows: %d, NNZ: %d\n", i, upper - lower, nnzCountPerSubMtx[i]);
	}

	for(i = 0; i < kway; ++i)
	{
		int lower = outdimMinus[i];
		int upper = outdimMinus[i + 1];

		for(j = 0; j < mtx->NNZ; ++j)
		{
			if(mtx->x[j] >= lower && mtx->x[j] < upper)
			{
				if(mtx->is_binary)
					fprintf(files[i], "%d %d\n", mtx->x[j] + 1 - outdimMinus[i], mtx->y[j] + 1);
				else
					fprintf(files[i], "%d %d %f\n", mtx->x[j] + 1 - outdimMinus[i], mtx->y[j] + 1, mtx->v[j]);
			}
		}
	}

/*
	for(j = 0; j < mtx->NNZ; ++j)
	{
		if(mtx->x[j] < outdim)

		int partId = partVec[mtx->x[j]];

		if(mtx->is_binary)
			fprintf(files[partId], "%d %d\n", mtx->x[j] + 1, mtx->y[j] + 1);
		else
			fprintf(files[partId], "%d %d %f\n", mtx->x[j] + 1, mtx->y[j] + 1, mtx->v[j]);
	}
*/

	// close streams
	for(i = 0; i < kway; ++i)
		fclose(files[i]);

	//for(i = 0; i < kway; ++i)
		//free(files[i]);
	free(files);

}
