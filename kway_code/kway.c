#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>
#include "patoh.h"
#define STRINGSIZE 1024
//#define ROWNET 1
//#define COLNET 2
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
extern void printSubMatrices(int* partVec, int partVecLength, int kway, struct mm_data* mtx, char* subMtxPathPrefix, int* outdimVec);
extern void printToConsole(int* arr, int length);
extern void applyOrdering(struct mm_data* mtx, int* ordering);

int ROWNET;

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
	
	
	
	if(ROWNET)
		sscanf(strtmp, "%d %d %d", &nc, &nr, &nnz); //printf("%d %d %d\n", nr, nc, nnz);	
	else
		sscanf(strtmp, "%d %d %d", &nr, &nc, &nnz); //printf("%d %d %d\n", nr, nc, nnz);
	
	mm->N=nr;mm->M=nc;mm->NNZ=nnz;
	mm->x = (int*)calloc(nnz, sizeof(int));	
	mm->y = (int*)calloc(nnz, sizeof(int));	
	mm->v = (double*)calloc(nnz, sizeof(double));	
	int i;
	char inpLine[STRINGSIZE];
	for(i=0; i<nnz; i++){
		int r, c;
		float v = 1;



		if(fgets(inpLine, STRINGSIZE, fp)==NULL){
			printf("Error in reading file\n");
			exit(0);
		}

		if (strlen(inpLine) > STRINGSIZE - 4) {
			printf("ERROR the buffer size for reading lines is too small\n");
			exit(0);
		}
		int readCnt = sscanf(inpLine, "%d %d %f\n", &r, &c, &v);
		if (readCnt < 3)
			v = 1;


		mm->x[i]=r-1;
		mm->y[i]=c-1;
		if(ROWNET){
			mm->x[i]=c-1;
			mm->y[i]=r-1;
		
		}
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
#include <sys/time.h>
double elapsed_time(struct timeval tbegin, struct timeval tend) {

        long sec, msec;

        sec = tend.tv_sec-tbegin.tv_sec;
        msec = tend.tv_usec-tbegin.tv_usec;


        if (msec < 0) {
                sec--;
                msec+=1000000;
        }

        if (msec >= 1000000) {
                sec += msec/1000000;
                msec = msec%1000000;
        }

        double result;

        //printf("%d\n", sec);
        //printf("%d\n", msec);

        result = sec + ((double)msec)/1000000;
        //printf("%.2lf\n", result);

        return result;
}
int main( int argc, char** argv )
{
	int i;
	if(argc < 12) {
		printf("Usage: ./kway mtx_path [ROWNET|COLNET] isSymmetric isBinary K vertexOrder vertexDim netOrder netDim vertexPartVec netPartVec subMtxPathPrefix\n");
		printf("\nQuick start:\nmake;./kway t5-8-11.mtx COLNET 0 0 2  roworder rowdim colorder coldim vertexpartvector netpartvector\n");
		return 1;
	}
	char *inputfilename = argv[1];
//	char* strmethod = argv[2];int method = parseinput(strmethod, "ROWNET|COLNET", ROWNET, COLNET);
	
	if(argv[2][0]=='R')
		ROWNET=1;
	else
		ROWNET=0;
	int isSymmetric = atoi(argv[3]);
	int isBinary = atoi(argv[4]);
	int kway = atoi(argv[5]);
	char* strvo = argv[6];
	char* strvd = argv[7];
	char* strno = argv[8];
	char* strnd = argv[9];
	char* strvp = argv[10];
	char* strnp = argv[11];
	char* subMtxPathPrefix = NULL;
	if(argc > 12)
		 subMtxPathPrefix = argv[12];


	char strtmp[1000];

	int *Lptr, *Lind;	
	double* Lval;	
	struct mm_data *mmL = calloc(1, sizeof(struct mm_data)); sprintf(strtmp, "%s", inputfilename);
	mmL->is_symmetric = isSymmetric;
	mmL->is_binary = isBinary;

	mtxread(mmL, strtmp);

	int nrow = mmL->N;
	int ncol = mmL->M;
	mm2csr(ncol, nrow, mmL->NNZ, mmL->y, mmL->x, mmL->v, &Lptr, &Lind, &Lval);

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
		cwghts[i]=0;
	}
	for(i=0; i<_n;i++){
		int k;
		for(k=xpins[i]; k<xpins[i+1]; k++){
			cwghts[pins[k]]++;
		}
	}
	for(i=0; i<_n;i++){
		nwghts[i]=1;
	}

	PaToH_Alloc(&args, _c, _n, nconst, cwghts, nwghts, xpins, pins);
    struct timeval tbegin, tend;
    gettimeofday(&tbegin, NULL);
	PaToH_Partition(&args, _c, _n, cwghts, nwghts, xpins, pins, partvec, partweights, &cutsize);
	gettimeofday(&tend, NULL);
	double patoh_part_time = elapsed_time(tbegin, tend);

	PaToH_Free();

	int max = 0;
	int sum = 0;
	for(i=0; i<kway; i++){
		if(partweights[i] > max){
			max = partweights[i];
		}
		sum += partweights[i];
	}
	double avg = sum/(kway*1.0f);
	double imb = ((max/avg)-1)*100;

	int* ordering = NULL;
	int* outputDim = NULL;
	convertPartVecToOrdering(partvec, &ordering, &outputDim, _c, kway);

	int ncutnets = 0;
	int* netpartvec = calloc(_n, sizeof(int));
	{ // @KA: find net part vector induced by vertex part vector
	int _k = kway;
	int i;
	//int* map = calloc(_k, sizeof(int));
	//int mark = 1;
	for(i=0; i<_n; i++){
		int j;
		int netpart = -1;
		for(j=xpins[i]; j<xpins[i+1]; j++){
			int p = pins[j];
			int part = partvec[p];
			if(netpart == -1) {
				netpart = part;
			} else if (netpart != part) {
				ncutnets++;
				netpart = _k;
				break;
			}
		}
		if(netpart == -1) { // assign empty nets to part 0
			netpart = 0;
		}
		netpartvec[i] = netpart;
	}
	}
	int* netordering = NULL;
	int* netoutputDim = NULL;
	convertPartVecToOrdering(netpartvec, &netordering, &netoutputDim, _n, kway+1);
	
	printToFile(partvec, _c, strvp);
	printToFile(netpartvec, _n, strnp);

	
	
	if(ROWNET){
		printToFile(netordering, _n, strvo);
		printToFile(netoutputDim, kway+1, strvd);

		printToFile(ordering, _c, strno);
		printToFile(outputDim, kway, strnd);
	}
	else{
		printToFile(netordering, _n, strno);
		printToFile(netoutputDim, kway+1, strnd);

		printToFile(ordering, _c, strvo);
		printToFile(outputDim, kway, strvd);
	}

	
	printf("%s\t", inputfilename);
	printf("%d\t", _c);
	printf("%d\t", _n);
	printf("%d\t", kway);
	printf("%d\t", cutsize);
	printf("%d\t", ncutnets);
	printf("%.0g\t", imb);
	printf("%.3g\t", patoh_part_time);
	printf("\n");

	applyOrdering(mmL, ordering);
	printSubMatrices(partvec, _c, kway, mmL, subMtxPathPrefix, outputDim);

	free(ordering);
	free(outputDim);

	free(Lptr);
	free(Lind);
	free(Lval);

	free(mmL->x);
	free(mmL->y);
	free(mmL->v);

	printf("done!!\n");

	return EXIT_SUCCESS;
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

	fprintf(f, "%d\n", length);

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
		mtx->x[i] = ordering[mtx->x[i]];
	}
}

void printSubMatrices(int* partVec, int partVecLength, int kway, struct mm_data* mtx, char* subMtxPathPrefix, int* outdimVec)
{
	if(subMtxPathPrefix == NULL)
		return;

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
		sprintf(subMtxPath, "%s_%d", subMtxPathPrefix, i);

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

		fprintf(files[i], "%s %s_%d\n", "%", subMtxPathPrefix, i);
		fprintf(files[i], "%d %d %d\n", upper - lower, mtx->M, nnzCountPerSubMtx[i]);

		printf("sub_mtx_%d => rows: %d, NNZ: %d\n", i, upper - lower, nnzCountPerSubMtx[i]);
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
