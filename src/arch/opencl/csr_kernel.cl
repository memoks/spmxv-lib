
// Kernel code for M x M for 1 element of output matrix
__kernel
void opencl_csr(
		job_batch_t* batchArrPerBlock,
		REAL* x,
		REAL* y)
{
	// Get global position for X direction
	int column = get_global_id(0);
	
	sub_mtx_dim_t* subMtx = &jb->data.subMtx;

	float sum = 0.0f;
	// Calculate result of one element of Matrix C
	for(int i = 0; i < widthA; ++i)
		sum += inputA[row * widthA + i] * inputB[widthB * i + column];

	outputC[row * widthB + column] = sum;
}
