#include <cuda.h>
#include <stdio.h>
#include "mex.h"

__global__ void addEle(float* A, float* B, float* C){

	int i;
	i = threadIdx.x;
	C[i] = A[i] + B[i];


	return;
}

void addVec(float* h_A, float* h_B, float* h_C, int size){
	int dataSize, i;
	float *d_A, *d_B, *d_C;

	dataSize = sizeof(float)*size;

   printf("Initializing data...\n");

		/* Don't need because I am allocating it in MATLAB
			but if I were only doing C I would be all over this

      printf("...allocating CPU memory.\n");
	  h_A     = (float *)malloc(DATA_SZ);
      h_B     = (float *)malloc(DATA_SZ);
      h_C_CPU = (float *)malloc(RESULT_SZ);
      h_C_GPU = (float *)malloc(RESULT_SZ);	*/

      printf("...allocating GPU memory.\n");
      cudaMalloc((void **)&d_A, dataSize);
      cudaMalloc((void **)&d_B, dataSize);
      cudaMalloc((void **)&d_C, dataSize);

      printf("...copying input data to GPU mem.\n");
		//Copy options data to GPU memory for further processing 
      // COPY FROM CPU to GPU
      cudaMemcpy(d_A, h_A, dataSize, cudaMemcpyHostToDevice);
      cudaMemcpy(d_B, h_B, dataSize, cudaMemcpyHostToDevice);

	printf("Data init done.\n");

	cudaThreadSynchronize();		// you can synchronize, not synchronize, demas?
	addEle<<<1, size>>>(d_A, d_B, d_C);
	cudaThreadSynchronize();

 	printf("Reading back GPU result...\n");
   	//Read back GPU results 
      cudaMemcpy(h_C, d_C, dataSize, cudaMemcpyDeviceToHost);
	
	for (i=0; i<size; i++){
		mexPrintf("%f %f %f\n", h_A[i], h_B[i], h_C[i]);
	}

	printf("Shutting down...\n");
      cudaFree(d_C);
      cudaFree(d_B);
      cudaFree(d_A);

	cudaThreadExit();


	return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *A, *B, *C;
	float *h_A, *h_B, *h_C;
	int size, i, dataSize;
	

	if (nrhs!=2){
 		mexErrMsgTxt("RATS! The number of input arguments must be 2.");
	}
	if (nlhs!=1){
		mexErrMsgTxt("RATS! The number of output arguments must be 1.");
	}

	size = mxGetN(prhs[0]);
	A = mxGetPr(prhs[0]);	
	B = mxGetPr(prhs[1]);
	plhs[0] = mxCreateDoubleMatrix(1, size, mxREAL);
	C = mxGetPr(plhs[0]);

	dataSize = sizeof(float)*size;
	h_A = (float *)malloc(dataSize);
	h_B = (float *)malloc(dataSize);
	h_C = (float *)malloc(dataSize);

	for (i=0; i<size; i++){
		h_A[i] = (float)A[i];
		h_B[i] = (float)B[i];
	} 

	addVec(h_A, h_B, h_C, size);

	for (i=0; i<size; i++){
		C[i] = (double)h_C[i];
	}

	return;
}

int main(){

	int i;
	int A[10], B[10], C[10];

	for (i=0; i<10; i++){
		A[i] = i;
		B[i] = i+3;
		C[i] = 0;
	}

	addVec(A, B, C, 10);

	for (i=0; i<10; i++){
		printf("%d %d %d \n", A[i], B[i], C[i]);
	}

	return 1;

}
