// jf-cuda.h
// convenience macro's for calling cuda
// Copyright 2009-07-03, Jeff Fessler, University of Michigan

#ifndef jf_cuda_h
#define jf_cuda_h

#include "defs-env.h"
#include <cuda.h>

#define jf_gpu_malloc(p, n) \
	{ \
	cudaError_t jf_cuda_error = \
	cudaMalloc( (void **) &p, n * sizeof(*p) ); \
	if (jf_cuda_error != cudaSuccess ) \
		Fail1("cudaMalloc %s", cudaGetErrorString(jf_cuda_error)); \
	}

#define jf_gpu_free(p) \
	{ \
	cudaError_t jf_cuda_error = cudaFree( (void *) p ); \
	if (jf_cuda_error != cudaSuccess ) \
		Fail1("cudaFree %s", cudaGetErrorString(jf_cuda_error)); \
	}

#define jf_gpu_memset(p, v, n) \
	{ \
	cudaError_t jf_cuda_error = cudaMemset(p, v, n * sizeof(*p)); \
	if (jf_cuda_error != cudaSuccess ) \
		Fail1("cudaMemset %s", cudaGetErrorString(jf_cuda_error)); \
	}

#define jf_gpu_put(to_gpu, from_cpu, n) \
	{ \
	cudaError_t jf_cuda_error = \
        cudaMemcpy((void *) to_gpu, (void *) from_cpu, n * sizeof(*from_cpu), \
		cudaMemcpyHostToDevice); \
	if (jf_cuda_error != cudaSuccess ) \
		Fail1("cudaMemcpy %s", cudaGetErrorString(jf_cuda_error)); \
	}
//	Note1("put %d bytes", n * sizeof(*from_cpu))

#define jf_gpu_get(to_cpu, from_gpu, n) \
	{ \
	cudaError_t jf_cuda_error = \
        cudaMemcpy((void *) to_cpu, (void *) from_gpu, n * sizeof(*to_cpu), \
		cudaMemcpyDeviceToHost); \
	if (jf_cuda_error != cudaSuccess ) \
		Fail1("cudaMemcpy %s", cudaGetErrorString(jf_cuda_error)); \
	}
//	Note1("get %d bytes", n * sizeof(*to_cpu))

#endif // jf_cuda_h
