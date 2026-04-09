// fdk-ts-t.cu
// Threaded versions of FDK back-projection
// For detector index (t,s).
// Copyright 2008-10-09, Jeff Fessler, University of Michigan

#include "jf-cuda.h"
#include "def,fdk.h"
#include "fdk-gpu.h"
#include <pthread.h>
#include <stdio.h>
#include "mex.h"

#define TEXMEMO 1
#define MULT_GPU 0
#define USE_MASK 0

#if MULT_GPU
#define NUM_THREADS 4
void *PrintHello(void *threadid){
	long tid;
	tid = (long)threadid;
	printf("Hello World! It's me, thread #%d! \n", tid);
	return NULL;
}
#endif


#ifdef fdk_gpu
	//Mew: use testure memory
	texture< float, 2, cudaReadModeElementType> texRef;
	//Mew: use constant memory
	__constant__ float dc_wx;
	__constant__ float dc_wy;
	__constant__ float dc_wz;
	__constant__ float dc_ws;
	__constant__ float dc_wt;
	__constant__ int dc_nx;
	__constant__ int dc_ny;
	__constant__ int dc_nz;
	__constant__ int dc_ns;
	__constant__ int dc_nt;
	__constant__ float dc_dx;
	__constant__ float dc_dy;
	__constant__ float dc_dz;
	__constant__ float dc_ds;
	__constant__ float dc_dt;
	__constant__ float dc_dso;
	__constant__ float dc_dsd;
	__constant__ float dc_dfs;
	__constant__ int dc_nxy;

//////////////////////////////////////////////////////////////////////////////////////
///         kernel_1
//////////////////////////////////////////////////////////////////////////////////////
#if USE_MASK
static __global__ void fdk_ts_back1_kernel( float *s_val, truf is_arc,float sinb,float cosb,cbyte *mask2) //source angle [radians]
#else 
static __global__ void fdk_ts_back1_kernel( float *s_val, truf is_arc,float sinb,float cosb)
#endif
{
	// index into image array
	// determine the index of x, y
	cint ix = blockIdx.x * blockDim.x + threadIdx.x;	
	cint iy = blockIdx.y * blockDim.y + threadIdx.y;

	// if index is out of bound
#if USE_MASK
	if (ix >= dc_nx || iy >= dc_ny || !mask2[ix + iy*dc_nx])
#else	
	if (ix >= dc_nx || iy >= dc_ny )
#endif
		return;
	cfloat yy = dc_dy * iy - dc_wy;
	cfloat xx = dc_dx * ix - dc_wx;
	cfloat ybetas = dc_dso - (-xx * sinb + yy * cosb);
	cfloat xbeta = xx * cosb + yy * sinb;
	cfloat mag = dc_dsd / ybetas;
	float ss_bin;
	float w2;
	if ( is_arc ){
	ss_bin = dc_dsd * atan2f(xbeta, ybetas) / dc_ds + dc_ws;
	w2 = Sqr(dc_dsd) / (Sqr(ybetas) + Sqr(xbeta));
	}
	else{
	ss_bin =  mag * xbeta / dc_ds + dc_ws;
	w2 = mag* mag ;
	}

	s_val[ix + iy*dc_nx] = (float)ss_bin;
	s_val[ix + iy*dc_nx + dc_nxy] = (float)w2;
	s_val[ix + iy*dc_nx + 2*dc_nxy] = (float)mag;

}

//////////////////////////////////////////////////////////////////////////////////////
///         kernel_2
//////////////////////////////////////////////////////////////////////////////////////
#if TEXMEMO
#if USE_MASK	
__global__ void fdk_ts_back1_kernel_2(float *s_val, float *image,cbyte *mask2)
#else
__global__ void fdk_ts_back1_kernel_2(float *s_val, float *image)
#endif
{
	cint ix = blockIdx.x % dc_nx;
	cint iy = blockIdx.y; 
	cint iz = threadIdx.x+512*(dc_nz/512);	

	//NEW 5: use shared memory
	__shared__ float ss_bin;
	__shared__ float w2;
	__shared__ float mag;
	if (iz==1){
	ss_bin = s_val[ix + iy*dc_nx];
	w2 = s_val[ix + iy*dc_nx + dc_nxy];
	mag = s_val[ix + iy*dc_nx + 2*dc_nxy];
	}
	__syncthreads();
	
	//New: use texture memory for bilinear interpolation
	float zz = dc_dz * iz - dc_wz;
	float tt_bin = mag * zz / dc_dt + dc_wt;
#if USE_MASK
	if (tt_bin < 0 || tt_bin > dc_nt || ss_bin < 0 || ss_bin > dc_ns|| ix>=dc_nx || iy>=dc_ny || (mask2[ix + iy*dc_nx] != 1))
#else 	
	if ( ix>=dc_nx || iy>=dc_ny )
#endif
		return;
 	image[(ix + iy * dc_nx) * dc_nz + iz] += w2 * tex2D(texRef, tt_bin , ss_bin );
}	
	
#else
#if USE_MASK	
__global__ void fdk_ts_back1_kernel_2(float *s_val, float *image, cfloat *proj,cbyte *mask2)
#else
__global__ void fdk_ts_back1_kernel_2(float *s_val, float *image, cfloat *proj)
#endif
{

	// NEW 5:
	// index into image array
	cint ix = blockIdx.x % dc_nx;
	cint iy = blockIdx.y; 
	// if more than 512, iz=threadIdx.x+512
	// else iz=threadIdx.x 
	cint iz = threadIdx.x+512*(dc_nz/512);	

	//NEW 5: use shared memory
	__shared__ float ss_bin;
	__shared__ float w2;
	__shared__ float mag;
	if (iz==1){
	ss_bin = s_val[ix + iy*dc_nx];
	w2 = s_val[ix + iy*dc_nx + dc_nxy];
	mag = s_val[ix + iy*dc_nx + 2*dc_nxy];
	}
	__syncthreads();
	// index of s is "is"
	// index of nearest neighbor in "s"
	cint is = floorf(ss_bin); 

	// Check if index out of bound
#if USE_MASK
	if (is < 0 || is >= dc_ns-1 || ix>=dc_nx || iy>=dc_ny || !(mask2[ix + iy*dc_nx])) // each thread does its part only
#else	
	if (is < 0 || is >= dc_ns-1 || ix>=dc_nx || iy>=dc_ny )
#endif	
		return;

	// horizontal bilinear
	cfloat wr = ss_bin - is;
	// interpolation factors 
	cfloat wl = 1. - wr; 
	//find the image point
	image += (ix + iy * dc_nx) * dc_nz + iz;
	
	cfloat *pp1 = proj + is * dc_nt;
	cfloat *pp2 = proj + (is+1) * dc_nt;
	
	//vertical biliner
	cfloat zz = dc_dz * iz - dc_wz;
	cfloat tt = mag * zz;
	cfloat tt_bin = tt / dc_dt + dc_wt;
	// z value is used to determine index of t "it" nearest nbr in "t"
	cint it = floorf(tt_bin); 
	if (it < 0 || it >= dc_nt-1) // out of FOV
		return;
	// reconstructing the image
	else { 
			cfloat wu = tt_bin - it;
			cfloat wd = 1. - wu;
			cfloat p1 = wl * pp1[it]
				+ wr * pp2[it]; // interpolate
			cfloat p2 = wl * pp1[it+1]
				+ wr * pp2[it+1]; // horizontal

			// final vertical interpolation:
			*image += w2 * (wu * p1 + wd * p2);
	}
}
#endif
#endif

typedef struct {
	float *image; // [nz nx ny] <- trick!
	const cbct_ig *ig; // image geometry
	const cbct_cg *cg; // cone-beam CT system geometry
	int na; // # of views
	float *proj; // [nt ns na] <- trick! projection views
	cdouble *beta; // [na] source angles [radians]
} fdk_ts_s;

#ifdef fdk_gpu
static int iDivUp(int a, int b) {
	return (a % b != 0) ? (a / b + 1) : (a / b);
}
#endif

//////////////////////////////////////////////////////////////////////////////////////
///         fdk_ts_back_init()
//////////////////////////////////////////////////////////////////////////////////////
static sof fdk_ts_back_init(void *in, cint id, cint nthread)
{
	fdk_ts_s *pa = (fdk_ts_s *) in;
	const cbct_ig *ig = pa->ig;
	const cbct_cg *cg = pa->cg;
	cint na = pa->na;
	float *proj = pa->proj;
	cdouble *beta = pa->beta;
	//calculate the data size
	cint nst = cg->ns * cg->nt;
	cint nxy = ig->nx * ig->ny;
	cint nxyz = ig->nx * ig->ny * ig->nz;
	//printf("Ox=%d  Oy=%d   Oz=%d   Os=%d   Ot=%d \n", ig->offset_x,ig->offset_y, ig->offset_z, cg->offset_s, cg->offset_t);
	//calculate the offset
	//note: wxyz are different with wst
	cfloat wx = ig->dx *( (ig->nx-1)/2. + ig->offset_x );
	cfloat wy = ig->dy *( (ig->ny-1)/2. + ig->offset_y );
	cfloat wz = ig->dz *( (ig->nz-1)/2. + ig->offset_z );
	cfloat ws = (cg->ns-1)/2. + cg->offset_s;
	cfloat wt = (cg->nt-1)/2. + cg->offset_t;
	//printf("wx=%d  wy=%d   wz=%d   ws=%d   wt=%d \n", wx, wy, wz, ws, wt);
#ifdef fdk_gpu

	////prepare for multigpus
#if MULT_GPU

	int dCnt = 0;
	int selectedCudaDeviceId = 0;
	cudaGetDeviceCount(&dCnt) ;
	printf("number of cuda gpu devices: %d\n", dCnt);
	if (dCnt > 0) {
		if (dCnt > 1) {
			int multiprocessor_cnt = 0;
			cudaDeviceProp prop;
			for (int deviceId=0; deviceId<dCnt; ++deviceId) {
				if (cudaSuccess == cudaGetDeviceProperties(&prop, deviceId)) {
					if (prop.multiProcessorCount > multiprocessor_cnt) {
						multiprocessor_cnt = prop.multiProcessorCount;
						selectedCudaDeviceId = deviceId;
					}
				}
			}
		} else {
			selectedCudaDeviceId = 0;
		}
		printf("selected device with most multiprocessors: %d\n", selectedCudaDeviceId);
		cudaSetDevice(selectedCudaDeviceId);
	}
	pthread_t threads[NUM_THREADS];
	int rc;
	long t;
	for(t=0 ; t< NUM_THREADS ; t++){
		printf("In main: creating thread %ld\n", t);
		rc = pthread_create(&threads[t], NULL, PrintHello, (void *)t);
      		if (rc){
         		printf("ERROR; return code from pthread_create() is %d\n", rc);
         		exit(-1);
     		 }
 	}
#endif 

	///// Load all this stuff into graphics memory	

	// image memory on device
	float *dev_img;
	jf_gpu_malloc(dev_img, nxyz) 
	// initialize device image to 0
	jf_gpu_memset(dev_img, 0, nxyz) 

#if USE_MASK	
	byte *dev_mask2;
	jf_gpu_malloc(dev_mask2, nxy) // 2D mask
	jf_gpu_put(dev_mask2, ig->mask2, nxy)
#endif

#if TEXMEMO
	//printf("Using texture memory \n");
	//use texture memory
	cudaArray *dev_proj;
	cudaChannelFormatDesc input_tex = cudaCreateChannelDesc<float>();
	//cudaChannelFormatDesc input_tex = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaMallocArray(&dev_proj, &input_tex, cg->nt, cg->ns);

	texRef.addressMode[0]	= cudaAddressModeWrap;
	texRef.addressMode[1]	= cudaAddressModeWrap;
	texRef.filterMode	= cudaFilterModeLinear;
	texRef.normalized	= 0;
	
#else
	//printf("Using global memory \n");
	// one projection view on device
	float *dev_proj;
	jf_gpu_malloc(dev_proj, nst) 

#endif
	// s values for each x,y pair on device
	float *dev_sval;			
	jf_gpu_malloc(dev_sval, nxy*3)
	// initialize values to 0
	jf_gpu_memset(dev_sval, 0, nxy)	

	//load the parameters to constant memory
	cudaMemcpyToSymbol( "dc_wx", &wx, sizeof(float) );
	cudaMemcpyToSymbol( "dc_wy", &wy, sizeof(float) );
	cudaMemcpyToSymbol( "dc_wz", &wz, sizeof(float) );
	cudaMemcpyToSymbol( "dc_ws", &ws, sizeof(float) );
	cudaMemcpyToSymbol( "dc_wt", &wt, sizeof(float) );
	cudaMemcpyToSymbol( "dc_nx", &(ig->nx), sizeof(int) );
	cudaMemcpyToSymbol( "dc_ny", &(ig->ny), sizeof(int) );
	cudaMemcpyToSymbol( "dc_nz", &(ig->nz), sizeof(int) );
	cudaMemcpyToSymbol( "dc_dx", &(ig->dx), sizeof(float) );
	cudaMemcpyToSymbol( "dc_dy", &(ig->dy), sizeof(float) );
	cudaMemcpyToSymbol( "dc_dz", &(ig->dz), sizeof(float) );
	cudaMemcpyToSymbol( "dc_ns", &(cg->ns), sizeof(int) );
	cudaMemcpyToSymbol( "dc_nt", &(cg->nt), sizeof(int) );
	cudaMemcpyToSymbol( "dc_ds", &(cg->ds), sizeof(float) );
	cudaMemcpyToSymbol( "dc_dt", &(cg->dt), sizeof(float) );
	cudaMemcpyToSymbol( "dc_dso", &(cg->dso), sizeof(float) );
	cudaMemcpyToSymbol( "dc_dsd", &(cg->dsd), sizeof(float) );
	cudaMemcpyToSymbol( "dc_dfs", &(cg->dfs), sizeof(float) );
	cudaMemcpyToSymbol( "dc_nxy", &(nxy), sizeof(int) );

	////decide the block and grid sturcture
	dim3 dimBlock(16, 16);  
 	dim3 dimGrid(iDivUp(ig->nx,dimBlock.x), iDivUp(ig->ny,dimBlock.y));
	// NEW 5: case where nz <=512
	// all the z's in the same block for a given (x,y)
	dim3 dimBlock2(ig->nz, 1, 1);		
	int numBlock_x = ((ig->nz/512)+1)*ig->nx;
   	dim3 dimGrid2(numBlock_x, ig->ny);

#endif
	//decide the shape of detector
	truf is_arc = 0;
	if (cg->dfs == 0)
		is_arc = 1;
	else if (!Isinf(cg->dfs))
		Warn("dfs not done - junk!")

#ifdef fdk_gpu

	for (int ia=0; ia < na; ia=ia+1, proj += nst) { 
	// copy this view to gpu
#if USE_MASK
#if TEXMEMO
	cudaMemcpyToArray(dev_proj, 0, 0, proj, nst*sizeof(float),cudaMemcpyHostToDevice);
	cudaBindTextureToArray(texRef, dev_proj,input_tex);
	fdk_ts_back1_kernel<<<dimGrid, dimBlock>>>(dev_sval,is_arc,sinf(beta[ia]),cosf(beta[ia]),dev_mask2);
	fdk_ts_back1_kernel_2<<<dimGrid2, dimBlock2>>>(	dev_sval,dev_img,dev_mask2);
#else
	jf_gpu_put(dev_proj, proj, nst);
	fdk_ts_back1_kernel<<<dimGrid, dimBlock>>>(dev_sval,is_arc,sinf(beta[ia]),cosf(beta[ia]),dev_mask2);
	fdk_ts_back1_kernel_2<<<dimGrid2, dimBlock2>>>(	dev_sval,dev_img,dev_proj,dev_mask2);
#endif	
#else 
#if TEXMEMO
	cudaMemcpyToArray(dev_proj, 0, 0, proj, nst*sizeof(float),cudaMemcpyHostToDevice);
	cudaBindTextureToArray(texRef, dev_proj,input_tex);
	fdk_ts_back1_kernel<<<dimGrid, dimBlock>>>(dev_sval,is_arc,sinf(beta[ia]),cosf(beta[ia]));
	fdk_ts_back1_kernel_2<<<dimGrid2, dimBlock2>>>(	dev_sval,dev_img);
#else
	jf_gpu_put(dev_proj, proj, nst);
	fdk_ts_back1_kernel<<<dimGrid, dimBlock>>>(dev_sval,is_arc,sinf(beta[ia]),cosf(beta[ia]));
	fdk_ts_back1_kernel_2<<<dimGrid2, dimBlock2>>>(	dev_sval,dev_img,dev_proj);
#endif
#endif	

	}

	//Note("Copying image to host")
	jf_gpu_get(pa->image, dev_img, nxyz) // caution: works only for 1 thread

	//Note("freeing memory\n")
	jf_gpu_free(dev_img)
	cudaFree(dev_proj);	
//	jf_gpu_free(dev_proj)
	jf_gpu_free(dev_sval)
#if USE_MASK
	jf_gpu_free(dev_mask2)
#endif

#else

	for (int ia=0; ia < na; ++ia, proj += nst) {

	float *dev_img = pa->image; // already zeroed
	cfloat *dev_proj = proj;
	//cbyte *dev_mask2 = ig->mask2;
	if (!fdk_ts_back(dev_img,	
			ig->nx, ig->ny, ig->nz,
			ig->dx, ig->dy, ig->dz,
			wx, wy , wz,
//			dev_mask2, id + 1, // each thread does some voxels only
			cg->dso, cg->dsd, cg->dfs,
			cg->ns, cg->nt,
			cg->ds, cg->dt, ws, wt,
			dev_proj, beta[ia]),is_arc)
			Fail("fdk_ts_back()")
	}	
		
#endif

	Ok
}

sof fdk_ts_back_t(
float *image, // [nz nx ny] <- trick!
const cbct_ig *ig,
const cbct_cg *cg,
cint na, // # of views
cfloat *proj, // [nt ns na] <- trick! projection views
cdouble *beta, // [na] source angles [radians]
cint nthread, // # of threads
cint chat)
{
	fdk_ts_s st;
#define put(arg) st.arg = arg;
	put(image)
	put(ig)
	put(cg)
	put(na)
	put(proj)
	put(beta)
#undef put

// TESTING
// printf("MY FILE COOLNESS! \n");

	Bzero(image, ig->nx * ig->ny * ig->nz) // initialize image volume to 0

//	Call(jf_thread1_top, (fdk_ts_back_init, NULL /* wrap up */, &st, nthread, Chat))
	fdk_ts_back_init(&st, 4, 0);
        Ok
}


