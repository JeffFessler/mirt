// fdk-ts-t.cu
// Threaded versions of FDK back-projection
// For detector index (t,s).
// Copyright 2008-10-09, Jeff Fessler, University of Michigan

// Multi-GPU VERSION 1
//Same as others 

#include "jf-cuda.h"
#include "def,fdk.h"
#include "fdk-gpu.h"
#include <pthread.h>
#include <stdio.h>
#include "mex.h"

typedef struct {
	float *image; // [nz nx ny] <- trick!
	const cbct_ig *ig; // image geometry
	const cbct_cg *cg; // cone-beam CT system geometry
	int na; // # of views
	float *proj; // [nt ns na] <- trick! projection views
	cdouble *beta; // [na] source angles [radians]
} fdk_ts_s;

#ifdef fdk_gpu
#define NUM_THREADS 4
///Mew: for multi-gpu, define the basic data structure for each thread
struct threadData{	
	float	*proj;
	cdouble	*beta;
	truf 	is_arc; 
	int 	nx; 
	int 	ny;
	int 	nz; 
	float 	dx; 
	float	dy; 
	float 	dz;
	float 	wx; 
	float 	wy;
	float	wz; 
	float 	dso; 
	float 	dsd; 
	int 	ns;
	int 	nt;
	float 	ds;
	float	dt;
	float 	ws; 
	float	wt;
	int 	nxy;
	int	na;
	float 	*image;
};

struct threadData thread_data;
///Using Condition Variables
pthread_mutex_t count_mutex;
pthread_cond_t count_threshold_cv[3];
int count = 100;

	//Mew: use testure memory
texture< float, 2, cudaReadModeElementType> texRef;

static int iDivUp(int a, int b) {
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

//////////////////////////////////////////////////////////////////////////////////////
///         kernel_1
//////////////////////////////////////////////////////////////////////////////////////
static __global__ void fdk_ts_back1_kernel( 
	float *s_val,  
	truf is_arc,  
	float sinb,  
	float cosb,  
	int nx,  
	int ny,  
	float dx,  
	float dy,  
	float wx,  
	float wy,  
	float dso,  
	float dsd,  
	float ds,  
	float ws,  
	int nxy)
{
	// index into image array
	// determine the index of x, y
	cint ix = blockIdx.x * blockDim.x + threadIdx.x;	
	cint iy = blockIdx.y * blockDim.y + threadIdx.y;

	// if index is out of bound

	if (ix >= nx || iy >= ny )
		return;
	cfloat yy = dy * iy - wy;
	cfloat xx = dx * ix - wx;
	cfloat ybetas = dso - (-xx * sinb + yy * cosb);
	cfloat xbeta = xx * cosb + yy * sinb;
	cfloat mag = dsd / ybetas;
	float ss_bin;
	float w2;
	if ( is_arc ){
	ss_bin = dsd * atan2f(xbeta, ybetas) / ds + ws;
	w2 = Sqr(dsd) / (Sqr(ybetas) + Sqr(xbeta));
	}
	else{
	ss_bin =  mag * xbeta / ds + ws;
	w2 = mag* mag ;
	}
	s_val[ix + iy*nx] = (float)ss_bin;
	s_val[ix + iy*nx + nxy] = (float)w2;
	s_val[ix + iy*nx + 2*nxy] = (float)mag;
}

//////////////////////////////////////////////////////////////////////////////////////
///         kernel_2
//////////////////////////////////////////////////////////////////////////////////////
__global__ void fdk_ts_back1_kernel_2( 
	float *s_val,  
	float *image,  
	int nx,  
	int nz,  
	float dz, 
	float dt,  
	int nt,  
	int ns,  
	int nxy,  
	float wz,  
	float wt)
{

	cint ix = blockIdx.x % nx;
	cint iy = blockIdx.y; 
	cint iz = threadIdx.x+512*(nz/512);	

	//NEW 5: use shared memory
	__shared__ float ss_bin;
	__shared__ float w2;
	__shared__ float mag;

	if (iz==1){
	ss_bin 	= s_val[ix + iy*nx];
	w2 	= s_val[ix + iy*nx + nxy];
	mag 	= s_val[ix + iy*nx + 2*nxy];
	}
	__syncthreads();
	
	//New: use texture memory for bilinear interpolation
	float zz = dz * iz - wz;
	float tt_bin = mag * zz / dt + wt;

	if ( tt_bin < 0 || tt_bin >= nt || ss_bin < 0 || ss_bin >= ns )
		return;
 	image[(ix + iy * nx) * nz + iz] += w2 * tex2D(texRef, tt_bin+1 , ss_bin+1 );
}

//////////////////////////////////////////////////////////////////////////////////////
///         kernel_3
//////////////////////////////////////////////////////////////////////////////////////
__global__ void fdk_ts_back1_kernel_3(int nx, int nz, float *temp_img, float *image)
{
	cint ix = blockIdx.x % nx;
	cint iy = blockIdx.y; 
	cint iz = threadIdx.x+512*(nz/512);
	image[(ix + iy * nx) * nz + iz] += temp_img[(ix + iy * nx) * nz + iz];
}
	
//////////////////////////////////////////////////////////////////////////////////////
///         thread_function
//////////////////////////////////////////////////////////////////////////////////////
void *fdk_ts_back_thread(void *arg)
{	
	//using taskid to set gpu device later
   	long taskid;
	//load data for each device
   	float *thread_proj;
	cdouble *thread_beta;
		
   	taskid = (long)arg;
	int nth = (int) thread_data.na / NUM_THREADS; 

	int nxyz = thread_data.nx * thread_data.ny * thread_data.nz;
	int nst = thread_data.ns * thread_data.nt;
	// image memory on device
	thread_proj = thread_data.proj + nst * nth * taskid;
	thread_beta = thread_data.beta + nth * taskid;

	////select which gpu device to use

	int dCnt = 0;
	int selectedCudaDeviceId = 0;
	cudaGetDeviceCount(&dCnt) ;
	if (dCnt > 0) {
		if (dCnt > 1) {
			int deviceId = taskid % dCnt; 			
			cudaDeviceProp prop;
			if (cudaSuccess == cudaGetDeviceProperties(&prop, deviceId)) {				selectedCudaDeviceId = deviceId;
			}
		}else {
			selectedCudaDeviceId = 0;
		}
		printf("Thread: %d  selected device: %d\n",taskid, selectedCudaDeviceId);
		cudaSetDevice(selectedCudaDeviceId);
	}
	else
		return(NULL);

	///// Load all this stuff into graphics memory	


	
	float *dev_img ;
	jf_gpu_malloc(dev_img, nxyz) 
	// initialize device image to 0
	jf_gpu_memset(dev_img, 0, nxyz) 
	
	//create a space to save image from host as a temp
	float *temp_img;
	jf_gpu_malloc(temp_img, nxyz) 
	jf_gpu_memset(temp_img, 0, nxyz)
	

	//use texture memory
	cudaArray *dev_proj;
	cudaChannelFormatDesc input_tex = cudaCreateChannelDesc<float>();
	cudaMallocArray(&dev_proj, &input_tex, thread_data.nt, thread_data.ns);

	texRef.addressMode[0]	= cudaAddressModeClamp;
	texRef.addressMode[1]	= cudaAddressModeClamp;
	texRef.filterMode	= cudaFilterModeLinear;
	texRef.normalized	= 0;
	
	// s values for each x,y pair on device
	float *dev_sval;			
	jf_gpu_malloc(dev_sval, thread_data.nxy*3)
	// initialize values to 0
	jf_gpu_memset(dev_sval, 0, thread_data.nxy*3)	


	////decide the block and grid sturcture

	dim3 dimBlock(16, 16);  
 	dim3 dimGrid(iDivUp(thread_data.nx,dimBlock.x), iDivUp(thread_data.ny,dimBlock.y));
	dim3 dimBlock2(thread_data.nz, 1, 1);		
	int numBlock_x = ((thread_data.nz/512)+1)*thread_data.nx;
   	dim3 dimGrid2(numBlock_x, thread_data.ny);
	
	int num_poj = nth;
	if (taskid == NUM_THREADS )
		num_poj = thread_data.na - nth * (NUM_THREADS-1);
	
	for (int ia=0; ia <  num_poj ; ia=ia+1, thread_proj += nst) { 
	// copy this view to gpu
	cudaMemcpyToArray(dev_proj, 0, 0, thread_proj, nst*sizeof(float),cudaMemcpyHostToDevice);
	cudaBindTextureToArray(texRef, dev_proj, input_tex);

	fdk_ts_back1_kernel<<<dimGrid, dimBlock>>>(
		dev_sval, 
		thread_data.is_arc, 
		sinf(thread_beta[ia]), 
		cosf(thread_beta[ia]),  
		thread_data.nx,  
		thread_data.ny,  
		thread_data.dx,  
		thread_data.dy,  
		thread_data.wx,  
		thread_data.wy,  
		thread_data.dso,  
		thread_data.dsd,  
		thread_data.ds,  
		thread_data.ws,  
		thread_data.nxy);

	fdk_ts_back1_kernel_2<<<dimGrid2, dimBlock2>>>(	 
		dev_sval,  
		dev_img,   
		thread_data.nx,   
		thread_data.nz,   
		thread_data.dz,   
		thread_data.dt,   
		thread_data.nt,   
		thread_data.ns,   
		thread_data.nxy,   
		thread_data.wz,   
		thread_data.wt);
	}

	///combine the images from each thread
	
	printf("thread %d finish calc \n", taskid);
	//stop other thread the wait
	pthread_mutex_lock(&count_mutex);
	if (taskid < 3 && count != taskid ){
		printf("thread %d is waiting \n", taskid);
		pthread_cond_wait(&count_threshold_cv[taskid], &count_mutex);
	}
	printf("here is %d\n",taskid );
	//copy the image from host to device
	jf_gpu_put(temp_img, thread_data.image, nxyz)
	//add together
	fdk_ts_back1_kernel_3<<<dimGrid2, dimBlock2>>>(thread_data.nx, thread_data.nz, temp_img, dev_img);
	//copy the image from device back to host
	jf_gpu_get(thread_data.image, dev_img, nxyz)
	printf("here is still %d\n",taskid );

	//allow other thread to run
	if(taskid > 0){
		printf("thread %d calls thread %d \n" ,taskid , taskid-1);
		pthread_cond_signal(&count_threshold_cv[taskid-1]);
		count = taskid-1;
	}
	pthread_mutex_unlock(&count_mutex);

	
	//Note("freeing memory\n")
	jf_gpu_free(dev_img)
	jf_gpu_free(temp_img)
	cudaFree(dev_proj);
	jf_gpu_free(dev_sval)
   	return(NULL);
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
	int nxy = ig->nx * ig->ny;


	//int nxyz = ig->nx * ig->ny * ig->nz;
	//int nst = cg->ns * cg->nt;
	//calculate the offset
	//note: wxyz are different with wst
	cfloat wx = ig->dx *( (ig->nx-1)/2. + ig->offset_x );
	cfloat wy = ig->dy *( (ig->ny-1)/2. + ig->offset_y );
	cfloat wz = ig->dz *( (ig->nz-1)/2. + ig->offset_z );
	cfloat ws = (cg->ns-1)/2. + cg->offset_s;
	cfloat wt = (cg->nt-1)/2. + cg->offset_t;

	//decide the shape of detector
	truf is_arc = 0;
	if (cg->dfs == 0)
		is_arc = 1;
	else if (!Isinf(cg->dfs))
		Warn("dfs not done - junk!")

	////prepare for multigpus
#ifdef fdk_gpu
	
	///GPU count
	int dCnt = 0;
	cudaGetDeviceCount(&dCnt) ;
	if (dCnt > 0)		printf("number of cuda gpu devices: %d\n", dCnt);
	else{
		printf("Hey dude, install a cuda device first. \n");
		Ok		
	}
	
	  /*
	  Lock mutex and wait for signal.  Note that the pthread_cond_wait 
	  routine will automatically and atomically unlock mutex while it waits. 
	  Also, note that if COUNT_LIMIT is reached before this routine is run by
	  the waiting thread, the loop will be skipped to prevent pthread_cond_wait
	  from never returning. 
	  */

	
	pthread_attr_t attr;
	pthread_attr_init(&attr);
 	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	pthread_mutex_init(&count_mutex, NULL);
	for(int init_cond = 0;init_cond <3 ;init_cond++ ){
		pthread_cond_init (&count_threshold_cv[init_cond], NULL);
	}

	//using a thread to control a gpu
	pthread_t threads[NUM_THREADS];
	
	int rc, t;
		
	//load to data to each thread
	thread_data.proj 	= proj;
	thread_data.beta	= beta;
	thread_data.is_arc	= is_arc; ; 
	thread_data.nx		= ig->nx; 
	thread_data.ny		= ig->ny;
	thread_data.nz		= ig->nz; 
	thread_data.dx		= ig->dx; 
	thread_data.dy		= ig->dy; 
	thread_data.dz		= ig->dz;
	thread_data.wx		= wx; 
	thread_data.wy		= wy; 
	thread_data.wz		= wz; 
	thread_data.dso		= cg->dso; 
	thread_data.dsd		= cg->dsd; 
	thread_data.ns		= cg->ns;
	thread_data.nt		= cg->nt;
	thread_data.ds		= cg->ds;
	thread_data.dt		= cg->dt;
	thread_data.ws		= ws; 
	thread_data.wt		= wt;
	thread_data.nxy		= nxy;
	thread_data.na		= na;
	thread_data.image	= pa->image;

	
	for(t=0;t<NUM_THREADS;t++) {
		//printf("Creating thread %d\n", t);
		rc = pthread_create(&threads[t],  &attr, fdk_ts_back_thread, (void *)t);
	  	if (rc) {
	    		printf("ERROR; return code from pthread_create() is %d\n", rc);
	    		exit(-1);
	    	}
	  }
	
	pthread_attr_destroy(&attr);

	///wait for all the thread finished	
	for (int i=0;i<NUM_THREADS;i++) {
		pthread_join(threads[i],NULL);
	}

	//destory the pthread arguments
	pthread_attr_destroy(&attr);
	pthread_mutex_destroy(&count_mutex);
	for(int dstr_cond = 0;dstr_cond <3 ;dstr_cond++ ){
		pthread_cond_destroy(&count_threshold_cv[dstr_cond]);
	}

#else 

	for (int ia=0; ia < na; ++ia, proj += nst) {

	float *dev_img = pa->image; // already zeroed
	cfloat *dev_proj = proj;
	//cbyte *dev_mask2 = ig->mask2;
	if (!fdk_ts_back(dev_img,	
			ig->nx, ig->ny, ig->nz,
			ig->dx, ig->dy, ig->dz,
			wx, wy , wz,
			//dev_mask2, id + 1, // each thread does some voxels only
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
	fdk_ts_back_init(&st, nthread, 0);
        Ok
}


