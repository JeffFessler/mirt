// fdk-ts-h.cu
#include <stdio.h>

void fdk_ts_help(void)
{
	printf("\n\
\n\
	image = function('fdk,ts,back', nx,ny,nz, dx,dy,dz, \n\
		offset_x, offset_y, offset_z, mask2, \n\
		dso, dsd, ds, dt, offset_s, offset_t, proj, beta, nthread)\n\
\n\
		image output is single [nz nx ny] <- trick!\n\
		nx,ny,nz: (int32) image size\n\
		dx,dy,dz: (double) voxel size\n\
		offset_x,_y,_z: (double) center offset in pixels (usually 0)\n\
		mask2: (uint8) [nx ny] 2D support mask\n\
		dso: (double) distance from source to isocenter\n\
		dsd: (double) distance from source to detector\n\
		dfs: (double) distance from focal point to source (0 or inf)\n\
		ds: (double) horizontal ray spacing\n\
		dt: (double) vertical ray spacing\n\
		offset_s: (double) channel offset [pixels]\n\
		offset_t: (double) vertical offset on detector [pixels]\n\
		nthread: (int32) # of processors\n\
		proj: (single) [nt ns na] (trick!) projection view for each beta\n\
		beta: (double) [na] source angle(s) [radians]\n\
	(CUDA version)\n\
\n");
}
