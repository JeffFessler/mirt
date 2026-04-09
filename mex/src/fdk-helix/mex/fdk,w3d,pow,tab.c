// fdk,w3d,pow,tab.c
// 2014-11-06, Donghwan Kim, University of Michigan
//
//
#include "defs-env.h"
#include "def,fdk.h"

sof fdk_w3d_pow_tab_init(
float *pow_tab,
cint nt,
cint nz,
cfloat dt,
cfloat dz,
cfloat dsd,
cint w3d,
cfloat pitch)
{
	// todo: offset t
	if (pitch)
		for (int it = 0; it < Ceilf(nt/2); it++, pow_tab++)
			*pow_tab = pow(it * dt / dsd, w3d * pitch); 
	else
	{
		for (int iz = 0; iz < Ceilf(nz/2); iz++)
		for (int it = 0; it < Ceilf(nt/2); it++, pow_tab++)
			*pow_tab = pow(it * dt / dsd, w3d * iz * dz);
	}

	Ok
}

sof fdk_w3d_pow_tab(
float *ww,
cfloat *pow_tab,
cint nt,
cfloat t_bin,
cfloat t_bin_c,
cfloat wt,
cint iz,
cfloat wz,
cfloat pitch)
{
	// todo: offset t, z
	cint id = pitch ? 0 : Abs(iz - wz);
	cint indx = Floorf(Min(Abs(t_bin - wt), nt/2)) + Ceilf(nt/2) * id;
	cint indx_c = Floorf(Min(Abs(t_bin_c - wt), nt/2)) + Ceilf(nt/2) * id;
	cfloat scale = pow_tab[indx];	
	cfloat scale_c = pow_tab[indx_c];

	cfloat eps = 0.001;
	*ww = 2 * Min(Max(scale_c / (scale_c + scale), eps), 1-eps);
	
	Ok
}


