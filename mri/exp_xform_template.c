/*
* exp_xform_template.c
* Template for exp_xform_*, to support both single and double precision
* Copyright 2005-7-28, Jeff Fessler, The University of Michigan
*/

/*
* y(m,l) = sum_n x(n,l) exp(-sum(u(:,n) .* v(:,m)))
*/

(
dtype *y_r,		/* [M,1] */
dtype *y_i,
const dtype *x_r,	/* [N,1] */
const dtype *x_i,
const dtype *un_r,	/* [D,N] */
const dtype *un_i,
const dtype *vm_r,	/* [D,M] */
const dtype *vm_i,
const int DD,		/* # of dimensions, e.g., 2 spatial, 1 temporal */
const int MM,
const int NN)
{
	int mm, nn;

	for (mm=0; mm < MM; ++mm) {
		const dtype *xr = x_r;
		const dtype *xi = x_i;
		const dtype *vr_m = vm_r+mm*DD;
		const dtype *vi_m = vm_i+mm*DD;
		double sumr = 0., sumi = 0.;

		for (nn=0; nn < NN; ++nn, ++xr, ++xi) {
			double ee, er, ei;
			const dtype *ur = un_r+nn*DD;
			const dtype *ui = un_i+nn*DD;

			/* dot product of v_m with u_n */
			const dtype *vr = vr_m;
			const dtype *vi = vi_m;
			double uv_r = 0.;
			double uv_i = 0.;
			int dd;
			for (dd=0; dd < DD; ++dd, ++ur, ++ui, ++vr, ++vi) {
				uv_r += *ur * *vr - *ui * *vi;
				uv_i += *ur * *vi + *ui * *vr;
			}

			/* real and imag part of exp(-u.v) */
			ee = exp(-uv_r);
			er = ee * cos(uv_i);
			ei = ee * sin(-uv_i);

			/* sum += xn * exp(-un.v) */
			sumr += *xr * er - *xi * ei;
			sumi += *xr * ei + *xi * er;
		}

		*y_r++ = sumr;
		*y_i++ = sumi;
	}
}
