/*
* mri_exp_mult_template.c
* Template for mri_exp_mult_* to support both double and float.
* Sangwoo Lee, The University of Michigan, Jun, 2004
* 2004-6-21 modified by JF
*/

/*
in:
	A	[N,L]	complex matrix
	u	[N,1]	vector
	v	[M,1]	vector
			one (and only one) of u and v must be complex!
out:
	D	[L,M]	complex vector, D = A' * exp(-u * v.')
			D_lm = sum_n A_nl^* exp(-u_n v_m)
*/

/*
* static bool mri_exp_mult_*
*/
(
dtype *D_r,		/* [L,M] */
dtype *D_i,
double *work_r,		/* [N,1] work space */
double *work_i,
const dtype *A_r,	/* [N,L] */
const dtype *A_i,
const dtype *un_r,	/* [N,1] */
const dtype *un_i,
const dtype *vm_r,	/* [M,1] */
const dtype *vm_i,
const int LL,
const int MM,
const int NN)
{
	int ll, mm, nn;

	for (mm=0; mm < MM; ++mm) {
		const dtype *Ar = A_r;
		const dtype *Ai = A_i;

		/* real and imag part of exp(-u(:) v(i)) */
		if (!vm_i) { /* vm is real */
			const double nvmr = -vm_r[mm];
			for (nn=0; nn < NN; ++nn) {
				const double ee = exp(nvmr * un_r[nn]);
				work_r[nn] = ee * cos(nvmr * un_i[nn]);
				work_i[nn] = ee * sin(nvmr * un_i[nn]);
			}
		}

		else if (!un_i) { /* un is real */
			const double nvmr = -vm_r[mm];
			const double nvmi = -vm_i[mm];
			for (nn=0; nn < NN; ++nn) {
				const double ee = exp(nvmr * un_r[nn]);
				work_r[nn] = ee * cos(nvmi * un_r[nn]);
				work_i[nn] = ee * sin(nvmi * un_r[nn]);
			}
		}

		else
			Fail("bug")

		for (ll=0; ll < LL; ++ll) {
			const double *wr = work_r;
			const double *wi = work_i;
			double sumr = 0., sumi = 0.;

			for (nn=0; nn < NN; ++nn, ++Ar, ++Ai, ++wr, ++wi) {
				sumr += *Ar * *wr + *Ai * *wi;
				sumi += *Ar * *wi - *Ai * *wr;
			}

			*D_r++ = sumr;
			*D_i++ = sumi;
		}
	}

	return 1;
}
