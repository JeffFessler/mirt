NUFFT - nonuniform FFT toolbox for Matlab
by Jeff Fessler, Brad Sutton, and Yingying Zhang

This directory contains a Matlab toolbox for computing nonuniform
fast Fourier transforms (NUFFT) for unequally spaced samples.
The interpolation methods have been optimized by a min-max criterion
that minimizes the worst-case error.

See "Nonuniform fast Fourier transforms using min-max interpolation"
by Jeffrey A. Fessler and Bradley P. Sutton,
IEEE Trans. Signal Processing, 51(2):560-74, Feb. 2003

This work was inspired by collaboration on MR imaging with Doug Noll.
For more details, look for MR-related papers on Fessler's web site:
	http://www.eecs.umich.edu/~fessler/

For a general NUFFT of any dimension (1,2,3,etc.), use the routine
nufft.m, which requires a structure that is initialized by nufft_init.m
as illustrated in the built-in test routine inside nufft.m itself.
Just run "nufft" without any arguments at the Matlab prompt and it
will run the test example and should print something like:
	minmax:kb       max%diff = 0.00... (the exact value depends on version)
It is quite a small error which is why these routines are useful!
For its adjoint, use "nufft_adj.m".

I am working on a replacement called newfft.m that uses real interpolation,
whereas the original code used complex interpolation unnecessarily.

For use in iterative algorithms, try the "Gnufft" object under "../systems"
as illustrated in ../example/mri_example.m, which is my recommended method
for reconstruction from nonuniform frequency-space samples.

Prior to developing nufft.m, the former principal routines were called
nufft1.m and nufft2.m, which I have retained for backward compatibility
in the "archive" subdirectory.

The routines in the "private" subdirectory are used by the routines
in this directory but are probably not useful otherwise.
Hiding them that way is a little-known but useful Matlab convention!

The routines in the "tsp2003figs" subdirectory are some of those used for the
figures in the 2003 T-SP paper.  If the one you want is not there, just ask.

The very first version of these files was written by Brad Sutton.
But they have evolved considerably since then.

A limitation of the original NUFFT code was that it precomputes all of
the J*M interpolation coefficients and stores those in a large sparse
matrix.  For large-scale problems this precomputing can use too much
memory, even when stored as a sparse matrix.  There is now a version
where the interpolator is tabulated only on a equispaced fine grid;
this 'table' method uses much less memory, at the price of somewhat
reduced accuracy.  See nufft_table_test.m for illustrations.

Before using this, you might want to run "fftn_fast test" in matlab to see if
fftn_fast.m is configured optimally for your machine.  See comments within it.

For completeness, I mention that Daniel Potts has a NUFFT C library:
	http://www.math.uni-luebeck.de/potts
