% example_nufft1_antenna.m
%
% This m-file is an example of applying the NUFFT method to compute
% the far-field beam pattern of a 2D array of nonuniformly spaced
% (point?) antenna elements.
% In the nomenclature of the 1999 SIAM J Sci. Comput. paper by Nguyen and Liu,
% this is "Problem 1."
%
% This provides an example of how to compute "the Fourier transform"
% of nonuniformly spaced data.
%
% The user should think very carefully about whether it is truly meaningful
% to compute "the FT" of unequally spaced data first.  It certainly is
% meaningful for attenna beam patterns, but may not always be what is really
% useful for all applications with nonuniformly sampled data.
%
% Copyright 2007, Jeff Fessler, University of Michigan

%
% antenna element locations in 2d plane
%
clf, pl = @(p) subplot(220+p);
for choice=1:2
	tmp = [0:40]'/41 * 2 * pi;
	yc = sin(choice*tmp); % funny bowtie pattern for illustration
	xc = cos(tmp); clear tmp
	pl(choice), plot(xc, yc, 'o'), axis square
	xlabel 'x', ylabel 'x', title 'antenna locations'

	% create NUFFT structure
	N = [1 1]*2^8;
	J = [5 5];	% interpolation neighborhood
	K = N*2;	% two-times oversampling
	om = [xc yc];	% 'frequencies' are locations here!

	% the following line probably needs to be changed
	% to get the proper scaling/dimensions in the pattern
	% but for now i just make it fill up [-pi/20,pi/20]
	% in hopes of getting a nice 'picture' of pattern
	om = (pi/20) * om / max(om(:));
	st = nufft_init(om, N, J, K, N/2, 'minmax:kb');

	weights = ones(size(xc)); % equal weights on each element; could change

	% call the *adjoint* NUFFT; this is what does "the FT of unequal data"
	pattern = nufft_adj(weights, st);

	pl(2+choice)
	im(pattern, 'pattern')
end
