  function dft_sym_check(xk, varargin)
%|function dft_sym_check(xk)
%| see if array xk (usually containing DFT coefficients)
%| obeys the conjugate symmetry conditions expected
%| for its (inverse) DFT to be real.

if nargin < 1, ir_usage, end
if nargin == 1 && streq(xk, 'test'), dft_sym_check_test, return, end

arg.tol = 1e-6;

if ndims(xk) == 2
	dft_sym_check2(xk, arg.tol)
else
	error 'only 2d done'
end


function dft_sym_check2(xk, tol)

[nx ny] = size(xk);

iix = 0:nx-1;
iiy = 0:ny-1;
[ix iy] = ndgrid(iix, iiy);

i1 = 1 + ix + iy * nx;
i2 = 1 + mod(nx-ix,nx) + mod(ny-iy,ny) * nx;
err = abs(xk(i1) - conj(xk(i2)));

err = err / max(max(abs(xk(:))), eps);
bad = find(err(:) > tol);
if length(bad)
	printm('#bad = %d, worst = %g%%', length(bad), max(err(:))*100)
	nn = min(5, length(bad));
	for ii=1:nn
		kk = bad(ii);
		printm('ix,iy=(%d,%d) err=%g%%', ix(kk), iy(kk), err(kk)*100)
	end
	im(iix, iiy, err), cbar
elseif any(err(:))
	printm('ok, worst = %g%%', max(err(:))*100)
end


function dft_sym_check_test
rng(0)
xk = fft2(rand(14,23));
dft_sym_check(xk)
xk(35) = 1.1 * xk(5);
dft_sym_check(xk)
