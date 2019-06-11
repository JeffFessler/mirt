 function X = dtft1(x, omega, n_shift)
%function X = dtft1(x, omega, n_shift)
%| Compute DTFT of 1D signal x at frequency locations wx
%| In
%|	x	[N,ncol]	signal values
%|	omega	[M,1]		frequency locations
%|	n_shift	[1,1]		use [0:(N-1)]-n_shift indices
%| Out
%|	X	[M,ncol]	DTFT values
%|
%| This needs enough memory to store M * N size matrices
%| Copyright 2000-1-9, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), dtft1_test, return, end
if nargin < 2, ir_usage(), end

N = size(x,1);
if ~isvar('n_shift') || isempty(n_shift), n_shift = 0; end
nn = [0:(N-1)] - n_shift;

X = exp(-1i*omega*nn) * x; % compute 1D DTFT


function dtft1_test()
N = 4;
x = [[1:N]', [1 1 2 2]'];	% two test signal
omega = 2*pi*[0:(N-1)]'/N;	% test with uniform frequency locations
Xd = dtft1(x, omega);
Xf = fft(x);
%disp([Xd Xf])
printm('max %% difference = %g', max_percent_diff(Xf,Xd))
