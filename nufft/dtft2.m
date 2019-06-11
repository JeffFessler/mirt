 function X = dtft2(x, omega, n_shift, useloop)
%function X = dtft2(x, omega, n_shift, useloop)
%| Compute 2D DTFT of 2D signal x at frequency locations (o1, o2)
%| in
%|	x	[N1,N2,L]	signal values
%|	omega	[M,2]		frequency locations (radians)
%|	n_shift [2,1]		use [0:N-1]-n_shift (default [0 0])
%|	useloop			1 to reduce memory use (slower)
%| out
%|	X	[M,L]		2D DTFT values
%|
%| Requires enough memory to store M * (N1*N2) size matrices (for testing)
%|
%| Copyright 2001-9-17, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), dtft2_test, return, end
if nargin < 2, ir_usage(), end

%warning 'dtft2 is obsolete.  use dtft'

if ~isvar('n_shift') || isempty(n_shift), n_shift = [0 0]; end
if ~isvar('useloop') || isempty(useloop), useloop = 0; end

[N1, N2] = size(x);
x = reshape(x, N1*N2, numel(x)/N1/N2);	% [N1*N2 L]

[nn1, nn2] = ndgrid([0:(N1-1)]-n_shift(1), [0:(N2-1)]-n_shift(2));

if useloop % loop way: slower but less memory

	M = length(omega);
	X = zeros(numel(x)/N1/N2,M); % [L M]
	t1 = (-1i) * nn1(:)';
	t2 = (-1i) * nn2(:)';
	for ii=1:M
		X(:,ii) = exp(omega(ii,1)*t1 + omega(ii,2)*t2) * x;
	end
	X = X.'; % [M L]

else % non-loop way

	X = exp(-1i*(omega(:,1)*nn1(:)' + omega(:,2)*nn2(:)')) * x;
end


function dtft2_test()
N1 = 4; N2 = 6;
n_shift = [2 1];
x = [[1:N1]'*ones(1,3), ones(N1,N2-3)];	% test signal
o1 = 2*pi*[0:(N1-1)]'/N1; % test with uniform frequency locations
o2 = 2*pi*[0:(N2-1)]'/N2; % test with uniform frequency locations
[o1 o2] = ndgrid(o1, o2);
om = [o1(:) o2(:)];
Xd = dtft2(x, om, n_shift);
Xl = dtft2(x, om, n_shift, 1);
printm('loop max %% difference = %g', max_percent_diff(Xl,Xd))
Xf = col(fft2(x));
Xf = Xf .* exp(1i * om * n_shift(:));
printm('fft2 max %% difference = %g', max_percent_diff(Xf,Xd))
