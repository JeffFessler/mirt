 function kern = nufft1_kernel(klist, N, J, K, alpha, beta)
%function kern = nufft1_kernel(klist, N, J, K, alpha, beta)
% compute equivalent interpolation kernel (real part)
% at points in klist

if nargin < 4
	help(mfilename)
	J = 7;	N = 100;	K = 2*N;
	k = linspace(-J/2-1,J/2+1,201)'; % fine sampling for display
	ker = nufft1_kernel(k, N, J, K, 1, 0);
	plot(k, ker), axis tight
	return
end

if ~isvar('alpha') || isempty(alpha), alpha = 1; end
if ~isvar('beta') || isempty(beta), beta = 0.5; end
[alpha beta] = nufft_alpha(N, J, K, alpha, beta);

n_shift1 = 0;

gam = 2*pi/K;

st = nufft_init(gam*klist, N, J, K, n_shift1, 'minmax:user', {alpha}, {beta});

kern = full(st.p(:,1));

%
% undo phase to make it real
%
om = st.om;
if 0
	koff = nufft_offset(om, J, K);		% [M,1]
	jj = 0 - nufft_offset(om, J, K);	% effective j
	oarg = om - gam * (koff + jj);
else
	oarg = gam * klist;
end

kern = kern .* exp(1i * oarg * (N-1)/2);
kern = reale(kern);
