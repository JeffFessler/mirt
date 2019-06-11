 function [out, eig1] = ir_nufft_dpswf1(varargin)
%function [out, eig1] = ir_nufft_dpswf1(varargin)
%|
%| Compute discrete prolate spheroidal wave function (DPSWF) in 1D
%|
%| required
%| 'J' # of filter taps
%| 'N' signal length
%|
%| optional
%| 'K' over-sampled FFT length (default: 2*N)
%| 'M' midpoint (default N if N is odd, N+1 if N is even)
%| 'chat' show pictures
%|
%| out
%| out [N] DPSWF
%| eig1 maximum eigenvalue
%|
%| 2016-11-02, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(varargin{1}, 'test'), ir_nufft_dpswf1_test, return, end

arg.J = []; % required
arg.N = []; % required
arg.K = []; % default 2 * N
arg.M = []; % default N or N+1, to make it odd
arg.chat = 0;

arg = vararg_pair(arg, varargin);

if isempty(arg.J) || isempty(arg.N)
	fail('J and N required')
end
if isempty(arg.K)
	arg.K = 2 * arg.N;
end
if isempty(arg.M)
	arg.M = arg.N + 1 - rem(arg.N, 2);
end

[out, eig1] = ir_nufft_dpswf1_do(arg.J, arg.M, arg.N, arg.K, arg.chat);


% ir_nufft_dpswf1_do()
function [out, eig1] = ir_nufft_dpswf1_do(J, M, N, K, chat)

n = 0:M-1;
m = 0:M-1;
[nn, mm] = ndgrid(n, m);

fun = @(n,m,J) ...
	J * sinc(J/K * (n - m)); % new, for exp
G = fun(nn, mm, J); % gram matrix
[v1, eig1] = eigs(G, 1);

out = v1 * sign(v1(imax(abs(v1)))); % sign flip if needed

if chat && im
	kb_m = 0;
	kb_alf = 2.34 * J;
	sn_kb = nufft_scale_kb(N, J, K, kb_alf, kb_m);

	im plc 2 3
	im(1, 'tick0', n, m, G, 'Gram matrix')
	xlabelf('$N=%d,\ M=%d$', N, M)
	ylabelf('$J=%d,\ K=%d$', J, K)

	im subplot 2
	plot(0:M-1, out, '.-')
	xlim([0 M-1]), xtick([0 M-1])
	xlabelf('$n$')
	ylabelf('$q[n]$')

	% show FT
	kap = linspace(-1,1,601)' * K/2;
	Cfun = @(kap, M) cos(2*pi/K * kap .* ([0:M-1] - (M-1)/2)); % [#kap M]
	kernel = Cfun(kap,M) * out;
	kernel = reale(kernel);

%	kernel_kb = Cfun(kap,M) * [sn_kb; 0];
	kernel_kb = Cfun(kap,N) * sn_kb;
	kernel_kb = reale(kernel_kb);
	kernel_diric = Cfun(kap,M) * ones(M,1);

	subplot(234)
clf
	plot(kap, kernel, '-', kap, kernel_kb, '-', kap, kernel_diric, '-')
	axis tight, xlim([-K/2 K/2])
	xtick([-K/2 -J/2 0 J/2 K/2])
	xlabelf('$\kappa$')
	ir_legend({'DSPWF1', 'KB', 'Diric'})
return

	subplot(235)
	plot(kap, kernel)
	axis tight, xlim([-1 1] * J * 2)
	xtick([-K/2 -J/2 0 J/2 K/2]), grid
	xlabelf('$\kappa$')

	subplot(233)
	plot(0:M-1, 1 ./ out, '.-', 0:N-1, kb, '-o')
	xlim([0 M-1]), xtick([0 M-1])
ylim([0 10])
	xlabelf('$n$')
	ylabelf('$1 / q[n]$')
	ir_legend({'DSPWF1', 'KB'})
end


function ir_nufft_dpswf1_test

J = 6;
N = 64;
M = N+1;
K = 2*N;

[out, eig1] = ir_nufft_dpswf1('J', J, 'M', M, 'N', N, 'K', K, 'chat', 1);
pr eig1
equivs(out, flipud(out))
