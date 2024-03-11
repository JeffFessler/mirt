 function sn = nufft_scale_kb(Nd, Jd, Kd, kb_alf, kb_m)
%function sn = nufft_scale_kb(Nd, Jd, Kd, kb_alf, kb_m)
%|
%| Compute KB scaling factors for NUFFT.
%| (No scaling along any dimension where J=1, e.g., for stack of spirals.)
%|
%| in
%|	N,J,K	[d]
%|	kb_alf	[d]
%|	kb_m	[d]
%|
%| out
%|	sn	[[Nd]]		scaling factors
%|
%| Copyright 2004-7-8, Jeff Fessler, University of Michigan
%| 2023-08-29, JF: handle J=1 cases

if nargin == 1 && streq(Nd, 'test'), nufft_scale_kb_test, return, end
if nargin < 5, ir_usage, end

if all(Jd == 1)
    sn = ones([Nd 1]);
    return
end

% scaling factors: "outer product" of 1D vectors
sn = 1;
dd = length(Nd);
for id = 1:dd
	if Jd(id) == 1
		tmp = ones(Nd(id),1);
	else
		Nmid = (Nd(id)-1)/2;
		nc = (0:Nd(id)-1)' - Nmid;
		tmp = 1 ./ kaiser_bessel_ft(nc/Kd(id), Jd(id), kb_alf(id), kb_m(id), 1);
	end
	sn = sn(:) * tmp';
end
if length(Nd) > 1
	sn = reshape(sn, Nd); % [(Nd)]
else
	sn = sn(:); % [N 1]
end


% nufft_scale_kb_test
function nufft_scale_kb_test

N = [8 4];
J = [1 1];
K = 2 * N;
kb_alf = 2.34 * J;
kb_m = zeros(2);
sn = nufft_scale_kb(N, J, K, kb_alf, kb_m);
assert(all(sn == ones(N), 'all'))

J = [6 1];
s2 = nufft_scale_kb(N, J, K, kb_alf, kb_m);
s1 = nufft_scale_kb(N(1), J(1), K(1), kb_alf(1), kb_m(1));
assert(all(s2 == s1(:) * ones(1,N(2)), 'all'))

N = 64; K = 2*N;
Jlist = 1:7;
kb_m = 0;
leg = {};
for jj=1:numel(Jlist)
	J = Jlist(jj);
	kb_alf = 2.34 * J;
	sn = nufft_scale_kb(N, J, K, kb_alf, kb_m);
	sn = reale(sn);
	plot(0:N-1, sn, '-o')
	if jj == 1, hold on, end
	leg{jj} = sprintf('$J = %d$', J);
end
hold off
xlim([0 N]), xtick([0 N/2 N-1])
ylim([0 1.2]), ytick([0 0.5 1])
xlabelf('$n$'), ylabelf('$s[n]$')
ir_legend(leg, 'location', 'best'), grid
