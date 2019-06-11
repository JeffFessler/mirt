 function st = nufft1_build(J, varargin)
%function st = nufft1_build(J, [option])
%| Build 1D LS-NUFFT interpolation coefficients by brute force.
%|
%| in
%|	J		neighborhood size: [-J/2,J/2]
%| option
%|	om	[M 1]	frequency locations, if not provided build fine table
%|	N		# of signal values (default: 2^8)
%|	K		# of DFT frequencies (default: 2*N)
%|	sn	[N 1]	scaling factors (default: uniform)
%|			or 'uniform', 'cos', 'gauss', 'kb'='kb:beatty', 'kb:mm'
%|	wn	[N 1]	weighting factors (default: uniform)
%|	type		'kb'='kb:beatty' or 'kb:mm' or default: 'minmax'
%| out
%|	st	strum with various data and methods
%|
%| methods:
%|	st.eon		[M N] exact exp(1i*omega*n)
%|	st.kap		[M J] "k" value associated with each interpolator value
%|	st.approx	[M N] approximation to exp(1i*omega*n)
%|	st.err		[M 1] (weighted) approximation error for each omega
%|	st.slow1(Fm)	[N 1] slow type1 "nufft" of [M,1] spectral data Fm
%| data:
%|	st.uj		[M J] complex interpolator
%|	st.core		[M 1] real interpolator core
%|
%| Copyright 2006-4-15, Jeff Fessler, University of Michigan

if nargin == 1 && streq(J, 'test'), nufft1_build_test, return, end
if nargin < 4, help(mfilename), error(mfilename), end

% defaults
st.J = J;
st.N = 2^8;
st.K = [];
st.sn = [];
st.wn = [];
st.om = [];
st.type = '';
st = vararg_pair(st, varargin);

if isempty(st.K), st.K = 2 * st.N; end
if isempty(st.wn), st.wn = ones(st.N,1); end
if isempty(st.type), st.type = 'minmax'; end
st.type = lower(st.type);
if streq(st.type, 'kb'), st.type = 'kb:beatty'; end
if isempty(st.om), error 'empty om not done', end

if isempty(st.sn)
	switch lower(st.type)
	case {'kb:beatty', 'kb:mm'}
		st.sn = st.type;
	case 'minmax'
		st.sn = ones(st.N,1);
	otherwise
		error('unknown type "%s", cannot infer sn', st.type)
	end
else
	if ischar(st.sn)
		st.sn = lower(st.sn);
		if streq(st.sn, 'kb'), st.sn = 'kb:beatty'; end
	end
end

if size(st.om,2) ~= 1, error 'om must be Mx1', end

st = nufft1_build_setup_init(st); % misc and sn
switch st.type
	case 'minmax'
		st = nufft1_build_setup_mm(st, st.om);
	case {'kb:beatty', 'kb:mm'}
		st = nufft1_build_setup_kb(st, st.om);
	otherwise
		error('unknown type "%s"', st.type)
end

meth = {
	'approx', @nufft1_build_approx, ...
	'eon', @nufft1_build_eon, ...
	'kap', @nufft1_build_kap, ...
	'err', @nufft1_build_err, ...
	'slow1', @nufft1_build_slow1, ...
};
st = strum(st, meth);


% nufft1_build_setup_init()
function st = nufft1_build_setup_init(st)
st.Nmid = (st.N-1)/2;
st.gam = 2*pi/st.K;
st.nn = [0:st.N-1]';
st.M = length(st.om);
st.koff = nufft_offset(st.om, st.J, st.K);
st = nufft1_build_sn(st);


% nufft1_build_setup_kb()
% Conventional KB interpolator, with reasonably optimized parameters.
function st = nufft1_build_setup_kb(st, om)

[kb_alf kb_m] = nufft1_build_kb_alf(st, st.type);

J = st.J;
for jj=1:J
	otmp = om - st.gam * (st.koff + jj);
	st.core(:,jj) = ...
		kaiser_bessel(otmp / st.gam, J, kb_alf, kb_m, st.K/st.N);
end

% include linear phase term of the ideal interpolator!
otmp = outer_sum(om-st.gam*st.koff, -st.gam*[1:J]); % [M,J]
lamjo = exp(-1i * st.Nmid * otmp); % [M,J]
%lamjo = exp(-1i * (st.N/2) * otmp); % [M,J] % todo: try wrong midpoint
%lamjo = 1 % exp(-1i * st.Nmid * otmp); % [M,J] % todo: try disregarding phase
st.uj = conj(lamjo) .* st.core;


% nufft1_build_setup_mm()
function st = nufft1_build_setup_mm(st, om)

J = st.J;
K = st.K;
N = st.N;
M = st.M;

% [J 1] <= [J N] * [N 1]
st.T = cos(st.gam*[0:J-1]'*(st.nn'-st.Nmid)) * (st.wn .* st.sn.^2);
st.T = toeplitz(st.T); % [J,J]

Joff = (J + rem(J,2)) / 2; % J/2 if even, (J+1)/2 if odd

st.ro = zeros(M,J);
for jj=1:J
	otmp = om - st.gam * (st.koff + jj);
	st.ro(:,jj) = cos(otmp * (st.nn'-st.Nmid)) * (st.wn .* st.sn); % [M 1]
end
st.core = st.ro * inv(st.T); % [M J]

otmp = outer_sum(om-st.gam*st.koff, -st.gam*[1:J]); % [M J]
lamjo = exp(-1i * st.Nmid * otmp); % [M J]
st.uj = conj(lamjo) .* st.core;


% nufft1_build_sn()
% handle default scaling factor types
function st = nufft1_build_sn(st)
if ~ischar(st.sn), return, end

tt = (st.nn - st.Nmid) / st.K;
K_N = st.K / st.N;

switch st.sn
case {'', 'uniform'}
	st.sn = ones(size(tt));

case {'cos', 'cosine'}
	st.sn = 1 ./ cos(pi * tt);

case {'gauss', 'gaussian'}
	[sig, ~, gauss_ft] = nufft_best_gauss(st.J, K_N, 'ft');
	sn_gauss = inline(gauss_ft, 't');
	st.sn = 1 ./ sn_gauss(tt);

case {'kb:beatty', 'kb:mm'}
	if streq(st.sn, st.type, 3) && ~streq(st.sn, st.type)
		warning 'kb type mismatch - are you sure you want this?'
	end
	[kb_alf kb_m] = nufft1_build_kb_alf(st, st.sn);
	st.sn = 1 ./ kaiser_bessel_ft(tt, st.J, kb_alf, kb_m, 1);

otherwise
	fail('unknown type "%s"', st.sn)
end


% nufft1_build_kb_alf()
% "optimal" alpha from beatty:05:rgr, roughly 2.34 * J for K/N=2 !
function [kb_alf, kb_m] = nufft1_build_kb_alf(st, type)
K_N = st.K / st.N;
kb_m = 0;
switch type
case 'kb:beatty'
	kb_alf = pi * sqrt( st.J^2 / K_N^2 * (K_N - 1/2)^2 - 0.8 );
case 'kb:mm'
	[kb kb_alf kb_m] = kaiser_bessel(0, st.J, 'best', 0, K_N);
otherwise
	error 'unknown kb type'
end
% pr kb_alf / J


% nufft1_build_kap()
function kap = nufft1_build_kap(st)
kap = outer_sum(st.om/st.gam-st.koff, -[1:st.J]); % [M J]


% nufft1_build_eon()
% return the exact exponentials: exp(1i * om * n) [M N]
function eon = nufft1_build_eon(st, varargin)
eon = exp(1i * st.om * st.nn');
eon = eon(varargin{:});


% nufft1_build_approx()
% return the NUFFT approximation to exp(1i * om * n) [M N]
function out = nufft1_build_approx(st, varargin)

out = zeros(st.M, st.N);
for jj=1:st.J
	uj = st.uj(:,jj);
	out = out + repmat(uj, [1 st.N]) ...
		.* exp(1i * st.gam * (st.koff + jj) * st.nn');
end
out = out .* repmat(st.sn', [st.M 1]); % [M N] apply scaling
out = out(varargin{:});


% nufft1_build_err()
% return the NUFFT approximation error for each frequency [M,1]
% note: this error is normalized by 1/N
function err = nufft1_build_err(st)
err = sqrt(sum(abs(st.eon - st.approx).^2, 2) / st.N);


% nufft1_build_slow1()
% a slow but easy to implement "type 1 nufft" of [M,L] data Fm
function fn = nufft1_build_slow1(st, Fm)
fn = st.approx.' * Fm; % [N,1] <= [N,M] * [M,L]


% nufft1_build_test
function nufft1_build_test
N = 2^5;
K = 2*N;
J = 5;
om = linspace(0, 2*pi/K, 101)';
types = {'Uniform', 'Cosine', 'Gaussian', 'KB:mm', 'KB:beatty'};
for it=1:length(types)
	type = types{it};
	st = nufft1_build(J, 'om', om, 'N', N, 'K', K, ...
		'sn', types{it}, 'type', 'minmax');
	err(:,it) = st.err;
end
if 1
	st = nufft1_build(J, 'om', om, 'N', N, 'K', K, 'type', 'KB:beatty');
	err(:,it+1) = st.err;
	types{it+1} = 'KB std be';
	st = nufft1_build(J, 'om', om, 'N', N, 'K', K, 'type', 'Kb:mm');
	err(:,it+2) = st.err;
	types{it+2} = 'KB std mm';
end
if im
	semilogy(om/st.gam, err), xlabel '\omega/\gamma', ylabel 'E'
	axisy(10.^[-J -2])
	axisy(10.^[-J 0]) % todo: for trying phase errors
	legend(types{:})
end
