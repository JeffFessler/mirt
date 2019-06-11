  function [B, C, hk, wk] = mri_fun_approx(ti, uj, vj, f1, f2, varargin)
%|function [B, C, hk, wk] = mri_fun_approx(ti, uj, vj, f1, f2, [options])
%|
%| UNDER DEVELOPMENT
%| For certain iterative MR image reconstruction problems, we have
%| signals of the form f1(uj ti) f2(vj ti), for some functions f1() and f2(),
%| and we want to approximate them by \sum_{l=1}^L b_{il} c_{lj},
%| generalizing "time segmentation" and "frequency segmentation" methods.
%| This is a key part of the Gmri object for field-corrected MR reconstruction.
%|
%| in
%|	ti	[M,1]	sample times, with reciprocal units to uj,vj
%|	uj	[N,1]	map 1 (typically a rate map)
%|	vj	[N,1]	map 2 (typically a field map)
%|	f1	handle	function, typically an exponential: exp(- * .)
%|	f2	handle	function, typically a complex exponential: exp(-1i * .)
%|			caution: be sure to include "2 pi" if uj is in Hz!
%|
%| options
%|	'type'	cell	What type of approximation (see choices below).
%|	'L'	int	number of components: 1 <= L << N.  (default: 5)
%|	'nrms'	dbl	Increase L until NRMS error <= nrms. (default: 1/eps)
%|			(Implemented only for some approximation types.)
%|	'Lmax'	int	Maximum 'L' when searching for RMS.  (default: 40)
%|	'nhist'	int	size of histogram: [Ku=#uj Kv=#vj] (default: [20 19])
%|
%|	'ctest'	1|0	return C as [L,nhist] instead of full [L,N], for testing
%|	'tol'	tol	tolerance for pinv(), see pinv_tol() below.
%|			use {'fro'} for tol = 1e-8 * norm(X,'fro') (for large L)
%|
%| out
%|	B	[M,L]	basis functions (orthogonal columns)
%|	C	[L,N]	coefficients
%|		B,C are such that B * C \approx f1(ti * uj.') .* f2(ti * vj.')
%|	hk,wk	[Ku,Kv]	histogram values and coordinates, if used
%|
%| type of approximation choices:
%|
%|	'hist,svd' (default)
%|		Nearly optimal approach that uses an SVD and a 2D histogram
%|		of the uj, vj values to form B.
%|
%|	'hist,time,unif'
%|		Time segmentation with uniform time samples and LS coefficients
%|		(using histogram).
%|		This approach works almost as well as the SVD method
%|		unless L is chosen so small that the error is large.
%|	'time,unif'
%|		Time segmentation with uniform time samples (LS coef).
%|		No histogram, so it is very slow.  not recommended.
%|
%|	'hist,fs,unif'
%|		Not recommended since it works poorly except for uniform distn.
%|	'hist,fs,prctile'
%|		Not recommended since it works quite poorly.
%|	'hist,fs,lbg'
%|		Frequency segmentation methods (exponential bases)
%|		with different choices for the nominal frequency components.
%|		in all cases the coefficients are chosen by LS (Man, MRM, 1997).
%|
%| Copyright 2006-9-8, Jeff Fessler, University of Michigan

if nargin == 1 && streq(ti, 'test'), mri_fun_approx_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

persistent warned
if isempty(warned), warned = 0; end

% defaults
o.ctest = false;
o.chat = false;
o.tol = {};
o.nhist = [20 19];
o.type = 'hist,svd';
o.L = 5;
o.Lmax = 40;
o.nrms = 1/eps;
o.jskip = 10;

o = vararg_pair(o, varargin);

[B, C, hk, wk] = mri_fun_approx_do(ti, uj, vj, f1, f2, ...
	o.type, o.nhist, o.L, o.Lmax, ...
	o.nrms, o.tol, o.jskip, o.ctest, o.chat);

% 
% mri_fun_approx()
%
function [B, C, hk, wk] = mri_fun_approx_do(ti, uj, vj, f1, f2, ...
	atype, nhist, LL, Lmax, nrms, tol, jskip, ctest, chat)

uj = uj(:); % [N,1]
vj = vj(:); % [N,1]
ti = ti(:); % [M,1]

% histogram the maps
if streq(atype, 'hist,', 5) || ctest
	if length(nhist) == 2
		[hk wk] = hist_equal([uj vj], nhist, 'ifsame', '1bin');
	else
		error 'bug'
	end

	[uk vk] = ndgrid(wk{:}); % each [Ku,Kv]
	Sh = f1(ti * uk(:).') .* f2(ti * vk(:).'); % [N,K=Ku*Kv] "signals"
	Sh_rms = sqrt(mean(abs(Sh).^2));
	minmax(Sh_rms)
end


%
% SVD approach (optimal, to within histogram approximation)
%
if streq(atype, 'hist,svd')
	Sw = Sh .* repmat(hk(:)', [nrow(Sh) 1]); % weighted signals
	[U S V] = svd(Sw, 0);

	LL = min(LL, ncol(U));
	if isinf(nrms)
		B = U(:,1:LL);	% [M,L] keep desired components

	else % search over LL until RMS criterion met
		rms = inf;
		while LL <= Lmax
%			printm('L=%d rms=%g Sh_rms=%g', LL, rms, Sh_rms)
			ticker(mfilename, LL, Lmax)

			B = U(:,1:LL);	% [M,L] keep desired components
			Ch = B' * Sh; % coefficients of the "histogram" signals

%			rms = sqrt(mean(abs(Sh - B * Ch).^2) * hc / sum(hc));
			rms = sqrt(mean(abs(Sh - B * Ch).^2)); % unweighted

			if all(rms/Sh_rms <= nrms), break, end
			LL = LL + 1;
		end

		if LL > Lmax, error 'max LL exceeded!?', end
	end


elseif 1
	error('currently only "hist,svd" done!')


%
% time segmentation approaches
% NOT DONE!
%
elseif 0 7& (streq(atype, 'time,unif') || streq(atype, 'hist,time,unif'))

	pn = jf_protected_names;

	while LL <= Lmax
		ticker(mfilename, LL, 0)

		% time sample locations [0 ... end]
		if LL == 1
			tl = mean(ti);
		else
			p = 100*[0:1:(LL-1)]/(LL-1);
			tl = col(pn.prctile(ti, p));
		end

		if rmsmax < 1 || ctest || streq(atype, 'hist,time,unif')
			Ch = exp(-tl * zk(:).'); % [L,K]
		end

		if LL > 9 && isempty(tol) && ~warned
			warning 'For large L, try tol=''fro'''
		end

		if streq(atype, 'time,unif')
			C = exp(-tl * zmap.'); % [L,N] - classic TS
			X = C.'; % [N,L]
			% X = X * inv(X'*X); % [N,L]
			X = pinv_tol(X, tol{:})'; % [N,L]
			X = complexify(X);
			B = mri_exp_mult_mex(X, complexify(zmap), ti).'; % [M,L]

		elseif streq(atype, 'hist,time,unif')
			W = spdiag(sqrt(hk), 'nowarn');
			P = pinv_tol(W * double(Ch.'), tol{:}) * W; % [L,K], weighted pinv
			if chat
				rP = rank(P);
				if rP ~= LL
					printf('Warn: rank=%d < L=%d', rP, LL)
				end
			end
			% B = exp(-ti*zk.') * P.'; % [M,K] * [K,L] = [M,L]
			P = complexify(P);
			B = mri_exp_mult_mex(P', complexify(zk(:)), ti).';

		else
			error 'bug'
		end

%		rms = sqrt(mean(abs(Sh - B * Ch).^2) * hk(:) / sum(hk(:)));
		rms = sqrt(mean(abs(Sh - B * Ch).^2)); % unweighted
		if rms / Sh_rms <= nrms, break, end
		LL = LL + 1;
	end

	if LL > Lmax, error 'max LL reached!?', end

	if ctest
		C = Ch;
	elseif ~streq(atype, 'time,unif')
		C = exp(-tl * zmap.'); % [L,N]
	end
return


%
% freq. segmentation approaches
%
elseif streq(atype, 'hist,fs,', 8)

	if streq(atype, 'hist,fs,unif') % uniform spacing
%		fl = linspace(0, max(abs(fmap)), LL); % not good
%		fl = linspace(min(fmap), max(fmap), LL); % not great
		fl = linspace(min(fmap), max(fmap), LL+2);
		fl = fl(2:end-1);
		rl = 0; % lazy: 'uniform' in 2D seems too arbitrary
			% this may stink if rmap is nonzero.

	elseif streq(atype, 'hist,fs,prctile') % histogram percentiles
		p = 100*[0:1:(LL+1)]/(LL+1); p = p(2:end-1);
		pn = jf_protected_names;
		fl = pn.prctile(fmap, p);
		rl = 0; % lazy again

	elseif streq(atype, 'hist,fs,lbg') % LBG quantization of histogram
		f0 = linspace(min(fmap), max(fmap), LL+2);
		f0 = f0(2:end-1); % initialize with uniform centers
		% this is my "complex enabled" lloyd-max quantizer designer
		if any(rmap)
			if min(fmap) ~= max(fmap)
				r0 = median(rmap);
				z0 = r0 + 2i*pi*f0;
			else
				r0 = linspace(min(rmap), max(rmap), LL+2);
				r0 = r0(2:end-1);
				z0 = r0 + 2i*pi*fmap(1);
			end
			zl = lloyd_max_hist(zmap, z0, nhist);
			rl = real(zl); fl = imag(zl) / (2i*pi);
		else
			fl = lloyd_max_hist(fmap, f0, nhist);
			rl = 0;
		end

	else
		error(sprintf('fs type "%s" unknown', atype))
	end

	zl = rl + (2i*pi) * fl;
	B = exp(-ti * zl(:).'); % [M,L] bases

else
	error('type "%s" unknown', atype)
end


%
% given basis, now compute the LS coefficients, i.e., (B'*B) \ B' * E
%
Binv = pinv_tol(B, tol{:}); % [L,M]
%Bpp = pinv_tol(B, tol{:})'; % [M,L]
%Bpp = complexify(Bpp);

if ctest
	C = Binv * Sh; % [L,K] coefficients
%	C = mri_exp_mult_mex(Bpp, ti, complexify(zk(:)));
else
	nj = length(uj);
	C = zeros(LL, nj);
	for ii=1:jskip:nj
		ticker(mfilename, ii, nj)
		jj = ii:min(ii+jskip,nj);
		C(:,jj) = Binv * (f1(ti * uj(jj).') .* f2(ti * vj(jj).'));
	end
%	C = mri_exp_mult_mex(Bpp, ti, complexify(zmap)); % [L,N] coefficients
end


%
% pinv_tol()
% pinv with tolerance based on frobenious norm
% todo: add diagnostics if cond # is too small!
%
function p = pinv_tol(x, tol)
if ~isvar('tol') || isempty(tol)
	p = pinv(x);
elseif streq(tol, 'fro')
	p = pinv(x, 1e-8 * norm(x, 'fro'));
elseif iscell(tol) && length(tol)==2 && streq(tol{1}, 'fro')
	p = pinv(x, tol{2} * norm(x, 'fro'));
elseif iscell(tol)
	p = pinv(x, tol{:});
else
	error 'unknown tolerance specifier for pinv_tol'
end


%
% self test / example
%
function mri_fun_approx_test
f.dt = 5e-6;
ti = [0:f.dt:25e-3]'; % 25ms readout with 5us sampling

% field map
fmap = zeros(64,64);
fmap(32+[-25:25],10:20) = 90;
fmap(16+[-10:10],30:60) = 30;
fmap(48+[-10:10],30:60) = 60;
fmap = 10 + conv2(fmap, ones(5)/5^2, 'same'); % Hz
mask = true(size(fmap));
im plc 2 3
im(1, fmap, 'Field map'), cbar

% gradient map
gmap = [0 0 18 23 0 20; 6 0 8 8 0 3];
ig = image_geom('nx', 64, 'ny', 64, 'dx', 1);
gmap = 5 * ellipse_im(ig, gmap, 'oversample', 3);
im(4, gmap, 'gradient map'), cbar
nhist = [20 19];
f1 = @ (u) exp(-2i * pi * u); % field inhomogeneity
f2 = @ (u) nufft_sinc(0.4 * u); % 4mm rect slice profile
%profile on
cpu etic
[B C hk wk] = mri_fun_approx(ti, fmap, gmap, f1, f2, 'L', 5, 'nrms', 5e-2, ...
...%		'tol', {'fro'}, ...
		'ctest', 1, 'chat', 0, 'type', 'hist,svd', 'nhist', nhist);
cpu etoc 'approx time'
%profile report
[uk vk] = ndgrid(wk{:});
Sh = f1(ti * uk(:).') .* f2(ti * vk(:).'); % [N,K] "signals"
Sp = B * C; % prediction
err = abs(Sh - Sp);
%mse = mean(err.^2);
%wrms = sqrt( (mse * hk(:)) / sum(hk(:)) );
nrms = sqrt(mean(abs(err).^2)) / sqrt(mean(abs(Sh).^2));
im pl 1 3
im subplot 2
subplot_stack(1000*ti, B, 'Basis components')
nhist = numel(hk);
ik = [1 floor(nhist/2) nhist-[0:4]];
im subplot 3
subplot_stack(1000*ti, {Sh(:,ik), Sp(:,ik)}, 'True and Approx.', ...
	{'g', 'r', 'b--', 'y--'})
printm('L=%d nrms=%g max=%g', ncol(B), nrms, max(err(:)))
