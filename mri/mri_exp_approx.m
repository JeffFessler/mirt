  function [B, C, hk, zk] = mri_exp_approx(ti, zmap, LL, varargin)
%|function [B, C, hk, zk] = mri_exp_approx(ti, zmap, LL, [options])
%|
%| Build approximations to exponentials for iterative MR image reconstruction,
%| generalizing "time segmentation" and "frequency segmentation" methods.
%| This is a key part of the Gmri object for field-corrected MR reconstruction.
%| in
%|	ti	[M 1]	sample times
%|	zmap	[N 1]	rate map: relax_map + 2i*pi*field_map
%|			*_map and ti should have reciprocal units.
%|			usually relax_map is 0, so zmap is purely imaginary!
%|	LL		number of components: 1 <= L << N
%|			Or, use {Linit, rmsmax}: give initial L to try,
%|			then increase L until RMS error <= rmsmax < 1.
%|			This is implemented only for 'hist,time,unif' type.
%|
%| options
%|	'ctest'	1|0	return C as [L,nhist] instead of full [L N], for testing
%|	'acorr'	1|0	autocorrelate the fmap histogram (for Toeplitz case)
%|	'tol'	tol	tolerance for pinv(), see pinv_tol() below.
%|			use {'fro'} for tol = 1e-8 * norm(X,'fro') (for large L)
%|	'type'	type	what type of approximation (see choices below)
%|
%| out
%|	B	[M L]	basis functions
%|	C	[L N]	coefficients, such that B * C \approx exp(-ti * zmap.')
%|	hk,zk	[Ko Kr]	histogram values and 'frequencies', if used
%|
%| type of approximation choices:
%|
%|	{'hist,time,unif', nhist}
%|		This recommended (and default) approach uses time segmentation
%|		with uniform time samples and LS coefficients (using histogram).
%|		This approach works almost as well as the SVD method
%|		unless L is chosen so small that the error is large.
%|	'time,unif'
%|		Time segmentation with uniform time samples (LS coef).
%|		No histogram, so it is very slow.  not recommended.
%|
%|	{'hist,svd', nhist}
%|		Nearly optimal approach that uses an SVD and a zmap histogram.
%|		nhist is the # of histogram bins; about 40 is recommended.
%|
%|	{'hist,fs,unif', nhist}
%|		Not recommended since it works poorly except for uniform distn.
%|	{'hist,fs,prctile', nhist}
%|		Not recommended since it works quite poorly.
%|	{'hist,fs,lbg', nhist}
%|		Frequency segmentation methods (exponential bases)
%|		with different choices for the nominal frequency components.
%|		in all cases the coefficients are chosen by LS (Man, MRM, 1997).
%|
%|	for relaxation cases, nhist should be [Ko=#omap Kr=#rmap]
%|
%| Copyright 2004-7-1, Jeff Fessler, University of Michigan

if nargin == 1 && streq(ti, 'test'), mri_exp_approx_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

persistent warned
if isempty(warned), warned = 0; end

% defaults
o.ctest = 0;
o.acorr = 0;
o.chat = 0;
o.tol = {};
% todo: should handle pure real zmap too
if any(real(zmap(:)))
	o.type = {'hist,time,unif', [40 10]};
else
	o.type = {'hist,time,unif', [40]};
end

o = vararg_pair(o, varargin);

atype = o.type{1};

% make doubles for mri_exp_mult_mex happy
zmap = double(zmap);
ti = double(ti);

if iscell(LL)
	rmsmax = LL{2};
	LL = LL{1};
else
	rmsmax = 1/eps;
end

zmap = zmap(:);	% [N 1]
ti = ti(:);	% [M 1]

rmap = real(zmap);
fmap = imag(zmap) / (2*pi);

% histogram the field map
if streq(atype, 'hist,', 5) || o.ctest
	nhist = o.type{2};
	if length(nhist) == 1 % fmap only

		if any(rmap(:)), error 'rmap requires length(nhist)=2', end
		[hk zc] = hist(imag(zmap(:)), nhist);
		zk = 0 + 1i * zc(:); % [K 1]
		if o.chat
			bar(imag(zk)/(2*pi), hk)
		prompt
		end

		if o.acorr
			hk = conv(hk, fliplr(hk)); % autocorr of fmap, formerly xcorr
			zk = [-(nhist-1):(nhist-1)]' * (zk(2)-zk(1));
			if o.chat
				bar(imag(zk)/(2*pi), hk)
				prompt
			end
		end

	else

		[hk zc] = hist_equal([imag(zmap(:)) real(zmap(:))], nhist); % 2d histogram
		zk = outer_sum(1i*zc{1}, zc{2}); % [K1 K2]

		if o.acorr % code by valur olafsson
			hk = conv2(hk, flipdim(hk,1)); % acorr fmap, aconv on r2smap
			zc{1} = [-(nhist(1)-1):(nhist(1)-1)]' * (zc{1}(2)-zc{1}(1));
			zc{2} = linspace(2*min(zc{2}), 2*max(zc{2}), 2*nhist(2) - 1);
			zk = outer_sum(1i*zc{1}, zc{2});
		end

	end

	hk = hk(:);

	Eh = exp(-ti * zk(:).'); % [N K]
end


%
% SVD approach (optimal, to within histogram approximation)
%
if streq(atype, 'hist,svd')
	Ew = Eh * spdiag(sqrt(hk));
	[U S V] = svd(Ew, 0);

	B = U(:,1:LL);	% [M L] keep desired components

%
% time segmentation approaches (recommended)
%
elseif streq(atype, 'time,unif') || streq(atype, 'hist,time,unif')

	pn = jf_protected_names;
	rms = Inf;
	while rms > rmsmax && LL < 40
		ticker(mfilename, LL, 0)

		% time sample locations [0 ... end]
		if LL == 1
			tl = mean(ti);
		else
			p = 100*[0:1:(LL-1)]/(LL-1);
			tl = col(pn.prctile(ti, p));
		end

		if rmsmax < 1 || o.ctest || streq(atype, 'hist,time,unif')
			Ch = exp(-tl * zk(:).'); % [L K]
		end

		if LL > 9 && isempty(o.tol) && ~warned
			warning 'For large L, try tol=''fro'''
		end

		if streq(atype, 'time,unif')
			C = exp(-tl * zmap.'); % [L N] - classic TS
			X = C.'; % [N L]
			% X = X * inv(X'*X); % [N L]
			X = pinv_tol(X, o.tol{:})'; % [N L]
			X = complexify(X);
			B = mri_exp_mult_mex(X, complexify(zmap), ti).'; % [M L]

		elseif streq(atype, 'hist,time,unif')
			W = spdiag(sqrt(hk), 'nowarn');
			P = pinv_tol(W * double(Ch.'), o.tol{:}) * W; % [L K], weighted pinv
			if o.chat
				rP = rank(P);
				if rP ~= LL
					printf('Warn: rank=%d < L=%d', rP, LL)
				end
			end
			% B = exp(-ti*zk.') * P.'; % [M K] * [K L] = [M L]
			P = complexify(P);
			B = mri_exp_mult_mex(P', complexify(zk(:)), ti).';

		else
			error 'bug'
		end

		if rmsmax < 1
			LL = LL + 1;
			rms = sqrt(mean(abs(Eh - B * Ch).^2) * hk(:) / sum(hk(:)));
		else
			break
		end
	end

	if LL == 40, error 'max LL reached!?', end

	if o.ctest
		C = Ch;
	elseif ~streq(atype, 'time,unif')
		C = exp(-tl * zmap.'); % [L N]
	end
return

%
% freq. segmentation approaches
%
elseif streq(atype, 'hist,fs,', 8)

	pn = jf_protected_names;

	if streq(atype, 'hist,fs,unif') % uniform spacing
%		fl = linspace(0, max(abs(fmap)), LL); % not good
%		fl = linspace(min(fmap), max(fmap), LL); % not great
		fl = linspace(min(fmap), max(fmap), LL+2);
		fl = fl(2:end-1);
		rl = 0; % lazy: 'uniform' in 2D seems too arbitrary
			% this may stink if rmap is nonzero.

	elseif streq(atype, 'hist,fs,prctile') % histogram percentiles
		p = 100*[0:1:(LL+1)]/(LL+1); p = p(2:end-1);
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
	B = exp(-ti * zl(:).'); % [M L] bases

else
	error(sprintf('type "%s" unknown', atype))
end


%
% given basis, now compute the LS coefficients, i.e., (B'*B) \ B' * E
%
Bpp = pinv_tol(B, o.tol{:})'; % [M L]
Bpp = complexify(Bpp);
if o.ctest
	C = mri_exp_mult_mex(Bpp, ti, complexify(zk(:))); % [L K] coefficients
else
	C = mri_exp_mult_mex(Bpp, ti, complexify(zmap)); % [L N] coefficients
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
function mri_exp_approx_test
f.dt = 5e-6;
ti = [0:f.dt:25e-3]'; % 25ms readout with 5us sampling
if 1
	fmap = zeros(64,64);
	fmap(32+[-25:25],10:20) = 90;
	fmap(16+[-10:10],30:60) = 30;
	fmap(48+[-10:10],30:60) = 60;
	fmap = 10 + conv2(fmap, ones(5)/5^2, 'same');
	mask = true(size(fmap));
	im clf, pl = 130;
	im(231, fmap, 'Field map'), cbar
end
%rmap = 0;
rmap = [0 0 18 23 0 20; 6 0 8 8 0 3];
ig = image_geom('nx', 64, 'ny', 64, 'dx', 1);
rmap = 1 * ellipse_im(ig, rmap, 'oversample', 3);
im(234, rmap, 'Relax map'), cbar
zmap = rmap + (2i*pi) * fmap;
if rmap == 0
	nhist = 40;
else
	nhist = [40 10];
end
acorr_arg = {'acorr', 1}; % test autocorrelation version
acorr_arg = {};
if 1
	ti = dsingle(ti);
	zmap = dsingle(zmap);
end
[B C hk zk] = mri_exp_approx(ti, zmap, {1, 5e-2}, ...
...%		'tol', {'fro'}, ...
		acorr_arg{:}, ...
		'ctest', 1, 'chat', 1, 'type', {'hist,time,unif', nhist});
Eh = exp(-ti * zk(:).');
Ep = B * C;
err = abs(Eh - Ep);
mse = mean(err.^2);
wrms = sqrt( (mse * hk(:)) / sum(hk(:)) );
subplot(pl+2)
subplot_stack(1000*ti, B, 'Basis components')
nhist = length(zk);
ik = [1 floor(nhist/2) nhist-[0:4]];
subplot(pl+3)
subplot_stack(1000*ti, {Eh(:,ik), Ep(:,ik)}, 'True and Approx.', ...
	{'g', 'r', 'b--', 'y--'})
printf('%s for L=%d rms=%g max=%g', mfilename, ncol(B), wrms, max(err(:)))
