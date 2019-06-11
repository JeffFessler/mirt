 function [smap, sinit] = ir_mri_sensemap_admm(ykj, varargin)
%function [smap, sinit] = ir_mri_sensemap_admm(ykj, [options])
%|
%| Regularized estimation of sensitivity maps for parallel MRI
%| using augmented Lagrangian methods with variable splitting.
%| See M.J. Allison et al., IEEE TMI, 32(3), 556-564, Mar. 2013.
%| (This code is less optimized/more general than the one used for the paper.)
%|
%| Signal model for "surface coil" image: y_kj = s_kj f_j + noise_j
%|	s_j	sensitivity map (relative to "bodycoil" reference image)
%|	f_j	unknown underlying object
%|	j: voxel index, k=1:#coils
%| Signal model for "body coil" image: z_j = f_j + noise_j.
%|
%| Cost function for each coil (done separately):
%|	argmin_s 1/2 ||z - Ds||_W^2 + lam/2 ||R s||^2
%|	subject to (optional) constraint that s is zero outside optional maskS
%| where:
%|	D = diag{y_kj} for one k at a time
%|	W = diag{w_j}
%|	R is finite differencing matrix with non-periodic boundary conditions.
%|	R = B * C where C'C is circulant and B has boundary effects
%|
%| Augmented Lagrangian for ADMM with u1=s and u0=C*s:
%|	L(u,s) = 1/2 ||h - A u||^2 + 1/2 ||u - Ts - eta||_V^2
%|	so u = [u1; u0], T = [I; C], h = [sqrt(W)*z; 0],
%|	A = [sqrt(W)*D 0 ; 0 sqrt(lam)*B], eta = [eta1; eta0],
%|	and V = [v1*I 0; 0 v0*I], where R = B*C and B has boundary effects
%|	subject to maskS support constraint on u1
%|	(support constraint not in original paper but is an easy extension)
%|
%| in
%|	ykj	[(N) ncoil]	noisy complex images for each coil
%|
%| options
%|	bodycoil [(N)]	body coil (reference) image
%|				(def: sum-of-squares with phase of 1st coil)
%|	init		type of initialization
%|				(see ir_mri_sensemap_init.m)
%|	l2b		log_2(beta), regularization parameter (def: -5)
%|	order		regularization order: 1|2 (default: 2)
%|	niter		# of iterations (default: 500)
%|	isave		which iterations to save (number
%|				list, 'last', 'all') (default: last)
%|	thresh		thresholding value for initialization (default: 0.05)
%|	maskD [(N)]	mask specifying pixels that include real (useful) data
%|				default: entire image, which is bad choice!
%|	maskS [(N)]	mask specifying pixels to estimate (smap support)
%|				default: entire image, as shown in paper
%|	etabtw		update eta between each alternating
%|					minimization step (def: 1)
%|	chat		verbosity (default: 0)
%|
%| out
%|	smap	[*N ncoil (isave)]	denoised sensitivity maps
%|	sinit	[*N ncoil]		initialization images
%|					(caller must reshape afterwards...)
%|
%| Copyright 2013, Michael Allison, University of Michigan
%| 2015-08-12 JF maskS constraints
%| 2017-10-26 JF generalized to support 3D (hopefully)

if nargin == 1 && streq(ykj, 'test')
	run_mfile_local('ir_mri_sensemap_admm_test')
return, end
if nargin < 2, ir_usage(), end

% defaults
arg.bodycoil = [];
arg.l2b = -5;
arg.order = 2;
arg.niter = 500;
arg.isave = [];
arg.thresh = 0.12;
arg.init = [];
arg.maskD = [];
arg.maskS = [];
arg.v0 = [];
arg.v1 = [];
arg.conds = [];
arg.condu0 = [];
arg.etabtw = 1;
arg.chat = 0;

arg = vararg_pair(arg, varargin);

% Try to infer surface coil image size from any of bodycoil, maskS, maskD.
% Otherwise default to using all but last dimensions of ykj.
% Also check that the sizes are consistent.
sdim = ir_mri_sensemap_admm_sdim(ykj, arg.bodycoil, arg.maskD, arg.maskS);

if isempty(arg.maskD)
	arg.maskD = true(sdim);
	warn('Are you sure the entire "bodycoil" reference is good???')
end

if isempty(arg.maskS)
	arg.maskS = true(sdim);
end


%% Check for common cases of misuse.
if ~islogical(arg.maskD)
	fail('arg.maskD must be logical datatype')
end

if ~isempty(arg.conds)
	if ~isempty(arg.v1)
		fail(['You are specifying the condition number for s, ' ...
			'so you cannot specify v1 directly.'])
	end
end
if ~isempty(arg.condu0)
	if ~isempty(arg.v0)
		fail(['You are specifying the condition number for u0, ' ...
			'so you cannot specify v0 directly.'])
	end
end

% create parameters
arg.isave = iter_saver(arg.isave, arg.niter);
%[nx ny ncoil] = size(ykj);
%sdim = [nx ny]; % pre 2017-10-26
ncoil = size(ykj, numel(sdim)+1); % might be 1 if single surface coil
%pr size(ykj)
%pr sdim
%pr ncoil
N = prod(sdim);
%M = sum(arg.maskD);

% initialize smaps and bodycoil image (reference) if not given
[sinit, arg.bodycoil] = ir_mri_sensemap_init(ykj, ...
	'bodycoil', arg.bodycoil, 'type', arg.init, 'thresh', arg.thresh);
sinit = reshape(sinit, [], ncoil); % [*N ncoil]
ykj = reshape(ykj, [], ncoil); % [*N ncoil]

%{
switch ndims(arg.bodycoil)
case 2
case 3
otherwise
	fail 'bug'
end
%}

% create regularizer
[C, B, phi] = ir_reg_diff_zeroed(sdim, 'order', arg.order, 'mask', arg.maskS);
B = B(:);

%% create output matix
smap = zeros(N, ncoil, numel(arg.isave), 'single');

%% determine v0/v1 from condition numbers if needed (untimed)
if ~isempty(arg.condu0)
	arg.v0 = 2^arg.l2b / (arg.condu0 - 1);
end
if ~isempty(arg.conds)
	arg.v1 = arg.v0 * (max(abs(phi(:))) - arg.conds * min(abs(phi(:)))) ...
		/ (arg.conds - 1);
end
condu0Act = (2^arg.l2b + arg.v0) / arg.v0;
condsAct = max(abs(phi(:)) + arg.v1/arg.v0) / ...
		(min(abs(phi(:))) + arg.v1/arg.v0);
condu1Act = (1 + arg.v1) / arg.v1;
printm(['Actual v0:' num2str(arg.v0) ' v1:' num2str(arg.v1)]);
printm(['Actual K(u0):' num2str(condu0Act) ' K(s):' num2str(condsAct) ...
	' K(u1):' num2str(condu1Act)]);

% init convergence params
v1ov0 = arg.v1 / arg.v0;
lamov0 = (2^arg.l2b) / arg.v0;

% precompute terms that dont change with coil or iteration
bodycoilM = arg.bodycoil(arg.maskD); % [nd] apply maskD
Ginv1 = abs(arg.bodycoil(:)).^2 .* arg.maskD(:);
Ginv1 = 1 ./ (Ginv1 + arg.v1); % [*N] for u1 update "D2" in Fig. 2
GinvS = 1 ./ (phi + v1ov0); % [(N)] for s update inv(Phi2) in Fig. 2
% unique(B(:)) % trick: B is all 0 and 1 so B^2 = B in Fig. 2
Ginv0 = 1 ./ (lamov0 * B + 1); % [*N] B2 in Fig 2, for u0 update


%% ADMM Algorithm
for ic = 1:ncoil % it is better to do this in parallel as they are independent.
	yj = ykj(arg.maskD,ic); % [nd] masked

	sk = sinit(:, ic); % [*N] initial guess
	sk = sk .* arg.maskS(:); % enforce support constraint on initial

	u0 = C * sk;
	u1 = sk;
	eta0 = zeros(size(u0));
	eta1 = zeros(N,1);

	DMz = zeros(N, 1, 'single');
	DMz(arg.maskD) = single(conj(bodycoilM) .* yj);

	if any(arg.isave == 0)
		smap(:,ic,arg.isave == 0) = sk;
	end

	% begin iterations
	for iter = 1:arg.niter

		% update sk
		sk = C' * (u0 - eta0) + v1ov0 * (u1 - eta1); % rhs of (12)
		sk = ifftn(GinvS .* fftn(reshape(sk, sdim)));
		sk = sk(:);

		Cs = C * sk; % used often

		if arg.etabtw
			eta0 = eta0 - (u0 - Cs);
			eta1 = eta1 - (u1 - sk);
		end

		% update u1 and u0
		u1 = Ginv1 .* (DMz + arg.v1 * (sk + eta1));
		u1 = u1 .* arg.maskS(:); % enforce support constraint
		u0 = Ginv0 .* (Cs + eta0);

		% update etas
		eta0 = eta0 - (u0 - Cs);
		eta1 = eta1 - (u1 - sk);

		if any(arg.isave == iter) % save subset of iterations
			smap(:,ic,arg.isave == iter) = sk;
		end

		if arg.chat % display progress counter
			ticker(mfilename, iter, arg.niter)
		end
	end % iter
end % ic

end % ir_mri_sensemap_admm()


% ir_mri_sensemap_admm_sdim()
function sdim = ir_mri_sensemap_admm_sdim(ykj, bodycoil, maskD, maskS)
if ~isempty(bodycoil)
	sdim = size(bodycoil);
elseif ~isempty(maskD)
	sdim = size(maskD);
elseif ~isempty(maskS)
	sdim = size(maskS);
else
	warn('assuming ncoil > 1')
	tmp = size(ykj);
	sdim = tmp(1:(end-1)); % assumes last dimension is ncoil
	pr sdim
end

% now check for dimension consistency
if ~isempty(bodycoil)
	if any(size(bodycoil) ~= sdim)
		pr sdim
		pr size(bodycoil)
		fail('bad bodycoil size')
	end
end

if ~isempty(maskD)
	if any(size(maskD) ~= sdim)
		pr sdim
		pr size(maskD)
		fail('bad maskD size')
	end
end

if ~isempty(maskS)
	if any(size(maskS) ~= sdim)
		pr sdim
		pr size(maskS)
		fail('bad maskS size')
	end
end

tmp = size(ykj);
if any(tmp(1:numel(sdim)) ~= sdim)
	pr sdim
	pr size(ykj)
	fail('bad ykj size')
end

end % ir_mri_sensemap_admm_sdim()
