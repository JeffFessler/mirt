 function [smap, sinit] = ir_mri_sensemap_admm_3D_linear(ykj, varargin)
%function [smap, sinit] = ir_mri_sensemap_admm_3D_linear(ykj, [options])
%|
%| Regularized estimation of sensitivity maps for parallel MRI
%| using augmented Lagrangian methods with variable splitting.
%| See M.J. Allison et al., IEEE TMI, 32(3), 556-564, Mar. 2013.
%| (This code is less optimized/more general than the one used for the paper.)
%|
%| Signal model for "surface coil" image: y_j = s_j f_j + noise_j
%|	s_j	sensitivity map (relative to "bodycoil" reference image)
%|	f_j	unknown underlying object
%|	j: voxel index
%| Signal model for "body coil" image: z_j = f_j + noise_j.
%|
%| Cost function:
%|	argmin 1/2 ||z - Ds||_W^2 + lam/2 ||R s||^2
%| where:
%|	D = diag{y_j}
%|	W = diag{w_j}
%|	R is finite differencing matrix with non-periodic boundary conditions.
%|
%| Augmented Lagrangian for ADMM with u1=s and u0=C*s:
%|	L(u,s) = 1/2 ||h - A u||^2 + 1/2 ||u - Ts - eta||_V^2
%|	so u = [u1; u0], T = [I; C], h = [sqrt(W)*z; 0],
%|	A = [sqrt(W)*D 0 ; 0 sqrt(lam)*B], eta = [eta1; eta0],
%|	and V = [v1*I 0; 0 v0*I], where R = B*C and B has boundary effects
%|
%| in
%|	ykj	[nx ny nz ncoil]	noisy complex images for each coil
%|
%| options
%|	bodycoil	[nx ny nz]	body coil (reference) image
%|				(def: sum-of-squares with phase of 1st coil)
%|	init			type of initialization
%|					 (see ir_mri_sensemap_init.m)
%|	l2b			log_2(beta), regularization parameter (def: -5)
%|	order			regularization order: 1|2 (default: 2)
%|	niter			# of iterations (default: 500)
%|	isave			which iterations to save (number
%|					list, 'last', 'all') (default: last)
%|	thresh			thresholding value for initialization (default: 0.05)
%|	quiet			hide iteration display if 1 (default:0)
%|	maskD			mask specifying pixels that include real data
%|	etabtw			update eta between each alternating
%|					minimization step (def: 1)
%| out
%|	smap	[nx*ny ncoil (isave)]	denoised sensitivity maps
%|	sinit	[nx*ny ncoil]		initialization images
%|
%| Copyright 2013, Michael Allison, University of Michigan
%| 2015-08 Lianli Liu  Extended to 3D using LADMM method

% todo; see fessler/l/tex/paper/submit/liu,bias/
% argsADMM = { 'bodycoil', bodycoil_s, 'maskD', maskObj_s,        'conds', ksADMM, 'condu0', ku0ADMM,'l2b', l2b,'niter',niter,...

if nargin == 1 && streq(ykj, 'test'), ir_mri_sensemap_admm_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% defaults
arg.bodycoil = [];
arg.l2b = -5;
arg.order = 2;
arg.niter = 500;
arg.isave = [];
arg.quiet = 0;
arg.thresh = 0.12;
arg.init = [];
arg.maskD = true(size(ykj(:,:,1)));
arg.v0 = [];
arg.v1 = [];
arg.conds = [];
arg.condu0 = [];
arg.etabtw = 1;

arg = vararg_pair(arg, varargin);

%% Check for common cases of misuse.
if ~islogical(arg.maskD)
	error('arg.maskD must be logical datatype')
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
[nx ny nz ncoil] = size(ykj);
N = nx * ny * nz;
M = sum(arg.maskD);

% default reference image
if isempty(arg.bodycoil)
	arg.bodycoil = sqrt(sum(abs(ykj).^2, 3)); % sum-of-squares
	tmp = angle(ykj(:,:,1)); % crude estimate of phase of image f_j
	arg.bodycoil = arg.bodycoil .* exp(1i * angle(ykj(:,:,1)));
end
arg.bodycoil = arg.bodycoil(:); % lexicographical

%% create regularizer
%tic;
if nz > 1
	[C,B,phi] = ir_circ_zeroed_reg_3D(arg.order, true([nx ny nz]), 0);
else
	[C,B,phi] = ir_circ_zeroed_reg(arg.order, true([nx ny nz]), 0);
end
%printm('regulizer initialization:');
%toc

%% create output matix
smap = zeros(N,ncoil,length(arg.isave));
sinit = zeros(N,ncoil);

%% determine v0/v1 from condition numbers if needed (untimed)
if ~isempty(arg.condu0)
	arg.v0 = 2^arg.l2b / (arg.condu0 - 1);
end
if ~isempty(arg.conds)
	arg.v1 = (arg.v0 * (max(abs(phi)) - arg.conds * min(abs(phi)))) ...
		/ (arg.conds - 1);
end
condu0Act = (2^arg.l2b + arg.v0) / arg.v0;
condsAct = (max(abs(phi)) + arg.v1/arg.v0) / ...
		(min(abs(phi)) + arg.v1/arg.v0);
condu1Act = (1 + arg.v1) / arg.v1;
printm(['Actual v0:' num2str(arg.v0) ' v1:' num2str(arg.v1)]);
printm(['Actual K(u0):' num2str(condu0Act) ' K(s):' num2str(condsAct) ...
	' K(u1):' num2str(condu1Act)]);

%% ADMM Algorithm
for i = 1:ncoil % it is better to do this in parallel as they are independent.
	%yj = ykj(:,:,i);
	yj = ykj(:);

	%% initialize AL algorithm

	% initialize smap
	sk = ir_mri_sensemap_init_3D(arg.init, yj, arg.bodycoil, [nx ny nz], ...
			arg.thresh, arg.maskD);
	sk = sk(:);
	sinit(:,i) = sk;

	u1 = sk;
	u0 = C * u1;
	eta0 = zeros(size(u0));
	eta1 = zeros(N,1);

	% apply maskD
	yj = yj(arg.maskD);
	arg.bodycoilM = arg.bodycoil(arg.maskD);

	% init convergence params
	v1 = arg.v1;
	v0 = arg.v0;
	v1ov0 = arg.v1 / arg.v0;
	lamov0 = (2^arg.l2b) / arg.v0;

	% precompute inverse terms that dont change with coil or iteration
	if ~exist('GinvS','var')
		Ginv1 = zeros(N,1);
		Ginv1(arg.maskD) = abs(arg.bodycoilM).^2;
		Ginv1 = Ginv1 + arg.v1;
		Ginv1 = 1 ./ Ginv1;

		GinvS = phi(:) + v1ov0;
		GinvS = 1 ./ GinvS;
		GinvS = reshape(GinvS, [nx ny nz]);

		Ginv0 = lamov0 * B + 1; % B is 0 or 1 and DIAGONAL
		Ginv0 = 1 ./ Ginv0;
	end

	DMz = zeros(N,1);
	DMz(arg.maskD) = conj(arg.bodycoilM) .* yj;

	if any(arg.isave == 0)
		smap(:,i,arg.isave == 0) = sinit(:,i);
	end
    
	printm('pre computation:');
	%toc

	% begin iterations
	for j = 1:arg.niter

		% update sk
% 		sk = ifftn(GinvS .* fftn(reshape(C' * (u0 - eta0) + v1ov0 * ...
% 			(u1 - eta1), [nx ny nz])));
        	rou = 1./49;
        	%sk = sk-rou.*(sk+v1ov0*(C'*C)*sk-C' * (u0 - eta0) - v1ov0 * (u1 - eta1));
        	temp = v1ov0*sk+(C'*C)*sk-v1ov0*(u1-eta1)-C'*(u0-eta0);
        	sk = sk-rou*temp;
        
		sk = sk(:);

		% this is used often
		Cs = C * sk;

		if arg.etabtw
			eta0 = eta0 - (u0 - Cs);
			eta1 = eta1 - (u1 - sk);
		end

		% update u1 and u0
		u1 = Ginv1 .* (DMz + arg.v1 * (sk + eta1));
		u0 = Ginv0 .* (Cs + eta0);

		% update etas
		eta0 = eta0 - (u0 - Cs);
		eta1 = eta1 - (u1 - sk);

		% save subset of iterations
		if any(arg.isave == j)
			smap(:,i,arg.isave == j) = sk; % eta1;
		end

		% display progress counter
		if (mod(j,500) == 1) && (~arg.quiet)
		%	printm('outputting eta1');
			printm([num2str(j) ' of ' num2str(arg.niter)])
		end
	end % for j
end % for i

end % ir_mri_sensemap_admm()


% create_mask()
function mask = create_mask(img, thresh, dial)
img = abs(img);
img = img ./ max(max(img));

mask = (img >= thresh);

if dial > 0
	SE = strel('disk', dial, 0);

	mask = imdilate(mask, SE);
	mask = logical(mask);
end

end % create_mask()


function ir_mri_sensemap_admm_test()
%| A test function for the toolbox version of ir_mri_sensemap_admm.
%| Note: This function is not optimized like the one in the paper
%| (uses Fatrix etc.) and thus the results of the test will not match those published.
%|
%| Copywrite, M Allison, U of Michigan, 2013.

%% set parameters
ns = [0.0275,0.005]; % noise level for [bcoil, surfcoil] images
threshD = 0.13;
sCoil = 1;
order = 2;
l2b = 5;
niter = 1000;
ksADMM = 650;
ku0ADMM = 255;

%% create senstivity data
% create data (based on M Allison file generateSENSEdata.m)
printm('Creating sensitivity data.')
f.dir = [path_find_dir('mri') '/../data/mri/'];
f.xtrue = [f.dir 'brainweb_t1.jpg'];
ftrue = double6(imread(f.xtrue)');
ftrue = ftrue(2:end-1,2:end-1); % make it 256^2

% add phase component to real image
[sX, sY] = size(ftrue);
tmp1 = [1:sX] ./ sX;
tmp2 = [1:sY] ./ sY;
[X,Y] = meshgrid(tmp1,tmp2);
ctrue = (Y .* X .* pi) - pi/2; ctrue = ctrue';
ftrue = ftrue.*exp(j*ctrue);
clear tmp1; clear tmp2; clear X; clear Y; clear ctrue

%normalize image
ftrue = ftrue ./ max(max(abs(ftrue)));

%create sense maps (uses mri_sensemap_sim.m from toolbox)
ncoil = 4;
smap = mri_sensemap_sim('nx', sX, 'ny', sY, 'dx', 192/sX, ...
		'ncoil', ncoil, 'orbit', 90*ncoil, 'rcoil', 100);
if ncoil == 4, smap = smap(:,:,[4 1 3 2]); end
smap = smap ./ max(max(max(abs(smap))));

bcoilI = ftrue + ns(1) * (randn(size(ftrue)) + 1i * randn(size(ftrue)));
lcoilI = smap .* repmat(ftrue, [1 1 ncoil]) + ns(2) * ...
	(randn(size(smap)) + 1i * randn(size(smap)));

% flip so it is in correct orientation for radiologists
ftrue = fliplr(ftrue);
bcoilI = fliplr(bcoilI);
lcoilI = flipdim(lcoilI,2);
smap = flipdim(smap,2);

%% create mask of real data
maskObj = create_mask(bcoilI, threshD, 0);

%% normalize coils and apply W mask outside of code (for Cholesky)
bcoilZ = bcoilI;
mB = max(max(abs(bcoilI)));
bcoilI = bcoilI ./ mB;
lcoilI = lcoilI ./ mB;
bcoilZ(~maskObj) = 0;

%% select a coil to work with
lcoilI = lcoilI(:,:,sCoil);
smap = smap(:,:,sCoil);
printm('-- Done')

%% Compute the Cholesky based estimate for this small problem
if 1
	printm('Starting Cholesky based method')
	printm('-- Computing C');

	[jnk,jnk,jnk,C] = ir_circ_zeroed_reg(order,true(size(bcoilI)),1); % same as above three lines

	% implement maskR
	bI = col(bcoilZ(:));
	sI = col(lcoilI(:));

	% compute using sparse Cholesky
	printm('-- Computing H');
	bI2 = conj(bI) .* bI;
	H = spdiag(bI2) + 2^(l2b) * C' * C;
	printm('-- Computing rhs');
	rhs = (conj(bI) .* sI);
	printm('-- Computing lhs');
	tic;
	lhs = H \ double(rhs);
	toc
	est_chol = lhs;
	est_chol = embed(est_chol,true(size(bcoilI)));
	printm('-- Done')
end

%% Run the ADMM Estimator
printm('Starting ADMM based method')
argsADMM = {'order', order, 'bodycoil', bcoilI, 'l2b', ...
		l2b, 'niter', niter, 'isave', 'last', 'thresh', threshD, ...
		'conds', ksADMM, 'condu0', ku0ADMM, 'maskD', maskObj, 'etabtw', ...
	1, 'init', 'avg'};
tADMM= tic;
[est_admm] = ir_mri_sensemap_admm(lcoilI, argsADMM{:});
toc(tADMM)
est_admm = reshape(est_admm,size(bcoilI));
printm('-- Done')

%% compare results to truth and Cholesky estimator
im pl 2 3
im(1, smap)
im(2, est_chol)
im(3, est_admm)
im(5, est_chol-smap)
im(6, est_admm-smap)

end % ir_mri_sensemap_admm_test()
