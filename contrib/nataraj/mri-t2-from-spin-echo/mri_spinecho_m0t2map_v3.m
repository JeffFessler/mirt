 function [M0_ml, M0_rls, T2_ml, T2_rls, ml_time, rls_time, cost] ...
	= mri_spinecho_m0t2map_v3(y, TE, varargin)
%function [M0_ml, M0_rls, T2_ml, T2_rls, ml_time, rls_time, cost] ...
%	= mri_spinecho_m0t2map_v3(y, TE, varargin)
%|
%| Penalized weighted least-squares T2 mapping from (single) spin-echo
%| image datasets, using alternating minimization.
%|
%| Can be used for multiple spin echo data as well if stimulated echo
%| contributions are neglected, though this can yield substantial bias.
%| A typical ad hoc workaround is to weight the shorter echo times less,
%| which is supported by the "weights" option.
%|
%| A more rigorous (and complicated) solution based on Bloch simulations:
%| Ben-Eliezer N., Sodickson, D. K., and Block, K. T.
%| Magn. Reson. Med. 2015 73(2):809-17. doi: 10.1002/mrm.25156.
%|
%| Signal model: y_jm = M0_j * exp(-TE_m / T2_j) + e_jm
%| M0_j	spin density at jth voxel
%| T2_j	spin-spin relaxation time at jth voxel
%| TE_m	Echo time of mth dataset
%|
%| Inputs
%| y	[(odims) M]	Complex coil-combined data (2D or 3D)
%| TE	[M]		Echo times (ms)
%|
%| Options
%| M0_init [(odims)]	Initial M0 guess	(def: ML estimate)
%| T2_init [(odims)]	Initial T2 guess (ms)	(def: ML estimate)
%| flip	[(odims)]	Excitation flip map	(def: pi/2 constant)
%| mask	[(odims)]	Data mask		(def: true(odims))
%| T2min [1]		Min T2 (large->fast)	(def: 10^0)
%| T2max [1]		Max T2 (small->fast)	(def: 10^3.5)
%| n_dict [1]		Num dictionary elements (def: 300)
%| weights [M]		Echo weights		(def: ones(M, 1))
%| beta_m [1]		M0 reg parameter	(def: M * 10^-9)
%| beta_2 [1]		T2 reg parameter	(def: M * 10^-8)
%| delta_m [1]		M0 reg curvature	(def: 10^-0.5)
%| delta_2 [1]		T2 reg curvature	(def: 10^0.5)
%| thresh [1]		Fraction of mag max
%|			used for median initial
%|			value in "background"	(def: 0.05)
%| n_outer [1]		Max num of outer iter	(def: 20)
%| n_innerM[1]		Max num of M0 iter	(def: 50)
%| n_inner2[1]		Max num of T2 iter	(def: 50)
%| tolM	[1]		M0 stop tolerance	(def: 10^-7)
%| tol2	[1]		T2 stop tolerance	(def: 10^-7)
%| norm	0|1		Normalize data on/off
%|			to make regularization
%|			scale-invariant		(def: 1)
%| is_mag 0|1		Use magnitude sig model (def: 1)
%| disp	0|1		Show image updates	(def: 1)
%|
%| Outputs
%| M0_ml [(odims)]	Max-likelihood M0 estimate
%| M0_rls [(odims)]	Regularized LS M0 estimate
%| T2_ml [(odims)]	Max-likelihood T2 estimate (ms)
%| T2_rls [(odims)]	Regularized LS T2 estimate (ms)
%| ml_time [1]		ML recon time (ms) (if used; NaN otherwise)
%| rls_time[1]		RLS recon time (ms)
%| cost	[2*n_outer+1]	Cost function (should monotonically descend!)
%|
%| Related work:
%| Nataraj G., Nielsen J. F., and Fessler, J. A.
%| Optimizing MR Scan Design for Model-Based T1,T2 Estimation
%| from Steady-State Sequences. In preparation.
%|
%| Version 1.1: 2014-10-14 Simple PGD routine w/ curvature bounded at min T2
%| Version 2.1: 2015-06-10 Flip angle compensation added, PGD corrected
%|	2.2: 2015-06-16 Added is_mag option to helper functions
%|	2.3: 2015-06-17 Curvature bound wrong: corrected with grid search
%|	2.4: 2015-08-04 Trick: normalize data before processing
%|					Outer iteration early exit added
%| Version 3.1: 2015-08-27 Optional ML initialization added
%|
%| Written By:		Gopal Nataraj
%|			University of Michigan
%|			Copyright 2015

% Dimensions
tmp = size(y);
odims = tmp(1:end-1);	% Object dimensions
M = tmp(end);		% Number of echo datasets
if (M ~= length(TE))
	error('Number of datasets does not match number of echoes.');
end

% Default values
arg.M0_init = [];
arg.T2_init = [];
arg.flip = pi/2 * ones(odims);
arg.mask = [];
arg.T2min = 10^0;
arg.T2max = 10^3.5;
arg.n_dict = 300;
arg.weights = ones(M, 1);
arg.beta_m = M * 10^-9;
arg.beta_2 = M * 10^-8;
arg.delta_m = M * 10^-0.5;
arg.delta_2 = M * 10^0.5;
arg.thresh = 0.05;
arg.n_outer = 20;
arg.n_innerM = 50;
arg.n_inner2 = 50;
arg.tolM = 10^-7;
arg.tol2 = 10^-7;
arg.norm = true;
arg.is_mag = true;
arg.disp = true;

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% If no mask specified, use all voxels
% Set N to number of voxels of interest
if isempty(arg.mask)
	arg.mask = true(odims);
	N = prod(odims);
else
	N = numel(arg.mask(arg.mask));
end

% Reshape and extract data and flip angles of interest
y = reshape(y(repmat(arg.mask, [ones(1, length(odims)) M])), [N M]);
flip = arg.flip(arg.mask);

% Create sparse weighting matrix [N*M N*M]
W = Gdiag(col(repmat(arg.weights', [N 1])));

% Set T2_init from ML estimate, if necessary
if isempty(arg.T2_init)
	% For dictionary, use median of excitation flips
	med_flip_ex = median(flip);

	% Create a dictionary of pathways
	T2dict = logspace(log10(arg.T2min), log10(arg.T2max), arg.n_dict)';
	D = NaN(M, arg.n_dict);
	for t2 = 1:length(T2dict)
		for m = 1:M
			D(m, t2) = f_se(1, T2dict(t2), med_flip_ex, TE(m), arg.is_mag);
		end
	end

	% ML Estimation
	% Dictionary-based estimation via variable-projection method
	Wdict = spdiags(arg.weights, 0, M, M); % Single-voxel weights
	ydict = y.';				% [M N]
	maxProd = zeros(N, 1);			% Max inner product
	t2_idx = ones(N, 1);						% Corresponding index
	tic;										% Start ML timer

	% Efficiently find maximum inner product for each voxel
	for k = 1:arg.n_dict
		% Compute kth inner product
		hess = abs(D(:,k)' * Wdict * D(:,k));
		ytild = (D(:,k)' * Wdict * ydict).' / sqrt(hess);
		newProd = abs(ytild).^2;

		% If the kth inner product is largest, save k
		update = newProd > maxProd;
		maxProd(update) = newProd(update);
		t2_idx(update) = k;
	end
	ml_time = toc;				% Stop ML timer

	% Extract indices for ML maps
	arg.T2_init = embed(T2dict(t2_idx), arg.mask);
	if (arg.disp)
		im(4, arg.T2_init, [0 300], 'T2 init'), cbar, drawnow
	end
else
	ml_time = NaN;
end

% Initialize T2 to T2_init
T2_ml = arg.T2_init;
T2 = T2_ml(arg.mask);

% Set M0_init from ML estimate, if necessary
if isempty(arg.M0_init)
	% Update M0 system matrix, given initial T2 guess
	tmp = cell(M, 1);
	for m = 1:M
		tmp{m} = Gdiag(df_se_M0(T2, flip, TE(m), arg.is_mag));
	end
	Am = vertcat(tmp{:});

	% Least-squares estimate
	arg.M0_init = embed(div_safe(((Am'*W) * col(y)), ((Am'*W*Am) * ...
		ones(size(Am,2), 1))), arg.mask);
	if (arg.disp)
		im(3, abs(arg.M0_init), 'M0 init'), cbar, drawnow
	end
end

% Initialize M0 to M0_init
M0_ml = arg.M0_init;
M0 = M0_ml(arg.mask);

% Construct regularizer objects for M0 and T2
% Potential functions are corner-rounded L1 (edge-preserving)
Rm = Reg1(arg.mask, 'pot_arg', {'hyper3', arg.delta_m},...
	'beta', arg.beta_m, 'type_penal', 'mat');
R2 = Reg1(arg.mask, 'pot_arg', {'hyper3', arg.delta_2},...
	'beta', arg.beta_2, 'type_penal', 'mat');

% Trick: normalize data by median of sum of non-background values
% so that the effective regularization beta is scale-invariant
if (arg.norm)
	tmp = sum(abs(y), 2);
	tmp = median(tmp(tmp > arg.thresh * max(tmp(:))));
	y = y / tmp;
end

% Cost function initialization
cost = NaN(2*arg.n_outer+1, 1);
cost_idx = 1;							% Counter
cost(cost_idx) = Rm.penal(Rm, M0) + R2.penal(R2, T2);
for m = 1:M
	cost(cost_idx) = cost(cost_idx) ...
		+ datacost(y(:,m), f_se(M0, T2, flip, TE(m), arg.is_mag));
end
cost_idx = cost_idx + 1;

%% Regularized LS via Alternating Minimization
tic;											% Start RLS timer
for outer = 1:arg.n_outer
	% Display outer iteration
	printm('Outer Iteration: %d of %d', outer, arg.n_outer);

	%% M0 update
	% If not the first iteration, update M0 system matrix
	if (outer ~= 1) || (~exist('Am'))
		tmp = cell(M, 1);
		for m = 1:M
			tmp{m} = Gdiag(df_se_M0(T2, flip, TE(m), arg.is_mag));
		end
		Am = vertcat(tmp{:});
	end

	% Update preconditioning matrix
	Dm = (Am' * Am) + Gdiag(Rm.denom(Rm, M0));
	Pm = Gdiag(1 ./ (Dm * ones(size(Dm,1), 1)));

	% Preconditioned conjugate gradient (equiv. to WLS for beta_m = 0)
	[M0, M0_info] = pwls_pcg1(M0, Am, W, col(y), Rm,...
		'niter', arg.n_innerM, 'precon', Pm,...
		'stop_diff_tol', arg.tolM, 'chat', arg.disp);

	% Cost function evaluation (after M0 update)
	cost(cost_idx) = Rm.penal(Rm, M0) + R2.penal(R2, T2);
	for m = 1:M
		cost(cost_idx) = cost(cost_idx) ...
			+ datacost(y(:,m), f_se(M0, T2, flip, TE(m), arg.is_mag));
	end
	cost_idx = cost_idx + 1;

	%% T2 update
	% After a few iterations, tighten bound of [T2min, T2max]
	if (outer > 5)
		if (min(T2) > arg.T2min)
			fprintf('T2min changed from %0.2f to %0.2f at iter %u.\n',...
				arg.T2min, min(T2), outer);
			arg.T2min = min(T2);
		end
		if (max(T2) < arg.T2max)
			fprintf('T2max changed from %0.2f to %0.2f at iter %u.\n',...
				arg.T2max, max(T2), outer);
			arg.T2max = max(T2);
		end
	end

	% Preconditioner update
	% Quick hack: evaluate over full T2 range and take max value
	T2eval = logspace(log10(arg.T2min), log10(arg.T2max), arg.n_dict);
	d2eval = zeros(N, arg.n_dict);
	for k = 1:arg.n_dict
		for m = 1:M
			% Sum derivatives of cost function at various TE values
			d2eval(:, k) = d2eval(:, k)...
				+ real(conj(f_se(M0, T2eval(k), flip, TE(m), arg.is_mag)...
				- y(:,m)) .* ddf_se_T2(M0, T2eval(k), flip, TE(m), arg.is_mag))...
				+ abs(df_se_T2(M0, T2eval(k), flip, TE(m), arg.is_mag)).^2;
		end
	end

	% Assign preconditioner using max value of Hessian over T2 test points
	d2 = R2.denom(R2, T2) + max(d2eval, [], 2);
	P2 = Gdiag(1 ./ d2);

	% Initialize Nesterov auxiliary variable and momentum parameter
	z2 = T2;
	t = 1;

	% PGD w/ Nesterov Momentum for T2 update
	for inner2 = 1:arg.n_inner2
		% Compute the next momentum parameter
		t_next = (1 + sqrt(1 + 4*t^2))/2;

		% Compute column gradient of cost function at z2
		grad = R2.cgrad(R2, z2);
		for m = 1:M
			grad = grad + real(conj(f_se(M0, z2, flip, TE(m), arg.is_mag) - y(:,m))...
				.* (arg.weights(m) * ones(N,1))...
				.* df_se_T2(M0, z2, flip, TE(m), arg.is_mag));
		end

		% PGD over z2 is stored in T2
		T2_prev = T2;
		T2 = z2 - (P2 * grad);

		% Project between [T2min, T2max]
		T2 = max(T2, arg.T2min);
		T2 = min(T2, arg.T2max);

		% Momentum over T2 is store in z2
		z2 = T2 + ((t-1)/t_next) * (T2 - T2_prev);

		% Update the momentum parameter
		t = t_next;

		% Display progress
		if (rem(inner2,5)==0)
			printm('T2 PGD: %d of %d', inner2, arg.n_inner2);
		end

		% Exit early if change in solution < tol2
		change = norm(T2 - T2_prev) / norm(T2);
		if (change < arg.tol2)
			printm('Exited at inner iteration %u of %u', inner2, arg.n_inner2);
			break;
		end
	end

	% Cost function evaluation (after T2 update)
	cost(cost_idx) = Rm.penal(Rm, M0) + R2.penal(R2, T2);
	for m = 1:M
		cost(cost_idx) = cost(cost_idx) ...
			+ datacost(y(:,m), f_se(M0, T2, flip, TE(m), arg.is_mag));
	end
	cost_idx = cost_idx + 1;

	% Observe images over outer iteration
	if (arg.disp)
		im(3, embed(abs(M0), arg.mask), '|M0|'), cbar
		im(4, embed(T2, arg.mask), [0 300], 'T2'), cbar, drawnow
	end

	%% If all inner updates exited after 1st iteration, exit outer iterations
	innerM = size(M0_info, 1);
	quit = (innerM <= 1) && (inner2 <= 2);
	if (quit)
		printm('Exited program at outer iteration %u of %u', outer, arg.n_outer);
		break;
	end
end
rls_time = toc;						% Stop RLS timer

% Extract M0 and T2 RLS estimates
M0_rls = embed(M0, arg.mask);
T2_rls = embed(T2, arg.mask);
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Helper functions
function cost_df = datacost(y, f)
	cost_df = (1/2) * norm(y - f).^2;
end

function quot = div_safe(dividend, divisor)
	undef = (dividend == 0) & (divisor == 0);
	quot = dividend ./ divisor;
	quot(undef) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin echo signal model
function s_se = f_se(M0, T2, a_ex, TE, is_mag)
	s_se = M0.*exp(-TE./T2).*sin(a_ex).*1i;
	if(is_mag), s_se = abs(s_se); end;
end

% Spin echo first derivative w.r.t. M0
function ds_se_M0 = df_se_M0(T2, a_ex, TE, is_mag)
	ds_se_M0 = exp(-TE./T2).*sin(a_ex).*1i;
	if(is_mag), ds_se_M0 = abs(ds_se_M0); end;
end

% Spin echo first derivative w.r.t. T2
function ds_se_T2 = df_se_T2(M0, T2, a_ex, TE, is_mag)
	ds_se_T2 = M0.*1.0./T2.^2.*TE.*exp(-TE./T2).*sin(a_ex).*1i;
	if(is_mag), ds_se_T2 = abs(ds_se_T2); end;
end

% Spin echo second derivative w.r.t. T2
function dds_se_T2 = ddf_se_T2(M0, T2, a_ex, TE, is_mag)
	dds_se_T2 = -M0.*1.0./T2.^4.*TE.*sin(a_ex).*(T2.*2.0i-TE.*1i);
	if(is_mag), dds_se_T2 = abs(dds_se_T2); end;
end
