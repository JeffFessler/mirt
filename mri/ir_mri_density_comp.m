 function wi = ir_mri_density_comp(kspace, dtype, varargin)
%function wi = ir_mri_density_comp(kspace, dtype, varargin)
%|
%| todo: THIS NEEDS A LOT OF WORK!
%|
%| Compute density compensation factors for the conjugate phase
%| method for image reconstruction from Fourier samples.
%|
%| in
%|	kspace	[M 1]	kspace sample locations, e.g., spiral
%|	dtype	char	which density compensation method (see below)
%|			'voronoi', 'jackson', 'pipe', 'qian'
%| options
%|	G	?
%|	fix_edge 0|1|2	for voronoi, (default: 2 - 2nd-order poly extrapolation)
%|	arg_pipe {}	options for Pipe&Menon method
%|
%| out
%|	wi	[M 1]	density compensation factors
%|
%| If voronoi, then "redundant" sampling at DC is corrected.
%| (But not if there are redundant samples at other locations in k-space.)
%|
%| Copyright 2003-7-29, Jeff Fessler, The University of Michigan
%| 2009-12-18, modified by Greg Lee to support pipe and jackson with table

% todo: src/matlab/alg/contrib/khalsa/Jop.m etc.

if ~nargin, ir_usage, end

if streq(kspace, 'test')
	if nargin == 1 || isempty(dtype)
%		dtype = 'jackson';
		dtype = 'pipe';
%		dtype = 'voronoi';
	end
	ir_mri_density_comp_test(dtype)
return
end

arg.G = [];
arg.fix_edge = 2;
arg.arg_pipe = {};
arg = vararg_pair(arg, varargin);

if nargin < 2, ir_usage, end

switch dtype
case 'voronoi'
	wi = ir_mri_dcf_voronoi0(kspace, arg.fix_edge);
case 'jackson'
	wi = ir_mri_dcf_jackson(kspace, arg.G);
case 'pipe'
	wi = ir_mri_dcf_pipe(kspace, arg.G, arg.arg_pipe{:});
%case 'qian'
%	wi = ir_mri_dcf_qian(kspace, arg.G);
otherwise
	fail('unknown DCF type %s', dtype)
end


% ir_mri_dcf_voronoi0()
% in radial imaging, k-space origin is sampled multiple times, and
% this non-uniqueness messes up matlab's voronoi routine.
% here we find those "redundant" zeros and remove all but one them
% for the voronoi call.  we then restore them with appropriate DCF.
%
function wi = ir_mri_dcf_voronoi0(kspace, fix_edge)
M = size(kspace, 1);
i0 = sum(abs(kspace), 2) == 0; % which points are at origin?
if sum(i0) > 1 % multiple DC points?
	i0f = find(i0);
	i0f = i0f(1); % keep the first zero point only
	i0(i0f) = false; % trick
	wi = zeros(M, 1);
	wi(~i0) = ir_mri_dcf_voronoi(kspace(~i0,:), fix_edge);
	i0(i0f) = true; % trick
	wi(i0) = wi(i0f) / sum(i0); % distribute dcf equally
else
	wi = ir_mri_dcf_voronoi(kspace, fix_edge);
end


% ir_mri_dcf_voronoi()
%
function wi = ir_mri_dcf_voronoi(kspace, fix_edge)
M = size(kspace, 1);

wi = zeros(M,1);
[v, c] = voronoin(double(kspace));
nbad = 0;
for mm=1:M
	ticker([mfilename ' (voronoi)'], mm, M)
	x = v(c{mm},:);
	if ~any(isinf(x))
		try
			[~, wi(mm)] = convhulln(x); % cell area
		catch
%			printm('bad %d', mm)
			nbad = nbad + 1;
		end
	end
end
if nbad
	printm('bad edge points %d of %d', nbad, M)
end

% points at the outer edges of k-space have infinite voronoi cell area
% so are assigned wi=0 above.  to improve on 0, here we extrapolate
% based on the points near the edge.
switch fix_edge
case 2
	rho = sum(kspace.^2, 2); % radial frequency coordinate
	igood = (rho > 0.6 * max(rho)) & (wi > 0);
	pp = polyfit(rho(igood), wi(igood), 2);
	wi(wi == 0) = polyval(pp, rho(wi == 0)); % extrapolate

% old way: look for points close to convex hull and use max of other points?
case 1
	printm('trying to fix %d zeros of %d', sum(wi==0), M)
	ii = logical(zeros(size(wi)));
	fac = 0.98;
	for id=1:ncol(kspace) % find cartesian edges of k-space
		k = kspace(:,id);
		ii = ii | (k > fac * max(k)) | (k < fac * min(k));
	end
	if ncol(kspace) >= 2
		k = sqrt(kspace(:,1).^2 + kspace(:,2).^2);
		ii = ii | (k > fac * max(k)); % cylindrical edge
	end
	if ncol(kspace) >= 3
		k = sqrt(kspace(:,1).^2 + kspace(:,2).^2 + kspace(:,3).^2);
		ii = ii | (k > fac * max(k)); % spherical edge
	end

	pn = jf_protected_names;
	wmax = 2 * pn.prctile(wi(~ii), 95); % fix: this is not working well
	wi = min(wi, wmax);
	wi(wi==0) = max(wi);

otherwise
	if ~isequal(fix_edge, 0), error('bad fix_edge argument'), end
end


% ir_mri_dcf_jackson()
%
function wi = ir_mri_dcf_jackson(kspace, G)
M = size(kspace, 1);

% todo: this is not *really* Jackson's method!  need to work on it!
if streq(G.arg.st.alpha, 'kaiser')
	kb_m = G.arg.st.kb_m(1);
	kb_alf = G.arg.st.kb_alf(1);
%else
	% ??
end

if isfield(G.arg.st, 'interp_table')
	tmp = G.arg.st.interp_table_adj(G.arg.st, ones(M,1));
	tmp = G.arg.st.interp_table(G.arg.st, tmp);
	wi = reale(1 ./ tmp, 'warn');
else
	P = G.arg.st.p;
	wi = reale(1 ./ (P * (P' * ones(M,1))), 'warn');
	% wi = reale(1 ./ (G * (G' * ones(M,1))), 'warn');
	% wi = w * G.arg.st.sn(end/2,end/2)^(-2);
	%	/ fov^2 / prod(G.arg.st.Kd) * N0^2;
end
    

% ir_mri_dcf_pipe()
% Pipe&Menon, based on equalty 1 = A A' w so w = w ./ (A A' w)
%
function wi = ir_mri_dcf_pipe(kspace, G, varargin)
arg.isave = [];
arg.niter = 20;
arg.thresh = 0.02;
arg.fov = 1;
arg.unitv = [];
arg.wi = []; % see below
arg = vararg_pair(arg, varargin);

P = G.arg.st.p;
goal = inf;
iter = 0;
%saver = zeros(arg.niter,1);

%{
removed 2017-04-11 in favor
scale = G.arg.st.sn(end/2,end/2)^(-2) / arg.fov^2 ...
	/ prod(G.arg.st.Kd) * prod(G.arg.st.Nd);
scale = reale(scale);
%}
if isempty(arg.wi)
	wi = ones(size(kspace,1), 1); % default initial
else
%	minmax(arg.wi)
%	wi = arg.wi / scale; % todo: why?
	wi = arg.wi;
end
wi_save = wi;
while(max(abs(goal-1)) > arg.thresh)
	iter = iter + 1;
	if isfield(G.arg.st,'interp_table')
		goal = G.arg.st.interp_table(G.arg.st, ...
			G.arg.st.interp_table_adj(G.arg.st, wi) );
	else
		goal = P * (P' * wi); % warn: complex results!?
	end
	wi = wi ./ real(goal);
	%wi = wi ./ abs(goal);
	if iter > arg.niter
		warn 'iteration stuck?'
		break
	end
%	saver(iter) = max(abs(goal-1));
	if ~isempty(arg.isave)
		wi_save = [wi_save, wi];
	end
end
printm('pipe ended at iteration %d with %g', iter, max(abs(goal-1)))
%plot(saver(2:end))

% scale
if isempty(arg.unitv)
	warn 'no unitv so no scaling - user must fix scaling'
else
	fail 'todo: this needs tested1'
	last = wi_save(:,1);
	psf_raw = G' * last;
	e0 = ig.unitv % todo
	scale = dot(e0, psf_raw) / norm(psf_raw(:)).^2;
	wi_save = scale * wi_save;
end

if ~isempty(arg.isave)
	wi = wi_save;
end

%wi = wi * scale; % todo: why?


% ir_mri_density_comp_test()
% self-test routine
function ir_mri_density_comp_test(dtype)

ig = image_geom_mri('nx', 2^5, 'fov', 256); % typical brain FOV
fov = ig.fov(1);
N0 = ig.nx;

if 0
	t = linspace(0, N0/2*2*pi, N0^2+3)'; % crude spiral:
	kspace = N0/2*(1/fov)*[cos(t) sin(t)] .* (t(:,[1 1]) / max(t));
else
	ktype = 'cartesian';
	ktype = 'radial';
	[kspace, om, wi_r] = mri_trajectory(ktype, {}, ig.dim, ig.fov);
end

im plc 2 3
if im
	im subplot 1
	plot(kspace(:,1), kspace(:,2), '.')
	axis(1.1*[-1 1 -1 1]*N0/2/fov), axis square
	xlabel 'k_1 [mm^{-1}]', ylabel 'k_2 [mm^{-1}]'
	titlef('%d k-space samples', size(kspace,1))
end

% create Gnufft object
omega = 2*pi*kspace*fov/N0;
if streq(dtype, 'pipe')
	N = ig.dim;
	G = Gnufft(ig.mask, {om, N, [6 6], 2*N, N/2, 'table', 2^10, 'minmax:kb'});
else
	%G = Gnufft({omega, ig.dim, [6 6], 2*ig.dim, 1*ig.dim/2, 'kaiser'});
	G = Gdsft(omega, ig.dim, 'n_shift', 1*ig.dim/2, 'mask', ig.mask);
end

% true object and analytical k-space data
obj = mri_objects('fov', ig.fov, 'rect2half');
xtrue = obj.image(ig.xg, ig.yg);
ytrue = obj.kspace(kspace(:,1), kspace(:,2));

if 0 % check forward direction (works, after fixing offsets with image_geom_mri)
	tmp = G * xtrue(ig.mask);
	tmp = tmp * abs(ig.dx * ig.dy);
	im plc 2 1
	im subplot 1
	plot([real(tmp) real(ytrue)]), axis tight
	im subplot 2
	plot([imag(tmp) imag(ytrue)]), axis tight
return
end

if 0 % testing
	P = G.arg.st.p;
	P = P.arg.G; % trick: because st.p is now a Gsparse
	tmp = P * P(:,:)';
	%tmp = conj(P) * P.';
	printm('minmax of real,imag of PP'' are:')
	disp(full(minmax(real(tmp(:))))')
	disp(full(minmax(imag(tmp(:))))')
end
if isempty(dtype)
	keyboard
end

switch dtype
case 'jackson'
	wi = ir_mri_density_comp(kspace, dtype, 'G', G);
case 'pipe'
	wi = ir_mri_density_comp(kspace, dtype, 'G', G, 'arg_pipe', ...
		{'fov', ig.fov});
otherwise
	wi = ir_mri_density_comp(kspace, dtype, 'fix_edge', 2);
end

rho = sqrt(sum(kspace.^2,2));
if 0 % compare voronoi vs analytical for radial trajectory
	im clf
	slope = (wi_r'* rho) / norm(rho)^2;
	slope = 0;
	plot(rho, wi - slope*rho, 'yo', rho, wi_r - slope*rho, 'c+'), axis tight
return
end

% CP recon
%wi = wi_r;
if im
	im subplot 2
	semilogy(rho, wi, '.'), titlef('DCF %s', dtype), axis tight
	xlabel '|k|', ylabel 'wi'
end

clim = [0 2];
im(4, ig.x, ig.y, xtrue, 'x true', clim), cbar

xcp = ig.embed(G' * (wi .* ytrue));
im(5, ig.x, ig.y, real(xcp), 'conj. phase'), cbar

%sum(xcp(:)) / sum(xtrue(:))
xlabelf('nrms %g\%%', 100*nrms(xcp(:), xtrue(:)))

if im
	subplot(133)
	plot(ig.x, xtrue(:,end/2), 'c-', ig.y, real(xcp(:,end/2)), 'y.-')
	axis tight
	legend('true', 'CP', 'location', 'east')
end
