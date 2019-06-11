 function ir_proj3_compare1(fun_proj, fun_im, varargin)
%function ir_proj3_compare1(fun_proj, fun_im, varargin)
%|
%| Compute a set of 2d line-integral projection views of one or more ellipsoids.
%| Works for both parallel-beam and cone-beam geometry.
%|
%| in
%|	fun_proj	@(cg)	make projections given ct_geom()
%|	fun_im		@(ig)	make 3d image given image_geom()
%|
%| options
%|	(see code below)
%|
%| out
%|
%| Copyright 2013-06-10, Jeff Fessler, University of Michigan

if nargin == 1 && streq(fun_proj, 'test'), ir_proj3_compare1_test, return, end
if nargin < 2, ir_usage, end

arg.downp = 30;
arg.downi = 4;
arg.dsd = [949 949 inf]; % cone-arc, cone-flat, parallel-beam
arg.dod = 408;
arg.dfs = [0 inf 0]; % must have same length as arg.dsd
arg.pitch = [0.5 0];
arg.nt = 64;
arg.na = 12;
arg.orbit_start = 0;
arg.orbit_start = 355; % todo
arg.source_z0 = -200; % stress test helix
arg.arg_image_geom = {};
arg.arg_ct_geom = {}; % caller can override all parameters with this
arg.chat = false;
arg = vararg_pair(arg, varargin);

if numel(arg.dsd) ~= numel(arg.dfs)
	fail('dsd and dfs must have same numel')
end

for ii=1:numel(arg.dsd)
	dsd = arg.dsd(ii);
	dfs = arg.dfs(ii);
	for pitch = arg.pitch
		if arg.chat
			pr '[dsd dfs pitch]'
		end
		ir_proj3_compare1_do(fun_proj, fun_im, ...
			arg.arg_image_geom, arg.arg_ct_geom, ...
			arg.downp, arg.downi, ...
			'orbit_start', arg.orbit_start, ...
			'na', arg.na, ...
			'nt', arg.nt, ...
			'pitch', pitch, 'source_z0', arg.source_z0, ...
			'dsd', dsd, 'dod', arg.dod, 'dfs', dfs);
	end
end

end % ir_proj3_compare1_test


% ir_proj3_compare1_do()
function ir_proj3_compare1_do(fun_proj, fun_im, ...
	arg_image_geom, arg_ct_geom, downp, downi, varargin)
cg = ct_geom('fan', 'ns', round(888/downp), ...
	'ds', 1.0*downp, 'dt', 1.1*downp, ...
	'down', 1, ... % only downsample s,t, beta
	'offset_s', 0.25, ... % quarter detector
	'offset_t', 0.0, ...
	varargin{:}, arg_ct_geom{:});

pa = fun_proj(cg); % analytical

ig = image_geom('nx', 512/downi, 'nz', 256/downi, ...
	'offset_x', -9, 'offset_y', 7, 'offset_z', -5, ... % stress test
	'dx', 1*downi, 'dz', 1*downi, 'down', 1, arg_image_geom{:});
x = fun_im(ig);

% clf; cg.plot3(ig); prompt; return % check geometry visually

if 1
	im plc 2 2
	im(x)
	if any(isnan(pa(:)))
		keyboard
	end
	im(pa), cbar
end

% compare with Gcone
if 1 && has_mex_jf, % printm 'compare to Gcone discrete projector'
	A = Gcone(cg, ig, 'type', 'sf1', 'use_hct2', 0);

	pc = A * x;

	if im
		im(pc), cbar, xlabel(A.arg.type)
		im(pc - pa, 'pc - pa'), cbar
	end
	max_percent_diff(pa, pc)

	if im, prompt, end
end

% compare with aspire
if 0 && has_mex_jf && isinf(cg.dfs) && ~isinf(cg.dsd) % only flat cone-beam
	printm 'compare to 3l@ discrete projector'
	systype = aspire_pair(cg, ig, 'system', '3l');
	A = Gtomo3(systype, ig.mask, ig.nx, ig.ny, ig.nz, ...
		'chat', 0, 'checkmask', false, 'permute213', true);

	p3 = A * x;

	if im
		im plc 2 2
		im(1, x)
		im(2, pa), cbar
		im(3, p3), cbar
		im(4, p3 - pa, 'p3 - pa'), cbar
	end

	max_percent_diff(pa, p3)

	if im, prompt, end

	if im && 0
		im_toggle(pa, p3, [0 max(pa(:))])
		clf, im(p3-pa)
		prompt
%		im clf, movie2(pa)
%	prompt
	end
end

end % ir_proj3_compare1_do()


% ir_proj3_compare1_test
function ir_proj3_compare1_test

ell = [20 0*50 -40 200 100 50 90 0 10;
	0 50 100 80 80 20 0 0 10];

fun_proj = @(cg) ellipsoid_proj(cg, ell, 'oversample', 2); % analytical
fun_im = @(ig) ellipsoid_im(ig, ell, 'oversample', 2, 'checkfov', true);

ir_proj3_compare1(fun_proj, fun_im, 'chat', 1);

end % ir_proj3_compare1_test()
