  function [grad out com] = ir_hct_dgrad(x, A, yi, wi, varargin)
%|function [grad out com] = ir_hct_dgrad(x, A, yi, wi, varargin)
%|
%| compute gradient of WLS data-fit term using hct binary
%| for UM testing only
%|
%| in
%| option
%|	(many - see below)
%|
%| Jeff Fessler, 2012-06-17

if nargin == 1 && streq(x, 'test'), ir_hct_dgrad_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.clean = true; % remove files after?
arg.dir = ''; % work directory
arg.file_x = 'ir-hct-dgrad-x.fld'; % work file for x
arg.file_grad = 'ir-hct-dgrad-grad.fld'; % work file for grad
arg.file_mask = 'ir-hct-dgrad-mask.fld'; % work file for mask
arg.file_yi = 'ir-hct-dgrad-yi.fld'; % work file for yi
arg.file_wi = 'ir-hct-dgrad-wi.fld'; % work file for wi
arg.nthread = jf('ncore');
arg.chat = 0;

arg = vararg_pair(arg, varargin);
if isempty(arg.dir)
	arg.dir = test_dir;
end

pn = jf_protected_names;
if ~pn.has_hct2
	fail 'need hct2'
end

arg.file_x = [arg.dir arg.file_x];
arg.file_grad = [arg.dir arg.file_grad];
arg.file_mask = [arg.dir arg.file_mask];
arg.file_wi = [arg.dir arg.file_wi];
arg.file_yi = [arg.dir arg.file_yi];

% write files
%tox = @(z) permute(z, [2 3 1]);
fld_write(arg.file_mask, A.imask, 'type', 'byte')
if isempty(x)
	arg.file_x = '-';
else
	fld_write(arg.file_x, x)
end
if isempty(yi)
	arg.file_yi = '-';
else
	fld_write(arg.file_yi, yi)
end
if isempty(wi)
	arg.file_wi = '-';
else
	fld_write(arg.file_wi, wi)
end

if exist(arg.file_grad, 'file')
	clean = ['/bin/rm ' arg.file_grad]
	os_run(clean)
end

tmp = pn.hct_arg(A.cg, A.ig);

%{ old way:
com = sprintf(['hct2 what iter chat %d sysmod %s nthread %d %s ' ...
	'file_denom_in %s ' ... % trick: to avoid computing denom!
	'file_yi %s file_wi %s ' ...
	'file_mask %s file_init %s file_dgrad %s'], ...
	arg.chat, A.arg.type, arg.nthread, tmp, ...
	arg.file_x, ... % trick
	arg.file_yi, arg.file_wi, ...
	arg.file_mask, arg.file_x, arg.file_grad);
%}

% new way:
com = sprintf(['hct2 what wls_grad chat %d sysmod %s nthread %d %s ' ...
	'file_yi %s file_wi %s ' ...
	'file_mask %s file_init %s file_dgrad %s'], ...
	arg.chat, A.arg.type, arg.nthread, tmp, ...
	arg.file_yi, arg.file_wi, ...
	arg.file_mask, arg.file_x, arg.file_grad);

out = os_run(com); % run hct2

grad = fld_read(arg.file_grad);

if arg.clean
	clean = ['/bin/rm ' arg.file_x ' ' arg.file_grad ' ' arg.file_mask ' ' arg.file_yi ' ' arg.file_wi ' '];
	clean = strrep(clean, ' - ', ' ')
	os_run(clean)
end


% ir_hct_dgrad_test()
function ir_hct_dgrad_test
rng(0)
ig = image_geom('nx', 512, 'nz', 32, 'fov', 500, 'down', 4);
if 0
	tmp = ig.mask;
	tmp([1 end],:,:) = false;
	tmp(:,[1 end],:) = false;
	ig.mask = tmp; % remove 1-pixel border
else
	ig.mask = ig.circ > 0;
end
cg = ct_geom('ge1', 'nt', 16, 'down', 4);
%cg.plot3(ig)

A = Gcone(cg, ig, 'type', 'sf2');

if 1
	xb = ellipsoid_im(ig, ''); % [0 2]
	yb = A * xb; % [0 490]
	yi = yb + 9 * randn(size(yb));
	wi = exp(-yb * (3 / max(yb(:))));
	minmax(wi)
	%sino = @(y) permute(y, [1 3 2]);

	x0 = xb + 0.3 * ig.mask .* randn(ig.dim);
	W = Gdiag(wi);
	mat3 = A' * W * (A * x0 - yi);

	printm 'calling hct (3)'
	cpu etic
	[hct3 out3 com3] = ir_hct_dgrad(x0, A, yi, wi);
	cpu etoc hct2:dgrad:time
	printm(['hct2 output: \n\n' out3])
	%printm(com3)

	err3 = hct3 - mat3;
	clim = [];
	im plc 1 3
	im('row', 4, 1, mat3, clim)
	im('row', 4, 2, hct3, clim)
	im('row', 4, 3, err3)
	equivs(hct3, mat3, 'thresh', 5e-6)
	prompt
end

%% now with no file i/o
if 1
	x0_fake = 1000 * ig.mask; % default init_value in hct2
	wi_fake = 1.2 * cg.ones; % hard-coded values in wls-grad.c
	yi_fake = 4e5 * cg.ones; % 400mm * 1000 HU in wls-grad.c
	W_fake = Gdiag(wi_fake);
	mat0 = A' * W_fake * (A * x0_fake - yi_fake);
	if 0
		tmp = A * x0_fake - yi_fake;
		minmax(tmp)
		s2 = @(y) squeeze(y(:,2,:));
		im plc 2 2, im(1, s2(A*x0_fake)), im(2, s2(yi_fake))
		im(3, s2(tmp)), im(4, mat0)
		cg.plot(ig)
	return
	end

	printm 'calling hct (0)'
	cpu etic
	[hct0 out0 com0] = ir_hct_dgrad([], A, [], []);
	cpu etoc hct2:dgrad:time
	printm(['hct2 output: \n\n' out0])
	%printm(com)

	err0 = hct0 - mat0;
	clim = [];
	im plc 1 3
	im('row', 4, 1, mat0, clim)
	im('row', 4, 2, hct0, clim)
	im('row', 4, 3, err0)
	equivs(hct0, mat0, 'thresh', 2e-6)
end
