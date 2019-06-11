  function [grad, out, com] = ir_hct_rgrad(x, varargin)
%|function [grad, out, com] = ir_hct_rgrad(x, varargin)
%|
%| compute 3D regularizer gradient using hct binary
%| for UM testing only
%|
%| in
%|	x	in xyz order
%| option
%|	(many - see below)
%|
%| Notes: because hct2 uses 'zxy' order, it is essential to use Reg1
%| with zxy order also, to ensure gradient matches over entire object.
%| If Reg1 uses the usual xyz order, then slices [1 end] do not match.
%|
%| Jeff Fessler, 2012-06-15

if nargin == 1 && streq(x, 'test'), ir_hct_rgrad_test, return, end
if nargin < 2, ir_usage(), end

arg.R = [];
arg.bij_center = '';
arg.bij_above = '';
arg.log2reg_x = 0; % 0 by default because beta built into R.beta
arg.log2reg_z = []; % ""
arg.kappa_min = 0;
arg.clean = true; % remove files after?
arg.dir = ''; % work directory
arg.file_x = 'ir-hct-rgrad-x.fld'; % work file for x
arg.file_grad = 'ir-hct-rgrad-grad.fld'; % work file for grad
arg.file_mask = 'ir-hct-rgrad-mask.fld'; % work file for mask
arg.file_kappa = 'ir-hct-rgrad-kappa.fld'; % work file for kappa
arg.chat = 0;

arg = vararg_pair(arg, varargin);
if isempty(arg.dir)
	arg.dir = test_dir;
end

R = arg.R;

if isempty(arg.log2reg_z)
	arg.log2reg_z = arg.log2reg_x;
end

arg.file_x = [arg.dir arg.file_x];
arg.file_grad = [arg.dir arg.file_grad];
arg.file_mask = [arg.dir arg.file_mask];
arg.file_kappa = [arg.dir arg.file_kappa];

pn = jf_protected_names;
if ~pn.has_hct2
	fail 'need hct2'
end

if ~streq(R.type_penal, 'zxy') || ~R.offsets_is_zxy
	fail 'zxy required to ensure match'
end

% write files
tox = @(z) permute(z, [2 3 1]);
fld_write(arg.file_x, x)
fld_write(arg.file_mask, tox(R.mask), 'type', 'byte')
fld_write(arg.file_kappa, tox(R.wt.kappa))

if exist(arg.file_grad, 'file')
	eval(['!/bin/rm ' arg.file_grad])
end

if ~isempty(R.pot_params)
	delta = R.pot_params(1);
	qgg_q = R.pot_params(2);
else
	delta = 0;
	qgg_q = 0;
end

if isempty(arg.bij_center)
	if isfield(R, 'beta') && numel(R.beta) > 1
		beta = R.beta;
		displace_zxy = penalty_displace(R.offsets, R.dim);
		oz = displace_zxy(:,1);
		ox = displace_zxy(:,2);
		oy = displace_zxy(:,3);
		bij = zeros(3,3,3);
		ii = sub2ind([3 3 3], 2+ox, 2+oy, 2+oz);
		bij(ii) = beta;
		ii = sub2ind([3 3 3], 2-ox, 2-oy, 2-oz);
		bij(ii) = beta;
		bij = bij / 2^arg.log2reg_x; % scale to match ge style in hct2
		if arg.log2reg_x ~= arg.log2reg_z, fail 'not done', end
		jf_equal(bij(:,:,1), bij(:,:,3)) % verify symmetry
		b = @(i,j) bij(1+i, 1+j, 1); % above
		arg.bij_above = sprintf('%g,%g,%g,%g,%g,%g,%g,%g,%g', ...
                        b(0,0), b(0,1), b(0,2), ...
                        b(1,0), b(1,1), b(1,2), ...
                        b(2,0), b(2,1), b(2,2));
		b = @(i,j) bij(1+i, 1+j, 2); % center
		arg.bij_center = sprintf('%g,%g,%g,%g,%g,%g,%g,%g,%g', ...
                        b(0,0), b(0,1), b(0,2), ...
                        b(1,0), b(1,1), b(1,2), ...
                        b(2,0), b(2,1), b(2,2));
	else
		arg.bij_center = '-';
		arg.bij_above = '-';
	end
end

% trick due to zxy order:
nx = R.dim(2);
ny = R.dim(3);
nz = R.dim(1);

com = sprintf(['hct2 what rgrad chat %d ' ...
	'ns 8 nt 16 na 12 ', ... % fake
	'nx %d ny %d nz %d ' ...
	'bij_center %s bij_above %s ' ...
	'distance_power %g ', ...
	'pot %s delta %g qgg_q %g use_reg_zxy 1 ' ...
	'log2reg_x %.6f log2reg_z %.6f kappa_min %g ' ...
	'file_kappa %s file_mask %s file_init %s file_xh %s'], ...
	...
	arg.chat, nx, ny, nz, ...
	arg.bij_center, arg.bij_above, ...
	R.distance_power, ...
	R.pot_type, delta, qgg_q, ...
	arg.log2reg_x, arg.log2reg_z, arg.kappa_min, ...
	arg.file_kappa, arg.file_mask, arg.file_x, arg.file_grad);

out = os_run(com); % run hct2

grad = fld_read(arg.file_grad);

if arg.clean
	clean = ['!/bin/rm ' arg.file_x ' ' arg.file_grad ' ' arg.file_mask ' ' arg.file_kappa];
	eval(clean)
end


% ir_hct_rgrad_test()
function ir_hct_rgrad_test

rng(0)
kappa = 2 + rand(18,20,10);
%kappa = ones(18,20,10);
kappa([1 end],:,:) = 0; % zero border required for zxy
kappa(:,[1 end],:) = 0;

toz = @(x) permute(x, [3 1 2]);
tox = @(z) permute(z, [2 3 1]);
R = Reg1(toz(kappa), 'beta', 4, 'nthread', jf('ncore'), ...
	'type_penal', 'zxy', ...
	'offsets', '3d:26', ...
	'type_wt', 'fly', ...
	'pot_arg', {'qgg2', 10, 1.2}, ...
	'distance_power', 0);
Rx = Reg1(kappa, 'beta', 4, 'nthread', jf('ncore'), ...
	'offsets', '3d:26', ...
	'type_wt', 'fly', ...
	'pot_arg', {'qgg2', 10, 1.2}, ...
	'distance_power', 0);
%x = zeros(size(kappa));
%x(end/2,end/2,end/2) = 0;
%x(end/2,end/4,1) = 1;
rng(0)
x = randn(size(kappa)) .* (kappa ~= 0);
if 0 % kludge to make penalty value match
	x(:,:,1) = 0;
	x(:,:,end) = 0;
end
mat = R.cgrad(R, toz(x));
mat = tox(mat);

xlim = [-1 1] * 3;
clim = [-1 1] * 1000;
im plc 2 2

im(1, x, xlim)
im(mat, clim)

printm 'calling hct'
%tic
[hct, out] = ir_hct_rgrad(x, 'R', R, 'log2reg_x', 4, 'chat', 99);
%toc
printm(['hct2 output: \n ' out])
pr Rx.penal(Rx, x);

im(hct, clim)
im(hct - mat)
equivs(hct, mat, 'thresh', 2e-6)
