 function [C, B, spectrum, R] = ir_reg_diff_zeroed(idim, varargin)
%function [C, B, spectrum, R] = ir_reg_diff_zeroed(idim, varargin)
%|
%| Regularizing finite-differences R = diag(B)*C where C'C is circulant
%| and B is a binary "mask" array that eliminates undesired boundary diffs.
%| (Useful in variable splitting methods like ADMM.)
%|
%| in
%|	idim	[1 ndims]	dimensions of image (1D, 2D or 3D)
%|
%| option
%|	'order'		1|2	finite difference order (default: 2)
%|	'class'			'fatrix2' (default) or 'Fatrix' or 'sparse'
%|	'mask'			supported only for 'sparse' version (todo)
%|
%| out
%|	C	[M (*N)]	finite-differencer (fatrix2 or sparse)
%|				M = (*N) * K differences where
%|				K is # of directions (1 1D, 4 2D, 3 3D)
%|	B	[(N) K]		binary mask to remove edge wrapping
%|				also does 'tight' option of Rweights if mask
%|	spectrum [(N)]		eigenvalues (spectrum) of C'C
%|					C'*C*x = ifftn(spectrum .* fftn(x))
%|	R	[M (*N)]	R = B*C in fatrix2 or sparse form
%|
%| Michael Allison (original version ir_circ_zeroed_reg.m)
%| 2015-08-09 Jeff Fessler revised to do 1D,2D,3D with built-in test

if nargin < 1, help(mfilename), error(mfilename), end
if streq(idim, 'test'), ir_reg_diff_zeroed_test, return, end

arg.order = 2;
arg.mask = [];
arg.class = 'fatrix2';
arg = vararg_pair(arg, varargin);

if isempty(arg.mask), arg.mask = true([idim 1]); end

% sparse special case
if streq(arg.class, 'sparse')
	printm('Using sparse matrix to create R (no C,B,spectrum returned)')

	R = Reg1(true([idim 1]), 'order', arg.order, ...
		'distance_power', 2, ...
		'type_diff', 'spmat', 'type_penal', 'mat');
	R = R.C;

	% get approprite 'B'
	[~, B] = ir_reg_diff_zeroed(idim, 'mask', arg.mask, 'order', arg.order);
	R = spdiag(B(:), 'nowarn') * R;

	% to avoid confusion, null these
	C = [];
	B = [];
	spectrum = [];
return
end

switch arg.order
case 2
	[C B] = ir_reg_diff_zeroed_order2(idim, arg.mask, arg.class);
otherwise
	fail('only arg.order = 2 done')
end

if nargout < 3, return, end

% compute spectrum when requested
tmp = zeros(idim);
tmp(1) = 1; % unit vector e_0
tmp = C' * (C * tmp); % [(N)]
spectrum = fftn(tmp);
spectrum = reale(spectrum);
spectrum = max(spectrum,0); % not essential for order 2

if nargout == 4
	Bd = Gdiag(B, 'class', arg.class);
	R = Bd * C; % slower than B .* where B is a vector?
end

end % ir_reg_diff_zeroed()


% ir_reg_diff_zeroed_order2()
% 2nd-order finite differences
function [C B] = ir_reg_diff_zeroed_order2(idim, mask, classC)

gtype = 'imfilter,circ';
switch numel(idim)
case 1 % 1D
	offsets = [1];
	psfs = {[-1; 2; -1]};

case 2 % 2D

	nx = idim(1);
	offsets = [1 nx nx+1 nx-1];
	psfs = {[-1; 2; -1];
		[-1 2 -1];
		[-1 0 0; 0 2 0; 0 0 -1] / sqrt(2);
		[0 0 -1; 0 2 0; -1 0 0] / sqrt(2); % for distance_power 2
		};

case 3 % 3D

	nx = idim(1);
	ny = idim(2);
	offsets = [1 nx nx*ny];
	z3 = zeros(3,3,3);
	psfs = {z3, z3, z3};
	psfs{1}(:,2,2) = [-1 2 -1];
	psfs{2}(2,:,2) = [-1 2 -1];
	psfs{3}(2,2,:) = [-1 2 -1];
	gtype = 'conv,per';

otherwise
	fail('only 1D 2D 3D done')
end

Cs = cell(numel(psfs), 1);
for ip=1:numel(psfs)
	% trick: use true() for mask here so C'C is circulant
	% trick: [idim 1] needed for 1D case
	Cs{ip} = Gblur(true([idim 1]), 'psf', psfs{ip}, ...
		'type', gtype, 'class', classC);
end
C = vertcat(Cs{:});

% here we use the actual mask
Rw = Rweights(mask, offsets, 'type_wt', 'array', ...
	'edge_type', 'tight', 'order', 2, 'distance_power', 0);

B = reshape(Rw, [idim numel(psfs)]);

end % ir_reg_diff_zeroed_order2_2d()


% ir_reg_diff_zeroed_test1()
% 1D test
function ir_reg_diff_zeroed_test1

idim = [6];
[C B spectrum R2] = ir_reg_diff_zeroed(idim);
[~, ~, ~, Rs] = ir_reg_diff_zeroed(idim, 'class', 'sparse');
im plc 2 2
im(1, full(C)', 'C 1D')
im(2, diag(B), 'B')
im(3, full(Rs)', 'R sparse')
im(4, full(R2)', 'R fatrix')
jf_equal(full(R2), full(Rs))

end % ir_reg_diff_zeroed_test1()


% ir_reg_diff_zeroed_test2()
% 2D test
function ir_reg_diff_zeroed_test2

idim = [6 5];
mask = true(idim); mask(1:2) = 0; % stress test mask
[C B spectrum R2] = ir_reg_diff_zeroed(idim, 'mask', mask);
[~, ~, ~, Rs] = ir_reg_diff_zeroed(idim, 'mask', mask, 'class', 'sparse');
R2f = full(R2);
Rsf = full(Rs);
rs = @(x) reshape(x, [prod(idim)*[1 1] 4]);
im plc 1 4
im('col', 1, 1, rs(full(C)'), 'C 2D')
im('col', 1, 2, rs(R2f'), 'R fatrix')
im('col', 1, 3, rs(Rsf'), 'R sparse')
%im('col', 1, 4, rs(Rsf' - R2f'))
im('col', 1, 4, B)
%jf_equal(R2f, Rsf) % fails due to sqrt(2)
equivs(R2f, Rsf)

if 0 % compare to old version
	[Co Bo Po Ro] = ir_circ_zeroed_reg(2, true(idim), 1);
%	im('col', 1, 5, rs(Ro'))
	jf_equal(Ro, Rs)
	[Co Bo Po Ro] = ir_circ_zeroed_reg(2, true(idim), 0);
	jf_equal(full(Ro), Rs) % ok
	jf_equal(full(Ro), full(Bo * Co)) % ok
end

end % ir_reg_diff_zeroed_test2()


% ir_reg_diff_zeroed_test3()
% 3D test
function ir_reg_diff_zeroed_test3

idim = [6 5 4];
[C B spectrum R2] = ir_reg_diff_zeroed(idim);
[~, ~, ~, Rs] = ir_reg_diff_zeroed(idim, 'class', 'sparse');
R2f = full(R2);
Rsf = full(Rs);
rs = @(x) reshape(x, [prod(idim)*[1 1] 3]);
im plc 1 4
im('col', 1, 1, rs(full(C)'), 'C 3D')
im('col', 1, 2, rs(R2f'), 'R fatrix')
im('col', 1, 3, rs(Rsf'), 'R sparse')
im('row', 3, 4, B)
%im('col', 1, 5, rs(Rsf' - R2f'))
jf_equal(R2f, Rsf)

if 0 % compare to old version
	[Co Bo Po Ro] = ir_circ_zeroed_reg_3d(2, true(idim), 1);
	jf_equal(Ro, Rs)
	[Co Bo Po Ro] = ir_circ_zeroed_reg_3d(2, true(idim), 0);
	jf_equal(full(Ro), Rs) % ok
	jf_equal(full(Ro), full(Gdiag(Bo) * Co)) % ok
end

end % ir_reg_diff_zeroed_test3()


% ir_reg_diff_zeroed_test()
function ir_reg_diff_zeroed_test
ir_reg_diff_zeroed_test1
prompt
ir_reg_diff_zeroed_test2
prompt
ir_reg_diff_zeroed_test3

end % ir_reg_diff_zeroed_test()
