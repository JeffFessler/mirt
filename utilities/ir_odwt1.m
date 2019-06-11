 function [coef codes] = ir_odwt1(x, varargin)
%function [coef codes] = ir_odwt1(x, varargin)
%|
%| 1D orthonormal discrete wavelet transform (ODWT) of (each column of) x
%|
%| in
%|	x	[N L]	L signals of length N
%|
%| option
%|	'level'	int	default: 1
%|	'wname'	char	default: 'haar'
%|	'ortho'	0|1	default: 1
%|	'adj'	0|1	default: 0 (if 1 then adjoint, which is inverse)
%|	'abs'	0|1	default: 0 (if 1 then take abs of filters)
%|
%| out
%|	coef	[N L]	L sets of 1D wavelet coefficients
%|	codes	[N 1]	uint2 level "codes" (only if not adj)
%|			code = 10 * level + (0 approx, 1 detail)
%|
%| 2012-05-17, Jeff Fessler, Univ. of Michigan

if nargin < 1, ir_usage, end
if streq(x, 'test'), ir_odwt1_test, return, end

arg.level = 1;
arg.wname = 'haar';
arg.ortho = true;
arg.adj = false;
arg.abs = false;
arg.usemat = false;
arg = vararg_pair(arg, varargin);

if arg.level == 0
	coef = x;
	codes = zeros(size(x,1), 'uint16');
return
end

% decomposition filters
[lo hi] = ir_dwt_filters(arg.wname, ...
	'abs', arg.abs, 'ortho', arg.ortho, 'usemat', arg.usemat);

if arg.adj % synthesis (reconstruction)
	coef = ir_odwt1_back(x, lo, hi, arg.level); % trick
else % analysis (decomposition)
	[coef codes] = ir_odwt1_forw(x, lo, hi, arg.level);
end


% ir_odwt1_forw()
function [coef codes] = ir_odwt1_forw(x, lo, hi, level)
N = size(x,1);
[coef codes] = ir_odwt1_do(x, lo, hi, []);
for ll=2:level
	ii = 1:(N / 2^(ll-1));
	[coef(ii,:) codes(ii)] = ir_odwt1_do(coef(ii,:), lo, hi, codes(ii));
end


% ir_odwt1_back()
function x = ir_odwt1_back(coef, lo, hi, level)
lo_r = flipud(lo);
hi_r = flipud(hi);
N = size(coef, 1);
for ll=level:-1:2
	i1 = 1:(N / 2^(ll-1));
	coef(i1,:) = ir_odwt1_adj(coef(i1,:), lo_r, hi_r);
end
x = ir_odwt1_adj(coef, lo_r, hi_r);


% ir_odwt1_do()
function [coef codes] = ir_odwt1_do(x, lo, hi, codes)

N = size(x,1);
if mod(N,2) % odd
	fail('N=%d not even', N)
end

if isempty(codes)
	codes = zeros(N, 1, 'uint16');
end
codes(1:N/2) = codes(1:N/2) + 10;
codes((N/2+1):end) = codes((N/2+1):end) + 11;

coef_lo = ir_conv(x, lo, 'per', true);
coef_hi = ir_conv(x, hi, 'per', true);
i1 = 1:2:N; % [N/2] 2015-04-12 change 2:2 to 1:2 after ir_conv change
coef_lo = coef_lo(i1,:); % [N/2 L]
coef_hi = coef_hi(i1,:); % [N/2 L]
coef = [coef_lo; coef_hi]; % [N L]


% ir_odwt1_adj()
% synthesis
function x = ir_odwt1_adj(coef, lo, hi)

N = size(coef,1);
if mod(N,2) % odd
	fail('N=%d not even', N)
end

tmp = zeros(size(coef), class(coef));
coef_lo = tmp;
coef_hi = tmp;
coef_lo(1:2:end,:) = coef(1:end/2,:); % 2015-04-12 changed 2:2 to 1:2
coef_hi(1:2:end,:) = coef((1+end/2):end,:); % 2015-04-12 ""

% trick: flipud back because of 'adj'
x = ir_conv(coef_lo, flipud(lo), 'per', true, 'adj', 1) ...
	+ ir_conv(coef_hi, flipud(hi), 'per', true, 'adj', 1); % synthesis


% ir_odwt1_test_wname()
function ir_odwt1_test_wname(wname)

warg = {'wname', wname};
has_wave = exist('wfilters', 'file') == 2;

x = eye(2^3);
coef = ir_odwt1(x, warg{:});

if has_wave
	tmp = ir_odwt1(x, warg{:}, 'usemat', 1);
	equivs(tmp, coef)
end

if has_wave
	dwtmode('per');
	mat1 = zeros(size(x));
	for ll=1:ncol(x)
		[ca cd] = dwt(x(:,ll), wname);
		mat1(:,ll) = [ca; cd];
	end
	equivs(mat1, coef, 'format', ' %6.4f')
%	im plc 1 3, im(1, coef), im(2, mat1), xlabel(wname)
end

y = ir_odwt1(coef, 'adj', 1, warg{:});

if has_wave
	tmp = ir_odwt1(coef, warg{:}, 'adj', 1, 'usemat', 1);
	equivs(tmp, y)
end

if has_wave
	mat2 = zeros(size(x));
	for ll=1:ncol(x)
		ca = coef(1:end/2,ll);
		cd = coef((end/2+1):end,ll);
		mat2(:,ll) = idwt(ca, cd, wname);
	end
	equivs(mat2, y, 'format', ' %6.4f')
%	im plc 1 3, im(1, mat2), im(2, y), xlabel(wname)
end

equivs(y, x, 'format', ' %+7.4f')

for level=1:3
	[coef codes] = ir_odwt1(x, 'level', level, 'ortho', 0, 'abs', 0, ...
		warg{:}, 'usemat', 0);
	if im
		pr level
		disp(num2str(coef', ' %2.2f'))
		disp(num2str(codes', ' %5.1f'))
	end
end

level = 3;
coef = ir_odwt1(x, 'level', level, warg{:});
y = ir_odwt1(coef, 'level', level, warg{:}, 'adj', 1);
equivs(y, x)

if has_wave
	tmp = ir_odwt1(x, 'level', level, warg{:}, 'usemat', 1);
	equivs(coef, tmp)

	tmp = ir_odwt1(coef, 'level', level, warg{:}, 'adj', 1, 'usemat', 1);
	equivs(y, tmp)

	mat2 = zeros(size(coef));
	for ll=1:ncol(x)
		mat2(:,ll) = wavedec(x(:,ll), level, wname);
	end
	jf2 = ir_odwt1(x, 'level', level, warg{:});
	equivs(mat2, jf2)
end


% ir_odwt1_test()
function ir_odwt1_test
list = {'haar', 'sym2'};
for ii=1:numel(list)
	wname = list{ii};
	pr wname
	ir_odwt1_test_wname(wname)
end
