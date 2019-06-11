 function [coef codes] = ir_odwt2(x, varargin)
%function [coef codes] = ir_odwt2(x, varargin)
%|
%| 2D orthonormal discrete wavelet transform (ODWT) of 2D image x.
%| Also provides the adjoint, which is the inverse ODWT.
%|
%| in
%|	x	[N M]	image
%|
%| option
%|	'level'	int	default: 1
%|	'wname'	char	default: 'haar'
%|	'ortho'	0|1	default: 1
%|	'adj'	0|1	default: 0 (if 1 then adjoint, which is inverse)
%|	'abs'	0|1	default: 0 (if 1 then take abs of filters)
%|
%| out
%|	coef	[N M]	2D wavelet coefficients
%|	codes	[N M]	uint2 level "codes" (only if not adj)
%|			code = 10 * level + (0 approx, 2 lohi, 4 hilo, 8 hihi)
%|
%| 2012-05-18, Jeff Fessler, Univ. of Michigan

if nargin < 1, ir_usage, end
if streq(x, 'test'), ir_odwt2_test, return, end

arg.level = 1;
arg.wname = 'haar';
arg.ortho = true;
arg.adj = false;
arg.abs = false;
arg.usemat = false;
arg = vararg_pair(arg, varargin);

if arg.level == 0
	coef = x;
	codes = zeros(size(x), 'uint16');
return
end

% decomposition filters
[lo hi] = ir_dwt_filters(arg.wname, ...
	'abs', arg.abs, 'ortho', arg.ortho, 'usemat', arg.usemat);

if arg.adj % synthesis (reconstruction)
	coef = ir_odwt2_back(x, lo, hi, arg.level); % trick
else % analysis (decomposition)
	[coef codes] = ir_odwt2_forw(x, lo, hi, arg.level);
end


% ir_odwt2_forw()
function [coef codes] = ir_odwt2_forw(x, lo, hi, level)

[N1 N2] = size(x);
[coef codes] = ir_odwt2_do(x, lo, hi, []);
for ll=2:level
	i1 = 1:(N1 / 2^(ll-1));
	i2 = 1:(N2 / 2^(ll-1));
	[coef(i1,i2) codes(i1,i2)] = ...
		ir_odwt2_do(coef(i1,i2), lo, hi, codes(i1,i2));
end


% ir_odwt2_back()
function x = ir_odwt2_back(coef, lo, hi, level)
lo_r = flipud(lo);
hi_r = flipud(hi);

[N1 N2] = size(coef);
for ll=level:-1:2
	i1 = 1:(N1 / 2^(ll-1));
	i2 = 1:(N2 / 2^(ll-1));
	coef(i1,i2) = ir_odwt2_adj(coef(i1,i2), lo_r, hi_r);
end
x = ir_odwt2_adj(coef, lo_r, hi_r);


% ir_odwt2_do()
function [coef codes] = ir_odwt2_do(x, lo, hi, codes)

[N M] = size(x);
if any(mod(size(x),2)) % odd
	fail('[%d %d] not even', N, M)
end

if isempty(codes)
	codes = zeros(N, M, 'uint16');
end
i1 = 1:N/2; i2 = 1:M/2;
codes(i1, i2) = codes(i1, i2) + 10;
codes(i1+N/2, i2) = codes(i1+N/2, i2) + 12;
codes(i1, i2+M/2) = codes(i1, i2+M/2) + 14;
codes(i1+N/2, i2+M/2) = codes(i1+N/2, i2+M/2) + 18;

coef_lo = ir_conv(x, lo, 'per', true);
coef_hi = ir_conv(x, hi, 'per', true);
i1 = (1:2:N); % [N/2] 2015-04-12 changed 2:2 to 1:2 per ir_conv change
coef_lo = coef_lo(i1,:); % [N/2 M]
coef_hi = coef_hi(i1,:); % [N/2 M]
coef_lolo = ir_conv(coef_lo, lo', 'per', true); % [N/2 M]
coef_lohi = ir_conv(coef_lo, hi', 'per', true);
coef_hilo = ir_conv(coef_hi, lo', 'per', true);
coef_hihi = ir_conv(coef_hi, hi', 'per', true);
i2 = (1:2:M); % [M/2] 2015-04-12 ""
coef = [coef_lolo(:,i2), coef_lohi(:,i2); ...
	coef_hilo(:,i2), coef_hihi(:,i2)]; % [N M]

if 1
	if numel(lo) > 2
		pad = (numel(lo) - 2) / 2; % matlab defaults to 'b' in 'wextend'
		x = [x((N-pad+1):N,:); x; x(1:pad,:)]; % 'b' periodic end conditions
		x = [x(:,(M-pad+1):M), x, x(:,1:pad)]; % 'b' periodic end conditions
	else
		pad = 0;
	end
	coef_lo = conv2(x, lo, 'same'); % [N+2*pad M+2*pad]
	coef_hi = conv2(x, hi, 'same');
	i1 = pad + (1:2:N); % [N/2]
	coef_lo = coef_lo(i1,:); % [N/2 M+2*pad]
	coef_hi = coef_hi(i1,:);
	coef_lolo = conv2(coef_lo, lo', 'same'); % [N/2 M+2*pad]
	coef_lohi = conv2(coef_lo, hi', 'same');
	coef_hilo = conv2(coef_hi, lo', 'same');
	coef_hihi = conv2(coef_hi, hi', 'same');
	i2 = pad + (1:2:M); % [M/2]
	xxcoef = [coef_lolo(:,i2), coef_lohi(:,i2); ...
		coef_hilo(:,i2), coef_hihi(:,i2)]; % [N M]
	equivs(xxcoef, coef)
end


% ir_odwt2_adj()
% synthesis
function x = ir_odwt2_adj(coef, lo, hi)

[N M] = size(coef);
if any(mod(size(coef),2)) % odd
	fail('[%d %d] not even', N, M)
end

tmp = zeros(N/2, M, class(coef));
coef_lolo = tmp;
coef_lohi = tmp;
coef_hilo = tmp;
coef_hihi = tmp;

i2 = 1:2:M; % 2015-04-12 ""
j1 = 1:N/2; j2 = 1:M/2;
coef_lolo(:,i2) = coef(j1, j2); % [N/2 M]
coef_lohi(:,i2) = coef(j1, j2+M/2);
coef_hilo(:,i2) = coef(j1+N/2, j2);
coef_hihi(:,i2) = coef(j1+N/2, j2+M/2);

% trick: flipud back because of 'adj'
coef_lo = ir_conv(coef_lolo, flipud(lo)', 'per', true, 'adj', 1) ...
	+ ir_conv(coef_lohi, flipud(hi)', 'per', true, 'adj', 1);
coef_hi = ir_conv(coef_hilo, flipud(lo)', 'per', true, 'adj', 1) ...
	+ ir_conv(coef_hihi, flipud(hi)', 'per', true, 'adj', 1);

i1 = 1:2:N; % 2015-04-12 ""
tmp = zeros(N, M, class(coef));
up_lo = tmp;
up_hi = tmp;
up_lo(i1,:) = coef_lo; % [N/2 M]
up_hi(i1,:) = coef_hi;

x = ir_conv(up_lo, flipud(lo), 'per', true, 'adj', 1) ...
	+ ir_conv(up_hi, flipud(hi), 'per', true, 'adj', 1);


% ir_odwt2_test_wname()
function ir_odwt2_test_wname(wname)

warg = {'wname', wname};

x = ellipse_im([160 192]); % stress test non-square
x(1) = 1; x(end) = 7; % stress boundary conditions
im plc 3 3
for level=1:3
	[coef code] = ir_odwt2(x, 'level', level, warg{:});
	im(0+level, coef)
	xlabelf('level %d', level)
	im(3+level, code)

	y = ir_odwt2(coef, 'level', level, 'adj', 1, warg{:});
	im(6+level, y)
%	clf, im(y-x)
	equivs(y, x)

	if exist('wavedec2', 'file') == 2
		[mat2 ss] = wavedec2(x, level, wname);
%		mat2 = reshape(mat2, size(x)); % does not work!
		sz = ss(1,:);
		nn = prod(sz);
		ii = 1:nn;
		lolo = reshape(mat2(ii+0*nn), sz);
		hilo = reshape(mat2(ii+1*nn), sz);
		lohi = reshape(mat2(ii+2*nn), sz);
		hihi = reshape(mat2(ii+3*nn), sz);
		mat2 = [lolo, lohi; hilo, hihi];
		im(6+level, mat2)
		if level == 1
			equivs(coef, mat2, 'fail', 0)
		end
		if streq(wname, 'haar')
			equivs(coef(1:(2*sz(1)), 1:(2*sz(2))), mat2)
		end
	end
end
xlabelf('wname = %s', wname)


% ir_odwt2_test()
function ir_odwt2_test
list = {'haar', 'sym2'};
for ii=1:numel(list)
	wname = list{ii};
	pr wname
	ir_odwt2_test_wname(wname)
end
