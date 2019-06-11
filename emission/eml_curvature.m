  function ni = eml_curvature(yi, ci, ri, li, yb, ctype)
%|function ni = eml_curvature(yi, ci, ri, li, yb, ctype)
%|
%| compute surrogate parabola curvatures for Poisson emission model:
%| yi ~ Poisson(ci * li + ri)
%|
%| in
%|	yi	count data
%|	ci	calibration factors
%|	ri	mean background counts (scatter, randoms, etc.)
%|	li	"line integral" of activity (ci*li is unitless, i.e., "counts")
%|	yb	(optional) at this point: default yb = ci * li + ri
%|	ctype
%|		oc	erdogan's optimal curvatures
%|		pc	precomputed curvatures, ala empl3
%|			(not guaranteed to provide a majorizer)
%|				fix: align with C version of empl2
%|		todo: add max curv ? (precomputed and a majorizer)
%|
%| out
%|	ni	curvature of negative log-likeihood (which is nonnegative)
%|
%| Jeff Fessler

if nargin == 1 && streq(yi, 'test'), eml_curvature_test, return, end
if nargin < 5, help(mfilename), error(mfilename), end

if ~isvar('ctype') || isempty(ctype)
	ctype = 'oc'
end

if streq(ctype, 'pc')
	if any(ci(:) < 0), error 'nonpositive ci', end
	li = (yi - ri) ./ ci;
	li = max(li, 0);
	ybi = ci .* li + ri;
	if any(ybi(:) <= 0 & yi(:) > 0), error 'model mismatch', end
%	ni = ci.^2 .* yi ./ max(ybi, eps*max([1; ybi(:)])).^2;
	% see eml_curv_pre
	ni = zeros(size(li));
	ii = li > 0;
	ni(ii) = ci(ii).^2 ./ max(yi(ii),1);
	ii = li <= 0 & ri > 0;
	ni(ii) = (ci(ii) ./ ri(ii)).^2 .* yi(ii);

elseif streq(ctype, 'oc')

	if ~isvar('yb') || isempty(yb)
		yb = ci .* li + ri;
	end

	ni_max = yi .* (ci ./ ri).^2;	% curvature at l=0
	ni = 0.001 * (ci ./ ri).^2;	% small nonzero value for yi=0

	tmp = log(yb ./ ri) - ci.*li ./ yb;
	iy = yi ~= 0; % non-zero data

	if 0
		il = li <= 0;
	else % trick in C program due to numerical precision issues
		il = li < 0.1 * ri ./ ci;
	end
	i = iy & il;
	ni(i) = ni_max(i);

	i = iy & ~il;
	ni(i) = 2 * yi(i) ./ li(i).^2 .* tmp(i);

	if any(ni <= 0), error 'nonpositive ni', end
	if any((ni > ni_max) & yi)
		plot([ni_max(:) ni(:) ni(:)>ni_max(:)])
		error 'large ni'
	end
%	range(ni,2)

else
	error 'unknown curvature'
end


% eml_curvature_test
% show figure of surrogate with optimal curvature
function eml_curvature_test
if 1
	l = linspace(0,2,1001)';
	l = logspace(-6,1,101)';
	n = eml_curvature(2+0*l, 1+0*l, 1+0*l, l, [], 'oc');
	if im
		semilogx(l, n, '-o')
%	prompt
	end
end

y = 3; c = 2; r = 1;
h = @(y, c, r, l) y * log(c*l+r) - (c*l+r);
l0list = [0 2 3];
l = linspace(-0.25,4,101);
hh = h(y,c,r,l);
args = {l, hh, '-'};
lines = {'y--', 'm--', 'g--'};
point = {'yo', 'mo', 'go'};
for ii=1:numel(l0list)
	l0 = l0list(ii);

	nc = eml_curvature(y, c, r, l0, [], 'oc');
	h0 = h(y,c,r,l0);
	derh = c * (y ./ (c*l0+r) - 1);
%	nc = 2/l0^2 * (h0 - h(y,c,r,0) - l0*derh);
%	nerh = c^2 * y / (c*l+r)^2;
	q = h0 + derh * (l-l0) - 0.5 * nc * (l-l0).^2;
	args = {args{:}, l, q, lines{ii}, l0, h0, point{ii}};
end

if im
	plot(args{:})
	axis([min(l) max(l) min(hh) max(hh)+0.1]), grid
	xlabel 'l', legend('h(l)', 'q(l;l_0)')
	title 'Illustration of parabola surrogate for emission case'
end
