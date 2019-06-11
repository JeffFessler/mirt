 function ni = trl_curvature(yi, bi, ri, li, ctype)
%function ni = trl_curvature(yi, bi, ri, li, ctype)
%|
%| Compute surrogate parabola curvatures for Poisson transmission model
%| ctype:
%|	oc	erdogan's optimal curvatures
%|	pc	precomputed curvatures, ala trpl3
%|	nc	newton curvatures
%|
%| fix: align with C version of trpl2,3
%| The minimum returned curvature will be zero.
%| It is the user's responsibility to impose additional bounds
%| if desired for certain algorithms.
%|
%| Copyright 2002-1-28, Jeff Fessler, University of Michigan

if nargin == 1 && streq(yi, 'test'), trl_curvature_test, return, end
if nargin < 4, help(mfilename), error(mfilename), end

[h dh] = trl_h_dh;

if ~isvar('ctype') || isempty(ctype)
	ctype = 'oc'
end

% Compute optimal surrogate parabola curvatures
% for Poisson transmission model based on Erdogan's formula.
switch ctype
case 'oc'

	% compute curvature at l=0
	ni_max = zeros(size(yi));
	if numel(bi) == 1 % scalar bi (must be positive!)
		ni_max = bi .* (1 - yi .* ri ./ (bi + ri).^2);
	else
		i0 = bi > 0;
		if numel(ri) == 1
			rii = 1;
		else
			rii = ri(i0);
		end
		ni_max(i0) = bi(i0) .* (1 - yi(i0) .* rii ./ (bi(i0) + rii).^2);
	end
	ni_max = max(ni_max, 0);
	ni = ni_max;

	if 0
		il0 = li <= 0;
	else % trick in C program due to numerical precision issues
		il0 = li < 0.1;
	end

	tmp = h(yi,bi,ri,li) - h(yi,bi,ri,0) - li .* dh(yi,bi,ri,li);
	i = ~il0;
	ni(i) = 2 ./ li(i).^2 .* max(tmp(i),0);

	if any(ni > ni_max)
	%	plot([ni_max(:) ni(:) ni(:)>ni_max(:)])
		warning 'large ni'
	end


% Precomputed approximating parabola curvatures
% for Poisson transmission model.
% The minimum returned curvature will be zero.
% This is compatible with trpl/trp_init_der02_sino() in aspire.
case 'pc'

	% ni = (yi-ri)^2 / yi, if yi > ri >= 0 and bi > 0
	ii = (yi > ri) & (ri >= 0) & (bi > 0); % good rays
	ni = zeros(size(yi));
	ni(ii) = (yi(ii) - ri(ii)).^2 ./ yi(ii);

% newton curvatures (current 2nd derivative)
case 'nc'
	bel = bi .* exp(-li);
	yb = bel + ri;
	ni = (1 - ri.*yi./yb.^2) .* bel;

otherwise
	fail 'bug'
end


% trl_h_dh()
% transmission Poisson likelihood function
function [h, dh] = trl_h_dh
h = @(y,b,r,l) y .* log(b.*exp(-l)+r) - (b.*exp(-l)+r);
dh = @(y,b,r,l) (1 - y ./ (b.*exp(-l)+r)) .* b.*exp(-l);


% trl_curvature_test()
% demonstrate an example surrogate parabola!
function trl_curvature_test
[h dh] = trl_h_dh;
if 0
	l = linspace(0,2,101)';
	n = trl_curvature(2+0*l, 3, 3, l, 'oc');
	plot(l, n, '-o'), xlabel l, ylabel n
else
%	y = 4; b = 3; r = 1; ln = 3.8;
%	y = 2; b = 3; r = 3; ln = 0.05;
	y = 2; b = 3; r = 1; ln = 1;
	n.oc = trl_curvature(y, b, r, ln, 'oc');
	n.pc = trl_curvature(y, b, r, ln, 'pc');
	l = linspace(-0.2,5,101)';
	hn = h(y,b,r,ln);
	dhn = dh(y, b, r, ln);%
	hl = h(y,b,r,l);
	q.oc = hn + dhn * (l-ln) - 0.5 * n.oc * (l-ln).^2;
	q.pc = hn + dhn * (l-ln) - 0.5 * n.pc * (l-ln).^2;
	if im
		clf, plot(l, hl, '-', l, q.oc, '--', l, q.pc, '-.', ln, hn, 'o')
		grid
		axis([minmax(l)' minmax(hl)'])
		legend('h(l)',	sprintf('q_o(l) n_i=%g', n.oc), ...
				sprintf('q_p(l) n_i=%g', n.pc))
	end
end
