 function [deriv, curv] = xct_poly1_dercurv(data, li)
%function [deriv, curv] = xct_poly1_dercurv(data, li)
%
% evaluate derivatives and curvatures for polyenergetic Poisson transmission
% negative log-likelihoods: f(x) = \sum_i h_i(l), h_i(l) = mi(l) - yi log mi(l)
% mi(l) = bi exp(-Fi(l)) + ri
% where Fi models beam hardening effects due to polyenergetic spectrum.
% for a single material type (e.g., water)
%
% in
%	data	{cell}		{yi, bi, ri, Fi}
%				yi is [nd,1].  bi, ri, same or scalar
%				Fi is beam-hardening functions and derivatives
%				see xct_poly1_?? todo
%	li	[nd,1]
%
% out
%	deriv	[nd,1]		\dot hi(l)
%	curv	[nd,1]		surrogate curvature for hi at li
%
% Copyright 2004-2-1, Jeff Fessler, The University of Michigan

if nargin < 2, help(mfilename), error args, end

yi = data{1};
bi = data{2};
ri = data{3};
Fi = data{4};

% transmission Poisson likelihood function
[fi di] = feval(Fi.fun_der, Fi, li);
ei = exp(-fi);
mi = bi .* ei + ri;
deriv = (1 - yi ./ mi) .* (ri - mi) .* di;

% approximate curvatures:
% good for small ri and approximately linear Fi with maximum slope at l=0
dFi0 = feval(Fi.der0, Fi, li);
curv = yi .* dFi0.^2;
