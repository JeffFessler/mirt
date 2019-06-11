 function x = eml_psca(x, G, yi, ci, ri, nsubiter)
%function x = eml_psca(x, G, yi, ci, ri, nsubiter)
% Runs one iteration of the ML-PSCA algorithm for emission Poisson problem
% (paraboloidal surrogates coordinate ascent)
% model: Y_i ~ Poisson(c_i [G x]_i + r_i)
% in
%	see em_fbp.m for model, G, yi, ci, ri
% out
%	x [np,1]	updated image
%
% This is for TESTING ONLY.
% IT IS NOT USEFUL because it is unregularized and slow
%
% Copyright Mar 2000, Jeff Fessler, The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

[nb na] = size(yi);

if ~isvar('ci') || isempty(ci)
	ci = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end
if ~isvar('nsubiter') || isempty(nsubiter)
	nsubiter = 1;
end

%
% compute surrogate curvatures
%
li0 = reshape(G * x(:), size(yi));		% l=G*x "line integrals"
yb = ci .* li0 + ri;				% predicted measurement means 

doth = ci .* (yi ./ yb - 1);			% first derivative
ni = eml_curvature(yi, ci, ri, li0, yb, 'oc');	% curvatures
ni = ni(:);
%disp(range(ni)')

% backproject curvatures for denominator
dj = (G').^2 * ni(:);

%
% now we are maximizing the quadratic surrogate function:
% Q(x) == P(G*x ; G*x0)
%	where P(l;l0) = doth' (l-l0) - 1/2 (l-l0)' D(n) (l-l0)
% d/dli P(l;l0) = dothi - ni (li-l0i)
% Q(x) == doth' G (x-x0) - 1/2 (x-x0)' G' D(n) G (x-x0)
% ??? d/dx_j Q(x) = \sumi \gij \dothi - \sumi \gij \ni [G(x-x0)]_i
% d/dx_j Q(x) = \sumi \gij [d/dli P(l;l0)]
% -d^2/dx_j^2 Q(x) = \sumi \gij^2 \ni
%
dqi = doth(:);	% initial surrogate derivatives

%xs = zeros([length(x) nsubiter]);	% save and return the subiterations too

%x0 = x;		% save initial image

like0 = eql_obj(x, G, yi(:), ci(:), ri(:));

%
% loop over subiterations of paraboloid
%
for is=1:nsubiter

	%
	% loop over pixels
	%
	for jj=1:numel(x)
		g = G(:,jj);
		x0j = x(jj);
		x(jj) = x0j + (dqi' * g) / dj(jj);
		x(jj) = max(x(jj),0);
		dqi = dqi - ni .* (g * (x(jj) - x0j));

	end
%	xs(:,is) = x;
end
