 function xs = wls_pscd(x, G, yi, dqi, wi, dj, niter)
% Runs one iteration of coordinate descent in the style suitable
% for the inner part of paraboloidal surrogates coordinate descent updates
% cost function: J(x) = ?
% Output
%	x [np,1]	updated image
%
% This is for TESTING ONLY.
% IT IS NOT USEFUL because it is unregularized and slow
%
% Copyright Aug 2000, Jeff Fessler, The University of Michigan

NOT DONE!

if nargin < 3, help(mfilename), error(mfilename), end

if ~isvar('wi') || isempty(wi)
	wi = ones(size(yi));
end
if ~isvar('niter') || isempty(niter)
	niter = 1;
end

% backproject curvatures for denominator
if ~isvar('dj') || isempty(dj)
	dj = wi(:)' * G.^2;
end

%
% now we are minimizing the quadratic function:
% J(x) = \sum_k q_i([Gx]_i)
%	where q_i(l) = doth' (l-l0) - 1/2 (l-l0)' D(n) (l-l0)
% ?
% d/dli P(l;l0) = dothi - ni (li-l0i)
% Q(x) == doth' G (x-x0) - 1/2 (x-x0)' G' D(n) G (x-x0)
% ??? d/dx_j Q(x) = \sumi \gij \dothi - \sumi \gij \ni [G(x-x0)]_i
% d/dx_j Q(x) = \sumi \gij [d/dli P(l;l0)]
% -d^2/dx_j^2 Q(x) = \sumi \gij^2 \ni
%
dqi = doth(:);	% initial surrogate derivatives

xs = zeros(length(x), nsubiter);
xs(:,1) = x;

%
% loop over subiterations of paraboloid
%
for it=1:niter

	%
	% loop over pixels
	%
	for jj=1:numel(x)
		g = G(:,jj);
		x0j = x(jj);
		x(jj) = x0j + (dqi' * g) / dj(jj);
		x(jj) = max(x(jj),0);
		dqi = dqi - wi .* (g * (x(jj) - x0j));

	end
	xs(:,it) = x;
end
