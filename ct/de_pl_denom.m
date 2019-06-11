 function [denom, wi] = de_pl_denom(G, ymi, masseff, gi)
%function [denom, wi] = de_pl_denom(G, ymi, masseff, gi)
%
%	Precomputed curvatures and denominator (see fessler::) for
%	the Dual-Energy PL-OS-SPS algorithm for transmission Poisson problem.
%	(ordered subsets, separable paraboloidal surrogates)
%	in:
%		G	[nb*na,np]	system
%		ymi	[nb,na,2]	raw polyenergetic measurements
%		masseff	[2,2]		effective mass atten coef's
%					2 energies by 2 materials
%		gi			sum(G')'
%	out:
%		denom	[np,2]		(precomputed denominator)
%		wi	[nb,na,2]	weights (precomputed curvatures)
%
%	Copyright 2002-1-30	Jeff Fessler	The University of Michigan

if nargin < 4, help(mfilename), error(mfilename), end

%
%	precompute denominator
%	a single denominator shared by all subsets (scaled)
%

denom = zeros(ncol(G), 2);
for ll=1:2
	% precomputed curvatures:
%	wi = trl_curvature_pre(yi, ?, ri);
	wi(:,:,ll)	= ymi(:,:,1) * masseff(1,ll) * sum(masseff(1,:)) ...
			+ ymi(:,:,2) * masseff(2,ll) * sum(masseff(2,:));

	denom(:,ll) = G' * col(wi(:,:,ll) .* gi);
end
