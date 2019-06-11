 function [x,denom] = qpwls_sps(x, G, W, vb, C, denom)
%function [x,denom] = qpwls_sps(x, G, W, vb, C, denom)
%	quadratically penalized weighted least squares reconstruction
%	using separable paraboloidal surrogates algorithm
%	cost(x) = x' G' W G x / 2 - b'x + x' C' C x /2
%	Input
%		x	[np,1]		old estimate
%		G	[nn,np]		system matrix, aij >= 0 required!
%		W	[nn,nn]		weighting matrix
%		vb	[nn,1]		"b" vector
%		C	[nc,np]		penalty matrix
%		denom	[np,1]		optional, otherwise computed
%	Output
%		x	[np,1]		new estimate
%
%	Copyright Apr 1999	Jeff Fessler, The University of Michigan

warning 'obsolete.  use pwls_sps_os'

if ~isvar('C')
	C = zeros(1,length(x));
end

%
%	denominator
%
if ~isvar('denom') || isempty(denom)
	denom = full(G' * (W * sum(G')'));
	denom = denom + abs(C)' * sum(abs(C)')';
	printm('denom min,max=%g,%g', min(denom(:)), max(denom(:)))
end

	x = x + (vb - G' * (W * (G * x)) - C' * (C * x)) ./ denom;
