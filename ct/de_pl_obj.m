 function [obj, like, penal] = de_pl_obj(xs, G, ymi, Im, rmi, ftab, R, mask)
%function [obj, like, penal] = de_pl_obj(xs, G, ymi, Im, rmi, ftab, R, mask)
%
%	compute quadratically penalized Poisson likelihood for each column of x
%	for dual-energy x-ray problem
%	in
%		xs	[np,2,niter]	iterates
%			[nx,ny,2,niter]	if mask provided
%		G	[nd,np]		system matrix
%		ymi	[nb,na,2]	raw polyenergetic measurements
%		Im	[2]		source intensities
%		rmi	[nb,na,2]	(or scalar) background/scatter
%		ftab			F table
%		R			penalty object
%		mask	[nx,ny]		support mask
%	out
%		obj	[niter]		objective function
%
%	Copyright 2002-1-29	Jeff Fessler	The University of Michigan

if nargin < 6, help(mfilename), error(mfilename), end
if ~isvar('rmi') | isempty(rmi)
	rmi = ones(size(rmi));
end
if ~isvar('R') | isempty(R)
	warning 'empty R means no penalty'
end

%
%	if a nonempty mask is passed in, it means the xs actually of size
%	[nx ny 2 niter] and it needs to be made [np 2 niter]
%
if isvar('mask') & ~isempty(mask)
	if size(xs,3) ~= 2, error bug, end
	xs = reshape(xs, [size(xs,1)*size(xs,2) 2 size(xs,4)]);	% [nxy 2 niter]
	xs = xs(mask(:), :, :);					% [np 2 niter]
end

penal = zeros(size(xs,3),1);
if ~isempty(R)
	for ii=1:length(penal)
		penal(ii)	= R{1}.penal(R{1}, xs(:,1,ii)) ...
				+ R{2}.penal(R{2}, xs(:,2,ii));
	end
end

like = de_p_like(xs, G, ymi, Im, rmi, ftab);	% likelihood

obj = like - penal;
%printf('%g\t',[like penal obj]')

%
%	likelihood for DE XCT model
%	with Poisson assumption
%
function like = de_p_like(xs, G, ymi, Im, rmi, ftab)

[nb na ns] = size(ymi);
if ns ~= 2, error ns2, end

niter = size(xs,3);
like = zeros(niter,1);
for ii=1:niter
	x = xs(:,:,ii);
	s1 = reshape(G * x(:,1), nb, na);
	s2 = reshape(G * x(:,2), nb, na);
	fh = ftab.feval(ftab, s1, s2);

	ybmi = zeros(nb, na, 2);
	ybmi(:,:,1) = Im(1) * exp(-fh(:,:,1)) + rmi(:,:,1);
	ybmi(:,:,2) = Im(2) * exp(-fh(:,:,2)) + rmi(:,:,2);

	like(ii) = sum(ymi(ymi>0) .* log(ybmi(ymi>0))) - sum(ybmi(:));
end
