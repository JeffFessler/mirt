 function [xs, ni] = tml_bitab(x, Gb, yi, bi, ri, niter, pixmax, relax)
%function [xs, ni] = tml_bitab(x, Gb, yi, bi, ri, niter, pixmax, relax)
% The T-ML-BITAB algorithm for transmission Poisson problem
% (block iterative transmission with AB lower/upper constraints)
% (proposed by Nayaryanan, Byrne, King)
% model: Y_i ~ Poisson(b_i exp(-[G x]_i) + r_i)
% in
%	Gb	Gblock object
%	yi	transmission sinogram
%	bi	blank scan factors
%		ri	background (randoms, scatter, crosstalk, etc)
%	bi,ri:	optional (can use empty matrices)
%	yi,bi,ri must have identical dimensions
%	niter	# of iterations
%	pixmax	upper constraint for pixel values (lower = 0)
%		can be scalar (e.g. 'inf') or an array the size of x
%	relax	ad hoc relaxation parameter
% out
%	x [np,niter]	updated image vectors each iteration
%
% Copyright Apr 2000, Jeff Fessler, University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end

nblock = block_ob(Gb, 'n');

if ~isvar('bi') || isempty(bi)
	bi = ones(size(yi));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi));
end
if ~isvar('niter'),	niter = 2; end
if ~isvar('pixmax'),	pixmax = 0.2; end
if ~isvar('relax'),	relax = 1; end

trl_check(yi, bi, ri);
if any(ri~=0), warning 'ri nonzero may be nonmonotone!', end

[nb na] = size(yi);
starts = subset_start(nblock);

xs = zeros(length(x), niter);
x = max(x,0);				% enforce nonnegativity
x = min(x,pixmax);			% enforce upper bound constraint
xs(:,1) = x;

% the relaxation parameter, assuming pixmin = 0
r = relax * 4 * max(pixmax - 0) / sum(sum(Gb'.^2)' .* bi(:))

%
% loop over iterations
%
for iter = 2:niter

	%
	% loop over subsets
	%
	for iset=1:nblock
		iblock = starts(iset);
		ia = iblock:nblock:na;

		li = Gb{iblock} * x;			% l=G*x "line integrals"
		li = reshape(li, nb, length(ia));
		yb = bi(:,ia) .* exp(-li) + ri(:,ia);	% predicted meas. means 

		dothi = bi(:,ia) .* (1 - yi(:,ia) ./ yb) .* exp(-li);

		A = pixmax - x;
		B = (x-0) .* exp(r * (Gb{iblock}' * dothi(:)));	% "minimizing"

		x = A./(A+B) * 0 + B./(A+B) .* pixmax;	% the update!
	end

	xs(:,iter) = x;
end
