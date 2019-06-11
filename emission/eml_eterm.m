 function eterm = eml_eterm(x, Gb, yi, ci, ri, iblock)
%function eterm = eml_eterm(x, Gb, yi, ci, ri, iblock)
%
% compute E-step (for a given block)
% yb = A*x + r; e = A' * (y ./ yb), where A = D(c)*G
% in:
%	x	[np,1]	image
%	yi,ci,ri	must be made 2d by "reshaper"
%
% Copyright 2005-2-17, Jeff Fessler, The University of Michigan

if nargin < 5, help(mfilename), error(mfilename), end

if nargin == 6 % os case, should be phased out!
	[nb na] = size(yi);
	nblock = block_ob(Gb, 'n');
	ia = iblock:nblock:na;
	eterm = eml_eterm_1(x, Gb{iblock}, ...
		col(yi(:,ia)), col(ci(:,ia)), col(ri(:,ia)));

elseif nargin == 5
	eterm = eml_eterm_1(x, Gb, yi, ci, ri);

else
	error 'not done'
end


function eterm = eml_eterm_1(x, Gb, yi, ci, ri)
li = Gb * x;
yb = ci .* li + ri; % predicted measurements
if any(yi & ~yb), warning 'model mismatch', end
yb(yb == 0) = inf;	% avoids /0 error
ratio = ci .* yi ./ yb;
eterm = Gb' * ratio;
