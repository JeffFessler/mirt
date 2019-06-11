 function [deriv, curv] = Reg1_com_dercurv(sr, C1x)
%function [deriv, curv] = Reg1_com_dercurv(sr, C1x)
%|
%| evaluate [\dpoti] and [\wpoti]
%|
%| because it has two output arguments, this requires feval()
%| due to strum limitation inherited from matlab subsref limitation.
%| output sizes match input size
%|
%| in
%|	sr			Reg1 object
%|	C1x	[(N) M]		C1 * x
%|
%| out
%|	deriv	[(N) M]		dpot(C1 * x) derivative of potential function
%|	curv	[(N) M]		wpot(C1 * x) Huber weighting for potential fun.
%|
%| Jeff Fessler, University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

siz = size(C1x);
C1x = reshape(C1x, [prod(sr.dim) sr.M]); % [*N M]
deriv = zeros(size(C1x));
curv = zeros(size(C1x));
for mm=1:sr.M
	cm = C1x(:,mm);
	wt = sr.wt.col(mm);
	pot = sr.pot{mm};
	deriv(:,mm) = wt .* pot.dpot(cm); % use wpot(t) * t ?
	curv(:,mm) = wt .* pot.wpot(cm);
end
deriv = reshape(deriv, siz);
curv = reshape(curv, siz);
