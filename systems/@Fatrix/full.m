  function out = full(ob)
% function out = full(ob)
%| full(ob) = ob(:,:)
%| Copyright 2002-2-20, Jeff Fessler, University of Michigan

np = ob.dim(2);
% out = ob(:, 1:np); % subsref seems not to work!?

% one column at a time by matrix multiplication
out = zeros(ob.dim(1),np);
for jj=1:np
	x = zeros(np,1);
	x(jj) = 1;
	out(:,jj) = ob * x;
end
