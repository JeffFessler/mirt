 function out = make_block(in, nblock)
%function out = make_block(in, nblock)
% construct blocks of sparse matrix for Gtomo2_sparse object
% called by Gblock.m
% Copyright 2002-2-19, Jeff Fessler, The University of Michigan

if nargin < 2
	error(mfilename)
end

out = in; % copy entire object

nb = in.nb;
na = in.na;
G = in.G;

for iblock=1:nblock
	ia = iblock:nblock:na;
	ii = outer_sum(1:nb, (ia-1)*nb);
	out.blocks{iblock} = G(ii(:),:);
end
