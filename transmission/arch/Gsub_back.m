 function back = Gsub_back(Gsub, yi)
%
%	backprojector for Gsub object
%	in:
%		Gsub		Gsubset object
%		yi	[?,na]		data
%	out:
%		back	[np,2]	backprojection G'y
%
%	Copyright 2002-1-30	Jeff Fessler	The University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

na = Gsub.na;
np = Gsub.np;
nsubset = Gsub.nsubset;

yi = reshape(yi, [numel(yi)/na na]);	% [?,na]

%
%	loop over subsets
%
back = zeros(np, 1);
for iset=1:nsubset
	ia = Gsub.starts(iset):nsubset:na;
%	ia = Gsub.ia{iset};
	back = back + Gsub.cell{iset} * col(yi(:,ia));
end
