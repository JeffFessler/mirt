 function lists = subset_lists(nsubset, nview, varargin)
%function lists = subset_lists(nsubset, nview [, options])
%
% in
%	nsubset		# of subsets (1 to nview)
%	nview		# of views
% option
%	(none now, hard-wired to "bit-reversal" ordering)
% out
%	lists {cell}	{nsubset} lists of view indices for each subset
%
% Copyright 2000-4-?, Jeff Fessler, The University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

starts = subset_start(nsubset);

for ii=1:nsubset
	lists{ii} = [starts(ii):nsubset:nview];
end
