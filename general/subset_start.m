 function [starts, nsubset] = subset_start(nsubset)
%function [starts, nsubset] = subset_start(nsubset)
%|
%| Compute array of subset starting indices "starts" for OS algorithms.
%| If input is an empty matrix, then 1 subset is used.
%| If input is a scalar power of 2 != 1,
%| then the "bit-reversal ordering" is used.
%| If input is an array, then it is simply checked for completeness
%| and returned.
%| (nsubset can be 1 to #views)
%|
%| See guan:94:apa doi 10.1088/0031-9155/39/11/013
%|
%| Copyright 2000-4-?, Jeff Fessler, University of Michigan

% if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(nsubset, 'test'), subset_start_test, return, end

if ~nargin || isempty(nsubset)
	nsubset = 1;
	starts = 1;
return
end

if length(nsubset) == 1
	starts = 1 + bit_reverse(nsubset);
end

nsubset = length(starts);

if any(sort(starts(:)) ~= [1:nsubset]')
	error 'missing subset?'
end


% bit_reverse()
% return (0:mm-1) in bit reversed order
% if mm is not a power of 2, then round up and discard extras
%
function ii = bit_reverse(mm)
nn = 2^ceil(log2(mm));
%if mm ~= nn, warning 'not power of 2', end
ii = bin2dec(fliplr(dec2bin([0:(nn-1)])));
ii = ii(ii < mm);


% subset_start_test
function subset_start_test
disp([subset_start(8)-1]')
