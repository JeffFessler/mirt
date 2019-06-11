  function ob = Gblock(base, nblock)
%|function ob = Gblock(base, nblock)
%|
%| Construct Gblock object, which is a 'meta' system object,
%| designed to work with ordered-subsets (aka block iterative) algorithms.
%| See Gblock_test.m for example usage.
%|
%| A Gblock system is created from another 'block-izable' base system object
%| 'base' by calling: Gb = Gblock(base, nblock, 1);
%|
%| Any projector object can be used if nblock=1, otherwise the base object
%| must have a "mtimes_block" method for block multiplications.
%|
%| The overloaded capabilities of a Gblock object include:
%|	Ab * x			full forward projection
%|	Ab' * y			full back projection
%|	Ab{i_block} * x		project one block (e.g., a set of views)
%|	Ab{i_block}' * x	back project one block
%|				i_block = 1,2,...,nubset
%|
%| Other properties are inherited from the base system object.
%|
%| Copyright 2002-2-18, Jeff Fessler, University of Michigan

if nargin < 2, nblock = 1; end

% create default object, as required by Mathworks
ob.base = [];
ob.nblock = nblock;
ob.i_block = 0;

if nargin == 0	% required by Mathworks
	ob = class(ob, 'Gblock');
return
end

% handle matrix as a special case
if isnumeric(base)
	error 'Gblock no longer supports numeric; call Gmatrix or Gsparse first'
%	r = input('do you wish to continue? [n|y]', 's');
%	if ~streq(r, 'y'), error '', end
%	base = Gsparse(base);
%	if nblock ~= 1, error 'call Gsparse first', end
end

% ensure that a "mtimes_block" method is available (if needed)
if nblock > 1
	try
		mtimes_block(base, 'exists');
	catch
		error 'this base object does not have a mtimes_block method'
	end
end

ob.base = base;
ob = class(ob, 'Gblock');
