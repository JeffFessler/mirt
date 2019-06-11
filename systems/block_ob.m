 function out = block_ob(varargin)
%function out = block_ob(Gb, arg)
%|
%| Extract information from a "block" object.
%|
%| There are (at least) three types of block objects possible in this toolbox.
%| 1. Cell arrays of matrices (or of some other matrix-like object).
%| 2. A Fatrix object with block capabilities (mtimes_block).
%| 3. The Gblock class (which is being phased out).
%|
%| This routine provides a common interface to all of them.
%|
%| arguments:
%|	'is'		is it a block object?
%|	'n'		nblock: # of blocks
%|	'ensure'	make it at least a 1-block block object (if needed)
%|
%| Copyright 2005-6-19, Jeff Fessler, University of Michigan

out = block_op(varargin{:}); % todo: retire block_ob and use block_op
return

if nargin < 2, help(mfilename), error(mfilename), end

switch arg
case 'is'
	if iscell(Gb)
		out = true;
	elseif isa(Gb, 'Fatrix')
		out = ~isempty(Gb.nblock);
	elseif isa(Gb, 'Gblock')
		out = true;
	else
		out = false;
	end

case 'n'
	if ~block_ob(Gb, 'is'), error 'not a block object', end
	if iscell(Gb)
		out = length(Gb);
	elseif isa(Gb, 'Fatrix')
		out = Gb.nblock;
	elseif isa(Gb, 'Gblock')
		out = Gb.nblock;
	end

case 'ensure'
	if block_ob(Gb, 'is')
		out = Gb;
		return
	end

	if iscell(Gb)
		error 'bug: cell already is a block object'
	elseif isnumeric(Gb)
		Gb = Gsparse(Gb);
		Gb = Gblock(Gb, 1);
	elseif isa(Gb, 'Fatrix') && ~isempty(Gb.nblock)
		out = Gb;
	else
		out = Gblock(Gb, 1);
	end

otherwise
	fail(['unknown argument: ' arg])
end
