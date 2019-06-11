 function out = block_op(Ab, varargin)
%function out = block_op(Ab, varargin)
%|
%| Various operations on "block" objects and block data.
%|
%| There are (at least) three types of block objects possible in this toolbox.
%| 1. Cell arrays of matrices (or of some other matrix-like object).
%| 2. A Fatrix or fatrix2 object with block capabilities (mtimes_block).
%| 3. The Gblock class (which is being phased out).
%|
%| This routine provides a common interface to all of them.
%|
%| arguments for object operations:
%|	'is'		is it a block object?
%|	'n'		nblock: # of blocks
%|	'ensure'	make it at least a 1-block block object (if needed)
%|
%| arguments for data operations:
%|	'ensure_block_data', data	rearrange arrays into cells
%|
%| Copyright 2005-6-19, Jeff Fessler, University of Michigan

if nargin < 2, ir_usage, end

switch length(varargin)
case 1
	out = block_op_ob(Ab, varargin{1});
case 2
	out = block_op_data(Ab, varargin{1}, varargin{2});
otherwise
	fail 'bad args'
end


% block_op_data()
function out = block_op_data(Ab, arg, data)

if streq(arg, 'ensure_block_data')
	if iscell(data)
		out = data; % ok, already cell

	elseif isa(Ab, 'Fatrix') || isa(Ab, 'fatrix2')
		if ~isempty(Ab.handle_blockify_data)
			out = blockify_data(Ab, data);
			return
		end

		% attempt to be smart for backward compat!
		try
			nblock = block_op(Ab, 'n');
			starts = subset_start(nblock);
			out = cell(1,nblock);

			if ndims(data) == 3 % [ns,nt,na]
				data = reshaper(data, '2d'); % [ns*nt,na]
			end

			na = size(data,2);
			if nblock > na, fail 'tried 2d but failed', end

			for iset=1:nblock
				istart = starts(iset);
				ia = istart:nblock:na;
				out{istart} = col(data(:,ia));
			end

		catch
			fail('failed to blockify data')
		end

	else
		fail('unknown blockify case "%s"', class(Ab))
	end

else
	fail 'bug'
end


% block_op_ob()
function out = block_op_ob(Ab, arg)

switch arg
case 'is'
	switch class(Ab)
	case 'cell'
		out = true;
	case {'Fatrix', 'fatrix2'}
		out = ~isempty(Ab.nblock);
	case 'Gblock'
		out = true;
	otherwise
		out = false;
	end

case 'n'
	if ~block_ob(Ab, 'is'), fail 'not a block object', end
	switch class(Ab)
	case 'cell'
		out = length(Ab);
	case {'Fatrix', 'fatrix2'}
		out = Ab.nblock;
	case 'Gblock'
		out = Ab.nblock;
	otherwise
		fail 'bug'
	end

case 'ensure'
	if block_ob(Ab, 'is')
		out = Ab;
		return
	end

	if iscell(Ab)
		fail 'bug: cell already is a block object'
	elseif isnumeric(Ab)
		Ab = Gsparse(Ab);
		Ab = Gblock(Ab, 1);
	elseif (isa(Ab, 'Fatrix') || isa(Ab, 'fatrix2')) && ~isempty(Ab.nblock)
		out = Ab;
	else
		out = Gblock(Ab, 1);
	end

otherwise
	fail('unknown argument: %s', arg)
end
