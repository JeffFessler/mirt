 function ob = fatrix2_setup_dim(ob)
%function ob = fatrix2_setup_dim(ob)
%|
%| Set up idim and odim etc.
%|
%| Copyright 2010-12-17, Jeff Fessler, University of Michigan

arg = ob.arg; % this shall not be modified

if isempty(ob.imask) && isfield(arg, 'mask')
	ob.imask = arg.mask;
end


% input mask (rhs)

if isempty(ob.idim)
	if ~isempty(ob.imask)
		ob.idim = size(ob.imask);
	elseif isfield(arg, 'idim') && ~isempty(arg.idim)
		ob.idim = arg.idim;
	else
		fail('idim or imask or arg.idim required')
	end
end

if ~isempty(ob.imask)
	if ~islogical(ob.imask)
		fail('imask must be logical')
	end
	fatrix2_check_dim1(ob.imask, ob.idim)
	if isfield(arg, 'mask')
		if ~isequal(ob.imask, arg.mask)
			warn('mask inconsistent!?')
		end
	end
end


% output mask (lhs)

if isempty(ob.odim)
	if ~isempty(ob.omask)
		ob.odim = size(ob.omask);
	elseif isfield(arg, 'odim') && ~isempty(arg.odim)
		ob.odim = arg.odim;
	else
		fail('odim or omask or arg.odim required')
	end
end

if ~isempty(ob.omask)
	if ~islogical(ob.omask)
		fail('omask must be logical')
	end
	fatrix2_check_dim1(ob.omask, ob.odim)
end


ob.size = [prod(ob.odim) prod(ob.idim)];
if ~isempty(ob.imask)
	ob.size(2) = sum(ob.imask(:));
end
if ~isempty(ob.omask)
	ob.size(1) = sum(ob.omask(:));
end


% 1D input fully masked with >1D output is ambiguous
if ~ob.accept1d
	if numel(ob.idim) == 1 && numel(ob.odim) > 1 && ob.size(2) == prod(ob.idim)
%		fail 'use 1d odim for 1d idim case with full mask'
		warn '1d idim with full imask and >1d odim may not work'
	end

	if numel(ob.odim) == 1 && numel(ob.idim) > 1 && ob.size(1) == prod(ob.odim)
		warn '1d odim full omask, >1d idim: transpose may not work'
	end
end
