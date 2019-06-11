 function out = fatrix2_subsref_brace(ob, iblock)
%function out = fatrix2_subsref_brace(ob, iblock)
%|
%| Manage subscript references of the form ob{iblock}
%| This will called from ../subsref with ob
%|
%| User can use subset_starts to select blocks in other orders
%| in which case this "iblock" is really "istart"
%|
%| ob = Do * A * Di => ob{k} = Do{k} * A{k} * Di
%|
%| Copyright 2010-12-17, Jeff Fessler, University of Michigan

out = ob;
if isempty(ob.nblock), fail('not a block object'), end

if iblock < 1 || iblock > ob.nblock
	fail('bad block index %d of %d', iblock, ob.nblock)
end

odim = ob.odim;
na = odim(end); % subset along last dimension of output
ia = iblock:ob.nblock:na;
out.odim = [odim(1:end-1) length(ia)];

if ~isempty(out.omask)
	fail 'block object with omask unsupported'
end
out.size = [prod(out.odim) ob.size(2)];

if isempty(ob.handle_forw_block)
	if isequal(ob.nblock, 1)
		% trick: keep usual handle_forw and handle_back links
	else
		fail 'bug: no forw_block() method for this object'
	end
else
	tmp = sprintf('@(arg, x) ob.handle_forw_block(arg, x, %d, %d)', ...
		iblock, ob.nblock);
%	out.handle_forw = @(arg, x) ...
%		ob.handle_forw_block(arg, x, iblock, ob.nblock);
	out.handle_forw = eval(tmp); % trick
end

if isempty(ob.handle_back_block)
	if isequal(ob.nblock, 1)
		% trick: keep usual handle_forw and handle_back links
	else
		fail 'bug: no back_block() method for this object'
	end
else
	tmp = sprintf('@(arg, x) ob.handle_back_block(arg, x, %d, %d)', ...
		iblock, ob.nblock);
%	out.handle_back = @(arg, x) ...
%		ob.handle_back_block(arg, x, iblock, ob.nblock);
	out.handle_back = eval(tmp); % trick
end

% out.idiag is unchanged, but out.odiag must be reduced:
if ~isempty(out.odiag)
	tmp = out.odiag;
	tmp = reshape(tmp, [], na);
	tmp = tmp(:, ia);
	tmp = reshape(tmp, out.odim);
	out.odiag = tmp;
end
