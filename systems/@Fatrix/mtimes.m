 function y = mtimes(ob, x)
%function y = mtimes(ob, x)
% y = M * x	or x = M' * y

if isa(x, 'Fatrix') % {object1 | scalar} * object2, tested in Gcascade
	y = Gcascade(ob, x);
return
end

if ~isa(ob, 'Fatrix')
	error 'only multiplication on right is done'
end

%
% partial multiplication?
%
if ob.is_subref
	error 'subref not done'
end

%
% block multiplication (if needed)
%
if ~isempty(ob.nblock) % block object
	if ~isempty(ob.handle_mtimes_block) % proper block object
		if ~isempty(ob.iblock) % Gb{?} * ?
			y = mtimes_block(ob, x, ob.iblock, ob.nblock);
			return
		% else: % Gb * ?
		end
	% else should be a 1-block object
	elseif ob.nblock > 1
		error 'nblock > 1 but no mtimes_block?'
	end
end

%
% ordinary multiplication
%
if ob.is_transpose
	y = Fatrix_mtimes_back(ob, x);
else
	y = Fatrix_mtimes_forw(ob, x);
end


%
% full forward multiplication ("projection")
% y = after * ob * before * x;
%
function y = Fatrix_mtimes_forw(ob, x);

x = do_cascade(ob.cascade_before, x, false, 0, 1, true);
y = ob.handle_forw(ob.arg, x);
y = do_cascade(ob.cascade_after, y, false, 0, 1, false);


%
% full transposed multiplication ("back-projection")
% x = before' * ob' * after' * y;
%
function x = Fatrix_mtimes_back(ob, y);

y = do_cascade(ob.cascade_after, y, true, 0, 1, false);
x = ob.handle_back(ob.arg, y);
x = do_cascade(ob.cascade_before, x, true, 0, 1, true);
