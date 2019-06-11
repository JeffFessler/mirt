 function out = subsref(ob, arg)
%function out = subsref(ob, arg)
% evaluate reference of the form ob.ref or ob(ref,:)

%
% ob.?
% get from Gblock object if possible, otherwise from base object
%

if arg(1).type == '.'

	if isfield(struct(ob), arg(1).subs)
		out = struct(ob);
		out = out.(arg(1).subs);
		if length(arg) > 1
			out = subsref(out, arg(2:end));
		end
		return
	end

	if isfield(struct(ob.base), arg(1).subs)
		out = struct(ob.base);
		out = out.(arg(1).subs);
		if length(arg) > 1
			out = subsref(out, arg(2:end));
		end
		return
	end

	% otherwise, error
	disp 'field names of base system'
	disp(fieldnames(ob.base))
	printf('subscript "%s" not a field of object', arg.subs)
	error bug

elseif length(arg) > 1
	disp '-------------'
	disp(arg)
%	disp(arg(:).type)
%	disp(arg(:).subs)
	error 'not done'

%
% ob{i_block}
% this is primary method 'handled' by this object!
%
elseif arg.type == '{}' && length(arg.subs) == 1

	out = ob;
	out.i_block = arg.subs{1};
	if out.i_block < 1 || out.i_block > out.nblock
		error 'bad block index'
	end

%
% ob(i_block)
% included for backward compatibility
%
elseif arg.type == '()' && length(arg.subs) == 1

	warning 'please use {} instead of () now for block objects'
	out = ob;
	out.i_block = arg.subs{1};
	if out.i_block < 1 || out.i_block > out.nblock
		error 'bad block index'
	end

%
% anything else gets passed to base object
% and we hope that it can handle it!
%
else
	tmp = subsref(ob.base, arg);
	if isa(tmp, class(ob.base))
		out = ob;
		out.base = tmp;
	else
		out = tmp;
	end
end
