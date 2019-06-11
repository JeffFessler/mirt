 function ob = subsasgn(ob, sub, arg)
%function ob = subsasgn(ob, sub, arg)
% method for "ob.sub = arg" etc.

% trick: for struct.sub = object
% e.g., ob.arg.field = ob;
if isstruct(ob)
%	ob = subsasgn(ob, sub, arg); % no, causes recursion
%	ob = builtin('subsasgn', ob, sub, arg);	% samit basu dislikes this

	% recurse a.b.c.d... to avoid builtin:
	if streq(sub(1).type, '.')
		if length(sub) > 1
%			printf('descending to %s', sub(1).subs)
			ob.(sub(1).subs) = subsasgn(ob.(sub(1).subs), sub(2:end), arg);
		else
			ob.(sub(1).subs) = arg;
		end
	else
		error 'not done'
	end

% trick: just convert it to a structure and then back!
else
	cl = class(ob);
	ob = struct(ob);
	name = sub.subs;
	if isfield(ob.data, name)
%		ob.data.(name) = arg;
		% 2011-03-12 changed to handle ig.mask(1)=0
		ob.data = subsasgn(ob.data, sub, arg);
	elseif isfield(ob.meth, name)
		ob.meth.(name) = arg;
	else
		error(['"' name '"is not in data or meth'])
	end
%	ob = subsasgn(ob, sub, arg);
	ob = class(ob, cl);
end
