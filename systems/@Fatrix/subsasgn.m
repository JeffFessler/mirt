 function ob = subsasgn(ob, sub, arg)
%function ob = subsasgn(ob, sub, arg)
% method for "ob.sub = arg" etc.

% trick: for struct.sub = object
% e.g., G.arg.field = G;
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
	ob = subsasgn(ob, sub, arg);
	ob = class(ob, cl);
end
