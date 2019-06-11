 function ob = subsasgn(ob, sub, arg)
%function ob = subsasgn(ob, sub, arg)
% method for "ob.sub = arg" etc.
% caution: reassigning dimensions (.idim .odim .imask .omask) is dangerous!

if ~streq(sub(1).type, '.')
	fail('only "ob.field = value" is supported')
end

% trick: for struct.sub = object
% e.g., A.arg.field = A;
if isstruct(ob)
%	ob = subsasgn(ob, sub, arg); % no, causes recursion
%	ob = builtin('subsasgn', ob, sub, arg);	% samit basu dislikes this

	% recurse a.b.c.d... to avoid builtin:
	if numel(sub) > 1
		ob.(sub(1).subs) = subsasgn(ob.(sub(1).subs), sub(2:end), arg);
	else
		ob.(sub(1).subs) = arg;
	end

% trick: just convert it to a structure and then back!
else
	cl = class(ob);
	ob = struct(ob);
	name = sub.subs;

	% note: ob.arg.mask = new_mask may cause errors.
	if numel(sub) > 1 % ob.arg.name = ... or ob.meth.name =
		ob.(sub(1).subs) = subsasgn(ob.(sub(1).subs), sub(2:end), arg);

	elseif isfield(ob, name) % reassigning private object data is dangerous
		fail('too dangerous to subsasgn "%s"', name)

%	elseif isfield(ob.meth, name) % user may reassign method handles
%		ob.meth = subsasgn(ob.meth, sub, arg);

%	elseif isfield(ob.arg, name) % user may reassign arg internals
%		ob.arg = subsasgn(ob.arg, sub, arg);

	else
		fail('%s is not a field', name)
	end

	ob = class(ob, cl);
end
