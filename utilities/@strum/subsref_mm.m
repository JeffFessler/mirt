 function varargout = subsref_mm(ob, args)
%function varargout = subsref_mm(ob, args)
% handle subscript references like ob.ref or ob.ref(arg[s])
% Copyright 2006-1-19, Jeff Fessler, University of Michigan
% 2011-10-25 modified by Madison McGaffin to support multiple outputs

st = struct(ob);
varargout = cell(nargout,1);

if args(1).type ~= '.'
	error 'only ob.data and ob.meth(args) done'
end

%
% ob.name
%
name = args(1).subs;
if isfield(st, name)
	if length(args) == 1
%		out = st.(name);
		varargout = wrap_vararg(nargout, 'st.(name)', st, ob, args);
	else
%		out = subsref(st.(name), args(2:end));
		varargout = wrap_vararg(nargout, ...
			'subsref(st.(args(1).subs), args(2:end));', ...
			st, ob, args );
	end

%
% ob.name, where name is in data field
%
elseif isfield(st.data, name)
%	out = st.data.(name);
	varargout = wrap_vararg(nargout, 'st.data.(args(1).subs);', st, ob, args);

	if length(args) > 1

		% ob.name{:} does not work properly because there seems
		% to be no way to return an arbitrary number of arguments
		% (in a list) to the caller.  this is unfortunate because
		% a.b = {1,2} and a.b{:} works for ordinary structures.
% todo: see this solution involving numel()
% http://www.mathworks.com/support/solutions/data/1-1ABOD.html?solution=1-1ABOD
% n.b., link broken; will need to rewrite this for multiple return arguments.

%		if iscell(out) && ...
		if length(varargout) > 0 && iscell(varargout{1}) && ...
		(isequal(args(2).subs, {':'}) || length(args(2).subs{1}) > 1)
			warn 'instead of strum.data{:} etc., use two steps:'
			warn 'tmp = strum.data; tmp{:} or tmp{1:end} etc.'
			fail 'strum.data{:} etc. unsupported thanks to matlab'
		end

		% fix: strum1.strum2.plot() does not work here
%		out = subsref(out, args(2:end));
		varargout = {subsref(varargout{1}, args(2:end))};
	end

%
% ob.name or ob.name(arg[s]), where name is in meth field
%
elseif isfield(st.meth, name)
	fun = st.meth.(name); % function handle for method
	if length(args) == 1
%		if isfreemat || nargout(fun) % freemat: does not like nargout()
		if isfreemat || nargout % freemat: does not like nargout()
%			out = fun(ob);
			varargout = wrap_vararg(nargout, 'fun(ob);', st, ob, ...
				args, fun);
		else
			fun(ob);
		end
	elseif length(args) == 2
%		if isfreemat || nargout(fun) % freemat: does not like nargout()
		if isfreemat || nargout % freemat: does not like nargout()
%			out = fun(ob, args(2).subs{:});
			varargout = wrap_vararg(nargout, 'fun(ob, args(2).subs{:});', ...
					st, ob, args, fun );
		else
			fun(ob, args(2).subs{:});
		end
	else
		printm('unknown subs')
		keyboard
	end

else
	disp(st)
	fail(['unknown field name: ' name])
end

end % subsref


% wrap_vararg()
% return cell array with each output
function out = wrap_vararg(nout, script, st, ob, args, fun)

switch nout
case 0
	eval(script);
	out = {};
case 1
	v = eval(script);
	out = {v};
case 2
	[v1, v2] = eval(script);
	out = {v1, v2};

otherwise % general case using eval
	head = sprintf('v%d, ', 1:(nout-1));
	full = sprintf('[ %s vlast ] = %s', head, script);
	out = cell(nout,1);
	for i=1:(nout-1)
		out{i} = eval(sprintf('v%d', i));
	end
	out{end} = vlast;
end

end % wrap_vararg
