 function out = fatrix2_subsref_dot(st, args)
%function out = fatrix2_subsref_dot(st, args)
% handle subscript references like ob.ref or ob.ref(arg[s])
% This will called from ../subsref with struct(ob)
% and it will have arg(1).type = '.'
% Based on @strum/subsref
% Copyright 2006-1-19, Jeff Fessler, University of Michigan

%if args(1).type ~= '.'
%	error 'only ob.data and ob.meth(args) done'
%end

name = args(1).subs;
%printm('here, name=%s', name)

switch name
case 'imask_array'
	out = fatrix2_imask_array(st);
	return

case 'omask_array'
	out = fatrix2_omask_array(st);
	return

end


% ob.name
if isfield(st, name)
	if length(args) == 1
		out = st.(name);
	else
		out = subsref(st.(name), args(2:end));
	end

% ob.name or ob.name(arg[s]), where name is in meth field
elseif isfield(st.meth, name)
	fun = st.meth.(name); % function handle for method
	if length(args) == 1
		out = fun(st.arg);
	elseif length(args) == 2
		out = fun(st.arg, args(2).subs{:});
	else
		printm('unknown subs')
		keyboard
	end

% ob.name, where name is in data field 'arg' (last resort)
elseif isfield(st.arg, name)
	out = st.arg.(name);
	if length(args) > 1

		% ob.name{:} does not work properly because there seems
		% to be no way to return an arbitrary number of arguments
		% (in a list) to the caller.
		% This is unfortunate because
		% a.b = {1,2} and a.b{:} works for ordinary structures.
% todo: see this solution involving numel()
% http://www.mathworks.com/support/solutions/data/1-1ABOD.html?solution=1-1ABOD
		if iscell(out) && ...
		(isequal(args(2).subs, {':'}) || length(args(2).subs{1}) > 1)
			warn 'instead of strum.data{:} etc., use two steps:'
			warn 'tmp = strum.data; tmp{:} or tmp{1:end} etc.'
			fail 'strum.data{:} etc. unsupported thanks to matlab'
		end

		% fix: strum1.strum2.plot() does not work here
		out = subsref(out, args(2:end));
	end

else
	disp(st)
	fail(['unknown field name: ' name])
end
