  function ob = strum(data, methods, varargin)
%|function ob = strum(data, methods, [options])
%|
%| Construct a "strum" object.  This acts like an ordinary matlab "struct"
%| except that it also can contain user-defined "methods."
%| For example, if x is a fieldname of the structure, then ob.x will return
%| the corresponding value, just like an ordinary structure.
%| But x might instead be a handle to a function, say @xfun,
%| and invoking ob.x will return xfun(st).
%| Furthermore, those functions can accept arguments, i.e., ob.x(arg[s]),
%| which will invoke xfun(st, arg[s]).
%| The invoked functions are passed the structure (really, the entire object)
%| so that the invoked functions can access all the "not so private" data
%| in the structure.
%|
%| in
%|	data	struct	initial structure, i.e., 'data'
%|	methods	cell	{'name1', @handle1, 'name2', @handle2, ...}
%|		or	{'name1', @handle1, 'doc1; ...} (documented version)
%|
%| option
%|	'base'	strum	add data / methods to an existing strum
%|
%| out
%|	ob		strum object
%|
%| Copyright 2006-1-19, Jeff Fessler, University of Michigan

if nargin == 1 && streq(data, 'test'), run_mfile_local strum_test; return, end

ob.caller = caller_name(1);
ob.meth = struct;
ob.data = struct;
ob.docs = {};

if nargin == 0	% required by Mathworks
	ob = class(ob, mfilename);
return
end

if nargin < 2, ir_usage, end

if nargin > 2 % base object provided
	arg.base = struct([]);
	arg = vararg_pair(arg, varargin);
	if isempty(arg.base), keyboard, fail('base required'), end

	ob = struct(arg.base); % initialize with base data/methods
	ob.caller = [ob.caller '/' caller_name(1)]; % append caller

	% add new data to base object
	if ~isstruct(data), fail('data must be a struct'), end
	names = fieldnames(data);
	for ii=1:length(names)
		name = names{ii};
		if isfield(ob.data, name) || isfield(ob.meth, name)
			warn(['name reuse: ', name])
		end
		ob.data.(name) = data.(name); % set or over-ride with new
	end

	method = ob.meth;
	has_base = true;

else
	ob.data = data;
	method = struct;%([]); % no methods
	has_base = false;
end


%if isempty(methods)
%	warning 'empty methods: why use strum instead of struct?'
%end

if 3 == size(methods,2) % name,handle,doc triples

	if has_base && isempty(ob.docs)
		error 'docs must be consistent between base and new'
	end

	for ii=1:size(methods,1)
		name = methods{ii,1};
		handle = methods{ii,2};
		doc = methods{ii,3};
		if isfield(ob.data, name) || isfield(method, name)
			warn('name reuse: %s', name)
		end
		method.(name) = handle;
		ob.docs{end+1} = doc;

		strum_check_nargout(handle)
	end


elseif 0 == rem(length(methods),2) % name,handle pairs

	if has_base && ~isempty(ob.docs)
		error 'docs must be consistent between base and new'
	end

	for ii=1:2:length(methods)
		name = methods{ii};
		handle = methods{ii+1};
		if isfield(ob.data, name) || isfield(method, name)
			warn('name reuse: %s', name)
		end
		method.(name) = handle;
		strum_check_nargout(handle)
	end

else
	error 'need name,value pairs, or name,value,doc; triples'
end

ob.meth = method;

ob = class(ob, mfilename);


% strum_check_nargout()
% todo: warn if not exactly one output?
function strum_check_nargout(handle)
nout = ir_nargout(handle);
if nout == -1
	return % in matlab, -1 means "variable number of arguments"
end

if nout ~= 1
	pr handle
	keyboard
end
