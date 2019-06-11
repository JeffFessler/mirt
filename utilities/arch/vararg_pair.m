 function [opt, extra] = vararg_pair(opt, varargs, varargin)
%function [opt, extra] = vararg_pair(opt, varargs, [options])
%|
%| Process name / value pairs, replacing the "default" field values
%| of the opt structure with the user-specified values.
%| This allows flexible argument order and "named arguments" somewhat like IDL.
%| This newer version allows option names to differ from the field names.
%| If two output arguments are given, then any "extra" name value pairs
%| that are not associated with the fields in structure opt will be
%| returned as name/value pairs in the "extra" output.
%|
%| in
%|	opt	struct	opt.name1 opt.name2 etc., containing default values
%|	varargs	{'name1', 'value1', ... } as a cell array
%|
%| options
%|	'subs'	Nx2 cell	{'arg_name1', 'struct_name1'; ...}
%|		short argument name followed by longer structure name
%|	'allow_new'	0|1	allow new structure elements? (default: 0)
%|
%| out
%|	opt		opt.(name) = value
%|	extra		{'name_i', 'value_i', ...} all "unused" pairs
%|
%| Copyright 2005-6-18, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(opt, 'test'), vararg_pair_test, clear, return, end

if isempty(varargs) && numel(varargin) == 0
	extra = {};
	return
end

base = [caller_name ' arg: ']; % for printing

local.subs = {};
local.allow_new = 0;
local = vararg_pair_do(local, varargin, {}, 0, 0);

allow_extra = nargout > 1;
[opt extra] = vararg_pair_do(opt, varargs, local.subs, local.allow_new, ...
	allow_extra, base);

% Structure elements of opt with value "NaN" are checked as being required.
% NO, this just adds overhead.  Removed 2007-2-9.
% vararg_pair_check(opt);

%
% vararg_pair_do()
%
function [opt, extra] = ...
	vararg_pair_do(opt, varargs, subs, allow_new, allow_extra, base)

if allow_extra && allow_new
	error 'only one of "new" or "extra" is sensible for vararg_pair'
end

npair = floor(length(varargs) / 2);
if 2*npair ~= length(varargs), fail('need names and values in pairs'), end
args = varargs(1:2:end);
vals = varargs(2:2:end);

if ~isempty(subs)
	if size(subs,2) ~= 2, fail('subs must be Nx2'), end
%	subs1 = strvcat(subs(:,1));
	subs1 = char(subs(:,1)); % strmatch doc says fastest for char array!
	subs2 = subs(:,2);
	% todo: apply 'subs' all at once here?
end

extra = {};
for ii=1:npair
	arg = args{ii};
	val = vals{ii};

	if ~ischar(arg)
		fail('unknown option of class %s', class(arg))
	end

	if ~isempty(subs)
		arg = vararg_pair_sub(arg, subs1, subs2);
	end

	[opt ex] = vararg_pair_put(opt, arg, val, allow_new, allow_extra);
	if allow_extra && ~isempty(ex)
		extra = [extra(:); ex(:)]'; % faster than {extra{:}, ex{:}};
	end

	if isfield(opt, 'chat') && opt.chat && isempty(ex)
		show_pair(arg, val, base)
	end
end


%
% vararg_pair_sub()
% substitute option name for field name, if a Nx2 cell array is provided.
%
function arg = vararg_pair_sub(arg, subs1, subs2)
ii = strmatch(arg, subs1, 'exact');
if isempty(ii)
	return
elseif length(ii) == 1
	arg = subs2{ii};
else
	printm([subs1(ii) ' ' subs2{ii}])
	error 'double match?  only one substitution allowed'
end


%
% vararg_pair_put()
% opt.(arg) = value
% with extra effort to handle opt.(arg1.arg2) = value
% because matlab does not handle that case by itself, except via subsagn
%
function [opt, extra] = vararg_pair_put(opt, arg, value, allow_new, allow_extra)

extra = {};
if isfield(opt, arg) % simple case of name / value pair
	opt.(arg) = value;
	return
end

idot = strfind(arg, '.');
if isempty(idot)
	if allow_new
		opt.(arg) = value;
		return
	elseif allow_extra
		extra = {arg, value};
		return
	else
		fail('unknown option name "%s"', arg)
	end
end

% tricky case of a.b.c
if length(idot) > 1, error 'a.b.c.d not done', end

arg1 = arg([1:(idot-1)]);
arg2 = arg([(idot+1):end]);
if ~isfield(opt, arg1)
	fail('unknown option1 name "%s"', arg1)
end
if ~isfield(opt.(arg1), arg2)
	fail('unknown option2 name %s', arg2)
end
s = struct('type', {'.', '.'}, 'subs', ...
	{ arg([1:(idot-1)]), arg([(idot+1):end]) });
try
	opt = subsasgn(opt, s, value);
catch
	printf(lasterr)
	fail('subsasgn? unknown option name %s', arg)
end


%
% vararg_pair_check()
% check that the user supplied the required arguments,
% which are signified by 'nan' values.
%
function vararg_pair_check(opt)
names = fieldnames(opt);
ok = 1;
for ii=1:length(names)
	t = opt.(names{ii});
	if isnumeric(t) && any(isnan(t(:)))
		printf('Required option "%s" missing', names{ii});
		ok = 0;
	end
end
if ~ok, error 'missing required "option(s)"', end


%
% show_pair()
%
function show_pair(name, value, base)

pri = @(varargin) printf([base varargin{1}], varargin{2:end});

if isnumeric(value) || islogical(value)
	if max(size(value)) == 1	% only print scalars
		if ~streq(name, 'chat')
			pri('%s = %g', name, value)
		end
	else	% otherwise print sum (e.g., for mask)
		pri('sum(%s) = %g', name, sum(value(:)))
	end

elseif ischar(value)
	pri('%s = %s', name, value)

elseif isa(value, 'cell')
	pri('%s {cell}', name)

elseif isa(value, 'function_handle')
	pri('%s = @%s', name, func2str(value))

elseif isstruct(value)	% show first structure element, if string
	fields = fieldnames(value);
%	name = getfield(arg, fields{1});
	name1 = value.(fields{1});
	if ischar(name1)
		pri('%s.%s = %s', name, fields{1}, name1)
	end

else
	warn('display type "%s" not implemented', class(value))
end


%
% test routine
%
function vararg_pair_test
args = {'a', 1, 'b', 2, 'c.d', 3, 'g', 9, 'l', 5, 'll', 6};%, 'c.f', 4};
opt.a = 0;
opt.b = 0;
opt.f = 0;
opt.long1 = 0;
opt.long2 = 0;
%opt.g = nan;
opt.c.d = 0;
opt.c.e = 0;
sub = {'g', 'f'; 'l', 'long1'; 'll', 'long2'};
opt = vararg_pair(opt, args, 'subs', sub);
jf_equal(opt.a, 1)
jf_equal(opt.b, 2)
jf_equal(opt.c.d, 3)
jf_equal(opt.c.e, 0)
jf_equal(opt.f, 9)
jf_equal(opt.long1, 5)
jf_equal(opt.long2, 6)

opt = vararg_pair(opt, {'z', 'new'}, 'allow_new', 1);
if ~isequal(opt.z, 'new'), fail('bug'), end
try % typo in name should yield error
	opt = vararg_pair(opt, {'verify_bad_test!', 'new'}, 'allow_new', 0)
catch
	return
end
error 'should not get here'
