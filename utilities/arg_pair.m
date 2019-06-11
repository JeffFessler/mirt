  function arg = arg_pair(varargin)
%|function arg = arg_pair(varargin)
%|
%| construct AsciiArg pairs from input arguments
%| for building arguments akin to .dsc file 
%| example: arg = arg_pair('system', 2, 'nx', 128, ...)
%|
%| Copyright May 2000, Jeff Fessler, University of Michigan

if nargin == 1 && streq(varargin{1}, 'test'), arg_pair_test, return, end

if nargin == 1 && isstruct(varargin{1})
	arg = arg_pair_struct(varargin{1});
return
end

if nargin < 2, ir_usage, end

% see if first argument is already a char array; if so, then augment it.
arg = varargin{1};
if ischar(arg) && size(arg,1) > 1
	varargin = {varargin{2:end}};
else
	arg = [];
end

if rem(length(varargin),2), error 'even # of arguments required', end

for ii=1:2:length(varargin)
	a = varargin{ii};
	b = varargin{ii+1};
	if ~ischar(a), error a, end
	if ischar(b)
		arg = strvcat(arg, [a ' ' b]);
	elseif isscalar(b)
		arg = strvcat(arg, sprintf('%s %g', a, b));
	else
		pr size(b)
		fail('nonscalar argument class "%s"', class(b))
	end
end
arg = remove_spaces(arg);


%
% arg_pair_struct()
%
function arg = arg_pair_struct(st)
names = fieldnames(st);
arg = [];
for ii=1:length(names)
	a = names{ii};
	b = st.(a);
	if ischar(b)
		arg = strvcat(arg, [a ' ' b]);
	else
		arg = strvcat(arg, sprintf('%s %g', a, b));
	end
end
arg = remove_spaces(arg);


function arg_pair_test
a1 = arg_pair('a1', 1, 'a2', 'v2');
a2 = arg_pair(a1, 'liii', 3.3);
a0 = remove_spaces(strvcat('a1 1', 'a2 v2', 'liii 3.3'));
jf_equal(a0, a2)
