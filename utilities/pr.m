 function pr(varargin)
%function pr(command)
%|
%| print a message with the calling routine's name,
%| the argument that is evaluated, and the value thereof.
%|
%| circa 2000, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end

arg = [varargin{:}]; % handle cases like 'pr x + y'

tmp = evalin('caller', arg);

% scalar
[name line] = caller_name;
if ~isempty(line) && line ~= 0
	name = [name sprintf(' %d', line)];
end

if isscalar(tmp) && isnumeric(tmp) && isreal(tmp)
	printf('%s: %s = %g', name, arg, tmp)

% short vector
elseif min(size(tmp)) == 1 && length(tmp) <= 3 && isnumeric(tmp) && isreal(tmp)
	printf('%s: %s = %s', name, arg, num2str(tmp(:)', ' %g'))

% string
elseif ischar(tmp)
	printf('%s: %s = %s', name, arg, tmp)

% struct
elseif isstruct(tmp)
	printf('%s: %s = [struct]', name, arg)
	ir_display_struct(tmp)

% fatrix2
elseif isa(tmp, 'fatrix2') || isa(tmp, 'Fatrix')
	printf('%s: %s = [%s]', name, arg, class(tmp))
	ir_display_struct(struct(tmp))

% array, etc.
else
	printf('%s: %s =', name, arg)
	disp(tmp)
end
