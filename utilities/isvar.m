 function tf = isvar(name, varargin)
%function tf = isvar(name, varargin)
%|
%| Cetermine if "name" is a variable in the caller's workspace.
%|
%| If argument is of the form 'name.field' or 'name.field1.field2' etc.
%| then this uses isfield(st, 'field') recursively as needed>
%|
%| To have isvar always return false when called from a script, use:
%|	isvar set false_if_script
%| To reset to usual mode, use:
%|	isvar set normal
%| To querty to current mode, use:
%|	isvar set query
%|
%| Copyright 2000-01-01, Jeff Fessler, University of Michigan
%| modified 2010-04-21 to use 'exist' and 'isfield'


if nargin < 1, ir_usage, end

if ~nargout && streq(name, 'test'), isvar_test, return, end

persistent State
if isempty(State)
	State = 'normal';
end

if numel(varargin) >= 1
	if nargout || ~streq(name, 'set') || numel(varargin) > 1
		fail('bad "isvar set ..." usage')
	end
	switch varargin{1}
	case 'false_if_script'
		State = 'false_if_script';
	case 'query'
		printm('isvar State = %s', State)
	case 'normal'
		State = 'normal';
	otherwise
		fail('isvar set "%s" unsupported')
	end
return
end

if ~numel(dbstack) && streq(State, 'false_if_script')
	tf = false;
end

dots = strfind(name, '.'); % look for any field references

if isempty(dots)
	base = name;
else
	base = name(1:dots(1)-1);
	tail = name((dots(1)+1):end);
end

str = sprintf('exist(''%s'', ''var'');', base);
tf = evalin('caller', str);

while tf && ~isempty(dots)
	if length(dots) == 1
		str = sprintf('isfield(%s, ''%s'');', base, tail);
		tf = tf & evalin('caller', str);
		return
	else
		dots = dots(2:end) - dots(1);
		next = tail(1:dots(1)-1);
		str = sprintf('isfield(%s, ''%s'');', base, next);
		tf = tf & evalin('caller', str);
		base = [base '.' next]; %#ok<AGROW>
		tail = tail((dots(1)+1):end);
	end
end


% isvar_test
function isvar_test
var1 = 0;
jf_equal(isvar('var1'), 1)
jf_equal(isvar('var2'), 0)

if im
	isvar set query
	isvar set false_if_script
	isvar set query
	jf_equal(isvar('var1'), 1)
	jf_equal(isvar('var2'), 0)
	isvar set normal
	isvar set query
end
