  function ob = Gembed(varargin)
%|function ob = Gembed([args])
%| Construct Gembed object that embes a short vector into a longer vector.
%| with dimensions [(Nd)].
%|
%| option1:
%|	'list'	[N]		list of array entries
%|	'odim'	ndims
%|
%| option2:
%|	'samp'	logical [(Nd)]	sampling pattern of array entries
%|
%| out
%|	ob	[M N]	fatrix2 object, where M = numel(samp), N = length(list)
%|
%| The two versions are related by list = find(samp), odim = size(samp).
%|
%| Copyright 2010-12-01, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), Gembed_test, return, end

arg.samp = [];
arg.list = [];
arg.odim = [];
arg = vararg_pair(arg, varargin);

if ~isempty(arg.samp) && isempty(arg.list) && isempty(arg.odim)
	ob = Gembed_samp_build(arg.samp);

elseif isempty(arg.samp) && ~isempty(arg.list) && ~isempty(arg.odim)
	ob = Gembed_list_build(arg.odim, arg.list);

else
	fail('must give exactly one of samp, or list and odim')
end


% Gembed_list_build()
function ob = Gembed_list_build(odim, list)
arg.odim = odim;
arg.list = list;

ob = fatrix2( ...
	'accept1d', true, ... % trick: OK because of forw(forw()) in mtimes2
	'idim', numel(arg.list), ...
	'odim', arg.odim, ...
	'arg', arg, ...
	'forw', @Gembed_list_forw, ...
	'back', @Gembed_list_back, ...
	'caller', [mfilename '(list)']);
%	'mask', true(numel(arg.list),1), ...


% Gembed_list_forw(): y = A * x
% in
%	x	[N]
% out
%	y	[(M)]
%
function y = Gembed_list_forw(arg, x)

y = zeros([arg.odim 1], class(x)); % [(M)]
y(arg.list) = x;


% Gembed_list_back(): x = A' * y
% in
%	y	[(M)]
% out
%	x	[N]
%
function x = Gembed_list_back(arg, y)

x = y(arg.list);


% Gembed_samp_build()
function ob = Gembed_samp_build(samp)
arg.samp = samp;
arg.odim = size(samp);
if arg.odim(end) == 1
	arg.odim = arg.odim(1:end-1); % remove trailing '1'
end
if ~islogical(arg.samp), error 'samp must be logical', end

ob = fatrix2( ...
	'idim', sum(arg.samp(:)), ...
	'odim', arg.odim, ...
	'arg', arg, ...
	'forw', @Gembed_samp_forw, ...
	'back', @Gembed_samp_back, ...
	'caller', [mfilename '(list)']);


% Gembed_samp_forw()
function y = Gembed_samp_forw(arg, x)
y = zeros([arg.odim 1], class(x)); % [(M)]
y(arg.samp) = x;


% Gembed_samp_back()
function x = Gembed_samp_back(arg, y)
x = y(arg.samp);


% Gembed_test
function Gembed_test

odim = [5 2];
list = [2 1 9 7];
odim = prod(odim); % required because of fatrix2 full 1d ambiguity

for ii=1:2
	if ii==1 % type list
		A = Gembed('odim', odim, 'list', list);
	else % type samp
		samp = false([odim 1]);
		samp(list) = 1;
		A = Gembed('samp', samp);
	end

	fatrix2_tests(A)
	test_adjoint(A);
end
