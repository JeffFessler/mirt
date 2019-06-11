 function type_diff = Cdiff1_tune(isize, varargin)
%function type_diff = Cdiff1_tune(isize, [options])
%|
%| "Auto-tune" function for Cdiff1() that selects the fastest overall method
%|
%| in
%|	isize	[]		vector of object dimensions (N), e.g., [64 64]
%|
%| option:
%|	'offset' [int]		offset (default: 1)
%|	'order' 1 or 2		1st- or 2nd-order differences.  (default: 1)
%|	'do_for1'		include slow 'for1' method (default: 0)
%|	'do_ind'		include slow 'ind' method (default: 0)
%|
%| out
%|	type_diff		fastest 'type_diff' option for Cdiff1
%|				(typically the result will be 'mex')
%|
%| Copyright 2013-05-08, Jeff Fessler, University of Michigan

if nargin == 1 && streq(isize, 'test'), Cdiff1_tune_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

arg.offset = [1];
arg.order = 1;
arg.do_for1 = false;
arg.do_ind = false;
arg = vararg_pair(arg, varargin); % parse options

% list in fastest to slowest order (on ire):
list = {'diff', 'convn', 'circshift', 'mex', 'imfilter', 'sparse'};
if arg.do_ind
	list{end+1} = 'ind'; % slow, but if user insists
end
if arg.do_for1
	list{end+1} = 'for1'; % slow, but if user insists
end

%{
switch numel(isize)
case 1
	ig = image_geom('nx', isize(1), 'ny', 1, 'dx', 1);
case 2
	ig = image_geom('nx', isize(1), 'ny', isize(2), 'dx', 1);
case 3
	ig = image_geom('nx', isize(1), 'ny', isize(2), 'nz', isize(3), ...
		'dx', 1, 'dz', 1);
end
%}

x = 7*ones([isize 1], 'single');

ntype = numel(list);
time_forw = inf(ntype, 1);
time_back = inf(ntype, 1);
time_warm = inf(ntype, 1);
lab = cell(ntype, 1);

for ii = 1:ntype
	ctype = list{ii};

	if streq(ctype, 'diff') && (arg.order ~= 1 || ~isequal(arg.offset, 1))
		lab{ii} = ['(' ctype ')'];
		continue % 2nd order and other offsets not done for 'diff'
	end

	C = Cdiff1(isize, 'order', arg.order, 'offset', arg.offset, ...
		'type_diff', ctype);

	if 1 % on ire, warm-up is pointless
		cpu etic
		tmp = C * x; % warm up
		time_warm(ii) = cpu('etoc');
	else
		time_warm(ii) = 0;
	end

	cpu etic
	tmp = C * x;
	time_forw(ii) = cpu('etoc');

	cpu etic
	C' * tmp;
	time_back(ii) = cpu('etoc');

	if isfield(C.arg, 'type_diff')
		lab{ii} = C.arg.type_diff;
	else
		lab{ii} = C.caller;
	end
end

time_both = time_forw + time_back;
[time_both ii] = sort(time_both);
time_forw = time_forw(ii);
time_back = time_back(ii);
time_warm = time_warm(ii);
lab = {lab{ii}};

if 1
	format = 'Cdiff1 %10s\t%5.3f\t%5.3f\t%5.3f\t%5.3f';
	for ii = 1:ntype
		printm(format, lab{ii}, ...
		time_warm(ii), time_forw(ii), time_back(ii), time_both(ii))
	end
end
type_diff = lab{1};


% Cdiff1_tune_test
function Cdiff1_tune_test
for order = 1:2
	isize = [2^8 2^8 2^7];

%	offset = [1]; % for 'diff'
	offset = [1 0 1]; % 3d offset

	printf ' ' % blank line
	printm('order=%d [%d %d %d] [%s] warm forw back both', ...
		order, isize(1), isize(2), isize(3), num2str(offset, '%d '))

	do_ind = ir_is_octave; % fast in octave, slow in matlab 2011a!
	out = Cdiff1_tune(isize, 'order', order, 'offset', offset, ...
		'do_ind', do_ind);
end
