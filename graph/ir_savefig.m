  function ir_savefig(varargin)
%|function ir_savefig([rootdir,] base, ['c'])
%|
%| save current figure to eps file rootdir/base.eps with prompting.
%|
%| options:
%|	'eps_c'	save to color .eps file without inverting!
%|	'cw'	color .eps file for use in white paper
%|	'c'	save both color and b/w files (for talks)
%|		name.epsc	color with black background for presentations
%|		name.eps	black text on white background for printing
%|	'-rmap'	reverse colormap for color figure (NOT DONE!)
%|	'tall'	orient tall
%|
%| note: there is another savefig version at mathworks that probably does
%| a better job of cropping unwanted white space from borders of figures:
%|
%| http://www.mathworks.com/matlabcentral/fileexchange/10889
%| http://www.mathworks.com/access/helpdesk/jhelp/techdoc/printing/printfi7.shtml
%|
%| Copyright 2002-2-1, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), savefig_test, return, end

if ir_is_octave
	arg_tight = {'-tight'};
	arg_tight = {}; % not working on ire 2014-01-22
else
	arg_tight = {};
end

if streq(varargin{1}, '-tight')
	if ir_is_octave
		arg_tight = {'-tight'};
	else
		; % ignore for matlab
	end
	varargin = {varargin{2:end}};
end

persistent Savefig
persistent Savefig_type
if isempty(Savefig), Savefig = true; end
if streq(varargin{1}, 'off'), Savefig = false; return, end
if streq(varargin{1}, 'on'), Savefig = true; return, end
if ~Savefig, disp 'savefig off', return, end

if isempty(Savefig_type), Savefig_type = 'eps'; end
if streq(varargin{1}, 'png'), Savefig_type = 'png'; return, end

rootdir = '.';
base = 'default';
do_color = false;
do_eps_c = false;
do_cw = false;

switch Savefig_type
case 'eps'
	suff = '.eps';
case 'png'
	suff = '.png';
otherwise
	fail('bad type %s', Savefig_type)
end

old_orient = orient;

while length(varargin)
	arg = varargin{1};

	if ischar(arg) && exist(arg, 'dir') && streq(base, 'default')
		rootdir = arg;

	elseif ischar(arg) && streq(arg, 'eps_c')
		do_eps_c = true;

	elseif ischar(arg) && streq(arg, 'c')
		do_color = true;

	elseif ischar(arg) && streq(arg, 'cw')
		do_cw = true;

	elseif ischar(arg) && streq(arg, 'tall')
		orient tall

	else
		if ~ischar('arg'), error string, end
		if ~streq(base, 'default')
			printm('WARN %s: ignoring bad names "%s" "%s"', ...
				[mfilename '.m'], base, arg)
			return
		end
		base = arg;
	end
	varargin = {varargin{2:end}};
end

base = [rootdir '/' base];

keepask = 1;
while (keepask)
	keepask = 0;
	name = [base suff];

	if streq(name, './tmp.eps')
		ans = 'y';
		warn('writing to "%s"', name)
	else
		if exist(name, 'file')
			prompt = sprintf('overwrite figure "%s"? y/a/d/[n]: ', name);
		else
			prompt = sprintf('print new figure "%s"? y/n: ', name);
		end
		ans = input(prompt, 's');
	end

	% alternate filename
	if streq(ans, 'a')
		base = input('enter alternate basename (no .eps): ', 's');
		keepask = 1;
	end
end

% difference between current plot and saved version
% fix: no-color only!
if streq(ans, 'd')
	printm('running test to compare to "%s"', name)
	print('t-savefig-test', '-deps')
	os_run(['diff ' 't-savefig-test.eps ' name])
	delete('t-savefig-test.eps')
return
end

if isempty(ans) || ~streq(ans, 'y'), return, end

if streq(Savefig_type, 'png')
	print(base, '-dpng')

elseif do_eps_c
	print(base, '-depsc')

elseif do_cw % color version but allowing 'inverthardcopy' for white paper
%	print([base '.eps'], '-depsc', arg_tight{:})
	print([base '.eps'], '-depsc2', arg_tight{:})

elseif ~do_color
	if ir_is_octave
		opt = {[base '.eps'], '-deps'};
		print(opt{:}, arg_tight{:})
	else
		print(base, '-deps')
	end

else
	% color version
	set(gcf, 'InvertHardCopy', 'off')
	print([base '.epsc'], '-depsc')

	% b/w version
	c = colormap;
	colormap(flipud(c))
	set(gcf, 'InvertHardCopy', 'on')
	print(base, '-deps')
	colormap(c)
end

printm('printed ok to "%s"', name)

orient(old_orient)

function savefig_test
jf plc 1 2
jf sub 1
im([rand(9), zeros(9)])
jf sub 2
plot(rand(3))
savefig cw test1
