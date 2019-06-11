 function out = jf(varargin)
%function out = jf(varargin)
%| various personalization routines
%| 'mv var_old var_new'
%| 'ncore'	# of cores this machine has
%|		add optional 2nd argument to over-ride, e.g., jf('ncore', 4)
%| 'whos'
%| 'clf'
%| 'nobar'	preclude toolbar and menubar from figures to render faster!
%| 'nomex'	remove mex files from path
%| 'isum'	1 if UM, 0 else
%| 'path'
%| 'off'
%| 'on'
%| 'plc'		set up subplot with clf first
%| {'pl', 'pl-tight'}	set up subplot, don't clear first
%| 'sub'
%| 'title_no_tex'
%| 'test'

if ~nargin, ir_usage, end

%
% handle states
%
persistent state
if ~isvar('state') || isempty(state)
	state = jf_reset;
end

switch varargin{1}

case 'mv' % rename a variable, unix style
	% todo: make sure not renaming a field of a structure!
	if ~isempty(findstr(varargin{2}, '.')) ...
	|| ~isempty(findstr(varargin{2}, '.'))
		fail('struct field renames not implemented yet, but could be')
	end
	tmp = [varargin{3} ' = ' varargin{2} ';'];
	evalin('base', tmp);
	tmp = ['clear ' varargin{2} ';'];
	evalin('base', tmp);

case 'ncore'
	persistent ncore
	if length(varargin) == 2
		ncore = varargin{2};
		if ischar(ncore) % for this syntax: jf ncore 8
			if streq(ncore, 'reset')
				ncore = [];
				ncore = jf('ncore');
			else
				ncore = str2num(ncore);
			end
		end
	end

	if isvar('ncore') && ~isempty(ncore)
		out = ncore; % return previously determined value
	return
	end

	if ismac % only try this on a mac!
		try
			tmp = os_run('/usr/sbin/sysctl hw.ncpu');
			out = sscanf(tmp, 'hw.ncpu: %d');
			ncore = out;
			return
		catch
			persistent warned1
			if ~isvar('warned1') || isempty(warned1)
				warned1 = 1;
				warn 'sysctl did not work'
			end
		end
	end

	try
		maxNumCompThreads('automatic');
		out = maxNumCompThreads;
		if (out > 40)
			warn 'limiting to 40 threads (over-ride with caution)'
			out = 40;
		end
		ncore = out;
		return
	catch
		persistent warned2
		if ~isvar('warned2') || isempty(warned2)
			warned2 = 1;
			warn 'maxNumCompThreads failed; reverting to 1 thread'
		end
		out = 1;
		ncore = out;
		return
	end

case 'whos'
	st = evalin('caller', 'whos');
	byte = {st.bytes};
	byte = cell2mat(byte);
	byte = sum(byte);
	printm('total bytes = %d = %g Gb', byte, byte/1024^3)

case 'clf'
	if state.display, clf, end

case 'isum'
	out = exist('dd_ge2_mex') == 3;

case 'nobar'
	set(0, 'DefaultFigureToolbar', 'none')
	set(0, 'DefaultFigureMenubar', 'none')

case 'notoolbox' % remove all matlab extra toolboxes from path for testing
	jf_rm_toolbox

case 'nomex'
	tmp = path;
	tmp = strsplit(tmp, pathsep);
	for ii=1:length(tmp)
		if ~isempty(strfind(tmp{ii}, 'mex/v7'))
			printm('removing from path: "%s"', tmp{ii})
			rmpath(tmp{ii})
		end
	end
	printm('done removing mex dirs from path')

case 'path'
	path_jf

case 'off'
	state.display = false;

case 'on'
	state.display = true;

case 'plc' % set up subplot with clf first
	jf clf
	jf('pl', varargin{2:end});

case {'pl', 'pl-tight'} % set up subplot, don't clear first
	if nargin == 3
		state.sub_m = ensure_num(varargin{2});
		state.sub_n = ensure_num(varargin{3});
	elseif nargin == 1
		state.sub_m = [];
		state.sub_n = [];
	else
		error 'bad pl usage'
	end
	state.pl_tight = streq(varargin{1}, 'pl-tight');

case 'sub'
	arg = ensure_num(varargin{2});
	if arg > 99 % reset
		state.sub_m = [];
		state.sub_n = [];
		end
	if isempty(state.sub_m)
		if arg >= 111
			subplot(arg)
		else
			printm('ignoring subplot %d', arg)
		end
	else
		jf_subplot(state, arg)
	end

case 'test'
	jf_test

case 'title_no_tex'
	set(get(gca, 'title'), 'interpreter', 'none')

otherwise
	fail('unknown arg %s', varargin{1})
end

% split long string using
function out = strsplit(in, sep)
in = [in sep]; % add : at end
isep = strfind(in, sep);
for ii=1:(length(isep)-1)
	out{ii} = in((1+isep(ii)):(isep(ii+1)-1));
end


%
% jf_subplot()
%
function jf_subplot(state, num)
if state.display
	if state.pl_tight
		num = num - 1;
		x = 1 / state.sub_n;
		y = 1 / state.sub_m;
		ny = floor(num / state.sub_n);
		nx = num - ny * state.sub_n;
		ny = state.sub_m - ny - 1;
		subplot('position', [nx*x ny*y x y])
	else
		subplot(state.sub_m, state.sub_n, num)
	end
end


function x = ensure_num(x)
if ischar(x), x = str2num(x); end

%
% jf_state()
%
function state = jf_reset
state.sub_m = []; % for subplot
state.sub_n = [];
state.pl_tight = false;
state.display = true;


function jf_test
jf plc 2 3
jf sub 1, plot(rand(3))
