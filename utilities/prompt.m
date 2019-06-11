 function out = prompt(arg)
%function out = prompt(arg)
%|
%| In MIRT example scripts, there are various possible modes of execution
%| that are controlled by this m-file, as follows:
%|
%|	pause			tell user to press enter key to continue
%|				(or offer other options at command line)
%|	return | stop		force user to re-run to continue (rare)
%|	run | continue		run to end of script without stopping
%|	draw			'drawnow' then continue
%|	[k]eyboard		dbstop after prompt (then dbquit, dbclear all) 
%|	msgbox			use msgbox (gui) rather than prompt
%|	mode			return current prompt mode as "out"
%|
%| If a script is running as Matlab "Live Script" then this is ignored.
%|
%| Copyright 2001-8-30, Jeff Fessler, University of Michigan

persistent Prompt % stores state
if ~isvar('Prompt') || isempty(Prompt)
	Prompt = 'pause';
end

% query mode
if ~nargin && nargout
	out = Prompt;
return
end

% set mode (or give help)
if nargin && ~nargout
	switch arg
	case 'help'
		help(mfilename)
	case 'mode'
		if nargout
			out = Prompt;
		else
			printm(Prompt)
		end
		return

	case 'test'
		ir_prompt_test
	otherwise
		Prompt = arg;
	end
return
end


if ir_is_live_script
	return
end

% fix: we really need a "return all" call here!
% dbquit?
switch Prompt
case {'return', 'stop'}
%	disp ' '
%	disp 'returning'
%	return	% does not work!
%	evalin('base', 'return')
%	evalin('caller', 'return')
	error 'fake error to effect a "return all"'

case 'msgbox'
	drawnow
	uiwait(msgbox('MIRT: click ok to continue'))

case 'draw'
	drawnow	% do nothing, just continue

case {'run', 'continue'}
%	disp ' ' % do nothing, just continue

case 'pause'
	[name, line] = caller_name;
	if isempty(name)
		preface = [];
	else
		preface = sprintf('%s %d: ', name, line);
	end

	what = 'enter to continue (or [r]un [d]raw [q]uit [n]odraw): ';
	ans = input([preface  what], 's');
	if streq(ans, 'r', 1)
		Prompt = 'run';
	elseif streq(ans, 'd', 1)
		Prompt = 'draw';
	elseif streq(ans, 'n', 1)
		Prompt = 'run';
		im off
		close all
	elseif streq(ans, 'q', 1)
		fail('quitting')
	elseif streq(ans, 'k', 1) % keyboard mode after prompt
	%	printm('experimental keyboard mode')
		pr name
		pr line
		tmp = sprintf('dbstop in %s at %d', name, line+1);
		eval(tmp)
	elseif ~isempty(ans)
		Prompt = ans;
	end

otherwise
	warn('unknown Prompt mode "%s"; resetting to ''pause''', Prompt)
	Prompt = 'pause';
	prompt
end


function ir_prompt_test
x = 5;
prompt
y = 7;
z = 8;
