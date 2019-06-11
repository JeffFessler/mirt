 function run_mfile_local(arg, varargin)
%function run_mfile_local(arg)
%|
%| run an mfile in a local environment so that workspace variables
%| are untouched.  useful for tests.
%| arg can be a cell with many mfiles to run
%| options
%|	'draw'	1|0	1 to draw after each test (default: 0)
%|	'pause'	1|0	1 to pause after each test (default: 0)
%|	'abort'	1|0	1 to abort if any test fails (default: 0)
%|
%| Copyright 2005-6-21, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end

opt.draw = false;
opt.pause = false;
opt.abort = false;
opt = vararg_pair(opt, varargin);

% track which tests open a figure even though they should not if im disabled
check_fig = ~figure_opened && ~im;
bad_fig = {};

test_bad = {};
test_good = {};
if iscell(arg)

	tmp = cpu('etic');

	for ii=1:length(arg)
		printf('\n\nTesting: %s\n\n', arg{ii})
		try
			time1 = cpu('etic');
			run_mfile_local(arg{ii})
			time1 = cpu(time1);
			test_good{end+1} = sprintf('%s (%.2f)', arg{ii}, time1);
		catch
			time1 = cpu(time1);
			test_bad{end+1} = sprintf('%s (%.2f)', arg{ii}, time1);
			if opt.abort
				fail('test "%s" failed', arg{ii})
			end
		end
		if opt.draw, drawnow, end
		if opt.pause, printm 'pausing', pause, end

		drawnow;
		if check_fig && figure_opened
			bad_fig{end+1} = arg{ii};
			close all
		end
	end

	if length(test_good)
		printf('\n----\n')
		printm('in "%s" the following test(s) passed:\n', caller_name)
		disp(char(test_good))
	end

	if length(bad_fig)
		warn('in %s: the following test(s) had figure issues:', caller_name)
		disp(char(bad_fig))
	end

	cpu(tmp, [caller_name ' time'])

	if length(test_bad)
		warn('in "%s" the following test(s) failed:', caller_name)
		disp(char(test_bad))
		error 'the above tests failed'
	else
		printm('in "%s" all %d tests passed!', caller_name, length(arg))
	end

else
	eval(arg)
end


% try to determine if any figure window is opened
function out = figure_opened
if isfreemat
	out = true; % freemat will not tell, so assume yes
else
	out = ~isempty(get(0, 'children')); % matlab knows
end
