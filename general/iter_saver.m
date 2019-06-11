  function isave = iter_saver(isave, niter)
%|function isave = iter_saver(isave, niter)
%| process 'isave' and 'niter' options of iterative algorithms
%| supporting char options 'all' and 'last' (default if empty)

persistent warned
if ~isvar('warned') || isempty(warned)
	warned = false;
end

if isempty(isave)
	if ~warned
		warn 'using default isave "last"'
		warned = 1;
	end
	isave = 'last';
end

if streq(isave, 'last')
	isave = niter;
elseif streq(isave, 'all')
	isave = 0:niter;
end

if any(isave < 0) || any(isave > niter)
	pr isave
	pr niter
	fail('bad isave')
end
