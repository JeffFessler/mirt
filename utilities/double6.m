 function y = double6(x)
%function y = double6(x)
% optional conversion from single to double.
% options:
%	'small'		(default)
%			'double' if matlab 6 or before, otherwise 'single'
%	'double'	force conversion to double always
%	'single'	force  conversion to single always
%	'state'		return current state

persistent how	% stores state
if ~isvar('how') || isempty(how)
	how = 'small';
end

if nargin < 1, ir_usage, end

% query mode
if streq(x, 'state')
	y = how;
return
end

% set options
if ischar(x)
	if streq(x, 'double')
		how = 'double';
	elseif streq(x, 'single')
		how = 'single';
	elseif streq(x, 'small')
		how = 'small';
	else
		error(sprintf('double6 option "%s" is unknown', x))
	end
return
end

% convert
if streq(how, 'double')
	y = double(x);
elseif streq(how, 'single')
	y = single(x);
elseif streq(how, 'small')
	if is_pre_v7
		y = double(x);
	else
		y = single(x);
	end
else
	error 'unknown'
end
