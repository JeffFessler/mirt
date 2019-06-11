 function plot_cell(xc, yc, types, varargin)
%function plot_cell(xc, yc, types, varargin)
% plot several "y" curves against x curve(s) 
% in:
%	xc	{N} or 1	cell array of x points
%				or a single vector if all x's are the same
%	yc	{N}		cell array of y points
%	types	{N}		cell array of line types
% options:
%	'use', function_handle

if nargin == 1 && streq(xc, 'test'), plot_cell_test, return, end

yc = {yc{:}};
if ~isvar('types') || isempty(types)
	types = cell(size(yc));
end

handle = @plot;
while length(varargin)
	arg = varargin{1};
	if streq(arg, 'use')
		handle = varargin{2};
		varargin = {varargin{3:end}};
	else
		error 'unknown option'
	end
end

% build arguments to plot()
arg = {};
for ii=1:length(yc)
	if iscell(xc)
		x = xc{ii}
	else
		x = xc;
	end
	if isempty(types{ii})
		types{ii} = '-';
	end
	arg = {arg{:}, x, yc{ii}, types{ii}};
end
feval(handle, arg{:})
title(sprintf('showing %d plots', length(yc)))


function plot_cell_test
x = linspace(0,3,501);
y = {sin(2*pi*x), cos(2*pi*x), x};
if im
	clf
	subplot(121),
	plot_cell(x, y, {'r-', 'g--', 'm:'})
	subplot(122),
	plot_cell(10.^x, y, {'r-', 'g--', 'm:'}, 'use', @semilogx)
end
