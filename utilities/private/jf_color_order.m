 function out = jf_color_order(pn, varargin)
%function out = jf_color_order(pn, varargin)
% set color order to start with yellow then cyan etc., for black background

if ~length(varargin)
	out = jf_color_order_bw;
return
end

out = [];

switch varargin{1}
case 'revert'

	set(0, 'DefaultAxesColorOrder', 'default')
	return

case 'setup'
	set(0, 'DefaultAxesColorOrder', jf_color_order_bw)
%	set(gca, 'NextPlot', 'replacechildren')
%	set(gca, 'colororder', jf_color_order_bw)
	return

otherwise % 'test'

tmp = jf_color_order_bw;
N = nrow(tmp) * 2;

n = 1:N;
x = [n; n];
y = repmat([0 1]', [1, N]);

cla
set(gca, 'NextPlot', 'replacechildren')
set(gca, 'colororder', tmp)
plot(x, y, 'linewidth', 9)
axis([0 N+1 0 2])
xlabel 'color'

end


function out = jf_color_order_bw

out = [
1 1 0; % bright yellow
0 1 1; % bright cyan
1 0 1; % bright magenta
0 1 0; % green
1 0 0; % red
0 0 1; % blue
%0 3/4 3/4;
%3/4 0 3/4;
% 3/4 3/4 0;
1 1/2 0; % orange
1 1 1; % white
%1/2 1/2 1/2; % gray
];
