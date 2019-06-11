 function subplot_stack(x, ys, str_title, colors)
%function subplot_stack(x, ys, str_title, colors)
% a tight stack of subplots to show L signal components
% in
%	x	[N,1]
%	ys	[N,L]		or ?

if nargin == 1 && streq(x, 'test'), subplot_stack_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end
if ~isvar('colors') || isempty(colors), colors = {'c', 'y'}; end
if ~isvar('str_title') || isempty(str_title)
	str_title = '';
end

if ~iscell(ys)
	ys = {ys};
end
L = size(ys{1},2);

apos = get(gca, 'position'); % current axes position
for ll=1:L
%	pos = [0.1 0.1+0.8/L*(L-ll) 0.8 0.8/L];
	pos = [apos(1) apos(2)+apos(4)/L*(L-ll) apos(3) apos(4)/L];
	subplot('position', pos) % l b w h
	for ip=1:length(ys)
		plot(	x, real(ys{ip}(:,ll)), colors{1+2*(ip-1)}, ...
			x, imag(ys{ip}(:,ll)), colors{2+2*(ip-1)})
		if ip == 1, hold on, end
	end
	hold off
	axis tight
	ytick(0), set(gca, 'yticklabel', '')
	fontsize = 10;
	fontweight = 'normal';
	texts(1.02, 0.8, sprintf('%d.', ll), ...
		'fontsize', fontsize, 'fontweight', fontweight)

	if ll==1, title(str_title), end
	if ll<L
		xtick off
	end
end

function subplot_stack_test
x = linspace(0,1,101)';
y = exp(2i*pi*x*[1:5]);
if im
	clf, subplot(121)
	subplot_stack(x, y)
end
