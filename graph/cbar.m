 function hh = cbar(varargin)
%function hh = cbar(varargin)
%|
%| colorbar with options
%|	'h' 'horiz'	horizontal
%|	'v' 'vert'	vertical
%|	'below'		horizontal colorbar below current plot (jf)
%|	'hide'		make room for it, but hide it (invisible)
%|	'fSIZE'		font size
%|	1d-array	ytick
%|	'notick'	disable tick marks
%|
%| Copyright 2007, Jeff Fessler, University of Michigan

if nargin == 1 && streq(varargin{1}, 'test'), ir_cbar_test, return, end

if isfreemat
	return % freemat 3.6 colorbar does not work with subplot
end

if ~im('ison')
%	disp 'im disabled'
return
end


% handle state of display or not
persistent Display
if ~isvar('Display') || isempty(Display)
	Display = true;
end

st.digits = 3;
st.dotick = 1;
st.ytick = [];
st.orient = [];
st.fontsize = ir_fontsize('tick');
st.new = 0;
st.label = '';

while length(varargin)
	arg = varargin{1};

	if streq(arg, 'on')
		Display = true;
		printm 'enabling cbar'
		return

	elseif streq(arg, 'off')
		Display = false;
		printm 'disabling cbar'
		return

	% new
	elseif streq(arg, 'new')
		st.new = 1;

	% notick
	elseif streq(arg, 'notick')
		st.dotick = 0;

	% ytick
	elseif isa(arg, 'double')
		st.ytick = arg;

	% 'h' or 'horiz' for horizontal
	elseif ischar(arg) && (streq(arg, 'h') || streq(arg, 'horiz'))
		if ir_is_octave
			st.orient = 'southoutside';
		else
			st.orient = 'horiz';
%			colorbar horiz; return % fixed
		end

	% 'v' or 'vert' for vertical
	elseif ischar(arg) && (streq(arg, 'v') || streq(arg, 'vert'))
		st.orient = [];

	% 'below'
	elseif ischar(arg) && streq(arg, 'below')
		st.orient = 'below';

	% 'wo' or 'westoutside'
	elseif ischar(arg) && (streq(arg, 'wo') || streq(arg, 'westoutside'))
		st.orient = 'westoutside';

	% 'hide'
	elseif ischar(arg) && streq(arg, 'hide')
		set(colorbar, 'ytick', [], 'visible', 'off')
		return

	% 'fSIZE'
	elseif ischar(arg) && streq(arg, 'f', 1)
		st.fontsize = sscanf(arg, 'f%d');

	else
		if ischar(arg) && isempty(st.label)
			st.label = arg;
		else
			error 'arg'
		end
	end

	varargin = {varargin{2:end}};
end
clear arg

if ~Display
	return
end

if isempty(get(gcf, 'children'))
	warn 'no figure children?'
	help(mfilename)
return
end

% explore new way
if st.new
	ha = gca;
%	get(ha)
	hi = get(ha, 'children');
	hi = hi(end); % for pre_v7
%	get(hi)
	dat = get(hi, 'cdata');
	clim = get(ha, 'clim');
	[nv nh] = size(dat);
	if streq(st.orient, 'below')
		error 'not done'
	else
		arg.npad = ceil(0.08*nh);
		arg.nramp = ceil(0.1*nh);
		arg.padv = 0;
		ramp = linspace(clim(1), clim(2), nv)';
		ramp = flipud(ramp);
		dat(:,end+[1:arg.npad]) = arg.padv;
		dat(:,end+[1:arg.nramp]) = repmat(ramp, 1, arg.nramp);
	end
	set(hi, 'cdata', dat)
%	get(hi)
	nh = size(dat,2);
	set(ha, 'xlim', [0.5 nh+0.5])
	xlims = get(ha, 'xlim');
	ylims = get(ha, 'ylim');
	text(1.05*xlims(2), ylims(2), sprintf('%g', clim(1)))
	text(1.05*xlims(2), ylims(1), sprintf('%g', clim(2)))
%	set(ha, 'xlim', [0.5 size(dat,2)+0.5+arg.npad+arg.nramp])
%	minmax(dat)
%	axis off

	if ~isempty(st.label)
		text(1.05*xlims(2), mean(ylims(1:2)), st.label)
	end
return
end

if isempty(st.orient)
	h = colorbar;
elseif ~streq(st.orient, 'below')
	h = colorbar(st.orient);
else
	h = ir_cbar_below;
	st.orient = 'horiz';
end

%return % todo

if ir_is_octave && streq(st.orient, 'southoutside')
%	xtick = get(h, 'xtick');
	xtick = get(gca, 'clim');
	set(h, 'xtick', xtick([1 end]))
return % todo: other options below for octave?
end


if streq(st.orient, 'horiz')
	xtick = st.ytick;
	if isempty(xtick)
%		xtick = get(gca, 'clim');
		xtick = ir_cbar_tick(get(gca, 'clim'), st.digits); % truncate
	end
else
	if isempty(st.ytick)
		st.ytick = get(gca, 'clim');
		if ~isequal(class(st.ytick), 'int32') % trick
			st.ytick = ir_cbar_tick(st.ytick, st.digits); % truncate
		end
	end
end

if st.dotick
	if streq(st.orient, 'horiz')
		set(h, 'xtick', xtick)
	else
		set(h, 'ytick', st.ytick)
	end
	% 2016-03-12 R2015b, eliminated because it messes up "x 10^5"
	if 0 && ~ir_is_octave
		tmp = get(h, 'ticklabels');
		tmp = remove_trailing_zeros(tmp);
		set(h, 'ticklabels', tmp);
	end

%{
	% disabled because not working
		% trick: for v7, move ticks in slightly
		yticks = num2str(st.ytick');
		ytick = st.ytick + [1 -1] * 0.005 * diff(st.ytick);
		if 0 && length(ytick) == 2 % kludge:
			tmp1 = get(h, 'yticklabel');
			tmp2 = strvcat(yticks, tmp1);
			yticks = tmp2([1 4:end-1 2], :);
			yticks(2:end-1,:) = ' ';
			set(h, 'yticklabel', yticks)

	% this way should work but has had problems:
		set(h, ...
			'YTickMode', 'manual', ...
			'Ytick', ytick, ...
			'YTickLabelMode', 'manual', ...
			'YtickLabel', yticks)
%}

else
	set(h, 'ytick', [])
end

if ~isempty(st.fontsize)
	set(h, 'fontsize', st.fontsize)
end

if ~isempty(st.label)
	xlims = get(h, 'xlim'); % [-0.5 1.5]
	ylims = get(h, 'ylim');
	% new for R2014b:
	hl = get(h, 'label');
	set(hl, 'string', st.label);
%{
	htmp = gca;
	axes(h) % old
	text(2.2, mean(ylims), st.label, ...
		'rotation', 90, ...
		'verticalalign', 'top', ...
		'horizontalalign', 'center')
%	ylabel(label)
%	[x y] = ginput(1)
	axes(htmp) % return current axis to image
%}
end

if nargout
	hh = h;
end


% remove_trailing_zeros()
% change 3.0200 to 3.02
function out = remove_trailing_zeros(in)
out = in
for ii=1:numel(in)
	tmp = in{ii};
	if ~isempty(strfind(tmp, '.'))
		while (tmp(end) == '0')
			tmp(end) = '';
		end
		out{ii} = tmp;
	end
end


% ir_cbar_below()
function h = ir_cbar_below(vfrac)
pos = get(gca, 'position');
clim = get(gca, 'clim');
h = pos(4);
pos(2) = pos(2) - 0.11 * h;
pos(4) = 0.1 * h;
axes('position', pos)
x = linspace(clim(1),clim(2),101);
y = 1:10;
im(x, y, x'*ones(size(y)), clim, ' ');
h = gca;
ytick('off')
axis normal


% ir_cbar_tick()
% truncate to "digits" of precision:
% trick with ceil and floor to avoid losing the label
% need to "round" towards the middle of the clim range
function tick = ir_cbar_tick(clim, digits)
if clim(2) > 100
	clim(2) = floor(clim(2));
end
type2 = 'floor'; if clim(2) < 0, type2 = 'ceil'; end
type1 = 'ceil'; if clim(1) < 0, type1 = 'floor'; end
tick(1) = truncate_precision(clim(1), digits, 'type', type1);
tick(2) = truncate_precision(clim(2), digits, 'type', type2);


% ir_cbar_test
function ir_cbar_test
im plc 2 5
clim = [5 20];
x = 10 * eye(9);
if im
	x1 = 5.083e4 * x; x1(1) = -1.724e5;
	tmp = sprintf('max=%.2g', max(x1(:)));
	im(1, x1, tmp)
	colorbar
	im(6, x1, tmp)
	cbar
	xlabelf '$x$'
end

if im
	im(2, x1, tmp), colorbar horiz
	xlabelf '$x$'
	im(7, x1, tmp)
	cbar h
	xlabelf '$x$'
end
if im
	im(8, x, clim)
	cbar %notick
%	prompt
end
if im
	im(9, x, clim)
	cbar new
%	prompt
end
if 0 % too slow
	im(1, x, clim)
	cbar 'label'
%	prompt
end
