  function [sino, pos, ang] = rect_sino(sg, rects, varargin)
%|function [sino, pos, ang] = rect_sino(sg, rects, [options])
%|
%| Create sinogram projection of one or more rectangles.
%| Works for both parallel-beam geometry and for fan-beam geometry.
%|
%| in
%|	sg			sinogram geometry object from sino_geom()
%|	rects	[ne 6]		[centx centy widthx widthy angle_degrees value]
%|
%| options
%|	'oversample'		oversampling factor for emulating "strips"
%|				default 1: just 1 ray per detector element
%|	'xscale'		use -1 to flip in x (not recommended)
%|	'yscale'		use -1 to flip in y (not recommended)
%|
%| out
%|	sino	[nb na]		sinogram
%|	pos	[nb 1]		radial sample positions (or [nb na] if mojette)
%|	ang	[na]		angular sample locations (in radians)
%|
%| 2008-08-04, Yong Long, adapted from ellipse_sino()
%| Copyright 2003-10-22, Jeff Fessler, University of Michigan

if ~nargin, help(mfilename), error(mfilename), end
if streq(sg, 'test')
	rect_sino_test % no test on moj
return
end

% defaults
arg.xscale = 1;
arg.yscale = 1;
arg.oversample = 1;
arg = vararg_pair(arg, varargin);

switch sg.type
case 'fan'
	[sino pos] = rect_sino_go(rects, ...
		sg.nb, sg.ds, sg.offset_s, ...
		sg.na, sg.ar, ...
		sg.dso, sg.dod, sg.dfs, ...
		sg.source_offset, ...
		arg.xscale, arg.yscale, arg.oversample, 0);

case 'par'
	[sino pos] = rect_sino_go(rects, ...
		sg.nb, sg.dr, sg.offset_r, ...
		sg.na, sg.ar, ...
		inf, 1, 0, ...
		sg.source_offset, ...
		arg.xscale, arg.yscale, arg.oversample, 0);

case 'moj'
	[sino pos] = rect_sino_go(rects, ...
		sg.nb, [], sg.offset_r, ...
		sg.na, sg.ar, ...
		inf, 1, 0, ...
		sg.source_offset, ...
		arg.xscale, arg.yscale, arg.oversample, sg.dx);

otherwise
	fail('sino geom "%s" not done', sg.type)
end

if nargout > 2
	ang = sg.ar;
	warn 'this is redundant with sg.ar'
end


% rect_sino_go()
function [sino, pos] = rect_sino_go(rects, ...
	nb, ds, offset_s, ...
	na, ang, ...
	dso, dod, dfs, ...
	source_offset, ...
	xscale, yscale, ...
	oversample, mojette);

[pos pos2] = rect_sino_pos([], nb, ds, offset_s, ...
		oversample, mojette, ang);

sino = rect_sino_do(rects, pos2, ang(:)', ...
	xscale, yscale, ...
	dso, dod, dfs, source_offset);

if oversample
	sino = downsample2(sino, [oversample 1]);
end


% rect_sino_pos()
% determine usual and fine "radial" sample positions
function [pos_coarse pos_fine] = ...
	rect_sino_pos(pos, nb, ds, offset_s, nover, mojette, ang)

if isempty(pos)
	if isempty(nb), warning 'nb required when no pos provided', end
	wb = (nb-1)/2 + offset_s;
else
	if ~isempty(nb), warning 'nb ignored when pos provided', end
end

if mojette ~= 0 % tricky mojette radial sampling
	% trick: ray_spacing aka ds comes from dx which is in mojette
	if ~isempty(ds), error 'ds must be empty for mojette case', end
	dt = abs(mojette) * max(abs(cos(ang)), abs(sin(ang)))'; % [1,na]

	if ~isempty(pos), error 'mojette requires empty "pos"', end

	na = length(ang);
	pos_coarse = ([0:nb-1]-wb)' * dt(:)'; % [nb na]
	if nover > 1
		pos_fine = zeros(nover*nb,na);
		for ia=1:na
			tmp = [-(nover-1):2:(nover-1)]' / (2*nover) * dt(ia);
			pos_fine(:,ia) = col(outer_sum(tmp, pos_coarse(:,ia)'));
		end
	else
		pos_fine = pos_coarse;
	end

else % ordinary non-mojette sampling

	if isempty(pos)
		pos_coarse = ([0:nb-1]'-wb) * ds; % [nb 1]
	else
		pos_coarse = pos(:); % [nb 1]
		ds = pos(2) - pos(1);
		if any(abs(diff(pos) / ds - 1) > 1e-10)
			error 'uniform spacing required'
		end
	end

	if nover > 1
		% determine fine sampling positions
		pos_fine = [-(nover-1):2:(nover-1)]' / (2*nover) * ds;
		pos_fine = outer_sum(pos_fine, pos_coarse(:)'); % [nover nb]
		pos_fine = pos_fine(:); % [nb*nover 1]
	else
		pos_fine = pos_coarse;
	end
end


% rect_sino_do()
% analytical line-integral projections of rectangle
function sino = rect_sino_do(rects, pos, ang, ...
	xscale, yscale, dso, dod, dfs, source_offset)

nb = size(pos,1);
na = length(ang);

% effective radial and angular sample locations in parallel geometry!
if isinf(dso)
	if size(pos,2) > 1 % mojette
		rads = pos;
		angs = repmat(ang, [nb 1]); % [nb na]
	else
		[rads angs] = ndgrid(pos, ang); % [nb na]
	end
else % fan
	if size(pos,2) > 1, error 'mojette fan not supported', end
	dis_src_det = dso + dod;

	if isinf(dfs)	% flat detector
		gam = atan(pos / dis_src_det); % "gamma"
	else			% arc detector
		dis_foc_det = dfs + dis_src_det;
		alf = pos / dis_foc_det;
		gam = atan2(	dis_foc_det * sin(alf), ...
				dis_foc_det * cos(alf) - dfs );	% gamma
	end

	rad = dso * sin(gam) + source_offset * cos(gam);	
	rads = rad(:,ones(1,na)); % [nb na]
	angs = outer_sum(gam, ang); % [nb na] gamma + beta
end
clear alf gam rad pos ang

sino = zeros(nb,na);

cangs = cos(angs);
sangs = sin(angs);

% loop over rects
ticker reset
ne = size(rects,1);

if size(rects, 2) ~= 6, error '6 parameters per rect', end
for ie = 1:ne
	ticker(mfilename, ie, ne)
	rect = rects(ie,:); % [1 6]

	cx = xscale * rect(1);	wx = rect(3);
	cy = yscale * rect(2);	wy = rect(4);
	if (wx <= 0 || wy <= 0)
		fail 'need positive rectangle sizes'
	end
	eang = deg2rad(rect(5));
	if yscale == -1
		eang = -eang;
	end
	if xscale == -1
		eang = pi - eang;
	end
	val = rect(6);

	cos_rot = cangs * cos(eang) + sangs * sin(eang);
	sin_rot = sangs * cos(eang) - cangs * sin(eang);
	rp = sqrt((wx * cos_rot).^2 + (wy * sin_rot).^2); % projected radius

	sp = cx .* cangs + cy .* sangs; % radial shift
	dis = (rads - sp) ./ rp; % scaled distance from center

	% projection angle after affine scaling and rotate
	abs_cos_ang_pi = wx * abs(cos_rot) ./ rp;
	abs_sin_ang_pi = wy * abs(sin_rot) ./ rp;

	% break points of the trapezoid
	len = 1 ./ max(abs_cos_ang_pi, abs_sin_ang_pi);
	dmax = (abs_cos_ang_pi + abs_sin_ang_pi) / 2;
	dbreak = abs(abs_cos_ang_pi - abs_sin_ang_pi) / 2;
	dmax_break = dmax - dbreak;
	scale = val * wx * wy ./ rp .* len;

	% line integrals at the left part
	leftpart = -dmax < dis & dis < -dbreak;
	tmp = div0(leftpart, dmax_break); % trick to avoid / 0
	sino = sino + scale .* (dis + dmax) .* tmp;
	% line integrals at the middle part
	midpart = abs(dis) <= dbreak;
	sino = sino + scale .* midpart;
	% line integrals at the right part
	rightpart = dbreak < dis & dis < dmax;
	tmp = div0(rightpart, dmax_break); % trick to avoid / 0
	sino = sino + scale .* (dmax - dis) .* tmp;
end


% rect_sino_test()
% internal test routine: standard sampling
function rect_sino_test

down = 1;
% when down=1 and the rectangle(s) line up with the grid,
% wtfmex outperforms DD obviously
ig = image_geom('nx', 512, 'ny', 504, 'fov', 500, 'down', down);
rect = [[60 10.5 38 27]*ig.dx 0 1];
[xtrue rect] = rect_im(ig, rect, 'oversample', 3);
%unique(xtrue(:))

sg = sino_geom('ge1', 'strip_width', 'd', 'down', down);
%	'dfs', 0, ...%inf, ... % flat fan
%sg.plot(ig)

if im
	im plc 2 3
	im(1, ig.x, ig.y, xtrue,'xtrue'), cbar
	axis([-1 1 -1 1]*max(ig.x)), axis square
end

tik = @() xytick(0, [0:45:360]);
ya = rect_sino(sg, rect, 'oversample', 10); % "analytical" sinogram
im(4, sg.s, sg.ad, ya, 'ya'), cbar, tik(); xlabel 's'

if 0 % too big unless we downsample
	As = Gtomo2_strip(sg, ig);
	ys = sg.shape(As * xtrue(ig.mask));
	im(2, sg.s, sg.ad, ys, 'ys'), cbar
	im(3, sg.s, sg.ad, abs(ys-ya)), cbar, tik();
	titlef('Gtomo2_strip err %g%%', max_percent_diff(ya, ys))
	nrms(ya, ys)
return
end

if has_mex_jf
	Adsc = Gtomo2_dscmex(sg, ig);
%	Adsc = Gtomo2_wtmex(sg, ig); % too big unless we downsample a lot

	yw = Adsc * xtrue;
	im(2, sg.s, sg.ad, yw, 'yw'), cbar, tik();
	im(3, sg.s, sg.ad, abs(yw-ya)), cbar, tik();
	titlef('wtfmex max err %4.2f%%', max_percent_diff(ya, yw))
	nrms(yw, ya)
end

if 2 == exist('Gtomo_dd') % DD, only if user has it
	Add = Gtomo_dd(sg, ig, 'nthread', 1);

	yd = Add * xtrue;

	im(5, sg.s, sg.ad, yd, 'yd'), cbar, tik();
	im(6, sg.s, sg.ad, abs(yd-ya)), cbar, tik();
	titlef('dd max err %4.2f%%', max_percent_diff(ya, yd))
	nrms(yd, ya)
end
