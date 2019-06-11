 function ob = Gcone(cg, ig, varargin)
%function ob = Gcone(cg, ig, [options])
%|
%| Construct Gcone object, for 3D axial cone-beam or helical cone-beam
%| forward and back-projection.  Works with flat detectors or 3rd-gen CT arcs.
%| Also supports 3D parallel-beam geometry, but no polar angle.
%|
%| in
%|	cg	ct_geom object (scanner geometry)
%|	ig	image_geom object (recon grid dimensions, FOV and mask)
%|
%| option
%|	'type'		'sf2' : recommended (default)
%|			'sf2,v:all' : thread-within-view (for tomosynthesis)
%|				(default if nthread > cg.na)
%|			'pd1' : basic pixel-driven (not very good)
%|			'nn1' : nearest-neighbor pixel-driven (very bad)
%|			'user' : user must provide 'mexfun'
%|			'rf1' : rect/rect footprint [todo: buggy]
%|			'dd1' : distance-driven (GE's patented method, UM only)
%|			'dd2' : new improved DD version from GE (UM only)
%|	'is_ns_nt'	1 (default) for [ns nt na] or 0 for [nt ns na]
%|	'nthread'	# of threads, default: jf('ncore')
%|	'mexfun'	only needed for 'user' option; see usage in code.
%|	'mexarg'	cell array of optional arguments to mexfun
%|				e.g., single(dl) : [0 fovz/2), for sf4,sf5 only
%|			(dividing line for rectangle and trapezoid approx.)
%|	'scale'		scale mex output by this. (default: [])
%|	'use_hct2'	call hct2 instead of mex (UM only, for testing)
%|	'zxy'		if 1, assume image in zxy order.  default: 0 (xyz)
%|			(note: ig and ig.mask is always in xyz order.)
%|
%| out
%|	ob [nd np]	np = sum(mask(:)), so it is already "masked"
%|			nd = [ns * nt * na] if is_ns_nt = 1
%|			nd = [nt * ns * na] if is_ns_nt = 0
%|
%| For more help, type 'jf_mex cbct' or 'dd_ge1_mex' or 'dd_ge2_mex'
%| and see Gcone_test.m
%|
%| Copyright 2005-5-19, Jeff Fessler, University of Michigan

if nargin == 1 && streq(cg, 'test'), run_mfile_local Gcone_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% option defaults
arg.type = ''; % see below
arg.chat = false;
arg.scale = [];
arg.is_ns_nt = true; % default: [ns nt na] (channel fastest, row, view slowest)
arg.nthread = jf('ncore');
arg.mexfun = [];
arg.mexarg = {single(0)};
arg.use_hct2 = false; % UM only
arg.zxy = false;
arg.class = 'fatrix2';

% options specified by name/value pairs
arg = vararg_pair(arg, varargin);

if isempty(arg.type) % sf2 is default system type
	% 'sf2,v:all' recommended for (e.g.) tomosynthesis
	if arg.nthread > cg.na
	%	arg.type = ['sf2,v:' sprintf('%03d', arg.nthread) tmp];
		arg.type = 'sf2,v:all'; % use thread-within-view
	else
		arg.type = 'sf2';
	end
	if arg.chat
		printm('using system type="%s"', arg.type)
	end
end

% initialize generic geometry
arg = Gcone_setup(cg, ig, arg);

arg.dd_flip_x = false;
arg.dd_flip_y = false;
arg.dd_flip_z = false;

arg.cg = cg;
arg.ig = ig;
arg.nd = cg.ns * cg.nt * cg.na;
arg.np = sum(ig.mask(:));
dim = [arg.nd arg.np]; % trick: make it masked by default!

arg.mask2 = uint8(sum(ig.mask, 3) > 0);
arg.back_needs_mask3 = ~isequal(repmat(arg.mask2, [1 1 ig.nz]), ig.mask);

% case specific setup
if streq(arg.type, 'sf1', 3) || streq(arg.type, 'sf2', 3) ...
	|| streq(arg.type, 'sf6', 3)
	tmp = arg.type(1:3); % trick to handle sf1,... and sf2,...
else
	tmp = arg.type;
end
switch tmp % arg.type
case {'nn1', 'pd1', 'rf1', 'sf1', 'sf2', 'sf3', 'sf4', 'sf5', 'sf6'}
	arg.mexfun = Gcone_which_mex(arg.chat);
	arg.mexstr_back = ['cbct,back,' arg.type];
	arg.mexstr_proj = ['cbct,proj,' arg.type];

case 'user'
	if isempty(arg.mexfun), fail('mexfun required'), end

case 'dd1' % 1st GE version of DD (UM only)
	arg.mexfun = @dd_ge1_mex;
	arg = Gcone_setup_dd1(cg, ig, arg);
	arg.dd_flip_x = arg.ig.dx < 0;
	arg.dd_flip_y = arg.ig.dy > 0;
	arg.dd_flip_z = arg.ig.dz < 0;
	arg.back_needs_mask3 = true; % codo: dd_ge1_mex does not take mask

case 'dd2' % 2nd GE version of DD (UM only)
	arg.mexfun = @dd_ge2_mex;
	arg = Gcone_setup_dd2(cg, ig, arg);
	arg.dd_flip_x = arg.ig.dx < 0;
	arg.dd_flip_y = arg.ig.dy > 0;
	arg.dd_flip_z = arg.ig.dz < 0;
	arg.back_needs_mask3 = true; % codo: dd_ge2_mex does not take mask

otherwise
	fail('type %s unknown', arg.type)
end

if arg.back_needs_mask3
	arg.mask3_zxy = permute(ig.mask, [3 1 2]);
end

if arg.use_hct2
	pn = jf_protected_names;
	arg.mexfun = [];
	arg.hct2.x = [test_dir 'gcone-hct-x.fld'];
	arg.hct2.y = [test_dir 'gcone-hct-y.fld'];
	arg.hct2.m = [test_dir 'gcone-hct-mask-xyz.fld'];
	arg.hct2.name = 'hct2';

	com = [arg.hct2.name ' chat 0 what %s ' ...
		sprintf(' nthread %d', arg.nthread) ...
		sprintf(' sysmod %s', arg.type) ...
		sprintf(' file_yi %s file_init %s', arg.hct2.y, arg.hct2.x) ...
		sprintf(' file_mask %s', arg.hct2.m) ...
		pn.hct_arg(arg.cg, arg.ig) ...
	];
	arg.hct2.com = com;

else % default "scale" is xy voxel_size for dd1 and dd2 mex version
	if streq(arg.type, 'dd1') || streq(arg.type, 'dd2')
		if isempty(arg.scale)
			arg.scale = arg.voxel_size(1);
		else
			warn 'user scale factor for dd1 and dd2?'
		end
	end
end

if arg.is_ns_nt
	arg.odim = [arg.cg.ns arg.cg.nt arg.cg.na];
else
	arg.odim = [arg.cg.nt arg.cg.ns arg.cg.na];
end

switch arg.class
case 'fatrix2'
	if arg.use_hct2
		forw = @(arg, x) Gcone_mtimes_block_hct2_array(arg, 0, x, 1, 1);
		back = @(arg, y) Gcone_mtimes_block_hct2_array(arg, 1, y, 1, 1);
%		block = @Gcone_mtimes_block_hct2_array;
		forw_block = @(arg, x, iblock, nblock) ...
			Gcone_mtimes_block_hct2_array(arg, 0, x, iblock, nblock);
		back_block = @(arg, y, iblock, nblock) ...
			Gcone_mtimes_block_hct2_array(arg, 1, y, iblock, nblock);
	else
		forw = @(arg, x) Gcone_block_forw_array(arg, x, 1, 1);
		back = @(arg, y) Gcone_block_back_array(arg, y, 1, 1);
		forw_block = @(arg, x, iblock, nblock) ...
			Gcone_block_forw_array(arg, x, iblock, nblock);
		back_block = @(arg, y, iblock, nblock) ...
			Gcone_block_back_array(arg, y, iblock, nblock);
%		block = @Gcone_mtimes_block_array;
	end
%	forw = @(arg, x) block(arg, 0, x, 1, 1);
%	back = @(arg, y) block(arg, 1, y, 1, 1);

	if arg.zxy
		mask3 = permute(ig.mask, [3 1 2]);
	else
		mask3 = ig.mask;
	end
	ob = fatrix2('mask', mask3, 'arg', arg, 'does_many', 1, ...
		'abs', @(ob) ob, ... % aij nonnegative
		'forw', forw, 'back', back, ...
		'forw_block', forw_block, 'back_block', back_block);
	if ~isempty(arg.scale)
		ob = arg.scale * ob;
		arg.scale = [];
	end

case 'Fatrix'
	if arg.zxy, fail 'zxy done only for fatrix2 version', end
	if arg.use_hct2
%		forw = @(arg, x) Gcone_mtimes_block_hct2(arg, 0, x, 1, 1);
%		back = @(arg, y) Gcone_mtimes_block_hct2(arg, 1, y, 1, 1);
		block = @Gcone_mtimes_block_hct2;
	else
%		forw = @(arg, x) Gcone_block_forw(arg, x, 1, 1);
%		back = @(arg, y) Gcone_block_back(arg, y, 1, 1);
		block = @Gcone_mtimes_block;
	end
	forw = @(arg, x) block(arg, 0, x, 1, 1);
	back = @(arg, y) block(arg, 1, y, 1, 1);

	ob = Fatrix(dim, arg, ...
		'abs', @(ob) ob, ... % aij nonnegative
		'forw', forw, 'back', back, 'mtimes_block', block);
otherwise
	fail('class')
end


%
% Gcone_which_mex()
% pick cbct_mex if available (for development) otherwise use jf_mex
%
function fun = Gcone_which_mex(chat)
if exist('cbct_mex') == 3
	fun = @cbct_mex;
	if chat, printm 'using cbct_mex', end
elseif exist('jf_mex') == 3
	fun = @jf_mex;
else
	fail('cannot find jf_mex or cbct_mex')
end


% Gcone_setup()
function arg = Gcone_setup(cg, ig, arg)

arg.voxel_size = abs([ig.dx ig.dy ig.dz]);
arg.img_offset = single([ig.offset_x ig.offset_y ig.offset_z]);
arg.angles = single(cg.ar);


% Gcone_setup_dd2()
function arg = Gcone_setup_dd2(cg, ig, arg)

arg = Gcone_setup_dd1(cg, ig, arg); % trick: recycle most of the setup of dd1

% caution: z normalization differs for dd2:
zscale = arg.voxel_size(3) / arg.voxel_size(1);
arg.zds = arg.zds * zscale;
arg.zshifts = arg.zshifts * zscale; % cg.source_zs / dx
arg.pos_source(3) = arg.pos_source(3) * zscale;

arg.img_offset(1) = -arg.img_offset(1); % found empirically for dx > 0
arg.img_offset(3) = -arg.img_offset(3) * zscale; % found empirically for dz > 0


% Gcone_setup_dd1()
% note: dd1 code uses voxel_size normalized units
function arg = Gcone_setup_dd1(cg, ig, arg)

if abs(ig.dx) ~= abs(ig.dy)
	error 'need |dx| = |dy|'
end

arg.pos_source = [0 cg.dso 0];

% detector sample locations
arg.ss = cg.s; % sample centers along arc length of detector
arg.zds = cg.t;
if ~isequal(arg.zds, sort(arg.zds))
	error 'z positions of detectors, i.e., cg.t, must be ascending'
end

if isinf(cg.dso)
	error 'parallel beam not done'

% cone-beam with flat detector
elseif isinf(cg.dfs)
	arg.xds = arg.ss;
	arg.yds = repmat(-cg.dod, [1 cg.ns]);

% 3rd generation multi-slice CT with focal point of arc at isocenter
else
	if cg.dfs ~= 0, warning 'dis_foc_src nonzero?', end
	dis_foc_det = cg.dfs + cg.dsd; % 3rd gen
	t = arg.ss / dis_foc_det; % angle in radians
	arg.xds = dis_foc_det * sin(t);
	arg.yds = cg.dso - dis_foc_det * cos(t);
end


% dd1 code expects in normalized pixel units
arg.xds = single(arg.xds / arg.voxel_size(1));
arg.yds = single(arg.yds / arg.voxel_size(2));
arg.zds = single(arg.zds / arg.voxel_size(3)); % z pixels
arg.zshifts = single(cg.source_zs / arg.voxel_size(3)); % z pixels
arg.pos_source = single(arg.pos_source ./ arg.voxel_size);
arg.dz_dx = single(arg.voxel_size(3) / arg.voxel_size(1));


% Gcone3_dd_xyz_flipper3()
% assumes "x" is in xyz order
function x = Gcone3_dd_xyz_flipper3(x, arg);
if arg.dd_flip_x < 0
	x = flipdim(x, 1);
end
if arg.dd_flip_y > 0
	x = flipdim(x, 2);
end
if arg.dd_flip_z < 0
	x = flipdim(x, 3);
end


% Gcone3_dd_zxy_flipper3()
% assumes "x" is in zxy order
function x = Gcone3_dd_zxy_flipper3(x, arg);
if arg.dd_flip_x < 0
	x = flipdim(x, 2);
end
if arg.dd_flip_y > 0
	x = flipdim(x, 3);
end
if arg.dd_flip_z < 0
	x = flipdim(x, 1);
end


%
% fatrix2 versions "_array"
%

%{
% Gcone_mtimes_block_array()
function y = Gcone_mtimes_block_array(arg, is_transpose, x, istart, nblock)

if is_transpose
	y = Gcone_block_back_array(arg, x, istart, nblock);
else
	y = Gcone_block_forw_array(arg, x, istart, nblock);
end
%}


% Gcone_block_forw_array()
function y = Gcone_block_forw_array(arg, x, istart, nblock)

if ~arg.zxy % trick: dd1|dd2|nn1|pd1 mex code wants zxy order (z dim first)!
	x = permute(x, [3 1 2 4]);
end

ia = istart:nblock:arg.cg.na;

switch arg.type
case {'dd1', 'dd2'}

	x = Gcone3_dd_zxy_flipper3(x, arg); % only DD might need this
	y = arg.mexfun('proj3', arg.pos_source, arg.xds, arg.yds, arg.zds, ...
		arg.dz_dx, arg.img_offset, int32(arg.nthread), ...
		arg.angles(ia), arg.zshifts(ia), ...
		single(x));

otherwise % nn1 pd1 rf1 sf1 sf2 sf3 sf4 sf5 sf6

	if ((x(1) ~= 0) && (arg.mask2(1) == 0))
		warn 'projecting image with nonzero x(1) but mask2(1)=0'
		% codo: x = x .* ig.mask;
	end

	% codo: could add is_ns_nt here to save permute at end
	s = @(x) single(x);
%	cpu etic
	angles = arg.angles(ia);
	y = arg.mexfun(arg.mexstr_proj, ...
		int32([arg.cg.ns arg.cg.nt length(ia)]), ...
		s([arg.ig.dx arg.ig.dy arg.ig.dz]), ...
		s([arg.ig.offset_x arg.ig.offset_y arg.ig.offset_z]), ...
		arg.mask2, ...
		s(arg.cg.dso), s(arg.cg.dsd), s(arg.cg.dfs), ...
		s([arg.cg.ds arg.cg.dt]), ...
		angles, ...
		s(arg.cg.source_zs(ia)), ...
		s(arg.cg.offset_s), ...
		s(arg.cg.offset_t), ...
		s(x), ...
		int32(arg.nthread), ...
		arg.mexarg{:});
%	cpu etoc Gcone

end

if arg.is_ns_nt
	y = permute(y, [2 1 3 4]); % dd1|dd2|nn1|pd1|sf123 code makes [nt ns na]
end


% Gcone_block_back_array()
function x = Gcone_block_back_array(arg, y, istart, nblock)

ia = istart:nblock:arg.cg.na;

if arg.is_ns_nt
	y = permute(y, [2 1 3 4]); % dd1|dd2|nn1|pd1|sf123 code wants [nt ns na]
end

switch arg.type
case {'dd1', 'dd2'}

	x = arg.mexfun('back3', arg.pos_source, arg.xds, arg.yds, arg.zds, ...
		arg.dz_dx, arg.img_offset, int32(arg.nthread), ...
		arg.angles(ia), arg.zshifts(ia), ...
		single(y), int32([arg.ig.nz arg.ig.nx arg.ig.ny]));
	x = Gcone3_dd_xyz_flipper3(x, arg); % only DD might need this

otherwise % nn1 pd1 rf1 sf1 sf2 sf3 sf4 sf5 sf6

	% codo: could add is_ns_nt here to save permute at end
	s = @(x) single(x);
	x = arg.mexfun(arg.mexstr_back, ...
		int32([arg.ig.nx arg.ig.ny arg.ig.nz]), ...
		s([arg.ig.dx arg.ig.dy arg.ig.dz]), ...
		s([arg.ig.offset_x arg.ig.offset_y arg.ig.offset_z]), ...
		arg.mask2, ...
		s(arg.cg.dso), s(arg.cg.dsd), s(arg.cg.dfs), ...
		s([arg.cg.ds arg.cg.dt]), ...
		arg.angles(ia), ...
		s(arg.cg.source_zs(ia)), ...
		s(arg.cg.offset_s), ...
		s(arg.cg.offset_t), ...
		s(y), ...
		int32(arg.nthread), ...
		arg.mexarg{:});

end

if arg.back_needs_mask3
%	x = x .* repmat(arg.ig.mask, [1 1 1 size(y,4)]); % trick: apply 3d mask
	x = x .* repmat(arg.mask3_zxy, [1 1 1 size(y,4)]);
end

if ~arg.zxy % trick: dd1|dd2|nn1|pd1|sf123 code puts z dim first!
	x = permute(x, [2 3 1 4]);
end


%
% Fatrix versions
%

% Gcone_mtimes_block()
function y = Gcone_mtimes_block(arg, is_transpose, x, istart, nblock)

if is_transpose
	y = Gcone_block_back(arg, x, istart, nblock);
else
	y = Gcone_block_forw(arg, x, istart, nblock);
end

% Gcone_block_forw()
function y = Gcone_block_forw(arg, x, istart, nblock)

[x ei] = embed_in(x, arg.ig.mask, arg.np);
y = Gcone_block_forw_array(arg, x, istart, nblock);
if ~isempty(arg.scale)
	y = y * arg.scale;
end
y = ei.shape(y);


% Gcone_block_back()
function x = Gcone_block_back(arg, y, istart, nblock)

ia = istart:nblock:arg.cg.na;
[y eo] = embed_out(y, [arg.odim(1:end-1) length(ia)]); % [(M) *L]
if ~isempty(arg.scale)
	y = y * arg.scale;
end
x = Gcone_block_back_array(arg, y, istart, nblock);
x = eo.shape(x, arg.ig.mask, arg.np);


%
% hct2 versions (UM only)
%

% Gcone_hct2_com()
% set up
function [com cg ia] = Gcone_hct2_com(arg, istart, nblock)

fld_write(arg.hct2.m, arg.ig.mask, 'check', 0); % write xyz mask

com = arg.hct2.com;
cg = arg.cg;
ia = istart:nblock:cg.na;
if istart ~= 1 || nblock ~= 1
	ad = cg.ad(ia);
	cg.orbit_start = ad(1);
	cg.source_z0 = cg.source_zs(istart);
	% caution: this causes imperfect match for block versions:
	cg.orbit = (ad(end) - ad(1)) + (ad(2) - ad(1));
	cg.na = length(ia);
	com = [com sprintf(' na %d orbit %g orbit_start %g', ...
		cg.na, cg.orbit, cg.orbit_start)]; % trick to over-ride
	if ~isempty(cg.user_source_zs)
		fail('not done')
		com = [com sprintf(' file_source_z %s', file_source_z)];
	else
		com = [com sprintf(' source_z0 %g', cg.source_z0)]; % trick
	end
end


% Gcone_hct2_forw_array()
function y = Gcone_hct2_forw_array(arg, x, com, cg)

if arg.zxy
	x = permute(x, [2 3 1 4]); % zxy to xyz
end
x = Gcone3_dd_xyz_flipper3(x, arg);

com = sprintf(com, 'proj');
LL = size(x,4);
y = zeros([cg.ns cg.nt cg.na LL], 'single');
for ll=1:LL
	if exist(arg.hct2.y, 'file'), delete(arg.hct2.y), end
	fld_write(arg.hct2.x, x(:,:,:,ll), 'check', 0)
	tmp = os_run(com);
	y(:,:,:,ll) = fld_read(arg.hct2.y);
end

if 0
	if istart ~= 1 || nblock ~= 1 % trick: simple way to do subsets
		ia = istart:nblock:cg.na;
		y = y(:,:,ia,:); % [ns nt na LL]
	end
end

if ~arg.is_ns_nt
	y = permute(y, [2 1 3 4]); % -> [nt ns na *L]
end


% Gcone_hct2_back_array()
function x = Gcone_hct2_back_array(arg, y, com)

if ~arg.is_ns_nt
	y = permute(y, [2 1 3 4]); % -> [ns nt na *L]
end

com = sprintf(com, 'back');
LL = size(y,4);
x = zeros([arg.ig.dim LL], 'single');
for ll=1:LL
	if exist(arg.hct2.x, 'file'), delete(arg.hct2.x), end
	fld_write(arg.hct2.y, y(:,:,:,ll), 'check', 0)
	tmp = os_run(com);
	x(:,:,:,ll) = fld_read(arg.hct2.x);
end

x = Gcone3_dd_xyz_flipper3(x, arg);
if arg.zxy
	x = permute(x, [3 1 2 4]); % xyz to zxy
end


% Gcone_mtimes_block_hct2_array()
% forward projector based on calling unix binary "hct2" (UM only)
function y = Gcone_mtimes_block_hct2_array(arg, is_transpose, x, istart, nblock)

[com cg] = Gcone_hct2_com(arg, istart, nblock);

if is_transpose % back
	y = Gcone_hct2_back_array(arg, x, com);
else % forw
	y = Gcone_hct2_forw_array(arg, x, com, cg);
end


% Gcone_mtimes_block_hct2()
% forward projector based on calling unix binary "hct2" (UM only)
function y = Gcone_mtimes_block_hct2(arg, is_transpose, x, istart, nblock)

[com cg ia] = Gcone_hct2_com(arg, istart, nblock);

if is_transpose % back
	y = x;
	[y eo] = embed_out(y, [arg.odim(1:end-1) length(ia)]); % [(M) *L]
	if ~isempty(arg.scale)
		y = y * conj(arg.scale);
	end
	x = Gcone_hct2_back_array(arg, y, com);
	x = eo.shape(x, arg.ig.mask, arg.np);
	y = x;

else % forw
	[x ei] = embed_in(x, arg.ig.mask, arg.np); % [nx ny nz L]
	y = Gcone_hct2_forw_array(arg, x, com, cg);
	if ~isempty(arg.scale)
		y = y * arg.scale;
	end
	y = ei.shape(y);
end
