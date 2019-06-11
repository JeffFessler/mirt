 function [out, varargout] = fbp2(varargin)
%function geom = fbp2(sg, ig, [setup_options])
%function [image, sino_filt] = fbp2(sino, geom, [recon_options])
%|
%| FBP 2D tomographic image reconstruction for parallel-beam or fan-beam cases,
%| with either flat or arc detector for fan-beam case.
%|
%| To use this, you first call it with the sinogram and image geometry.
%| The routine returns the initialized "geom" structure.  Thereafter, to
%| to perform FBP reconstruction, you call this routine with that structure
%| (perhaps numerous times for the same geometry).
%| See fbp2_example.m for examples.
%|
%| in (for setup)
%|	sg			sino_geom() structure
%|	ig			image_geom() structure
%|
%| options (for setup)
%|	'type':			type of reconstruction (default: '')
%|				(not all are implemented for fan-beam)
%|		'', 'std:mex'	call mex file for fast backprojection
%|		'std:mat'	slower matlab backprojector
%|		'dsc'		backproject using Gtomo2_dsc with system9
%|		'mojette'	use mojette rebinning and Gtomo2_table
%|		'table'		use Gtomo2_table with tabulated linear interp
%|		'df,pull'	direct Fourier reconstruction with "pull" interp
%|	'window' [npad]		'' or 'hann', or array (default: '' for ramp)
%|				if array, then use samples [-K/2, K/2)
%|	'nthread'		default: jf('ncore')
%|
%| out (for setup)
%|	geom	(struct)	initialized structure
%|
%|
%| in (for recon)
%|	sino	[nb na *]	sinogram(s) (line integrals)
%|	geom			structure from first call
%|
%| options (for recon)
%|	'window'		e.g. 'hann'  (default: '', i.e., just ramp)
%|				see fbp2_window.m or "op fbp" options
%|
%| out (for recon)
%|	image		[nx ny *]	reconstructed image(s)
%|	sino_filt	[nb na *]	filtered sinogram(s)
%|
%| Copyright 2005-12-21, Jeff Fessler, University of Michigan

if nargin == 1 && streq(varargin{1}, 'test')
	run_mfile_local fbp2_example
	return
end
if nargin < 2, help(mfilename), error(mfilename), end

if ~isnumeric(varargin{1}) % setup
	out = fbp2_setup(varargin{1}, varargin{2}, varargin{3:end});
else % recon
	[out, sino_filt] = ...
		fbp2_recon(varargin{1}, varargin{2}, varargin{3:end});
	if nargout
		varargout{1} = sino_filt;
	end
end


%
% fbp2_setup()
%
function geom = fbp2_setup(sg, ig, varargin);

% defaults
arg.type = '';
arg.extra = {}; % extra arguments for individual recon methods
arg.window = '';
arg.nthread = jf('ncore');

arg = vararg_pair(arg, varargin);

if isempty(arg.type) || streq(arg.type, 'std')
	if has_mex_jf % choose default based on mex availability
		arg.type = 'std:mex';
	else
		arg.type = 'std:mat';
	end
end

switch arg.type
	case {'std:mat', 'std:mex'}
			geom = fbp2_setup_std(sg, ig, arg);
	case 'dsc',	geom = fbp2_setup_dsc(sg, ig, arg);
	case 'df,pull',	geom = fbp2_setup_df_pull(sg, ig, arg);
	case 'mojette',	geom = fbp2_setup_moj(sg, ig, arg);
	case 'table',	geom = fbp2_setup_tab(sg, ig, arg);
	otherwise	fail('unknown type: "%s"', arg.type)
end

geom.arg_save = arg;
geom.sg = sg;
geom.ig = ig;


%
% fbp2_setup_df_pull()
% DF recon with "pulling" interpolation from polar into cartesian
%
function st = fbp2_setup_df_pull(sg, ig, arg);

if ~streq(sg.type, 'par'), error 'only parallel-beam done for DF', end

% defaults
st.over = 2; % over-sampling
st.nu = [];
st.nv = [];
st.npad = [];
%st.interp = '*nearest';
%st.interp = '*spline';
%st.interp = '*linear';
st.interp = '*cubic';

st = vararg_pair(st, arg.extra);

nx = ig.nx;
ny = ig.ny;
if isempty(st.nu), st.nu = st.over * nx; end
if isempty(st.nv), st.nv = st.over * ny; end
if isempty(st.npad), st.npad = 2^ceil(log2(2*sg.nb-1)); end

nu = st.nu;
nv = st.nv;
npad = st.npad;

if ig.dx < 0 || ig.dy ~= ig.dx, error 'flip not done', end
% trick: go to n/2 for simplicity later
u = [-nu/2:nu/2]/nu/ig.dx;
v = [0:nv/2]/nv/ig.dy;
[uu vv] = ndgrid(u, v);

st.rho = [-npad/2:npad/2-1]'/npad / sg.dr;
st.phi = sg.ar;
st.qrho = sqrt(uu.^2 + vv.^2);
st.qphi = atan2(vv, uu);

% trick: add another row at pi by flipping, to help interpolator
if sg.orbit ~= 180 || sg.orbit_start ~= 0, error 'need 0,180', end
st.phi = [st.phi(:)' pi];

% clf, plot(qrho, qphi, '.', max(rho), pi/2, 'x'), return % see samples
% minmax(qphi)

wr = (sg.nb-1)/2 + sg.offset_r;
% trick: because wr is non-integer, we must fftshift the phase as follows:
st.phase1 = sg.dr * fftshift(exp(2i*pi*[-npad/2:npad/2-1]'/npad * wr));

% phase corresponding to image-domain half-pixel shift 
wx = (nx-1)/2 + ig.offset_x;
wy = (ny-1)/2 + ig.offset_y;
phase2 = outer_sum([-nu/2:nu/2]/nu * (wx-nx/2), [0:nv/2]/nv * (wy-ny/2));
st.phase2 = exp(-2i*pi*phase2);


%
% fbp2_setup_dsc()
%
function geom = fbp2_setup_dsc(sg, ig, arg);

if streq(sg.type, 'par')
	% just backproject into whatever support the user specified
	geom.G = Gtomo2_dscmex(sg, ig, 'nthread', arg.nthread, ...
			'pairs', {'system', 9});
	geom.scale = sg.dr; % trick: overcome 1/dr in system9

%elseif arg.dfs == 0 % arc fan-beam
%	error 'fix: todo'
%
%elseif isinf(arg.dfs) % flat fan-beam
%	error 'fix: todo'

else
	error 'unknown geometry'
end


%
% fbp2_setup_moj()
%
function geom = fbp2_setup_moj(sg, ig, arg);

nr = 2 * ceil(sqrt(2)/2 * sg.nb); % fix: kludge
switch sg.type
case 'par'
	geom.moj.sg = sino_geom('moj', 'nb', nr, 'na', sg.na, ...
		'orbit', sg.orbit, 'orbit_start', sg.orbit_start, ...
		'dx', ig.dx, 'offset_r', sg.offset_r);
case 'fan'
	error 'todo: moj fan'
%	geom.moj.phi_start = 0;
%	geom.moj.phi_orbit = 360;
%	geom.moj.offset_r = sg.offset_s;
%	geom.moj.dr = sg.ds;
otherwise
	error 'bad sino type'
end

geom.moj.ob_rebin = rebin_sino([], sg, geom.moj.sg, 'ob', 1);
%	's_interp', {'order', 3, 'ending', 'zero'}, ...
%	'beta_interp', {'order', 3, 'ending', 'periodic'}, ...

geom.moj.G = Gtomo2_table(geom.moj.sg, ig, {'mojette,back1'}, ...
		'nthread', arg.nthread);

geom.moj.H = fbp2_make_sino_filter_moj(nr, geom.moj.sg.na, geom.moj.sg.dx, ...
	geom.moj.sg.orbit, geom.moj.sg.orbit_start, arg.window);


% ir_fbp2_par_parker_wt()
function wt = ir_fbp2_par_parker_wt(sg)
orbit = abs(sg.orbit);
na = sg.na;
ad = abs(sg.ad - sg.orbit_start);
nb = sg.nb;

if orbit < 180
	warn('orbit %g < 180', orbit)
	wt = 1;
return
end

if orbit > 360
	fail 'only 180 <= orbit <= 360 supported for Parker weighting'
end
extra = orbit - 180; % extra beyond 180

wt = ones(1, na);
ii = ad < extra;
wt(ii) = sin(ad(ii) / extra * pi / 2).^2;
ii = ad >= 180;
wt(ii) = sin((orbit - ad(ii)) / extra * pi / 2).^2;
wt = wt * orbit / 180; % trick because of the back-projector normalization
wt = repmat(wt, [sg.nb 1]); % [nb na] sinogram sized


%
% fbp2_setup_std()
%
function geom = fbp2_setup_std(sg, ig, arg);

switch sg.type
case 'par'
	if abs(sg.orbit) ~= 180 && abs(sg.orbit) ~= 360
		geom.parallel_beam_parker_weight = ir_fbp2_par_parker_wt(sg);
	end
	switch arg.type
	case 'std:mex'
		geom.arg_back2 = {uint8(ig.mask), ...
			ig.dx, ig.dy, ...
			ig.offset_x, ...
			sign(ig.dy) * ig.offset_y, ... % trick: old backproject
			sg.dr, sg.offset, sg.orbit, sg.orbit_start, ...
			int32(arg.nthread)};
	case 'std:mat'
		geom.arg_back2 = {sg, ig};
	otherwise
		error 'bug'
	end

case 'fan'
	if sg.orbit ~= 360
		warn 'short-scan fan-beam Parker weighting not done'
	end
	d = @(x) double(x);
	geom.arg_fbp_fan = {uint8(ig.mask), ...
		d(ig.dx), d(ig.dy), ...
		d(ig.offset_x), d(ig.offset_y), ...
		d(sg.ds), d(sg.offset), d(sg.orbit), d(sg.orbit_start), ...
		d(sg.dsd), d(sg.dso), d(sg.dfs)};
	if ~isempty(arg.window), fail 'todo: window ignored', end

case 'moj'
	if sg.dx == abs(ig.dx)
		geom.moj.G = Gtomo2_table(sg, ig, {'mojette,back1'}, ...
			'nthread', arg.nthread);
	else
		warn('mojette sinogram with d=%g vs image with dx=%g', ...
			sg.dx, ig.dx)
	end

	geom.moj.H = fbp2_make_sino_filter_moj(sg.nb, sg.na, sg.dx, ...
		sg.orbit, sg.orbit_start, arg.window);

otherwise
	error 'bad sino type'
end

% todo: setup filter too


%
% fbp2_setup_tab()
%
function geom = fbp2_setup_tab(nb, na, arg);
if isinf(arg.dsd) % parallel-beam
	geom.tab.nr = nb;
else
	geom.tab.nr = 2 * nb; % fix: kludge
error 'not done'
end

geom.tab.phi_orbit = arg.orbit; % fix:
geom.tab.phi_start = arg.orbit_start; % fix:
geom.tab.offset_r = arg.offset_s;

scale = 1 / arg.dx; % trick: overcome dx scaling
'todo: use sg ig'
geom.tab.G = Gtomo2_table(ig.mask, {'linear'}, ...
	'dx', arg.flip_x * arg.dx * scale, ...
	'dy', -arg.flip_y * arg.dx * scale, ... % trick
	'nb', geom.tab.nr, ...
	'na', na, ...
	'offset_x', arg.offset_x, ...
	'offset_y', arg.offset_y, ...
	'offset_r', geom.tab.offset_r, ...
	'dr', arg.ds * scale, ...
	'orbit', geom.tab.phi_orbit, ...
	'orbit_start', geom.tab.phi_start, ...
	'nthread', arg.nthread);


%
% fbp2_recon()
%
function [image, sino_filt] = fbp2_recon(sino, geom, varargin);

if geom.sg.nb ~= size(sino, 1) || geom.sg.na ~= size(sino, 2)
	error 'bad sino size'
end
%dims = size(sino);
%sino = reshape(sino, sg.nb, sg.na, []);
%out = reshape(out, [size(out,1), size(out,2), dims(3:end)]);

% defaults for recon
opt.window = '';
opt = vararg_pair(opt, varargin);

% conventional FBP reconstruction
switch geom.arg_save.type
case {'', 'dsc', 'std:mex', 'std:mat'}
	[image sino_filt] = fbp2_recon_std(sino, geom, opt);

case 'df,pull'
	if ~isempty(opt.window), error 'no window for DF', end
	image = fbp2_recon_df_pull(sino, geom, opt);
	sino_filt = [];

case 'mojette'
	if ~isempty(opt.window), error 'window only at setup for mojette', end
	image = fbp2_recon_moj(sino, geom.moj);
	sino_filt = [];

case 'table'
	error 'not done'
otherwise
	error 'type bug'
end


%
% fbp2_recon_std()
%
function [image, sino] = fbp2_recon_std(sino, geom, opt);

arg = geom.arg_save;

switch geom.sg.type

case 'par' % parallel-beam
	if isvar('geom.parallel_beam_parker_weight')
		sino = bsxfun(@times, sino, geom.parallel_beam_parker_weight);
	end
	sino = fbp2_sino_filter('flat', sino, ...
			'ds', geom.sg.dr, 'window', opt.window);

	switch arg.type
	case 'dsc'
		image = geom.G' * sino;
		image =	geom.scale * image;
	case 'std:mex'
		image = jf_mex('back2', geom.arg_back2{:}, single(sino));
	case 'std:mat'
		image = fbp2_back(geom.arg_back2{:}, single(sino));
	otherwise
		error 'bug'
	end

case 'fan' % fan-beam
	if geom.sg.dfs ~= 0 && ~isinf(geom.sg.dfs) % arc or flat fan-beam
		error 'only arc or flat fan done'
	end

	switch arg.type
	case 'std:mex'

		image = jf_mex('fbp,fan', geom.arg_fbp_fan{:}, opt.window, ...
			single(sino));
		image = image .* repmat(geom.ig.mask, [1 1 size(image,3)]);

	case 'std:mat'

		if isinf(geom.sg.dfs)
			dtype = 'flat';
		elseif geom.sg.dfs == 0
			dtype = 'arc';
		else
			fail('bad detector dfs %g', geom.sg.dfs)
		end
		sino = fbp2_sino_weight(geom.sg, sino);
		sino = fbp2_sino_filter(dtype, sino, ...
			'ds', geom.sg.ds, 'dsd', geom.sg.dsd, ...
			'window', opt.window);
		image = fbp2_back_fan(geom.sg, geom.ig, sino);
%		im(image), cbar

	otherwise
		fail('fan not done for type "%s"', arg.type)
	end


case 'moj' % mojette, from mojette-rebinned sino
	% filter mojette sinogram
	sino = fbp2_apply_sino_filter_moj(sino, geom.moj.H);

	if geom.sg.dx == abs(geom.ig.dx)
		image = geom.moj.G' * sino; % backproject
		image = image * (pi / geom.sg.na); % account for "dphi" in integral
	else % revert to conventional pixel driven
		ig = geom.ig;
		sg = geom.sg;
		arg1 = {uint8(ig.mask), ig.dx, ig.dy, ig.offset_x, ...
			sign(ig.dy) * ig.offset_y}; % trick: old backproject
		arg2 = {sg.d(1:sg.na), sg.offset, sg.orbit, sg.orbit_start};
		image = jf_mex('back2', arg1{:}, arg2{:}, ...
				int32(arg.nthread), single(sino));
		image = image .* geom.ig.mask;
	end
otherwise
	error 'not done'
end


%
% fbp2_recon_df_pull()
%
function image = fbp2_recon_df_pull(sino, st, arg);

sg = st.sg;
ig = st.ig;

%na = size(sino, 2);
f_sino = fft(sino, st.npad, 1);
f_sino = repmat(st.phase1, [1 sg.na]) .* f_sino;
f_sino = fftshift(f_sino, 1);

% trick: add another row at pi by flipping, to help interpolator
f_sino = [f_sino, [conj(f_sino(1,1)); flipud(f_sino(2:end,1))]];

% [-nu/2,nu/2] x [0,nv/2]:
f_xd = interp2(st.phi, st.rho, f_sino, st.qphi, st.qrho, st.interp);
if any(isnan(f_xd)), error 'nan', end
% f_xd(isnan(f_xd)) = 0; % fix:

f_xd = st.phase2 .* f_xd;

% edges of spectrum must be real for iDFT to be real
nu = st.nu; nv = st.nv;
f_xd(1+[0:2]*nu/2, nv/2+1) = real(f_xd(1+[0:2]*nu/2, nv/2+1));
f_xd(1+[0:2]*nu/2, 1) = real(f_xd(1+[0:2]*nu/2, 1));
f_xd(1:nu/2,nv/2+1) = conj(flipud(f_xd(nu/2+2:end,nv/2+1)));

% form entire spectrum using conjugate symmetry
t1 = fliplr(flipud(f_xd(2:end-1, 2:end))); % [-nu/2+1,nu/2-1] x [1,nv/2]
t1 = [fliplr(f_xd(1, 2:end)); t1]; % [-nu/2] x [1,nv/2] added
t2 = f_xd(1:end-1, 1:end-1); % [-nu/2,nu/2-1] x [0,nv/2-1]
f_xd = [conj(t1), t2];

% clf, dft_sym_check(fftshift(f_xd)), return

dx = ig.dx;
nx = ig.nx;
ny = ig.ny;

image = fftshift(ifft2(fftshift(f_xd))) / dx^2;
image = reale(image, 'warn');
if (st.nu > nx)
	image = image([0:nx-1]+1+(st.nu-nx)/2, [0:ny-1]+1+(st.nv-ny)/2);
elseif st.nu ~= nx
	error 'only 1x and 2x done'
end


%
% fbp2_recon_moj()
%
function image = fbp2_recon_moj(sino, moj);

% bin sinogram to mojette sampling, whether it is fan or parallel
msino = moj.ob_rebin * sino;

msino = fbp2_apply_sino_filter_moj(msino, moj.H); % filter mojette sinogram
image = moj.G' * msino; % backproject
image = image * (pi / size(sino,2)); % account for "dphi" in integral


%
% fbp2_make_sino_filter_moj()
% make filter for each sinogram row.  trick: dr varies by row
%
function H = fbp2_make_sino_filter_moj(nb, na, dx, orbit, orbit_start, window)
ang = deg2rad(orbit_start + [0:na-1]/na * orbit);
npad = 2^ceil(log2(2*nb-1)); % padded size

dr = dx * max(abs(cos(ang)), abs(sin(ang)));
if 1
	[junk H] = fbp2_sino_filter('flat', ones(nb,1), 'ds', 1, ...
		'window', window, 'decon1', 0);
	% trick: ramp filter usually has 1/dr^2 in it, but convolution by sum
	% requires "dr" so we need dr / dr^2 = 1 / dr
	H = H * (1 ./ dr); % [npad,na]

else % trick: make sure cutoff frequencies match
	u0 = 1/2/dx; % standard cutoff for coarsest sampling
	r = [-(npad/2):(npad/2-1)]' * dr; % [nb na]
	h = u0^2 * (2 * nufft_sinc(2*u0*r) - nufft_sinc(u0*r).^2);
	h = h .* repmat(dr, [npad 1]); % extra dr for discrete-space convolution
%	clf, plot(r, h, '.'), keyboard
	H = reale(fft(fftshift(h,1), [], 1));
%	H = fbp_apodize(H, nb, window);
	if ~isempty(window), error 'window not done yet due to dr', end
end


%
% fbp2_apply_sino_filter_moj()
% fix: could allow optional window here
%
function sino = fbp2_apply_sino_filter_moj(sino, H)
nb = size(sino,1);
npad = 2^ceil(log2(2*nb-1)); % padded size
sinopad = [sino; zeros(npad-nb,size(sino,2))]; % padded sinogram
sino = ifft_sym(fft(sinopad, [], 1) .* H, [], 1);
sino = sino(1:nb,:);
