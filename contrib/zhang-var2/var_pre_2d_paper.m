% source code to produce analytical predictions
% for 2d fast variance predictions paper for IEEE.
%
% Empirical and FFT_based predictions are computed and saved
%
% Copyright by Yingying Zhang-O'Connor, University of Michigan 2006

clear
clc

% -------------------------------
% set up parameters
% -------------------------------
printm 'set up parameters'
% calibration for single integral
% true to compute on the fly; false to load the saved table
t.calibrate = true;

t.ellipse = false;
t.zubal = false;
t.ncat = true;
t.fov = 'ellipse'; % 'ellipse' or 'disc'

t.const = false; % standard: \til{R}_phi = kappa_{center}
t.kappa = true; % certainty-based \til{R}_phi = kappa_j

if t.zubal
	img = 1;
	nx = 128/img;
	ny = nx;
	nb = ceil(888/4/img);
	na = ceil(984/4/img);
	f.ray_spacing = 0.9765625*4*img; %mm
	f.strip_width = 0.9765625*4*img;
	f.pixel_size = 0.9765625*4*img;
elseif t.ncat
	img = 2;
	nx = 512/img;
	ny = nx;
	nb = ceil(888/img);
	na = ceil(984/img);
	f.ray_spacing = 0.9765625*img; %mm
	f.strip_width = 0.9765625*img;
	f.pixel_size = 0.9765625*img;
	f.channel_offset = 0;
end

% fan-beam parameters
f.orbit = 360;
f.orbit_start = 0;
f.dis_src_det = 949.075;
f.dis_iso_det = 408.075;

% phantom parameters
if t.ellipse
	m.water = 0.2; % cm^2/g
	m.bone = 0.1;
	d.ellipse = 0.7*m.water/10; % g's g/cc; pixel size is mm, so convert into mm^(-1)
	d.bone = (2.5-1)*m.bone/10;
	d.soft = 1*m.water/10;
	c_x = (nx+1)/2;
	c_y = (ny+1)/2;
	bx_center = 0;
	by_center = 0;
	bx_radius = (nx - 38)/2;
	by_radius = 25;
	disc_radius = 10;
end

% ------------------------------
% generate the true phantom
% ------------------------------
printm 'set up xtrue'
xtrue = zeros(nx,ny);
if t.ellipse
	x_center = [bx_center; bx_center+20];
	y_center = [by_center; by_center];
	x_radius = [bx_radius; disc_radius];
	y_radius = [by_radius; disc_radius];
	ampl = [d.ellipse; d.bone];
	params = [x_center y_center x_radius y_radius zeros(2,1) ampl];
	params(:,1:4) = params(:,1:4)*f.pixel_size;
	xtrue = zeros(nx,ny);
	xtrue = ellipses(nx, ny, params, f.pixel_size, f.pixel_size);

elseif t.zubal
	dir = '~fessler/l/src/matlab/alg/data';
	xtrue = read_zubal_attn(dir) * 2.1; % mulist = [0.002 0.0096 0.0120] *2.1 to make CT value;

	if nx < 128
		xtrue = downsample2(xtrue,img);
	end

elseif t.ncat
%	xtrue = fld_read('~yzz/data/phantoms/ncat/ncat,256,slice,140.fld')/2;
%	ddir = path_find_dir([filesep 'transmission']);
%	ddir = strrep(ddir, 'transmission', 'data');
%	xtrue = fld_read([ddir '/ncat,256,slice,140,ct,x100.fld']) / 2;
	ddir = '~fessler/l/dat/phantom,digital/unc/ncat1.1';
	xtrue = fld_read([ddir '/ncat,256,slice,140.fld']) / 2;

	if (nx < 1024) && size(xtrue,2) == 1024
		img = 1024/nx;
		xtrue = downsample2(xtrue,img);
	end
end


% -------------------------------
% generate mask
% -------------------------------
printm 'generate support'
if t.ellipse
	x_mask = bx_radius+5;
	y_mask = by_radius+5;
elseif t.zubal || t.ncat
	if streq(t.fov, 'ellipse')
		printm('ellipse FOV')
		x_mask = (nx - nx/16 + 10)/2;
		y_mask = (ny - (ny/4 + nx/16) + 50)/2;
	elseif streq(t.fov, 'disc')
		printm('circular FOV')
		ss = ((-(nb-1)/2:(nb-1)/2).' - f.channel_offset) * f.ray_spacing;
		smax = max(abs(ss));
		rmax = sin(smax/f.dis_src_det)*(f.dis_src_det - f.dis_iso_det);
		x_mask = floor((rmax-2)/2);
		y_mask = x_mask;
	end
end
f.support = sprintf('ellipse 0 0 %g %g', x_mask, y_mask);

% new code by jf showing how to generate system matrix
sg = sino_geom('fan', ...
	'nb', nb, 'na', na, ...
	'dsd', f.dis_src_det, ...
	'dod', f.dis_iso_det, ...
	'ds', f.ray_spacing, ...
	'strip_width', f.strip_width, ...
	'orbit', f.orbit, ...
	'orbit_start', f.orbit_start, ...
	'offset_s', 0);

ig = image_geom('nx', nx, 'ny', ny, 'dx', f.pixel_size);
ig.mask = ig.circ(max(x_mask, y_mask));
if x_mask ~= y_mask
	warn 'need to put in ellipse using ellipse_im'
end

G = Gtomo2_wtmex(sg, ig);

if 0 % this is an old way - do not use it
	sys = arg_pair('system', 14, ...
	'nx', nx, 'ny', ny, ...
	'nb', nb, 'na', na, ...
	'support',	f.support, ...
	'pixel_size',	f.pixel_size, ...
	'ray_spacing',	f.ray_spacing, ...
	'strip_width',	f.strip_width, ...
	'orbit',	f.orbit, ...
	'orbit_start',	f.orbit_start, ...
	'src_det_dis',	f.dis_src_det, ...
	'obj2det_x',	f.dis_iso_det, ...
	'obj2det_y',	f.dis_iso_det, ...
	'flip_y', 1, ...
	'scale', 0, ...
	'channel_offset', f.channel_offset, ...
	'source_offset', 0);
	Gdsc = Gtomo2_dscmex(sys, 'nthread', 2); % use it to generate the mask
	mask = Gdsc.arg.mask;
end

% -------------------------------
% system matrix
% -------------------------------
printm 'set up fan-beam system using DD projector'

if 0 % this is an old way - do not use it
	G = Gtomo2_dd(mask, nb, na, ...
	'pixel_size', f.pixel_size, ...
	'ray_spacing', f.ray_spacing, ...
	'strip_width', f.strip_width, ...
	'dis_src_det', f.dis_src_det, ...
	'dis_iso_det', f.dis_iso_det, ...
	'orbit', f.orbit, ...
	'orbit_start', f.orbit_start, ...
	'xscale', 1, ...
	'yscale', 1, ...
	'channel_offset', 0);
end

% ---------------------------------
% compute statistical weighting wi
% ---------------------------------
printm 'compute statistical weighting wi'
lb = reshape(G*xtrue(mask(:)), [nb na]);
bi = 10^6;
yb = bi*exp(-lb);
yi = poissrnd(yb);
wi = yb; % this is the true case. In reality no yb. Use yi plug-in method
lhat = -log(yi./bi);

W = diag_sp(wi(:));

% ---------------------------------
% choose \alpha for all penalty
% ---------------------------------
% regularization parameter
if nx == 128
	f.l2b = 12; % for 128x128
elseif nx == 32
	f.l2b = 14;
elseif (nx == 512) || (nx == 256)
	f.l2b = 11;
else
	f.l2b = 17;
end
f.beta = 2^f.l2b;

% ----------------------------------------
% Compute the H for all penalty cases
% ----------------------------------------
dx = f.pixel_size;
dy = dx;
np = G.arg.np; % number of pixels inside mask
ds = f.ray_spacing;
dp = 1; % distance power
dbeta = 2*pi/na;

printm('Polar frequency coordinates')
% decide rhomax
rhomax = 1 / (2 * dx);
drho = 1/nx/dx;
rho1 = (-nx/2:nx/2-1) / (nx * dx);
rho = rho1(nx/2 + 1:end);
rho = rho(:);
nrho = length(rho);
nPhi = na/4;
Phi = ((0:nPhi-1)')/nPhi * 2 * pi;
dPhi = abs(Phi(2)-Phi(1));

printm('compute w_0(Phi): angular dependent weighting')
nw = na/4;
phiw = ((0:nw-1)')/nw * pi; % to avoid = pi/2 since tan(pi/2) = inf
dphiw = abs(phiw(2)-phiw(1));
% w_0(Phi)
[wknot, tw] = adw_fan(G.arg, wi, mask, Phi);
sprintf('w_0(Phi) pre-compute: %g seconds', tw)

%
% pre-compute G(\rho, \Phi) in H(\rho, \Phi) at image center
%
% Jacobian: J(sj) \approx J(s_center) at image center: Ds0/Dsd
Jj = G.arg.dis_src_iso/G.arg.dis_src_det;

if streq(t.fov, 'ellipse')
	if 1 % d_0(\phi) at image center
		fov = x_mask .* y_mask ./ sqrt((x_mask .* cos(phiw)).^2 ...
			+ (y_mask .* sin(phiw)).^2);
		fov = 2 .* fov(:) .* dx;
	end

elseif streq(t.fov, 'disc')
	fov = sum(mask(:,nx/2)) .* dx;
	fov = repmat(fov, [length(phiw) 1]);
end

tt = cputime;
G_center = zeros(nrho, nPhi);
for ia = 1:length(phiw)
	G_center = G_center + fov(ia).* sinc(rho * (fov(ia) .* sin(Phi-phiw(ia))).').^2;
end
sprintf('G_0(rho, Phi) at center pre-compute: %g seconds', cputime - tt)
% A_0(\rho,\Phi) at image center
A = dx .* sinc(dx .* (rho * cos(Phi).')) .* ...
	dy .* sinc(dx .* (rho * sin(Phi).')) .* ...
	repmat((sinc(ds .* Jj .* rho)), [1 length(Phi)]); % Nrho x Phi

Ga = G_center .* (A.^2) .* dphiw ./ abs(ds .* dx .* dy .* dbeta);

% -----------------------------------------
% Compute the approx standard deviation map
% -----------------------------------------
if t.const || t.kappa

	kappan = G' * wi(:);
	kappad = G' * (ones(size(wi(:))));
	kappa = sqrt( kappan ./ kappad );

	f.data_dir = 'data/standard_penal/256image/';

	if t.const
	printm('standard quadratic penalty')
	printm('kappa is a const for all pixels')
	kappa_r = embed(kappa, mask);
	kappa_cen = kappa_r(nx/2+1, ny/2+1) .* ones(size(kappa));
	kappa_v = kappa; % varies for pixel

	his = f.l2b + log2(kappa_v.^2) - log2(kappa_cen.^2);
	his = hist(his, (6:0.5:14).'); % histogram
	clear kappa_r

	kappa = kappa_cen;
	kappa_r = embed(kappa, mask);

	% empirical and fft_based for comparison
	f.iot = [f.data_dir 'std_iot,iter71,ntime250,l2b11,256,nblock41,const.fld'];
	std_iot = fld_read([f.data_dir 'std_iot,iter71,ntime250,l2b11,256,nblock41,const.fld']);
	std_fft = fld_read([f.data_dir 'std_vh,dbsized_mask,iter71,ntime250,l2b11,256,nblock41,const.fld']);

	elseif t.kappa
	printm('certainty-base quadratic penalty')
	printm('kappa varies for pixels')
	kappa_r = embed(kappa, mask);

	% empirical and fft_based for comparison
%	std_iot = fld_read([f.data_dir 'std_iot,iter71,ntime250,l2b11,256,nblock41,kappa.fld']);
	std_iot = zeros(size(kappa_r)); % fake jf
%	std_fft = fld_read([f.data_dir 'std_vh,dbsized_mask,iter71,ntime250,l2b11,256,nblock41,kappa.fld']);
	std_fft = zeros(size(kappa_r)); % fake jf
	end

	r_sum_con = 1 + sqrt(2).^(2-dp);
	nl = [1 0 1 1]; ml = [0 -1 -1 1];
	phil = atan2(ml,nl);

	uu = rho * cos(Phi).';
	vv = rho * sin(Phi).';
	R1 = 0;
	for il = 1:4
	Rl = 4 .* sin((uu .* nl(il) + vv .* ml(il))*2*pi*dx/2).^2 ./ (sqrt(nl(il)^2 + ml(il)^2)).^dp;
	R1 = R1 + Rl; % rho x Phi
	end

	printm('double summation')
	tt = cputime;
	for ia = 1:size(wknot,2)
	ticker(mfilename, ia, size(wknot,2))
	Hd = wknot(:,ia) * Ga(:,ia).'; % Np x Nrho
	Hu = wknot(:,ia) * (Ga(:,ia) .* rho).'; % Np x Nrho

	R = kappa.^2 * R1(:,ia).'; % Np x Nrho
	int_rho(:,ia) = rhomax .* mean(Hu ./ ((Hd + f.beta .* R).^2), 2);
	end
	var_pre = mean(int_rho,2) .* (2*pi) .* (dx .* dy);
	std_pre = embed(sqrt(var_pre), mask);
	sprintf('db-integral prediction: %g second', cputime - tt)

	if t.calibrate
	% use extended mask for FFT in calibration
	maskcal = true([nx ny]);
	Gcal = Gtomo2_dd(maskcal, nb, na, ...
	'pixel_size', f.pixel_size, ...
	'ray_spacing', f.ray_spacing, ...
	'strip_width', f.strip_width, ...
	'dis_src_det', f.dis_src_det, ...
	'dis_iso_det', f.dis_iso_det, ...
	'orbit', f.orbit, ...
	'orbit_start', f.orbit_start, ...
	'xscale', 1, ...
	'yscale', 1, ...
	'channel_offset', 0);

	c.l2b = [2:0.15:20].'; % column
	c.beta = 2.^c.l2b;
	for indx = 1:length(c.beta)
	c.R0 = Robject(maskcal, 'edge_type', 'tight', 'beta', c.beta(indx), ...
	'type_denom', 'matlab', 'potential', 'quad', 'distance_power', dp);
	% mask size matters for smaller \beta
	[psf0, c.var0(indx), fwhm_t(indx)] = qpwls_psf(Gcal, c.R0, 1, maskcal, [], [0 0]);
	end
	c.var0 = c.var0(:); % column
	fwhm_t = fwhm_t(:);

	ej = zeros([nx ny]);
	ej(nx/2+1,ny/2+1) = 1;
	center_in = find(ej(mask(:))==1);
	[wknotcal, tw] = adw_fan(G.arg, ones(size(wi)), mask, Phi);
	CC = (dx .* dy) ./ (ds);
	KK = CC/(rhomax)^3;
	for indx = 1:length(c.beta)
	c.var_pre(indx) = mean(((dx .* dy) .* (2*pi) ./ 3) ./ ...
	(KK .* wknotcal(center_in,:) .* 2 + (2*pi*dx)^2*c.beta(indx)*r_sum_con));
	end
	c.var_pre = c.var_pre(:);

	ratio = sqrt(c.var0(:) ./ c.var_pre);

%	fld_write('data/cal,var_fft,fwhm,qpuls,l2b2to20,inl2b,2w0.fld', [c.var0 fwhm_t]);
%	fld_write('data/cal,ratio,in_std,qpuls,l2b2to20,inl2b,2w0.fld', ratio)
	else
		% load the calibration scale factor needed for single integral
		% Calibration is done for a range of fwhms
		% about (1.13)^2 in our case
		tmp = fld_read('data/cal,var_fft,fwhm,qpuls,l2b2to20,inl2b,2w0.fld');
		fwhm_t = tmp(:,2);
		ratio = fld_read('data/cal,ratio,in_std,qpuls,l2b2to20,inl2b,2w0.fld');
	end
	Rr = Robject(kappa_r, 'edge_type', 'tight', 'beta', f.beta, ...
		'type_denom', 'matlab', 'potential', 'quad', 'distance_power', dp);
	[psf, var, fwhm] = qpwls_psf(G, Rr, 1, mask, W, [0 0]);
	f_loc = find(abs(fwhm - fwhm_t) == min(abs(fwhm - fwhm_t)));

	CC = (dx .* dy) ./ (ds);
	KK = CC/(rhomax)^3;
	tt = cputime;
	var_pre1 = mean(1 ./ (KK .* wknot * 2 + (2*pi*dx)^2*f.beta*r_sum_con .* repmat(kappa.^2,[1 length(Phi)])), 2);
	std_pre1 = embed(sqrt(var_pre1 .* (dx .* dy) .* (2*pi) ./ 3), mask);
	std_pre1 = std_pre1 .* ratio(f_loc);
	sprintf('one-integral prediction: %g second', cputime - tt)
end

% jf: display results
im pl 2 2
im(1, std_pre), cbar
im(2, std_pre1), cbar
