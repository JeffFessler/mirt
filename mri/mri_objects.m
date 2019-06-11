  function st = mri_objects(varargin)
%|function st = mri_objects([options], 'type1', params1, 'type2', params2, ...)
%| Generate strum that describes image-domain objects and Fourier domain spectra
%| of simple structures such as rectangles, disks, superpositions thereof.
%| These functions are useful for simple "idealized" MRI simulations
%| where the data is modeled as analytical Fourier samples,
%| i.e., no field inhomogeneity and no relaxation effects.
%|
%| in
%|	(type, params)	e.g. 'rect', params, 'gauss', params, ...
%|
%| options
%|	'fov'	[1] | [Nd]	field of view (needed only for some tests)
%|
%| params:
%|	'dirac2'	[N 3]	[xcent ycent value]
%|	'dirac3'	[N 4]	[xcent ycent zcent value]
%|	'rect2'		[N 5]	[xcent ycent xwidth ywidth value]
%|	'rect3'		[N 7]	[xcent ycent zcent xwidth ywidth zwidth value]
%|	'gauss2'	[N 5]	[xcent ycent xwidth ywidth value]
%|	'gauss3'	[N 7]	[xcent ycent zcent xwidth ywidth zwidth value]
%|	'circ2'		[N 5]	[xcent ycent rad value]
%|	'cyl3'		[N 6]	[xcent ycent zcent xrad zwidth value]
%|
%| out
%|	st	strum object with methods:
%|	st.image(x,y)		returns 2D image-domain (sampled) picture
%|	st.image(x,y,z)		same but 3D
%|	st.kspace(u,v)		returns 2D Fourier-space samples
%|	st.kspace(u,v,w)	same but 3D
%|
%| Copyright 2007-6-28, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error args, end
if nargin == 1 && streq(varargin{1}, 'test'), mri_objects_test; return, end

fov = [];

while length(varargin)
	switch varargin{1}
	case 'fov'
		fov = varargin{2};
		varargin = {varargin{3:end}};
	case 'rect2half'
		st = mri_objects('rect2', [0 0 fov/2 fov/2 1]);
		return
	case 'rect3half'
		st = mri_objects('rect3', [0 0 0 fov/2 1]);
		return
	case 'case1'
		st = mri_objects_case1(fov, varargin{2:end});
		return
	case 'test4'
		tmp = mri_objects_test4(fov, varargin{2:end});
		st = mri_objects(tmp{:});
		return
	otherwise
		break
	end
end

if rem(length(varargin), 2), fail('bad args'), end
st.types = {varargin{1:2:end}};
st.params = {varargin{2:2:end}};

is3 = ~isempty(strfind(col(char(st.types)')', '3')); % any 3D?

if is3
	st = strum(st, { ...
		'image', @mri_objects_image3, '(x,y,z)';
		'kspace', @mri_objects_kspace3, '(u,v,w)';
		});
else
	st = strum(st, { ...
		'image', @mri_objects_image2, '(x,y)';
		'kspace', @mri_objects_kspace2, '(u,v)';
		});
end


%
% mri_objects_image2()
% x,y,z must all be the same size
%
function out = mri_objects_image2(st, x, y)

st_params = st.params; % fix: subsref limitation
st_types = st.types; % fix: subsref limitation

out = zeros(size(x));

for ii=1:length(st_types)
	params = st_params{ii};
	switch st_types{ii}
	case 'circ2'
		out = out + mri_objects_image_circ2(params, x, y);
	case 'dirac2'
		out = out + mri_objects_image_dirac2(params, x, y);
	case 'gauss2'
		out = out + mri_objects_image_gauss2(params, x, y);
	case 'rect2'
		out = out + mri_objects_image_rect2(params, x, y);
	otherwise
		fail('unknown object type %s', st_types{ii})
	end
end


%
% mri_objects_kspace2()
% x,y,z must all be the same size
%
function out = mri_objects_kspace2(st, u, v)

st_params = st.params; % fix: subsref limitation
st_types = st.types; % fix: subsref limitation

out = zeros(size(u));

for ii=1:length(st_types)
	params = st_params{ii};
	switch st_types{ii}
	case 'circ2'
		out = out + mri_objects_kspace_circ2(params, u, v);
	case 'dirac2'
		out = out + mri_objects_kspace_dirac2(params, u, v);
	case 'gauss2'
		out = out + mri_objects_kspace_gauss2(params, u, v);
	case 'rect2'
		out = out + mri_objects_kspace_rect2(params, u, v);
	otherwise
		fail('unknown object type %s', st_types{ii})
	end
end



%
% mri_objects_image3()
% x,y,z must all be the same size
%
function out = mri_objects_image3(st, x, y, z)

st_params = st.params; % fix: subsref limitation
st_types = st.types; % fix: subsref limitation

out = zeros(size(x));

for ii=1:length(st_types)
	params = st_params{ii};
	switch st_types{ii}
	case 'cyl3'
		out = out + mri_objects_image_cyl3(params, x, y, z);
	case 'dirac3'
		out = out + mri_objects_image_dirac3(params, x, y, z);
	case 'gauss3'
		out = out + mri_objects_image_gauss3(params, x, y, z);
	case 'rect3'
		out = out + mri_objects_image_rect3(params, x, y, z);
	otherwise
		fail('unknown object type %s', st_types{ii})
	end
end


%
% mri_objects_kspace3()
% x,y,z must all be the same size
%
function out = mri_objects_kspace3(st, u, v, w)

st_params = st.params; % fix: subsref limitation
st_types = st.types; % fix: subsref limitation

out = zeros(size(u));

for ii=1:length(st_types)
	params = st_params{ii};
	switch st_types{ii}
	case 'cyl3'
		out = out + mri_objects_kspace_cyl3(params, u, v, w);
	case 'dirac3'
		out = out + mri_objects_kspace_dirac3(params, u, v, w);
	case 'gauss3'
		out = out + mri_objects_kspace_gauss3(params, u, v, w);
	case 'rect3'
		out = out + mri_objects_kspace_rect3(params, u, v, w);
	otherwise
		fail('unknown object type %s', st_types{ii})
	end
end


%
% mri_objects_image_dirac2()
%
function out = mri_objects_image_dirac2(params, x,y)
if ncol(params) ~= 3, fail('dirac2 requires 3 parameters'), end
params = [params(:,1:2) zeros(nrow(params),1) params(:,3)];
out = mri_objects_image_dirac3(params, x,y,0);

%
% mri_objects_image_dirac3()
% param: [N,4] [xcent ycent zcent value]
%
function out = mri_objects_image_dirac3(params, x,y,z)
out = 0;
if ncol(params) ~= 4, fail('dirac3 requires 4 parameters'), end
for ii=1:nrow(params)
	par = params(ii,:);
	out = out + par(4) * (x == par(1)) .* (y == par(2)) .* (z == par(3));
end
out(out~=0) = inf;
warn 'image of Dirac impulses is invalid'


%
% mri_objects_kspace_dirac2()
%
function out = mri_objects_kspace_dirac2(params, u,v)
if ncol(params) ~= 3, fail('dirac2 requires 3 parameters'), end
params = [params(:,1:2) zeros(nrow(params),1) params(:,3)];
out = mri_objects_kspace_dirac3(params, u,v,0);

%
% mri_objects_kspace_dirac3()
% param: [N,4] [xcent ycent zcent value]
%
function out = mri_objects_kspace_dirac3(params, u,v,w)
out = 0;
if ncol(params) ~= 4, fail('dirac3 requires 4 parameters'), end
for ii=1:nrow(params)
	par = params(ii,:);
	out = out + par(4) * exp(-2i*pi*(u*par(1)+v*par(2)+w*par(3)));
end


%
% mri_objects_image_circ2()
%
function out = mri_objects_image_circ2(params, x,y)
if ncol(params) ~= 4, fail('circ2 requires 4 parameters'), end
z = zeros(nrow(params),1);
params = [params(:,1:2) z params(:,3) 1+z params(:,4)];
out = mri_objects_image_cyl3(params, x,y,0);


%
% mri_objects_image_cyl3()
% param: [N,6] [xcent ycent zcent rad zlen value]
%
function out = mri_objects_image_cyl3(params, x,y,z)
out = 0;
if ncol(params) ~= 6, fail('cyl3 requires 6 parameters'), end
for ii=1:nrow(params)
	par = params(ii,:);
	xc = par(1);
	yc = par(2);
	zc = par(3);
	rad = par(4);
	len = par(5);
	out = out + par(6) * ((x-xc).^2+(y-yc).^2 < rad^2) .* rect((z-zc)/len);
end


%
% mri_objects_kspace_circ2()
%
function out = mri_objects_kspace_circ2(params, u,v)
if ncol(params) ~= 4, fail('circ2 requires 4 parameters'), end
z = zeros(nrow(params),1);
params = [params(:,1:2) z params(:,3) 1+z params(:,4)];
out = mri_objects_kspace_cyl3(params, u,v,0);

%
% mri_objects_kspace_cyl3()
% param: [N,6] [xcent ycent zcent rad zlen value]
% note: circ(r) = rect(r/2) <=> 4*jinc(2*q)
%
function out = mri_objects_kspace_cyl3(params, u,v,w)
out = 0;
if ncol(params) ~= 6, fail('cyl3 requires 6 parameters'), end
for ii=1:nrow(params)
	par = params(ii,:);
	xc = par(1);
	yc = par(2);
	zc = par(3);
	rad = par(4);
	len = par(5);
	out = out + par(6) * rad^2 * 4*jinc(2*sqrt(u.^2+v.^2)*rad) ...
		.* (len * nufft_sinc(w*len)) ...
		.* exp(-2i*pi*(u*xc + v*yc + w*zc));
end


%
% mri_objects_image_gauss2()
%
function out = mri_objects_image_gauss2(params, x,y)
if ncol(params) ~= 5, fail('gauss2 requires 5 parameters'), end
z = zeros(nrow(params),1);
zw = inf(size(z)); % trick
params = [params(:,1:2) z params(:,3:4) zw params(:,5)];
out = mri_objects_image_gauss3(params, x,y,0);

%
% mri_objects_image_gauss3()
% param: [N,7] [xcent ycent zcent xfwhm yfwhm zfwhm value]
%
function out = mri_objects_image_gauss3(params, x,y,z)
out = 0;
if ncol(params) ~= 7, fail('gauss3 requires 7 parameters'), end
for ii=1:nrow(params)
	par = params(ii,:);
	xc = par(1);
	yc = par(2);
	zc = par(3);
	xw = par(4) / sqrt(log(256)) * sqrt(2*pi);
	yw = par(5) / sqrt(log(256)) * sqrt(2*pi);
	zw = par(6) / sqrt(log(256)) * sqrt(2*pi);
	out = out + par(7) ...
		.* exp(-pi * ((x-xc)/xw).^2) ...
		.* exp(-pi * ((y-yc)/yw).^2) ...
		.* exp(-pi * ((z-zc)/zw).^2);
end


%
% mri_objects_kspace_gauss2()
%
function out = mri_objects_kspace_gauss2(params, u,v)
if ncol(params) ~= 5, fail('gauss2 requires 5 parameters'), end
z = zeros(nrow(params),1);
zw = (1+z) * sqrt(log(256)) / sqrt(2*pi); % trick
params = [params(:,1:2) z params(:,3:4) zw params(:,5)];
out = mri_objects_kspace_gauss3(params, u,v,0);


%
% mri_objects_kspace_gauss3()
% param: [N,7] [xcent ycent zcent xfwhm yfwhm zfwhm value]
%
function out = mri_objects_kspace_gauss3(params, u,v,w)
out = 0;
if ncol(params) ~= 7, fail('gauss3 requires 7 parameters'), end
for ii=1:nrow(params)
	par = params(ii,:);
	xc = par(1);
	yc = par(2);
	zc = par(3);
	xw = par(4) / sqrt(log(256)) * sqrt(2*pi);
	yw = par(5) / sqrt(log(256)) * sqrt(2*pi);
	zw = par(6) / sqrt(log(256)) * sqrt(2*pi);
	out = out + par(7) * xw*yw*zw ...
		.* exp(-pi * (u*xw).^2) ...
		.* exp(-pi * (v*yw).^2) ...
		.* exp(-pi * (w*zw).^2) ...
		.* exp(-2i*pi*(u*xc + v*yc + w*zc));
end


%
% mri_objects_image_rect2()
%
function out = mri_objects_image_rect2(params, x,y)
if ncol(params) ~= 5, fail('rect2 requires 5 parameters'), end
z = zeros(nrow(params),1);
params = [params(:,1:2) z params(:,3:4) 1+z params(:,5)];
out = mri_objects_image_rect3(params, x,y,0);


%
% mri_objects_image_rect3()
% param: [N,7] [xcent ycent zcent xw yw zw value]
%
function out = mri_objects_image_rect3(params, x,y,z)
out = 0;
if ncol(params) ~= 7, fail('rect3 requires 7 parameters'), end
for ii=1:nrow(params)
	par = params(ii,:);
	xc = par(1);
	yc = par(2);
	zc = par(3);
	xw = par(4);
	yw = par(5);
	zw = par(6);
	out = out + par(7) ...
		.* rect((x-xc)/xw) ...
		.* rect((y-yc)/yw) ...
		.* rect((z-zc)/zw);
end


%
% mri_objects_kspace_rect2()
%
function out = mri_objects_kspace_rect2(params, u,v)
if ncol(params) ~= 5, fail('rect2 requires 5 parameters'), end
z = zeros(nrow(params),1);
params = [params(:,1:2) z params(:,3:4) 1+z params(:,5)];
out = mri_objects_kspace_rect3(params, u,v,0);


%
% mri_objects_kspace_rect3()
% param: [N,7] [xcent ycent zcent xw yw zw value]
%
function out = mri_objects_kspace_rect3(params, u,v,w)
out = 0;
if ncol(params) ~= 7, fail('rect3 requires 7 parameters'), end
for ii=1:nrow(params)
	par = params(ii,:);
	xc = par(1);
	yc = par(2);
	zc = par(3);
	xw = par(4);
	yw = par(5);
	zw = par(6);
	out = out + par(7) * xw*yw*zw ...
		.* nufft_sinc(u*xw) ...
		.* nufft_sinc(v*yw) ...
		.* nufft_sinc(w*zw) ...
		.* exp(-2i*pi*(u*xc + v*yc + w*zc));
end


%
% mri_objects_case1()
%
function out = mri_objects_case1(fov, arg)

rp = [ ...	% rect2
	0 0		200 200	1;
	-50 -50		40 40	1;
	50 -50		20 20	1;
	0 50		50 50	1;
];

gp = [ ...	% gauss2 bumps
	-70 0		1 1 1;
	-60 0		2 2 1;
	-50 0		3 3 1;
	-40 0		4 4 1;
	-20 0		5 5 1;
	00 0		6 6 1;
	20 0		7 7 1;
	50 0		8 8 1;
];

if isvar('arg') && streq(arg, 'cm')
	rp(:,1:6) = rp(:,1:6) / 10;
	gp(:,1:6) = gp(:,1:6) / 10;
end

out = mri_objects('rect2', rp, 'gauss2', gp);


%
% mri_objects_test4()
%
function out = mri_objects_test4(fov, arg)

cp = [0 0 0 fov(1)*0.4 fov(3)*1.0 1]; % cyl3

rp = [ ...	% rect3
	-50 -50   0	40 40 40	1;
	 50 -50  40	20 20 50	1;
	  0  50 -40	30 30 60	1;
];
rp(:,1:3) = rp(:,1:3)/256 .* repmat(fov, 3, 1);
rp(:,4:6) = rp(:,4:6)/256 .* repmat(fov, 3, 1);

gp = [ ...	% gauss3 bumps
	-70 0 0		1 1 1	1;
	-60 0 0		2 2 2	1;
	-50 0 0		3 3 3	1;
	-40 0 0		4 4 4	1;
	-20 0 0		5 5 5	1;
	00 0 0		6 6 6	1;
	20 0 0		7 7 7	1;
	50 0 0		8 8 8	1;
];
gp(:,1:3) = gp(:,1:3)/256 .* repmat(fov, 8, 1);
gp(:,4:6) = gp(:,4:6)/256 .* repmat(fov, 8, 1);

if isvar('arg') && streq(arg, 'cm')
	rp(:,1:6) = rp(:,1:6) / 10;
	gp(:,1:6) = gp(:,1:6) / 10;
end

out = {'cyl3', cp, 'rect3', rp, 'gauss3', gp};


%
% mri_objects_test2()
% 2d tests
%
function mri_objects_test2
ig = image_geom('nx', 2^7, 'offsets', 'dsp', 'dx', 4);

if 1 % test transforms of each object type
	shift = [0.1 0.2] .* ig.fovs;
	sizes = [0.15 0.1] .* ig.fovs;
	tests = {'circ2', [shift sizes(1) 2];
		'gauss2', [shift sizes 2];
		'rect2', [shift sizes 2];
		};
	im plc 3 3
	for ii=1:nrow(tests)
		otype = tests{ii,1};
		param = tests{ii,2};
		st = mri_objects(otype, param);
		i2 = st.image(ig.xg, ig.yg);
		s2 = fftshift(fftn(fftshift(i2)));
		s2 = s2 * abs(ig.dx * ig.dy);
		fg = ig.fg;
		f2 = st.kspace(fg{:});
		max_percent_diff(f2, s2)

		im(1+3*(ii-1), i2, otype), cbar
		im(2+3*(ii-1), abs(s2), 'fft'), cbar
		im(3+3*(ii-1), abs(f2), 'kspace'), cbar
	end
prompt
end

st = mri_objects('case1');
xt = st.image(ig.xg, ig.yg);
im clf, im(xt), cbar
prompt



%
% mri_objects_test3()
%
function mri_objects_test3
ig = image_geom('nx', 2^7, 'nz', 2^5, 'offsets', 'dsp', 'dx', 4, 'dz', 3);

if 1 % test transforms of each object type
	shift = [0.1 0.2 0.3] .* ig.fovs;
	sizes = [0.15 0.1 0.2] .* ig.fovs;
	tests = {'cyl3', [shift sizes([1 3]) 2];
		'gauss3', [shift sizes 2];
		'rect3', [shift sizes 2];
		};

	im plc 3 3
	for ii=1:nrow(tests)
		otype = tests{ii,1};
		param = tests{ii,2};
		st = mri_objects(otype, param);
		i3 = st.image(ig.xg, ig.yg, ig.zg);
		s3 = fftshift(fftn(fftshift(i3)));
		s3 = s3 * abs(ig.dx * ig.dy * ig.dz);
		fg = ig.fg;
		f3 = st.kspace(fg{:});
		max_percent_diff(f3, s3)
		im(1+3*(ii-1), i3), cbar
		im(2+3*(ii-1), abs(s3), 'fft'), cbar
		im(3+3*(ii-1), abs(f3), 'kspace'), cbar
	end
prompt
end

st = mri_objects('fov', ig.fovs, 'test4');
xt = st.image(ig.xg, ig.yg, ig.zg);
im clf, im(xt), cbar


%
% mri_objects_test()
%
function mri_objects_test
mri_objects_test2
mri_objects_test3
