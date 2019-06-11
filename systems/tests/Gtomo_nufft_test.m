% Gtomo_nufft_test.m
% Test the Gtomo_nufft object

% create Gtomo_nufft class object
if ~isvar('An'), printm 'setup'
	ig = image_geom('nx', 64, 'ny', 60, 'fov', 480);
	ig.mask = ig.circ > 0;
	if 1
		sg = sino_geom('par', 'nb', 70, 'na', 80, 'dr', 7, ...
			'orbit_start', -15, 'strip_width', 'd');
		fan_arg = {}; % parallel beam
	else
		sg = sino_geom('ge1', 'orbit_start', -15, 'down', 8, ...
			'strip_width', 'd');
		fan_arg = { % fan beam
			'offset_s', sg.offset_s, ...
			'dis_src_det', sg.dsd, ...
			'dis_iso_det', sg.dod, ...
		};
	end

	% test non-rect image basis function and beam shape (S. Matej)
	if 0
		basis.type  = 'KB';
		basis.diam  = 4;
		basis.shape = 10.4;
		basis.m     = 2;
		basis.dim   = 2;
		basis.kernel=[];

		% KB: FWHM~1.0 : J=6, alpha=40, m=2
		beam.type  = 'KB';
		beam.diam  = 6;
		beam.shape = 40.;
		beam.m     = 2;
	end

	cpu etic

%{
	arg_geom = { % for old version
		fan_arg{:}, ...
		'orbit', sg.orbit, ...
		'orbit_start', sg.orbit_start, ...
		'dx', ig.dx, ...
		'ds', sg.d, ...
		'yscale', -ig.dy/ig.dx, ...
		'strip_width', sg.strip_width, ...
	};
%}

	arg = {
		'interp', {'table', 2^11, 'minmax:kb'}, ...
...%		'basis', basis, ...
...%		'beam', beam, ...
		'is.complex', 0, ...
...%		'kaiser', ...	% use KB interpolator
...%		'uniform', ...
...%		'bilin', ...
...%		'is.test', 1, ...
		'class', 'fatrix2', ...
	};

%	Ao = Gtomo_nufft(ig.mask, [sg.nb sg.na], arg_geom{:}, arg{:});
	An = Gtomo_nufft_new(sg, ig, arg{:});
%	jf_equal(struct(Ao), struct(An)), return

	cpu etoc 'An pre time:'
end

%
% make test image
%
if 0
%	x = ig.zeros;
%	x(12, 8) = 1;
%	x((nx/4+1):(3*nx/4),(ny/4+1):(3*ny/4)) = 1;
	x = ig.unitv(ig.nx/4,ig.ny/3);	% point source
elseif 1
	x = ellipse_im(ig, [0 0 [1 1]*ig.ny*ig.dx*0.9/2 0 1], 'oversample', 2);
	x(end/4:end/2,end/4:end/2) = 1.25;
else
	ix = [-(nx-1)/2:(nx-1)/2]/nx;
	iy = [-(ny-1)/2:(ny-1)/2]/ny;
	[ix iy] = ndgrid(ix, iy);
	x = exp(-(ix.^2 + iy.^2) / 10^2);
end

% check real / imaginary
if 1
	y = An * x(ig.mask);
	y = sg.shape(y);
	im plc 3 3
	im(1, x, 'x')
	im(2, real(y), 'real(An*x)'), cbar
	im(3, imag(y), 'imag(An*x)'), cbar
prompt
end

% DSFT version for comparison
if ~isvar('Ad'), printm 'Ad'
%	Ao = Gtomo_nufft(ig.mask, [sg.nb sg.na], arg_geom{:}, arg{:}, ...
%		'is.dsft2', 1, 'is.dsft1', 1);
	Ad = Gtomo_nufft_new(sg, ig, arg{:}, 'is.dsft2', 1, 'is.dsft1', 1);
%	jf_equal(struct(Ao), struct(Ad)), return
end

% create a "strip integral" system object
if ~isvar('As'), printm 'As'
%	Ao = aspire_pair(sg, ig, 'strip_width', An.arg.strip_width, ...
%		'support', 'all');
%		'support', ig.mask);
%	Ao = Gtomo2_dscmex(Ao, 'nthread', 2); % multiprocessor!
	As = Gtomo2_dscmex(sg, ig, 'nthread', 2); % multiprocessor!
%	jf_equal(struct(Ao), struct(As)), return
end

% compare back projectors
if 0
	y = sg.ones;
	y = sg.unitv(sg.nb/2+7, 9);
%	y = yd / 100;
	xn = An' * y;
	xd = Ad' * y;
	xs = As' * y;
	im(4, xn, 'An''*y'), cbar
	im(5, xd, 'Ad''*y'), cbar
	im(6, xs, 'As''*y'), cbar
	im(8, xn-xd, 'An''*y - Ad''*y'), cbar
	im(9, xn-xs, 'An''*y - As''*y'), cbar
%	printf('back nrms = %g%%', nrms(x2(:), x1(:)) * 100)
prompt
end

% compare projectors
if 1
	yn = An * x;
	yd = Ad * x;
	ys = As * x;
	im(4, yd, 'Ad*x'), cbar
	im(5, yn-yd, 'An*x-Ad*x'), cbar
	im(7, ys, 'As*x'), cbar
	im(8, yn-ys, 'An*x-As*x'), cbar
%	printf('forward nrms = %g%%', nrms(y2(:), y1(:)) * 100)
prompt
end

if 1, printm 'single ray'
	y = sg.unitv(sg.nb/2+9,1);
	im pl 2 3
	im(1, As'*y, 'back ray strip'), cbar h
	im(2, An'*y, 'back ray nufft'), cbar h
	im(3, Ad'*y, 'back ray dsft'), cbar h
	im(4, An'*y - Ad'*y, 'nufft - dsft'), cbar h
	im(5, An'*y - As'*y, 'nufft - strip'), cbar h
	im(6, Ad'*y - As'*y, 'dsft - strip'), cbar h
return
end
