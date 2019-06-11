% Gtomo_nufft_3d_par_test.m
% Fourier-based forward projection for 3D parallel beam case
% Yong Long, 2007-07-27

if ~isvar('ig'), printm 'ig'
	down = 8; % for timing
	ig = image_geom('nx', 512, 'ny', 496, 'nz', 32*10, 'fov', 500, ...
		'down', down);
end

if ~isvar('cg'), printm 'cg'
	cg = ct_geom('fan', 'ns', 888, 'nt', ig.nz * ig.down, 'na', 984, ...
		'ds', 1, 'dt', 0.625 * 949/541, ...
		'dsd', inf,...%'dsd', 949, 'dso', 541, ...%
		'offset_s', 1.25, 'offset_t', 0, ...
		'down', 8, 'dfs', inf);
end

if ~isvar('xtrue'), printm 'xtrue'
	ell = [ ...
%		[0 0 0  200 150 100  0 0 10];
%		[60 1*40 2  50 40 6  0 0 10];
        [60 1*40 20  50 40 60  0 0 10];
%		[0 0 0  50 50 6  0 0 10];
%		[0 80 2.5  20 20 5  0 0 10];
		%[100 0 0 50 50 50 0 0 10];
		%[0 0 80 40 40 40 0 0 10];
		%[0 70 0 30 30 30 0 0 10];
		];
    f.over = 2;
	xtrue = ellipsoid_im(ig, ell, 'oversample', f.over);

    im plc 2 2
    t = sprintf('x true, z=%g to %g', ig.z(1), ig.z(end));
    im(1, ig.x, ig.y, xtrue, t), cbar
end

if ~isvar('proj'), printm ' analytical proj'
    if down == 1, f.over = 1; else, f.over = 2; end
    cpu tic
    proj = ellipsoid_proj(cg, ell, 'oversample', f.over);
    cpu toc 'analytical proj time:'
    
    tosino = inline('permute(x, [1 3 2])', 'x'); % sinogram display
    im(2, tosino(proj), 'true projections'), cbar
end

%
% create Gtomo_nufft class object
%
if ~isvar('Gn'), printm 'setup NUFFT_3d_test'
	cpu tic
	arg = {
		'dis_src_det', cg.dsd,...
        'dis_iso_det', cg.dod,...
        'offset_s', cg.offset_s,...
        'offset_t', cg.offset_t,...
		'orbit_a', cg.orbit, ...%azimuthal angle coverage 
		'orbit_a_start', cg.orbit_start, ...%first azimuthal angle [0 degrees]
        'orbit_p', 30, ...     %polar angle coverage [-orbit_p orbit_p]
		'dx', ig.dx, ...
        'dz', ig.dz, ...
		'ds', cg.ds, ...
        'dt', cg.dt, ...
		'rect_s', cg.ds,... %corresponding to 1D strip_width
        'rect_t', cg.dt,... %corresponding to 1D strip_width
		'interp', {'table', 2^11, 'minmax:kb'}, ...
		'is.complex', 0, ...
		'yscale', -ig.dy/ig.dx, ...
	};
    
   % npo = cg.nt; %number of polar angles
    npo = 1; %now only one polar angle(theta=0) is considered for test 
    %fuctions ct_geom and ellipsoid_proj set all polar angles equal to 0. 
    %Add parameters here to make Gtomo_nufft_3d_par handle 3D parallel beam 
    %case, not 2.5D.
    
    Gn = Gtomo_nufft_3d_par(ig.mask, [cg.ns cg.nt cg.na npo], arg{:});

	cpu toc 'Gn pre time:'
end

if 1  
    x = xtrue(ig.mask(:));
    cpu tic
    %profile on
	t = Gn * x;
    %profile viewer
    cpu toc
    %compare with analytical result(theta=0)
    im(3, tosino(t), 't'), cbar
    max_percent_diff(t, proj)
    100*nrms(t(:), proj(:))
end




