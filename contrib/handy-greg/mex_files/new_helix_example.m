% helix_example.m
% Example of how to use feldkamp.m for cone-beam CT reconstruction
% Copyright 2010-7-21, Gregory Handy and Jeff Fessler, 
% University of Michigan

if ~isvar('cg'), printm 'cg: cone-beam CT geometry'
	% see book chapter section (ask me) on cone-beam CT recon for notation
	down = 4; % down sample a lot to save time
	% default is arc detector; but allow flat for cone_beam_ct_example.m
	if ~isvar('dfs'), dfs = 0; end
	cg = ct_geom('fan', 'ns', 888, 'nt', 32, 'na', 17688, ...
		'ds', 1, 'dt', 2, ...
		'down', down, ...
		'dsd', 949, 'dod', 408, 'dfs', 0, 'pitch', .5, 'source_z0', -100, ...
        'orbit', 5760, 'orbit_start', 109.25);
	printm('fov rmax=%g', cg.rmax)
	clear dfs
end

if ~isvar('ig'), printm 'ig: image geometry'
    ig = image_geom('nx', 256, 'ny', 256, 'nz', 160, ...
	 	'down', down, 'dx', 2, 'dz', .625);
	mask2 = true([ig.nx ig.ny]);
	mask2(end) = 0; % trick: test it
	ig.mask = repmat(mask2, [1 1 ig.nz]);
	clear mask2
     nz = ig.nz;
end

if ~isvar('ell'), printm 'ell: ellipsoid object'
	ell = [ ...
		[20 10 10	150 150 380	0 0 0.01];
		[80 10 10	50 50 30	0 0 0.01];
		[-10 -40 75	40 40 40	0 0 0.01];
		[-10 80 -20	30 30 30	0 0 0.01];
	];
end

if ~isvar('xtrue'), printm 'xtrue: true image volume'
	xtrue = ellipsoid_im(ig, ell);
    prompt
end

if ~isvar('proj'), printm 'proj: analytical ellipsoid projection views'
	proj = ellipsoid_proj(cg, ell);
   li_hat = proj;
end

% Helical Cone-beam reconstruction
[sino z_orbit_start] = helix_rebin(cg, ig, li_hat);

  down = cg.down;
   img = zeros(ig.nx,ig.ny,nz);
   
 ig2 = image_geom('nx', ig.nx*down, 'ny', ig.ny*down, 'dx', ...
        ig.dx/down, 'down', down);
	mask2 = true([ig.nx ig.ny]);
	mask2(end) = 0; % trick: test it
	ig.mask = repmat(mask2, [1 1]);
	clear mask2
    
    printm('The rebinning step is done. Performing the filter back projection.')
    % filter back project for each z-slice
    for iz = 1:nz  
        sinost = sino_geom('fan', 'ns', cg.ns*down, 'na', size(sino,2)*down, ...
                    'ds', cg.ds/down, 'down', down, 'orbit', 'short', ...
                    'orbit_start', double(z_orbit_start(iz)), ...
                    'offset_s', cg.offset_s, ...
                    'dsd', cg.dsd, 'dod', cg.dod, 'dfs', cg.dfs); 
        geomsino = fbp2(sinost, ig2);
        img(:,:,iz) = fbp2(sino(:,:,iz), geomsino);
    end