% ir_fdk_helix_example.m
% Example of how to use feldkamp_helix.m for cone-beam axial and helical CT reconstruction
% Copyright 2014-11-11, Donghwan Kim and Jeff Fessler, University of Michigan

f.helical = 1;
f.down = 2;
f.pht_ell = 0; 
f.pht_cyl = 0;
f.xcat = 1;
if f.xcat
	f.down = 1;
end
f.fdk_type = 3; % FDK type: (1: SSRB, 2: Native, 3: Cone-parallel)
f.use_mex = 1;
if f.helical
	f.w3d = 10;  % not fully tuned
else
	f.w3d = 1; %1; % not fully tuned % see tang:05:atd
end
%f.w3d = 0;


if ~isvar('cgARC'), printm 'cg: cone-beam CT geometry'
	if f.helical && f.xcat
		f.nturn = 3; 
		f.pitch = 63/64;
		f.source_z0 = -60; 
		cgARC = ct_geom('ge2', 'pitch', f.pitch, ...
			'source_z0',f.source_z0, ...
			'na', 984*f.nturn, 'orbit', 360*f.nturn, ...
			'down', f.down);
	elseif f.helical
		f.nturn = 12;
		f.pitch = 0.5;
		f.source_z0 = -100; 
		cgARC = ct_geom('ge2', 'pitch', f.pitch, ...
			'source_z0',f.source_z0, ...
			'na', 984*f.nturn, 'orbit', 360*f.nturn, ...
			'down', f.down);
	elseif f.xcat
		cgARC = ct_geom('ge2', 'nt', 128, 'down', f.down);
	else
		cgARC = ct_geom('ge2', 'down', f.down);
	end
end

if ~isvar('ig'), printm 'ig: image geometry'
	ig = image_geom('nx', 512, 'ny', 512, 'nz', 154, ...
			'dx', 0.9766, 'dz', 0.625, 'down', f.down);
end

if ~isvar('ell') && f.pht_ell, printm 'ell: ellipsoid object'
	z0 = ig.offset_z * ig.dz;
	ell = [...
		[0 0 -z0        [0.3 0.1]*ig.fov 0.6*ig.zfov    0 0 1000];
		[80 10 -z0      30 30 10        0 0 1000];
		[-10 -40 75     40 40 40        0 0 1000];
		[-10 80 -20     30 30 30        0 0 1000];
		[0 0 -36        20 20 5 0 0 1000];
		[0 0 -12        20 20 5 0 0 1000];
		[0 0 12 20 20 5 0 0 1000];
		[0 0 36 20 20 5 0 0 1000];
	];
	clear z0
end

if ~isvar('cyl') && f.pht_cyl, printm 'cyl: cylinder object'
	cyl = [...
		[10 -30 0	[0.3 0.3]*ig.fov inf	0 1000];
		[40 30 10	40 40 inf		0 1000];
	];
end 

if ~isvar('xtrue'), printm 'xtrue: true image volume'
	if f.xcat
		xtrue2 = fld_read('XCAT.fld');
		xtrue = zeros(size(xtrue2,1)/2, size(xtrue2,2)/2, size(xtrue2,3));
		for iz=1:size(xtrue,3)
			xtrue(:,:,iz) = downsample2(xtrue2(:,:,iz), 2);
		end
	elseif f.pht_ell && f.pht_cyl
		xtrue = ellipsoid_im(ig, ell, 'oversample', 2) + ...
			cylinder_im(ig, cyl, 'oversample', 2);
	elseif f.pht_ell
		xtrue = ellipsoid_im(ig, ell, 'oversample', 2);
	else
		xtrue = cylinder_im(ig, cyl, 'oversample', 2);
	end

	clim = [800 1200]; 
	if im
		figure(1), im plc 2 3
		im(1, ig.x, ig.y, xtrue, clim), cbar
		titlef('x true, z=%g to %g', ig.z(1), ig.z(end))
	end
prompt
end

if ~isvar('projARC'), printm 'proj: analytical ellipsoid projection views'
	if f.xcat
		nz_add = 100; 
		ig2 = image_geom('nx', ig.nx*2, 'ny', ig.ny*2, 'nz', ig.nz + nz_add*2, ...
			'dx', ig.dx/2, 'dz', ig.dz);
		xtmp2 = zeros(ig2.nx, ig2.ny, ig2.nz);
		xtmp2(:,:,[1:nz_add]) = repmat(xtrue2(:,:,1), [1 1 nz_add]); 
		xtmp2(:,:,[1:nz_add] + ig.nz + nz_add) = repmat(xtrue2(:,:,end), [1 1 nz_add]);
		xtmp2(:,:,[1:ig.nz] + nz_add) = xtrue2;

		A = Gcone(cgARC, ig2, 'type', 'sf1', 'use_hct2', 0);
		projARC = A * xtmp2;
	elseif f.pht_ell && f.pht_cyl
		projARC = ellipsoid_proj(cgARC, ell) + cylinder_proj(cgARC, cyl);
	elseif f.pht_ell
		projARC = ellipsoid_proj(cgARC, ell);
	else
		projARC = cylinder_proj(cgARC, cyl);
	end
	figure(1), im(4, 'row', 1, permute(projARC, [1 3 2]), 'true helix sinograms'), cbar
prompt
end


if ~isvar('xfdk'), printm 'fdk'
	f.extrapolate_t = ceil(1.3 * cgARC.nt/2); 
	if f.use_mex == 0
		f.extrapolate_t = 1; % default
	end

	switch f.fdk_type
	case 1 % Single slice rebinning 
		f.short = 1; %1;
		f.itype = 'linear';
		[sino orbits used] = rebin_helix(cgARC, ig, projARC, ...
			'type', f.itype, 'short', f.short, 'collapse', 0);
		[xfdk tmp] = fbp_helix_stack(cgARC, ig, sino, orbits, ...
			'short', f.short, 'window', 'hanning,1.0');
	
	case 2 % Native cone-beam geometry
		xfdk = feldkamp_helix(cgARC, ig, projARC, ...
			'window', 'hann', ...
			'use_mex', f.use_mex, ...
			'extrapolate_t', f.extrapolate_t);

	case 3 % Cone-parallel geometry
		if ~isvar('cgPAR'), printm 'cg: cone-beam CT cone-parallel geometry'
			f.dsChange = cgARC.dsd / cgARC.dso;
			f.nr = 4 * ceil(cgARC.ns * f.dsChange / 4);
			f.dsChange = f.nr / cgARC.ns;
			f.dr = cgARC.ds / f.dsChange;

			% should be the same to cgARC except dfs ns ds
			if f.helical
				cgPAR = ct_geom('ge2', 'pitch', f.pitch, ...
					'source_z0', f.source_z0, ...
					'na', 984*f.nturn, 'orbit', 360*f.nturn, ...
					'ns', f.nr*f.down, ...
					'ds', f.dr/f.down, ...
					'down', f.down', 'dfs', -inf);
			else
				cgPAR = ct_geom('ge2', 'ns', f.nr*f.down, ...
					'ds', f.dr/f.down, ...
					'down', f.down, 'dfs', -inf);
			end
		end

		if ~isvar('projPAR'), printm 'sino: cone-parallel rebinning'
			% fan beam sinogram
			sf = sino_geom('fan', 'ns', cgARC.ns, 'na', cgARC.na, ...
				'ds', cgARC.ds, ...
				'offset_s', cgARC.offset_s, ... % quarter detector
				'strip_width', cgARC.ds, ...
				'dsd', cgARC.dsd, 'dod', cgARC.dod, 'dfs', cgARC.dfs, ...
				'orbit', cgARC.orbit);
			% parallel beam sinogram
			sp = sino_geom('par', 'nb', cgPAR.ns, 'na', cgPAR.na, ...
				'dr', cgPAR.ds, ...
				'offset_r', cgPAR.offset_s, ... % quarter detector
				'strip_width', cgPAR.ds, ...
				'orbit', cgPAR.orbit);
			% rebinning 
			tmp = permute(projARC, [1 3 2]); % sta -> sat
			fansproj = rebin_fan2par_helix(tmp, sf, sp, 'helical', f.helical); % DK
			projPAR = permute(fansproj, [1 3 2]); % sat -> sta
			projPAR = max(projPAR, 0); % DK

			figure(1); 
			im(5, 'row', 1, permute(projPAR, [1 3 2]), 'rebinned helix sinograms'), cbar
		end

		xfdk = feldkamp_helix(cgPAR, ig, projPAR, ...
			'window', 'hann', ...
			'use_mex', f.use_mex, ...
			'extrapolate_t', f.extrapolate_t, ...
			'cone_par', 1, ...
			'w3d', f.w3d);
	otherwise
		fail 'bug'
	end

	xfdk = xfdk .* (ig.circ(cgARC.rmax) > 0); % mask
	figure(1), im(2, xfdk, 'FDK recon', clim), cbar;
	figure(1), im(3, xfdk - xtrue, 'error', [-200 200]), cbar; 

	figure(2), im('mid3', xfdk, 'FDK recon', clim), cbar;
	figure(3), im('mid3', xfdk - xtrue, 'error', [-200 200]), cbar;
prompt
end

if im % profile, for checking HU accuracy
	figure(1), im subplot 6
	iz = 1:ig.nz;
	ix = ig.nx/2+1; iy = ig.ny/2+1;
	plot(   ig.z, squeeze(xtrue(ix,iy,iz)), 'o-', ...
		ig.z, squeeze(xfdk(ix,iy,iz)), '.--')
	axis([minmax(ig.z)' -200 2100])
	titlef('profile at (ix,iy)=(%g,%g)', ix,iy), xlabel z
end
