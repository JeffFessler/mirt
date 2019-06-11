cg = ct_geom('fan', 'ns', 888, 'nt', 32, 'na', 3625, ...
		'ds', 1.0239, 'dt', 1.096436, ...
		'dsd', 949.075012, 'dso', 541, 'dfs', 0, 'pitch', .53125, ...
        'source_z0', -50.5494160, 'offset_s', 1.25, ...
        'orbit', 1326.21951219512195, 'orbit_start', 109.12139);
	printm('fov rmax=%g', cg.rmax)
    
    ig = image_geom('nx', 320, 'ny', 320, 'nz', 61, 'dx', 2.191162, ...
	'dz', .6250, 'offset_z', 56.0);
	mask2 = true([ig.nx ig.ny]);
	mask2(end) = 0; % trick: test it
	ig.mask = repmat(mask2, [1 1 ig.nz]);
	clear mask2

sino = fld_read('/net/ir51/y/fessler/handy/helix-2795/yi-tsa.fld');
sino = flipdim(sino,1);
sino = permute(sino, [2 1 3]); % make it [ns nt na]
sino = sino * (cg.dsd/cg.dod);

[xfdk, sino2] = helix(cg, ig, sino, 0);

for i = 1:size(ig.x)
    for j = 1:size(ig.y)
        if sqrt(ig.y(j)^2+ig.x(i)^2) > 250
            xfdk(i,j,:) = 0;
        end
    end
end
    
