cg = ct_geom('fan', 'ns', 888, 'nt', 32, 'na', 11526, ...
		'ds', 1.0239, 'dt', 1.096436, ...
		'dsd', 949.075012, 'dso', 541, 'dfs', 0, 'pitch', .53125, ...
        'source_z0', -103.032, 'offset_s', 1.25, ...
        'orbit', 4216.83, 'orbit_start', 1309.33);
	printm('fov rmax=%g', cg.rmax)
    
    ig = image_geom('nx', 924, 'ny', 924, 'nz', 176, 'dx', .390625, ...
	'dz', .6250, 'offset_z', 65.3);
	mask2 = true([ig.nx ig.ny]);
	mask2(end) = 0; % trick: test it
	ig.mask = repmat(mask2, [1 1 ig.nz]);
	clear mask2

sino = fld_read('/net/ir51/y/fessler/handy/2009-11-11-sf1-v03/yi-tsa.fld');
% sino = flipdim(sino,1);
sino = permute(sino, [2 1 3]); % make it [ns nt na]
sino = sino * (cg.dsd/cg.dod);

[bigReal, sino] = helix(cg, ig, sino, 1);

% for i = 1:size(ig.x)
%     for j = 1:size(ig.y)
%         if sqrt(ig.y(j)^2+ig.x(i)^2) > 250
%             bigReal(i,j,:) = 0;
%         end
%     end
% end