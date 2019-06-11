 function img = helix(cg, ig, proj, short, varargin)
%|function img = helix(cg, ig, proj, varargin)
%|
%| A single slice rebinning method for cone-beam tomography data 
%| collected with a helical source trajectory 
%|
%| in
%|	cg			ct_geom()
%|	ig			image_geom()
%|	proj	[ns nt na]	cone-beam projection views (line integrals)
%|  short    1 for short-scan fanbeam; 0 for 360 fanbeam 
%|
%|
%| out
%|	img	[nx ny nz]	reconstructed image
%|
%| See helix_example, for information on how to call
%|
%| Equations used are taken from 
%| Noo F, Defrise M, Clackdoyle R and Kudo H 1999 Single-slice 
%| rebinning for helical cone-beam CT Phys. Med. Biol. 44 561?70
%|
%| Copyright 2010-7-21, Gregory Handy and Jeff Fessler, 
%| University of Michigan

img = helix_do(proj, cg.down, ig.z, ig.dz, ig.dx, ig.nx, ig.ny, ig.nz, ...
                         cg.na, cg.ns, cg.ds, cg.s, cg.t, cg.nt, cg.dt, cg.dsd, ...
                         cg.dso, cg.dod, cg.dfs, cg.offset_s, cg.wt,...
                         cg.orbit, cg.orbit_start, cg.pitch, cg.rmax, ...
                         cg.source_zs, short);
end % helix()

%
% helix_do()
%
function img = helix_do(proj, down, z, dz, dx, nx, ny, nz, na, ...
                         ns, ds, spoints, tpoints, nt, dt, dsd, dso, ... 
                         dod, dfs, offset_s, wt, orbit, orbit_start, ... 
                         pitch, rmax, source_zs, short)

% step 1: fix z-sampling
num_turns = orbit/360;
betas = mod(deg2rad(orbit_start + orbit*[0:(na/num_turns-1)]/na),2*pi);

num_betas = size(betas',1);

% a check to see if the data is used to its "full use"
if dz > (dso/dsd * dt)*1.1
    printm('Warning: Full use of CB data is not achieved')
end

% z locations of the z-slices
z_slices = z;

myPitch = pitch * nt * dso / dsd * dt;
delta = asin(rmax/dso);

% dist = d in the paper, and allows for either a short scan,
% or for a full 360 scan for each z-slice
if short == 1
    dist = .5 * myPitch * (pi + 2 * delta)/(2*pi); 
else
    dist = .5 * myPitch;
    
    max_t = (spoints(1)+dsd^2) / (dso*dsd) * dist;
    
    if max_t >= -tpoints(1)
       error('CT geometry does not allow for 360 degree rebinning'); 
    end
end

%step 2: rebin the cone beam data to fit a fan beam projection
fan_beam_proj = zeros(ns, num_betas, nz);

% loop over the different view angles
printm('Rebinning step beginning')
is = ndgrid(1:ns,1);
for ia=0:na-1
    % currentZ = source_zs(ia+1);
    currentZ = source_zs(1) + ia * myPitch/na*num_turns;
    
    upper_limit = currentZ + dist;
    lower_limit = currentZ - dist;
    
    zindex = 1;
   
    % enter into the acceptable range for the z-slices
    while zindex <= nz && z_slices(zindex) < lower_limit 
        zindex = zindex + 1;
    end
    
    % loop over the acceptable z-slices for the current view angle
    while zindex <= nz && z_slices(zindex) < upper_limit
        deltaZ = z_slices(zindex)-currentZ;
        
        % Calculate the values of t to be used for each s
        tpoints = (((spoints).^2+dsd.^2)./ (dso.*dsd)).* deltaZ;   
        t_index = (tpoints./dt + wt); 
        x0 = floor(t_index);
        x1 = 1 + x0;
        
        % Scaling factor due to different ray lenghts
        scale = sqrt(spoints.^2+dsd.^2)./sqrt(spoints.^2+tpoints.^2+dsd.^2);
     
        alpha = t_index - x0;
        
        for i = 1:ns
            y1 = proj(is(i), x1(i), ia+1);
            y0 = proj(is(i), x0(i),ia+1);
            fan_beam_proj(is(i), mod(ia,num_betas)+1, zindex) = scale(i).*...
                (alpha(i).*y1+(1-alpha(i)).*y0);
        end
       
        zindex = zindex + 1;
    end
  
end

if short == 1
    printm('The rebinning step is done. Preparing sinogram for FB2.')
    
    fb_orbit_start = zeros(nz,1);
    % decide on beginning orbit for each z-slice
    new_orbit = deg2rad(orbit_start);
    zloc = source_zs(1);
    change_betas = betas(2)-betas(1);
    for z = 1: nz
        orbitFound = 0;
        while orbitFound == 0
            if z_slices(z) < zloc + dist
                fb_orbit_start(z) = new_orbit;
                orbitFound = 1;
            else
                new_orbit = new_orbit + change_betas;
                zloc = zloc + (source_zs(2) - source_zs(1));
            end
        end
    end

    % calculate the number of angles for each z-slice has 
    newNA = floor(2*dist/myPitch*num_betas);
    new_fan_beam_proj = zeros(ns, newNA, nz);

    % deletes the empty rows of the sinograms
    for iz = 1 : nz
        angle = round((fb_orbit_start(iz)-deg2rad(orbit_start))...
                / change_betas);
        for ia = 0: newNA-1
            nextAngle = angle + ia;
            new_fan_beam_proj(:, ia+1, iz) = fan_beam_proj(:, mod(nextAngle, num_betas)+1, iz);
        end
    end
    
    % applies the Parker Weight to each sinogram
    for iz = 1:nz
         sinost = sino_geom('fan', 'ns', ns*down, 'na', newNA*down, ...
            'ds', ds/down, 'down', down, 'orbit', 'short', ...
            'orbit_start', rad2deg(fb_orbit_start(iz)), ...
            'offset_s', offset_s, ...
            'dsd', dsd, 'dod', dod, 'dfs', dfs);
         [wt, scale180] = fbp_fan_short_wt(sinost);

         for is = 1:ns
            for ia = 1:newNA   
                new_fan_beam_proj(is, ia, iz) = wt(is,ia)*new_fan_beam_proj(is, ia, iz);
            end
        end
    end
    
    % same image geometry as before, except for the z direction
    ig = image_geom('nx', nx*down, 'ny', ny*down, 'dx', ...
        dx/down, 'down', down);
    mask2 = true([nx ny]);
    mask2(end) = 0; % trick: test it
    ig.mask = repmat(mask2, [1 1]);
    clear mask2

    img = zeros(nx,ny,nz);
    
    printm('Performing the filter back projection.')
    % filter back project for each z-slice
    for iz = 1:nz  
        sinost = sino_geom('fan', 'ns', ns*down, 'na', newNA*down, ...
            'ds', ds/(down), 'down', down, 'orbit', 'short', ...
            'orbit_start', rad2deg(fb_orbit_start(iz)), ...
            'offset_s', offset_s, ...
            'dsd', dsd, 'dod', dod, 'dfs', dfs);
   
        geomsino = fbp2(sinost, ig);
        img(:,:,iz) = fbp2(scale180*new_fan_beam_proj(:,:,iz), geomsino);
    end
else
    newNA = floor(na/num_turns);
    
    ig = image_geom('nx', nx*down, 'ny', ny*down, 'dx', ...
        dx/down, 'down', down);
	mask2 = true([nx ny]);
	mask2(end) = 0; % trick: test it
	ig.mask = repmat(mask2, [1 1]);
	clear mask2

    img = zeros(nx,ny,nz);
    
    printm('The rebinning step is done. Performing the filter back projection.')
    % filter back project for each z-slice
    for iz = 1:nz  
        sinost = sino_geom('fan', 'ns', ns*down, 'na', newNA*down, ...
                    'ds', ds/down, 'down', down, 'orbit', 360, ...
                    'orbit_start', orbit_start, ...
                    'offset_s', offset_s, ...
                    'dsd', dsd, 'dod', dod, 'dfs', dfs); 
        geomsino = fbp2(sinost, ig);
        img(:,:,iz) = fbp2(fan_beam_proj(:,:,iz), geomsino);
    end
end
end

