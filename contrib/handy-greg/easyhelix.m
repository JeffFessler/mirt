  function img = easyhelix(cg, ig, proj, varargin)
%|
%| FBP reconstruction of cone-beam tomography data collected with
%| a circular source trajectory.
%| See feldkamp_example.m for example.
%|
%| in
%|	cg			ct_geom()
%|	ig			image_geom()
%|	proj	[ns nt na]	cone-beam projection views (line integrals)
%|
%|
%| out
%|	img	[nx ny nz]	reconstructed image
%|

% defaults
img = helix_do(proj, cg, ig, cg.na, cg.dt, cg.dsd, cg.dso, cg.orbit, cg.orbit_start);
end % feldkamp()

%
% feldkamp_do()
%
function img = helix_do(proj, cg, ig, na, dt, dsd, dso, orbit, orbit_start)

% step 1: fix z-sampling; for each CB source point, determine multifan
% betas = deg2rad(orbit_start + orbit * [0:na-1] / na); % [na] source angles

% assumes that na is even
num_turns = orbit/360;
phis = mod(deg2rad(orbit_start + orbit * [0:(na/num_turns-1)]/na), 2*pi);
numPhis = size(phis',1);

if ig.dz < .9*(dso/dsd * dt)
    min_spacing = ig.dz;
else
    min_spacing = .9*(dso/dsd * dt);
end

num_zsamp = ig.nz;
% num_zsamp = ceil((ig.z(size(ig.z,1))-ig.z(1))/min_spacing);
% zSamples = ig.z(1) + min_spacing * [0:num_zsamp-1];
zSamples = ig.z;
% -cg.source_zs(1)+cg.source_zs(2)
% myPitch = (-cg.source_zs(1)+cg.source_zs(2))*(na/2-1)
myPitch = cg.pitch * cg.nt * dso / dsd * cg.dt;

delta = asin(cg.rmax/dso);
dist = .5 * myPitch * (pi + 2 * delta)/(2*pi);

%step 2: rebin the cone beam data to fit a fan beam projection
[ns, nt, na] = size(proj);
fanBeam = zeros(ns, numPhis, num_zsamp);
for i=0:na-1
    currentZ = cg.source_zs(i+1);
    upper = currentZ + dist;
    lower = currentZ - dist;

    lcount = 1;

    while lcount < num_zsamp+1 && zSamples(lcount) < lower
        lcount = lcount + 1;
    end

    while lcount < num_zsamp+1 && zSamples(lcount) < upper
        deltaZ = zSamples(lcount)-currentZ;
        shortscan = (pi/2+delta)*(1-deltaZ/dist);

        is = ndgrid(1:cg.ns,1);
        spoints = cg.s;

        % tpoints = ((spoints).^2+dsd*dsd)./ (dso*dsd).* deltaZ;
        tpoints = ones(cg.ns,1).* deltaZ;

        scale = sqrt(spoints.^2+dsd^2)./sqrt(spoints.^2+tpoints.^2+dsd^2);

        t_in = (tpoints./cg.dt + cg.wt);

        x0 = floor(t_in);
        x1 = 1 + x0;

        alpha = t_in - x0;

        for j = 1:size(x0,1)
            weight = myParker(spoints(j),shortscan,delta,dsd);
            fanBeam(is(j), mod(i,numPhis)+1, lcount) = weight*proj(is(j), nt/2, i+1);
        end

        lcount = lcount + 1;
    end

end

% Everything below here goes really fast
disp('we are here')
    orbitStart = zeros(num_zsamp,1);
    % decide on beginning orbit
    orbit = deg2rad(cg.orbit_start);;
    zloc = cg.source_z0;
    change_phis = phis(2)-phis(1);
    for z = 1: num_zsamp
        orbitFound = 0;
        while orbitFound == 0
            if zSamples(z) < zloc + dist
                orbitStart(z) = orbit;
                orbitFound = 1;
            else
                orbit = orbit + change_phis;
                zloc = zloc + (cg.source_zs(2) - cg.source_zs(1));
            end
        end
    end

    newNA = floor(2*dist/myPitch*numPhis);
    newFanBeam = zeros(ns, newNA, num_zsamp);

    for z = 1 : num_zsamp
        angle = round((orbitStart(z)-deg2rad(cg.orbit_start)) / change_phis);
        for o = 0: newNA-1
            nextAngle = angle + o;
            newFanBeam(:, o+1, z) = fanBeam(:, mod(nextAngle, numPhis)+1, z);
        end
    end

    down = cg.down;
    wantedZ=ceil(ig.nz/2);
     im(newFanBeam(:,:,wantedZ), 'test'), cbar

    for z = 1:num_zsamp
        sinost = sino_geom('fan', 'ns', cg.ns*down, 'na', newNA*down, ...
            'ds', cg.ds/down, 'down', down, 'orbit', 'short', ...
            'orbit_start', rad2deg(orbitStart(z)), ...
            'offset_s', cg.offset_s, ...
            'dsd', cg.dsd, 'dod', cg.dod, 'dfs', cg.dfs);
        wt = fbp_fan_short_wt(sinost);

        for i = 1:ns
            for j = 1:newNA
                newFanBeam(i, j, z) = wt(i,j)*newFanBeam(i, j, z);
            end
        end
    end

    im(newFanBeam(:,:,:), 'test'), cbar

   % same image geometry as before, except for the z direction
    ig = image_geom('nx', ig.nx*ig.down, 'ny', ig.ny*ig.down, 'dx', ...
        ig.dx/ig.down, 'down', ig.down);
    mask2 = true([ig.nx ig.ny]);
    mask2(end) = 0; % trick: test it
    ig.mask = repmat(mask2, [1 1]);
    clear mask2

    img = zeros(ig.nx,ig.ny,num_zsamp);

    for i = 1:num_zsamp
        sinost = sino_geom('fan', 'ns', cg.ns*down, 'na', newNA*down, ...
            'ds', cg.ds/(down), 'down', down, 'orbit', 'short', ...
            'orbit_start', rad2deg(orbitStart(i)), ...
            'offset_s', cg.offset_s, ...
            'dsd', cg.dsd, 'dod', cg.dod, 'dfs', cg.dfs);
        geomsino = fbp2(sinost, ig);
        img(:,:,i) = fbp2(newFanBeam(:,:,i), geomsino);
    end

    im(img(:,:,round(num_zsamp/2)+1)), cbar;

end
