function img = cbct_back_mat_par(interpRect,rectGrid,proj, ns, nt, na, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	source_zs, ...
	mask, nz, dx, dy, dz, offset_xyz, ia_skip, scale_dang, extrapolate_t)

if any(source_zs ~= 0)
	warn('helix not yet tested')
end

[nx ny] = size(mask);
betas = deg2rad(orbit_start + orbit * [0:na-1] / na); % [na] source angles

% precompute as much as possible
wx = (nx-1)/2;
wy = (ny-1)/2;
wz = (nz-1)/2;
[xc yc] = ndgrid(([0:nx-1] - wx) * dx, ([0:ny-1] - wy) * dy);
zc = ([0:nz-1] - wz) * dz;

clear wx wy wz rr smax rmax

xc = xc(mask); % [np] pixels within mask
yc = yc(mask);

ws = (ns+1)/2; % trick: +1 because matlab starts from 1
wt = (nt+1)/2;

% loop over slices
img = zeros([size(mask) nz]);
sdim = [ns+3 nt+3]; % trick: extra zeros saves indexing in loop
proj1 = zeros(sdim);
ticker reset
for iz=1:nz

	ia_min = 1;
	ia_max = na;

	% loop over each projection angle
	img2 = 0;
	for ia=ia_min:ia_skip:ia_max
		ticker(mfilename, [iz ia], [nz na])
		beta = betas(ia);

		x_beta = +xc * cos(beta) + yc * sin(beta);
		y_betas = (-xc * sin(beta) + yc * cos(beta));
        
        if interpRect
            mag = sqrt(dso^2-(x_beta).^2) ./ (sqrt(dso^2-(x_beta).^2)+y_betas);     
        else
            mag = dsd ./ (sqrt(dso^2-(x_beta).^2)+y_betas);
        end

        sprime = x_beta;
		
        tprime = mag .* (zc(iz)-source_zs(ia));
        
		bs = sprime / ds + ws;
           
        if interpRect
            newdt = rectGrid(2)-rectGrid(1);
            bt = tprime / newdt + wt;
        else
            bt = tprime / dt + wt;
        end
        
       
		% bi-linear interpolation:
		is = floor(bs); % left bin
		it = floor(bt);
       
		wr = bs - is;	% left weight
		wl = 1 - wr;	% right weight
		wu = bt - it;	% upper weight
		wd = 1 - wu;	% lower weight

		ibad = (is < 0) | (is > ns) | (it < 0) | (it > nt);
		is(ibad) = ns+1; % trick! point at harmless zeros
		it(ibad) = nt+1;
       
		proj1(1+[1:ns],1+[1:nt]) = proj(:,:,ia); % trick: left side
        
		p1 =	wl .* proj1(sub2ind(sdim, is+1,it+1)) + ...
			wr .* proj1(sub2ind(sdim, is+2,it+1));
		p2 =	wl .* proj1(sub2ind(sdim, is+1,it+2)) + ...
			wr .* proj1(sub2ind(sdim, is+2,it+2));

		p0 = wu .* p1 + wd .* p2; % vertical interpolation

        %no backprojection weighting needed
		img2 = img2 + p0;
	end % ia

	img(:,:,iz) = embed(img2, mask);
end % iz

if scale_dang % final "\der angle" scale:
	img = (0.5 * deg2rad(abs(orbit)) / (na/ia_skip)) * img;
end

end % cbct_back_mat_par()

