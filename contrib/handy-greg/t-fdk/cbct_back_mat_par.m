% 2012-fall by Greg Handy, based on cbct_back.m
% 2013-04-04 refined by Rebecca Malinas
% 2013-04-07 refined by Jeff Fessler

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

ws = (ns+1)/2 + offset_s; % trick: +1 because matlab starts from 1
wt = (nt+1)/2 + offset_t; % offset_t should be zero for cone-par and tent geometry

% loop over slices
img = zeros([size(mask) nz]);
sdim = [ns+1 nt]; % trick: extra zeros saves indexing in loop
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
			mag = sqrt(dso^2-(x_beta).^2) ./ (sqrt(dso^2-(x_beta).^2)+y_betas); %T-FDK?
		else
			mag = dsd ./ (sqrt(dso^2-(x_beta).^2)+y_betas); %P-FDK?
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


		bt = max(bt, 1);
		bt = min(bt, nt);

		% bi-linear interpolation:
		is = floor(bs); % left bin
		it = floor(bt);

		itbad = (it==nt);
		it(itbad)=nt-1;

		isbad = (bs < 1) | (bs > ns);
		is(isbad) = ns;

		wr = bs - is;	% left weight
		wr(isbad)=1;

		wl = 1 - wr;	% right weight
		wu = bt - it;	% upper weight
		wd = 1 - wu;	% lower weight

%		ibad = (is < 0) | (is > ns) | (it < 0) | (it > nt);
%		is(ibad) = ns+1; % trick! point at harmless zeros
%		it(ibad) = nt+1;

 %		proj1(1+[1:ns],1+[1:nt]) = proj(:,:,ia); % trick: left side


%		p1 =	wl .* proj1(sub2ind(sdim, is+1,it+1)) + ...
%			wr .* proj1(sub2ind(sdim, is+2,it+1));
%		p2 =	wl .* proj1(sub2ind(sdim, is+1,it+2)) + ...
%			wr .* proj1(sub2ind(sdim, is+2,it+2));

									proj1(1:ns,1:nt) = proj(:,:,ia);%now only have extra row of zeros at ns+1 % trick: left side

		p1 =	wl .* proj1(sub2ind(sdim, is, it)) + ...
			wr .* proj1(sub2ind(sdim, is+1,it));

		p2 =	wl .* proj1(sub2ind(sdim, is, it+1)) + ...
			wr .* proj1(sub2ind(sdim, is+1,it+1));
		p0 = wd .* p1 + wu .* p2; % vertical interpolation

		% no backprojection weighting needed
		img2 = img2 + p0;
	end % ia

	img(:,:,iz) = embed(img2, mask);
end % iz

if scale_dang % final "\der angle" scale:
	img = (0.5 * deg2rad(abs(orbit)) / (na/ia_skip)) * img;
end

end % cbct_back_mat_par()

