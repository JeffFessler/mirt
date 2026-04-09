  function back = cbct_back(proj, cg, ig, varargin)
% function back = cbct_back(proj, cg, ig, varargin)
%|
%| cone-beam backprojector for feldkamp.m
%| in
%|	proj	[ns nt na]	cone-beam projection views
%| option
%|	'use_mex' 1|2|3		default: 1 mex with loop in mex
%|					2 mex with loop in matlab
%|					3 mex st
%|					0 no mex; use matlab
%|	'ia_skip'		default: 1
%|	'offset_source'		default: 0
%|	'nthread'		default: jf('ncore')
%| out
%|	back	[nx ny nz]	back projection result
%|
%| Copyright 2004-8-28 Jeff Fessler, University of Michigan

if nargin == 1 && streq(proj, 'test'), cbct_back_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

arg.ia_skip = 1;
arg.use_mex = 1;
arg.offset_source = 0;
arg.nthread = jf('ncore');
arg.back_call = @fdk_mex_call;
arg = vararg_pair(arg, varargin);

back = cbct_back_do(proj, cg.ns, cg.nt, cg.na, ...
	cg.ds, cg.dt, cg.offset_s, cg.offset_t, arg.offset_source, ...
	cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ...
	ig.mask_or, ig.nz, ig.dx, ig.dy, ig.dz, ...
	[ig.offset_x ig.offset_y ig.offset_z], ...
	arg.ia_skip, arg.use_mex, arg.nthread, arg.back_call);


function img = cbct_back_do(proj, ns, nt, na, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	mask, nz, dx, dy, dz, offset_xyz, ia_skip, use_mex, nthread, back_call)

[nx ny] = size(mask);

betas = deg2rad(orbit_start + orbit * [0:na-1] / na); % [na] source angles

% mex backprojector
if arg.use_mex
	proj = single(proj);
	nthread = int32(arg.nthread);

	if arg.use_mex == 1 % loop in mex
		proj = permute(proj, [2 1 3]); % ts
		tmp = back_call('fdk,ts,back', ...
				int32([nx ny nz]), [dx dy dz], ...
				offset_xyz, uint8(mask), ...
				dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
				proj, betas, nthread);
		img = double6(tmp);
		img = permute(img, [2 3 1]); % zxy -> xyz

	elseif arg.use_mex == 2 % loop in matlab (for testing)
		proj = permute(proj, [2 1 3]); % ts
		img = 0;
		for ia=1:na
			ticker(mfilename, ia, na)
			% note: 2006-5-30: replaced -dy with dy
			tmp = back_call('fdk,ts,back', ...
				int32([nx ny nz]), [dx dy dz], ...
				offset_xyz, uint8(mask), ...
				dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
				proj(:,:,ia), betas(ia), nthread);
			tmp = double6(tmp);
			img = img + tmp;
		end
		img = permute(img, [2 3 1]); % zxy -> xyz

	elseif arg.use_mex == 3 % fdk,st (for testing only - slower!)
		tmp = back_call('fdk,st,back', ...
				int32([nx ny nz]), [dx dy dz], ...
				offset_xyz, uint8(mask), ...
				dso, dsd, dfs, [ds dt], [offset_s offset_t], ...
				proj, betas, nthread);
		img = double6(tmp);

	else
		fail 'bug'
	end

	% final "\der angle" scale:
	img = (0.5 * deg2rad(abs(orbit)) / (na/arg.ia_skip)) * img;
return
end

% matlab backprojector (slower)

% precompute as much as possible
wx = (nx-1)/2 + offset_xyz(1);
wy = (ny-1)/2 + offset_xyz(2);
wz = (nz-1)/2 + offset_xyz(3);
[xc yc] = ndgrid(([0:nx-1] - wx) * dx, ([0:ny-1] - wy) * dy);
zc = ([0:nz-1] - wz) * dz;

if 0 % limit back-projection to FOV?  removed 2008-10-9
	rr = sqrt(xc.^2 + yc.^2); % [nx,ny]
	smax = ((ns-1)/2-abs(offset_s)) * ds; % maximum detector s coordinate

	if isinf(dfs)
		gamma_max = atan(smax/dsd);
	elseif dfs == 0
		gamma_max = smax / dsd;
	end

	rmax = dso * sin(gamma_max);
	mask = mask & (rr < rmax);
end
clear wx wy wz rr smax rmax

xc = xc(mask); % [np] pixels within mask
yc = yc(mask);

% loop over slices
img = zeros([size(mask) nz]);
ticker reset
for iz=1:nz

	% loop over each projection angle
	img2 = 0;
	for ia=1:na
		ticker(mfilename, [iz ia], [nz na])
		beta = betas(ia);

		x_beta = +xc * cos(beta) + yc * sin(beta);
		y_beta = -xc * sin(beta) + yc * cos(beta);

		% detector indices
		mag = dsd ./ (dso - y_beta);

		if isinf(dfs)
			sprime = mag .* x_beta;
		elseif (dfs == 0)
			d_loop = dso - y_beta;
			r_loop = x_beta - offset_source;
			sprime = dsd * atan2(r_loop, d_loop);
		end

		tprime = mag * zc(iz);
		bh = sprime / ds + (ns+1)/2 + offset_s;
		bv = tprime / dt + (nt+1)/2 + offset_t;

		% bi-linear interpolation:
		ih = floor(bh); % left bin
		iv = floor(bv);
		igood = 1<=ih & ih<ns & 1<= iv & iv<nt;

		ih(~igood) = ns+1; % trick! point at harmless zeros
		iv(~igood) = 1;

		wr = bh - ih;	% left weight
		wl = 1 - wr;	% right weight
		wu = bv - iv;	% upper weight
		wd = 1 - wu;	% lower weight

		sdim = size(proj);
		ia1 = ia * ones(size(ih));

		p1 =	wl .* proj(sub2ind(sdim, ih,iv,ia1)) + ...
			wr .* proj(sub2ind(sdim, ih+1,iv,ia1));
		p2 =	wl .* proj(sub2ind(sdim, ih,iv+1,ia1)) + ...
			wr .* proj(sub2ind(sdim, ih+1,iv+1,ia1));

		p0 = wu .* p1 + wd .* p2; % vertical interpolation

		if isinf(dfs)
			p0 = p0 .* mag.^2; % back-projection weighting for flat
		elseif dfs == 0
			p0 = p0 .* (dsd.^2) ./ (r_loop.^2 + d_loop.^2);
		end

		img2 = img2 + p0;
	end % ia

	img(:,:,iz) = embed(img2, mask);
end % iz

img = (0.5 * deg2rad(orbit) / (na/arg.ia_skip)) * img; % final "\der angle" scale
end % cbct_back_do()

end % cbct_back()

%
% fdk_mex_call()
%
function out = fdk_mex_call(varargin)
if exist('fdk_mex') == 3
	printm 'using fdk_mex'
	out = fdk_mex(varargin{:});
elseif exist('jf_mex') == 3
	out = jf_mex(varargin{:});
else
	fail('bug: neither fdk_mex nor jf_mex found')
end
end % fdk_mex_call()


function cbct_back_test

down = 4;
cg = ct_geom('ge1', 'nt', 800, 'down', down);
cg.dt = -cg.dt;
ig = image_geom('nx', 512, 'ny', 496, 'nz', 480, 'fov', 500, ...
	'mask', 'all-but-edge', ...
	'down', down);

pr [ig.nx ig.ny ig.nz]
pr [cg.ns cg.nt cg.na]
printm('image Mbytes %d', ig.nx*ig.ny*ig.nz*4 / 1024^2)
printm('proj Mbytes %d', cg.ns*cg.nt*cg.na*4 / 1024^2)

proj = cg.zeros;
if 1 % sprinkle of nonzero values
	for ia=1:cg.na
		proj(mod(ia,cg.ns)+1, mod(ia,cg.nt)+1, ia) = mod(ia,3)+1;
	end
else
	ia = round(cg.na/3);
	proj(:,:,ia) = 1; % one view nonzero
	%proj(round(end/2), 3, ia) = 1;
	proj = cg.ones;
end

use_mex = 1;

if exist('fdk_mex') == 3, printm 'testing cuda version'

	blist = {@jf_mex, @fdk_mex};
	for ii=1:length(blist)
		cpu etic
		back{ii} = cbct_back(proj, cg, ig, 'use_mex', use_mex, ...
			'back_call', blist{ii});
%			'nthread', 1, ...
		times(ii) = cpu('etoc');
		printm('time %s = %g', func2str(blist{ii}), times(ii))
	end
	im(back{2})
	max_percent_diff(back{1}, back{2})
%	equivs(back{1}, back{2}, 'thresh', 1e-5)
else
	cpu etic
	back = cbct_back(proj, cg, ig, 'use_mex', use_mex);
	cpu etoc 'back time: '
	im clf, im(back)
end

end % cbct_back_test()
