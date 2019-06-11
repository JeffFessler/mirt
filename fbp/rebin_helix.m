 function [sino orbits used] = rebin_helix(cg, ig, proj, varargin)
%function [sino orbits used] = rebin_helix(cg, ig, proj, varargin)
%|
%| A single-slice rebinning (SSRB) method for cone-beam tomography data
%| collected with a helical source trajectory.
%|
%| in
%|	cg			ct_geom()
%|	ig			image_geom()
%|	proj	[ns nt na]	cone-beam projection views (line integrals)
%|
%| option
%|	'type'		''	'linear' (default) or 'nearest' (too crude)
%|	'short'		1|0	1 for short-scan fan beam (default); 0 for 360
%|	'collapse'	1|0	collapse each view to one row first
%|	'chat'		1|0	verbosity
%|
%| out
%|	sino		[ns na_rebin nz] fan-beam sinograms
%|	orbits		[nz 2]	[orbit orbit_start] for each slice (degrees)
%|	used		[nt na]	show which views were used
%|
%| See rebin_helix_example, for information on how to call
%|
%| Equations used are taken from
%| Noo F, Defrise M, Clackdoyle R and Kudo H; Phys. Med. Biol. 44:561-70, 1999
%| "Single-slice rebinning for helical cone-beam CT"
%| @u doi 10.1088/0031-9155/44/2/019
%|
%| Copyright 2010-07-21, Gregory Handy and Jeff Fessler, University of Michigan

if nargin == 1 && streq(cg, 'test')
	run_mfile_local helix_example
return
end

if nargin < 3, help(mfilename), error(mfilename), end

arg.type = 'linear';
arg.short = true;
arg.chat = false;
arg.collapse = false;
arg.use_mex = []; % see below
arg = vararg_pair(arg, varargin);

if isempty(arg.use_mex)
	if arg.short
		arg.use_mex = exist('helix_rebin_mex') == 3; % use mex if available
	else
		arg.use_mex = 0;
	end
end

if arg.use_mex && ~arg.short
	warn 'mex only for short scan; reverting to matlab'
	arg.use_mex = 0;
end

if arg.collapse % collapse each view to a single row by averaging
	proj = mean(proj, 2); % [ns 1 na]
	cg.dt = cg.dt * cg.nt;
	if cg.offset_t, fail 'offset_t needs adjusted for collapse', end
	cg.nt = 1;
end

if arg.use_mex
	nthread = 1; % for now only 1 thread supported
	i = @(x) int32(x);
	s = @(x) single(x);

	[sino orbits] = helix_rebin_mex('helix,rebin,ssrb', ...
		i(ig.nz), s(ig.dz), s(ig.offset_z), ...
		s(cg.dso), s(cg.dsd), s(cg.dfs), ...
		s([cg.ds cg.dt]), s([cg.offset_s cg.offset_t]), ...
		s(cg.pitch), s(cg.source_z0), ...
		s(cg.orbit), s(cg.orbit_start), ...
		s(proj), i(nthread));
	used = [];

else

	[sino orbits used] = rebin_helix_do(proj, ig.mask_or, ...
		ig.z, ig.dz, ig.dx, ig.nx, ig.ny, ig.nz, ...
		cg.na, cg.ns, cg.ds, cg.s, ...
		cg.nt, cg.dt, ... % cg.t, ...
		cg.dsd, cg.dso, cg.dod, cg.dfs, cg.offset_s, cg.wt, ...
		cg.orbit, cg.orbit_start, cg.pitch, cg.rmax, ...
		cg.source_zs, cg.gamma, ...
		arg.type, arg.short, arg.chat);
end

end % rebin_helix()


% rebin_helix_do()
function [sino orbits used] = rebin_helix_do(proj, mask2, ...
		zslice, dz, dx, nx, ny, nz, ...
		na, ns, ds, ss, ...
		nt, dt, ... % tpoints, ...
		dsd, dso, dod, dfs, offset_s, wt, ...
		orbit, orbit_start, pitch, rmax, ...
		source_zs, gamma, ...
		type, short, chat)

if dfs ~= 0 && dfs ~= inf
	fail 'only arc and flat done'
end

na1 = na/orbit*360; % # of views in one turn - should be integer!
if abs(round(na1) - na1) > 1e-4
	fail 'bad na/orbit'
end

orbits = zeros(nz, 2);

if short
	orbit_short_ideal = 180 + 2 * rad2deg(max(abs(gamma)));
	na1 = 2 * ceil(orbit_short_ideal / (orbit / na) / 2); % make even
	orbit_short = na1 * (orbit / na); % actual
	na1_half = na1 / 2; % because even
	orbits(:,1) = orbit_short;
else
	orbits(:,1) = 360;
	na1_half = ceil(na1 / 2);
end

sino = zeros(ns, na1, nz);
used = zeros(nt, na);
%used = zeros(nz, na);

for iz=1:nz % loop over slices (to make one fan-beam sinogram per slice)
	ticker(mfilename, iz, nz)

	% determine which helix views to use
	zmid = zslice(iz); % z coordinate of center of this slice
	ia1_middle = imin(abs(zmid - source_zs)); % "1" because matlab indexing
	ia1_list = ia1_middle - na1_half + [0:na1-1];

	% see if view extrapolation needed
	if ia1_list(1) < 1
		ia1_list = 1:na1; % resort to first views
	elseif ia1_list(end) > na
		ia1_list = [1:na1] - na1 + na; % resort to last views
	end
	orbits(iz,2) = orbit_start + (ia1_list(1)-1) / na * orbit;

	for i1=1:na1
		ia1 = ia1_list(i1);

		zdiff = zmid - source_zs(ia1); % between fan-beam plane and src

		% which rows (t coordinate) of helical projection views
		% based on point where CB ray hits fan-beam plane
		switch dfs
		case inf % flat noo:99:ssr (2)
			tt = zdiff * (dsd^2 + ss.^2) / (dsd * dso); % [ns]
		case 0 % arc
			tt = zdiff * (dsd / dso) ./ cos(gamma); % [ns]
		otherwise
			fail 'bug'
		end

		itr = tt/dt + wt; % [ns] float
		itr = max(itr, 0); % extrapolate using first row
		itr = min(itr, nt-1); % extrapolate using last row

		switch type

		case 'nearest' % nearest neighbor
			itn = round(itr); % [ns] nearest neighbor
			it1 = 1 + itn; % matlab indexing
			it1 = max(it1, 1); % extrapolate using 1st row
			it1 = min(it1, nt); % extrapolate using last row
			tt = ((it1-1) - wt) * dt; % trick: actual t value!

			% todo: scaling may not be quite right if extrapolated!
			scale = rebin_helix_scale(dfs, dsd, ss, tt); % [ns]

			tmp = [1:ns]' + (it1-1)*ns + (ia1-1)*ns*nt; % [ns]
			view = proj(tmp); % [ns]
%			view = proj(:, it1, ia1); % [ns]

			used(it1, ia1) = 1;

		case 'linear'
			tt = (itr - wt) * dt; % trick: actual t value!
			it0 = floor(itr); % [ns] nearest neighbor
			it0 = min(it0, nt-2); % so that it0+1 is ok
			frac = itr - it0; % [ns] for linear interpolation

			scale = rebin_helix_scale(dfs, dsd, ss, tt); % [ns]

			tmp = [1:ns]' + it0*ns + (ia1-1)*ns*nt; % trick
			tmp0 = proj(tmp); % [ns]
			tmp1 = proj(tmp + ns); % [ns]
%			view = (1 - frac) .* tmp0 + frac .* tmp1;
			view = tmp0 + frac .* (tmp1 - tmp0);

			used([it0+1; it0+2], ia1) = 1;

		otherwise
			fail 'bad type'
		end

		sino(:, i1, iz) = scale .* view;
	end
%	used(iz, ia1_list) = 1;

end

end % rebin_helix_do()


% rebin_helix_scale()
% scaling factor due to different ray lengths
function scale = rebin_helix_scale(dfs, dsd, ss, tt)

switch dfs
case inf % flat noo:99:ssr (1)
	scale = sqrt(ss.^2 + tt.^2 + dsd^2) ./ sqrt(tt.^2 + dsd^2);
case 0 % arc
	scale = dsd ./ sqrt(tt.^2 + dsd^2);
otherwise
	fail 'bad dfs'
end

end % rebin_helix_scale()
