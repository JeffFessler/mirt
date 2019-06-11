% parallel_feldkamp_example.m
% Example of how to use the fan to parallel rebinning feldkamp code 
% for cone-beam CT reconstruction
% Written by Gregory Handy

%% Preliminaries: create the CT geom, image geom, true image, and proj
if ~isvar('cg'), printm 'cg: cone-beam CT geometry'
	if ~isvar('down'), down = 4; end % down sample a lot to save time
	if ~isvar('dfs'), dfs = 0; end
	cg = ct_geom_par('fan', 'ns', 256, 'nt', 240, 'na', 288, ...
		'ds', 1024/256, 'dt', 1024/256, ...
		'down', down, ...
		'offset_s', 1.25, ... % quarter detector
		'offset_t', 0.0, ...
		'dsd', 949, 'dod', 408, 'dfs', dfs,'rebinned',0);
	printm('fov rmax=%g', cg.rmax)
	clear dfs
end  

if ~isvar('ig'), printm 'ig: image geometry'
	ig = image_geom('nx', 256, 'ny', 240, 'nz', 200, 'fov', 500, ...
		'down', down);
	mask2 = true([ig.nx ig.ny]);
	mask2(end) = 0; % trick: test it
	ig.mask = repmat(mask2, [1 1 ig.nz]);
	clear mask2
end
if ~isvar('ell'), printm 'ell: ellipsoid object'
	ell = [ ...
		[20 10 10	150 150 380	0 0 0.01]; % 30cm diam "cylinder
		[80 10 10	50 50 30	0 0 0.01]; % bone-like inserts
		[-10 -40 75	40 40 40	0 0 0.01];
		[-10 80 -20	30 30 30	0 0 0.01];
	];
end
if ~isvar('xtrue'), printm 'xtrue: true image volume'
	xtrue = ellipsoid_im(ig, ell);
	t = sprintf('x true, z=%g to %g', ig.z(1), ig.z(end));
	im(xtrue, t)
prompt
end

if ~isvar('proj'), printm 'proj: analytical ellipsoid projection views'
	proj = ellipsoid_proj(cg, ell);
	im(proj)
prompt
end

%% Rebinning Stage (first define 2D fan and par sinograms)
if ~isvar('sf'), printm 'sf: sinogram fan CT geometry'
	sf = sino_geom('fan', 'ns', cg.ns, 'na', cg.na, ...
		'ds', cg.ds,...
		'offset_s', cg.offset_s, ... % quarter detector
        'strip_width', cg.ds,...
		'dsd', cg.dsd, 'dod', cg.dod, 'dfs', cg.dfs);
	clear dfs
end
%be sure to keep cg.ns*cg.ds constant
dsChange = 10;
if ~isvar('sp'), printm 'sp: sinogram parallel CT geometry'
	sp = sino_geom('par', 'nb', cg.ns*dsChange, 'na', cg.na, ...
		'dr', cg.ds/dsChange,...
		'offset_r', .25, ... % quarter detector...does not work if 1.25, can be 0
        'strip_width', cg.ds,...
        'orbit', 360);
    clear dfs
end 
       
%rearranges the projections indices for the rebinning step 
%rebin_fan2par only accepts sinograms that are 2D
temp = permute(proj, [1 3 2]);
fansproj = zeros(cg.ns*dsChange,cg.na,cg.nt);
%REBINS THE FAN TO PAR
for i = 1:cg.nt
    sino2D = temp(:,:,i);
    %input is fsino, sinoFan geom, sinoPar geom)
    fansproj(:,:,i) = rebin_fan2par(sino2D, sf, sp);
end

%flips the projection indices back
finalParProj = permute(fansproj, [1 3 2]);

%interpolation/extrapolation step
if 1 || ~isvar('tfdk') %performs T-FDK

%IMPORTANT: redefine the CG geom to account for changes in the ds.
%if ~isvar('rebinnedCG'), printm 'cg: cone-beam CT geometry'
	rebinnedCG = ct_geom_par('fan', 'ns', cg.ns*dsChange, ...
		'na', cg.na, ...
		'nt', cg.nt, ... % old way
...%		'nt', numel(rectGrid), ... % todo
		'ds', cg.ds/dsChange, ...
		'dt', cg.dt, ... % caution: not actually used!
		'offset_s', sp.offset_r, ... % quarter detector
		'offset_t', 0.0, ... % todo: check!
		'dsd', cg.dsd, 'dod', cg.dod, 'dfs', cg.dfs,...
		'orbit', 360,'rebinned',true);  
%end

%if ~isvar('rectGrid')
	[interpParProj rectGrid] = gridding(rebinnedCG, finalParProj);
	printm 'Interpolation/Extrapolation Completed'
%end

	% first input is whether or not you want to use the interp/extrap
	tfdk = feldkamp_par(1,rectGrid,rebinnedCG, ig, interpParProj,'use_mex',0);
end

if ~isvar('pfdk') %performs P-FDK
	pfdk = feldkamp_par(0,rectGrid,rebinnedCG, ig, finalParProj,'use_mex',0);
end

if ~isvar('xfdk') %performs FDK
	xfdk = feldkamp(cg,ig,proj,'use_mex',0);
end

%%prints out a figure of the first and middle slice for comparison 
im plc 3 4
im(1,xtrue(:,:,1),[0,.02],'True')
im(2,tfdk(:,:,1),[0,.02],'T-FDK')
im(3,pfdk(:,:,1),[0,.02],'P-FDK')
im(4,xfdk(:,:,1),[0,.02],'FDK')

im(5,xtrue(:,:,25),[0,.02], 'True')
im(6,tfdk(:,:,25),[0,.02], 'T-FDK')
im(7,pfdk(:,:,25),[0,.02],'P-FDK')
im(8,xfdk(:,:,25),[0,.02],'FDK')


im(9, xtrue(:,end/2,:),[0,.02], 'True')
im(10, tfdk(:,end/2,:),[0,.02], 'T-FDK')
im(11, pfdk(:,end/2,:),[0,.02], 'P-FDK')
im(12, xfdk(:,end/2,:),[0,.02], 'FDK')

