%% Test Script for Regularized M0/T2 Estimation from 2D Single SE Data at Variable TE
%       05.15.2015 -- Pulse sequences prepared
%       06.11.2015 -- NIST Phantom Data acquired
%       08.27.2015 -- varargin options tested
%
% Collected at two optimized TE times
% Written by: Gopal Nataraj

%% IRT Setup
%{
if (~exist('irtdir', 'var'))
	curdir = pwd; cd ../../../irt; setup(); cd(curdir);
end
%}

%% Raw data extraction
if ~isvar('y_im'), printm 'load the (possibly complex) image data'

	% addpath('../../NIST_Phantom_10,09,14/SE/');
	% load('ims_ge_se.mat');

%	datroot = 'http://web.eecs.umich.edu/~fessler/irt/data/';
%	daturl = [datroot 'nataraj-t2/2014_10_09_nist_phantom_ge3t_se5.mat'];
	urlsuff = 'nataraj-t2/2014_10_09_nist_phantom_ge3t_se5.mat';
	tmp = ir_webread(urlsuff);
	ims = tmp.ims;
	TE_ms = tmp.TE_ms;

	% code below wants magnitude data by default
	% If 2D: y_im should be [nx ny nTE]
	% If 3D: y_im should be [nx ny nz nTE]

	y_im = abs(permute(flipdim(ims, 1), [2 1 3])); 

	% Imaging parameters
	%TE = [15 30 60 120 240]';           % ms
	TE = TE_ms;
	clear TE_ms tmp ims

	im plc 2 2
	im(1, y_im)
end


% Create tight and loose masks
if ~isvar('loose_mask'), printm 'mask'
	tmp = squeeze(abs(y_im(:,:,ceil(1))) > ...
		0.1*max(col(abs(y_im(:,:,ceil(1))))));
	im(4, tmp)
%	tight_mask = imfill(tmp);
	tight_mask = ir_imfill1(ir_imfill1(tmp)')';
	im(3, tight_mask)
	tmp = strel('disk', 3);
	if ir_is_octave
		tmp = struct(tmp);
		tmp = tmp.nhood;
		loose_mask = conv2(tight_mask, tmp, 'same') > 0;
	else
		loose_mask = imdilate(tight_mask, tmp);
	end
	tmp = sum(y_im, ndims(y_im));
	tmp = tmp / max(tmp(:));
	im(2, -loose_mask + 1.8*tmp);
drawnow
end


if ~isvar('T2_rls'), printm 'T2'

	% Modifications to default parameters
	n_outer = 5;                % Less outer iterations
	is_mag = 1;                 % Use magnitude data (recommended)
	wghts = [0 1 1 1 1]';       % Drop the first echo, for instance

	% Regularized (M0, T2) estimation
	[M0_ml, M0_rls, T2_ml, T2_rls, ml_time, rls_time, cost] ...
		= mri_spinecho_m0t2map_v3(y_im, TE, 'n_outer', n_outer, ...
		'mask', loose_mask, 'is_mag', is_mag, 'weights', wghts);
end

	t2lim = [0 700];
	im(3, T2_ml, t2lim, 'ML'), cbar
	im(4, T2_rls, t2lim, 'RLS'), cbar
