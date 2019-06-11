% [footblur, radii] = blob_li_blur(R, dR, alpha, kb_m, blur, du, dy, Ny, chat)
%   Calculate line integrals of a radially symmetric Kaiser-Bessel blob
%   and blur them with each one of a set of 2D blur PSF's
%
%   R:		Radius of blob support [mm] (default = 2*dy)
%   dR:		Sample spacing [mm] in radial direction
%   		of precomputed blob line integrals (default = dy/10)
%   alpha:	Blob shape parameter (default = 10.4)
%   kb_m:	Order parameter (smoothness of blob,
%		controls rate of decay of FT, default = 2)
%
%   blur [Nx_psf Nz_psf Npsf]:	Blur PSF's (default is no blur) 
%				or use 'uhe' to generate PSF's for UHE collimator and I-131
%				or use 'tc' to generate PSF's for Tc studies
%   du:		Bin spacing [mm] of imaging system (default = dy)
%   dy:		Pixel size [mm] of imaging system (default = 4)
%   Ny:		Matrix size of imaging system (default = 64)

function [footblur, radii] = blob_li_blur(R, dR, alpha, kb_m, blur, du, dy, Ny, chat)

J = 2*R;			% Blob diameter

if ~isvar('Ny') | isempty(Ny), Ny = 64; end
if ~isvar('dy') | isempty(dy), dy = 4; end
if ~isvar('du') | isempty(du), du = dy; end
if ~isvar('dR') | isempty(dR), dR = dy/10; end
if ~isvar('R') | isempty(R), R = 2*dy; end
if ~isvar('alpha') | isempty(alpha), alpha = 10.4; end
if ~isvar('kb_m') | isempty(kb_m), kb_m = 2; end

if ~isvar('blur') | isempty(blur), sb = '';
elseif ischar(blur), sb = [',' blur]; else, sb = ',blur'; end
projfile = sprintf('blobli%s,R%g,dR%g,alpha%g,m%d.fld', sb, R, dR, alpha, kb_m);
radfile = sprintf('radii%s,R%g,dR%g,alpha%g,m%d.fld', sb, R, dR, alpha, kb_m);

if ~isvar('chat') | isempty(chat), chat = true; end

if exist(projfile) == 2
	if chat
		disp 'Loading projection tables from file... ';
		disp(projfile);
	end

	footblur = double(fld_read(projfile));
	radii = double(fld_read(radfile));
else
	if chat, disp 'Generating projection tables... ', end

	tic

	if ~isvar('blur') | isempty(blur)		% Default: No blur
		blur = [1];
	elseif strcmp(blur, 'uhe')
		if chat, disp 'Generating UHE Gaussian PSF''s for I studies...', end

		obj_det = 230;		% Object-to-detector distance (mm)
		blur = psfs_uhe([], obj_det, 2, Ny, dR, dy);

		if chat, disp '... done.', end
	elseif strcmp(blur, 'tc')
		if chat, disp 'Generating Gaussian PSF''s for Tc studies...', end

		obj_det = 230;		% Object-to-detector distance (mm)
		blur = psfs_tc([], obj_det, 2, Ny, dR, dy);

		if chat, disp '... done.', end
	end

	if chat, disp 'Generating line integrals...', end

	Nr = 2*ceil(R/dR) + 1;

	r = [-(Nr-1)/2 : (Nr-1)/2]' * dR;
	[rx, rz] = ndgrid(r);
	rad = sqrt(rx.^2 + rz.^2);			% Distance of current point from center of blob

	foot = kaiser_bessel_li(rad, J, alpha, kb_m);	% Line integrals of Kaiser-Bessel blob

	if chat, disp '... done.', end
	if chat, disp 'Blurring line integrals...', end

	[Nx_psf, Nz_psf, Npsf] = size(blur);		% Number of PSF's
	Nrx = Nx_psf + Nr - 1;
	Nrz = Nz_psf + Nr - 1;
	Nx_pad = 2^nextpow2(Nrx);
	Nz_pad = 2^nextpow2(Nrz);

	footpad = zpad0(foot, Nx_pad, Nz_pad);
	footft = fft2(fftshift(footpad));
	blurpad = zpad0(blur, Nx_pad, Nz_pad, Npsf);

	for k = 1:Npsf
		if chat, k, end

		if sum(abs(col(blur(:, :, k)))) > 0
			%blurft = reale(fft2(fftshift(blurpad(:, :, k))));
			blurft = fft2(fftshift(blurpad(:, :, k)));
			tmp = reale(fftshift(ifft2(blurft .* footft)));
			norm = sum(tmp(:)) * dR^2 / du^2;
			tmp = tmp / norm;	% Normalize to get unity sum
			tmp(find(tmp < 0.000001)) = 0.0;
			footblur(:, :, k) = tmp(ceil((Nx_pad - Nrx)/2) + (1:Nrx), ...
						ceil((Nz_pad - Nrz)/2) + (1:Nrz));
		else	% No blur
			tmp = zpad0(foot, Nrx, Nrz);
			norm = sum(tmp(:)) * dR^2 / du^2;
			footblur(:, :, k) = tmp / norm;	% Normalize to get unity sum
		end
	end

	%plot(r, footblur(:, :, end)), title 'Blob blurred by the widest blur kernel', prompt

	if chat, disp '... done.', end

	% Save only one quadrant (assume symmetry)
	footblur = footblur(ceil(Nrx/2):end, ceil(Nrz/2):end, :);

	% Crop surrounding zero values
	cropx = min(sum(footblur(:, :, end) == 0, 1));
	cropz = min(sum(footblur(:, :, end) == 0, 2));
	footblur = footblur(1:end-cropx, 1:end-cropz, :);

	% Save blurred footprints
	fld_write(projfile, footblur);
	toc

	% Save radii of blurred footprints in mm
	radx = max(sum(footblur > 0, 1));
	radz = max(sum(footblur > 0, 2));
	radii = max(radx(:), radz(:));		% In number of samples
	radii = (radii-1) * dR;			% In mm
	fld_write(radfile, radii);
end

