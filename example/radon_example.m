% radon_example
% This m-file was based on one written by Gianni Schena at Univ. Trieste.
% It is a natural place to start because it directly compares matlab's
% radon and iradon routines with an algorithm provided in this suite.

has_radon = exist('radon.m') == 2;

if ~isvar('xtrue'), printm 'xtrue and xmat'
	ig = image_geom('nx', 64, 'ny', 64, 'dx', 1);
	xtrue = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
	xtrue = max(xtrue,0);
	im plc 2 4
	im(1, xtrue, 'xtrue'), cbar

	% project xtrue using matlab way
	na = 180/2;
	angle = [0:(na-1)]/na*180;
	if has_radon
		[yi_mat, rad] = radon(xtrue', angle); % trick: transpose!
		printm('ray spacing = %g,%g', min(diff(rad)), max(diff(rad)))
		im(2, yi_mat, '"radon" sino'), cbar
		nb = size(yi_mat,1);
		dr = rad(2) - rad(1);
	else
		warn(['\nThis uses radon() in image processing toolbox\n' ...
		'but it shows that you do not really need to buy it\n' ...
		'thanks to the tools in this free toolbox!'])
		nb = 90;
		dr = 1;
	end

	sg = sino_geom('par', 'nb', nb, 'na', na, 'dr', dr);

	% matlab reconstruction from sinogram
	if has_radon
		xmat = iradon(yi_mat, sg.ad, 'linear', 'Ram-Lak', 1, ig.nx)'; % transpose!
		im(3, xmat, '"iradon" recon'), cbar
		im(4, abs(xmat-xtrue), '|error|'), cbar
	end
prompt
end


% make system model: use Gtomo2_wtmex if possible,
% otherwise fall back to Gtomo2_strip
if ~isvar('A'),	printm 'A'

	if has_mex_jf
		A = Gtomo2_wtmex(sg, ig, 'pairs', {'strip_width', sg.dr'});
	else
		printm 'no mex, so reverting to slow Gtomo2_strip: be patient!'
		A = Gtomo2_strip(sg, ig, 'strip_width', sg.dr);
	end
prompt
end

if ~isvar('yi'), printm 'yi' % make sinogram
	yi = A * xtrue; % forward projector, aka DRR
	im(7, yi, 'sinogram'), cbar
	if has_radon
		im(8, abs(yi-yi_mat), 'sino |diff|'), cbar
	end
prompt
end

warn('todo: im toggle shows some rotation (?) between yi and yi_mat')
if 0 && has_radon
	clf, im(yi-yi_mat, 'sino diff'), cbar
	im_toggle(yi, yi_mat)
end


% reconstruct using em_fbp FBP algorthm
if ~isvar('xfbp1'), printm 'xfbp1'
	fbp_kernel = [1]; % noiseless data, so no filtering
	xfbp1 = em_fbp(sg, ig, max(yi,0), [], [], 'kernel', fbp_kernel);
	im(5, xfbp1, 'em\_fbp'), cbar
	im(6, abs(xfbp1-xtrue), 'error'), cbar
prompt
end

if ~isvar('Rq'), printm 'Rq'
	l2b = 2; % log_2(beta)
	Rq = Reg1(ig.mask, 'beta', 2^l2b);

	if 1
		psf = qpwls_psf(A, Rq, 1, ig.mask);
		printm('expected FWHM = %g', fwhm2(psf))
	end
prompt
end


% CG iterative reconstruction
if ~isvar('xcg'), printm 'xcg'
	W = 1;
	niter = 10;
	xcg = qpwls_pcg(xfbp1(ig.mask), A, W, yi(:), 0, ...
		Rq.C, 1, niter, ig.mask, 1);
	xcg = ig.embed(xcg(:,end));
	im(7, xcg, 'CG'), cbar
	im(8, abs(xcg-xtrue), 'error'), cbar
prompt
end


% run the block-iterative SPS-OS algorithm
if ~isvar('xsps'), printm 'xsps'
	nblock = 10;
	Ab = Gblock(A, nblock);

	if ~isvar('R'),	printm 'R'
		R = Reg1(ig.mask, 'type_denom', 'matlab', ...
			'pot_arg', {'hyper3', 0.1}, 'beta', 2^l2b);
	end

	xsps = pwls_sps_os(xfbp1(ig.mask), yi, [], Ab, R, ...
		niter, inf, [], [], 1, 1);
	xsps = ig.embed(xsps(:,end));
	im(7, xsps, 'SPS-OS'), cbar
	im(8, abs(xsps-xtrue), 'error'), cbar
end


% compare NRMSE
pr nrms(xfbp1(:), xtrue(:))
pr nrms(xcg(:), xtrue(:))
pr nrms(xsps(:), xtrue(:))
if has_radon
	pr nrms(xmat(:), xtrue(:))
end
