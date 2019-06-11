% Simulation for time segmentation code.
% for mri object, must use mexfile mvmult.c
% to mex, simply type: mex mvmult.c at command prompt
%
% This file will form the mri structures and compare iterative
% reconstruction times for field-corrected imaging
%
% Brad Sutton, University of Michigan
%	June 2002
%
% Code contained in this has appeared in 2 articles:
%	B.P. Sutton and J.A. Fessler and D.C. Noll.
%		Fast, iterative, field-corrected image reconstruction for MRI.
%		IEEE Trans. Med. Im., Feb. 2003
%
%	J.A. Fessler and B.P. Sutton.
%		Nonuniform fast fourier transforms using min-max interpolation.
%		IEEE Trans. Signal Processing, Feb. 2003
%
%
%SIMULATION VERSION

%path(path,'./NUFFT/')


%Setup scan parameters
FOV = 22;	%in cm
N = 64;	%reconstruction size
gamp = 2.2;
gslew = 180;
TE = 25e-03;	%TE in seconds
tau = 2e-03;	%dTE between acquisitions for field map
nl = 1;	%number of interleaves
datseg = 3770;

warn 'data for this routine is no longer distributed.'
warn 'use Gmri instead (see readme).'
warn 'ask jf if you *really* need it for some reason'

%Load in simulation object and field map

%load('obj256')
ig = image_geom('nx', 256, 'dx', 1);
ell = [0 0 83 102 0 1;
	0 51 19 19 0 0.2;
	-32 -19 19 52 0 -0.3;
	32 -19 19 52 0 -0.3];
objh = ellipse_im(ig, ell, 'oversample', 3);

%load('obj64')
obj = downsample2(objh, 4);

%load('we256,140')
weh = fld_read('hide/we256,140.fld');

%load('we64,140')
we = downsample2(weh, 4);


% Load in k-space trajectory
sprintf('Loading k-space trajectory \n')
%load ktraj
ktraj = fld_read('hide/ktraj.fld');
kx = ktraj(:,1);
ky = ktraj(:,2);
ndat = length(kx)/nl;	%num data pts for one spiral
datseg = ndat;		%set num data pts to one spiral


% Load density compensation weighting
% load ww
ww = fld_read('hide/ww.fld');

% Set up timing of slice acquisition
tt = [(0:(ndat-1))*5e-06+TE]'; %Timing of non-delayed acquisition
%Shifting timing vector
%Put back into timing vector for all interleaves, tt_ext
tt_setup = 'tt_ext = [tt';
for mm = 1:nl-1
	tt_setup = sprintf('%s; tt',tt_setup);
end
tt_setup = sprintf('%s];',tt_setup);
eval(tt_setup);

%SPACE COORDS
%********************************************************
% Set up coordinates of space domain
outsz = N;
[x,y] = meshgrid([-outsz/2:outsz/2-1]./outsz);
xvals = x(:)+(1/(2*N));
yvals = y(:)+(1/(2*N));


%Set up mask
mask = true(outsz,outsz);
% Now mask out values of xval, yval, we
l = find(mask(:)>0);
xval = xvals(l);
yval = yvals(l);
npm = length(xval);	% gives number of points masked

% DATA SIGNAL - SIMULATION: simulate at 256 x 256
if ~isvar('dat2')
	sprintf('Forming data from 256x256 object by slow "DFT."\nMay take time; please be patient!')

	[x1,y1] = meshgrid([-256/2:256/2-1]./256);
	xvals1 = x1(:) +(1/(2*256));
	yvals1 = y1(:) +(1/(2*256));
	A2 = mridft(kx, ky, tt_ext, weh(:), xvals1, yvals1, FOV, 256);
	objm = objh(:);
	sigma = 25;
	nn = sigma/2*(randn(nl*datseg,1)+i*randn(nl*datseg,1));
	dat2 = A2*objm;
	SNR = sqrt(real(dat2'*dat2)/real(nn'*nn))
	apSNR = sqrt(real(dat2'*dat2)/(datseg*nl/2*sigma^2))
	dat2 = dat2+nn;
	clear xvals1 yvals1 weh objh objm x1 y1
end

%CGM
%******************************************************
% Set up parameters for CGM
niter = 10;
xinit = zeros(npm,1);


sprintf('Forming fast MRI object \n')
if 1
	A_fast = fast_mr(kx,ky,FOV,N,2*N,4,tt_ext,we(:),1,5,1);
	A_fastn = fast_mr(A_fast,'flgswt',0);
else	% USE this if memory is a problem or if runs slow the other way
	[bin_vals, bin_cens] = hist(we(:),256);
	we_histo = [bin_cens', bin_vals'];
	A_fast = fast_mr(kx,ky,FOV,N,2*N,5,tt_ext,we(:),1,5,2,we_histo);
	A_fastn = fast_mr(A_fast,'flgswt',0);
end

sprintf('Forming DFT object \n')
A_slow = mridft(kx,ky,tt_ext, we(:), xval, yval, FOV, N);

% PENALTY weighting for penalized CGM
C = C2sparse('leak', mask, 8, 0);
C = sqrt(2^(9)) * C;	% 2^(10)

if ~isvar('it_slow')
	sprintf('Now for slow iterative reconstruction \n')
	tic
	it_slow = qpwls_pcg(A_fastn'*(ww.*dat2)./(N*N),A_slow,1,dat2,0,C,1,niter+1,mask,0);
	t_slow = toc;
	sprintf('Slow iterative reconstruction done \n')
end

	sprintf('Fast iterative reconstruction starting \n')
	tic
	it_fast = qpwls_pcg(A_fastn'*(ww.*dat2)./(N*N),A_fast,1,dat2,0,C,1,niter+1,mask,0);
	t_fast = toc;
	sprintf('Fast iterative reconstruction done \n')

	cp_fast = A_fastn'*(ww.*dat2)./(N*N);


sprintf('Forming images for display \n')
	subplot(221),imagesc(abs(obj.'),[0 1.4])
	title('Simulation Object')
	axis image
	colorbar
	colormap gray

	subplot(222),imagesc(we.'./(2*pi))
	colorbar
	axis image
	title('Simulation Field Map in Hz')
	colormap gray

	subplot(234),imagesc(abs(reshape(cp_fast,64,64)./16).',[0 1.4])
	title('Fast Conjugate Phase')
	colormap gray
	axis image

	subplot(235),imagesc(abs(reshape(it_fast(:,end),64,64)./16).',[0 1.4])
	title('Fast Reconstruction')
	axis image
	colormap gray

	subplot(236),imagesc(abs((reshape(it_slow(:,end)-it_fast(:,end),64,64)./16).'))
	axis image
	colorbar
	title('Slow - Fast')
	colormap gray


sprintf('Slow iterative recon time was %f s \n fast iterative recon time was %f s \n',t_slow,t_fast)
