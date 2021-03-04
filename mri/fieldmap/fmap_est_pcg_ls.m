 function [out,cost,time] = fmap_est_pcg_ls(w, y, delta, smap, varargin)
%function [out,cost,time] = fmap_est_pcg_ls(w, y, delta, smap, varargin)
%|
%| Phase unwrapping of multiple data sets (w denotes frequency)
%| using a NCG with monotonic line search
%| 3D, multi-coil version
%|
%| cost(w) = sum(i=0 to n-1) sum(j=0 to n-1)
%|		|yi*yj| wj (1 - cos(w*(dj-di) + \angle(yi) - \angle(yj)) + R(w)
%|
%| in
%|	w	[np 1]		initial estimate
%|	y	[np nc n]	n sets of measurements for nc coils
%|	delta	[1 n]	row vector of n echo time offsets
%|  smap [np nc]    coil maps
%|
%| options
%|	stepper {'qs',# its}	monotonic line search parameters (default: {})
%|	niter			# of iterations (def: 30)
%|	maskR	[(np)]	logical reconstruction mask (required!)
%|	order			order of the finite diff matrix (def: 2)
%|	l2b             regularization parameter (2^) (def: -6)
%|	gammaType		CG direction: PR = Polak-Ribiere (default) or FR = Fletcher-Reeves
%|	precon			Preconditioner: 'diag', 'chol', 'ichol' (def: 'ichol')
%|	reset			# of iterations before resetting direction (def: inf)
%|  df              delta f value in water-fat imaging (def: 0)
%|  relamp          relative amplitude in multipeak water-fat  (def: 1)
%|  tol             tolerance for ichol (def: 1e-3)
%|
%| out
%|	out.ws	[np niter+1]	iterates
%|  out.xw / out.wf [np 1] water / fat images if arg.df~=0
%|  cost   [niter+1 1]   (nonconvex) cost for each iteration
%|	time   [niter+1 1]	time for each iteration
%|
%| This algorithm is based on the paper:
%| C Y Lin, J A Fessler, 
%| "Efficient Regularized Field Map Estimation in 3D MRI", TCI 2020

if nargin >= 1 && streq(w, 'test')
	out = fmap_est_pcg_ls_test(w, varargin{:});
	if ~nargout, clear out, end
return
end

if nargin < 4, help(mfilename), error(mfilename), end

arg.stepper = {'qs', 3}; % quad surr with this # of subiterations
arg.niter = 30;
arg.maskR = [];
arg.order = 2;
arg.l2b = -6;
arg.gammaType = 'PR';
arg.precon = 'ichol';
arg.reset = inf;
arg.df = 0;
arg.relamp = 1;
arg.tol = 1e-3; % factor of tol, for ichol
arg = vararg_pair(arg, varargin);

%% prepare variables and finite differencing matrices
w = double(w);
y = double(y);

if isempty(arg.maskR)
	fail('Mask required!')
end

% create the sparse regularization matrix
R = Reg1(arg.maskR, 'beta', 2^arg.l2b, 'order', arg.order, ...
	'distance_power', 2, 'type_diff', 'spmat', 'type_penal', ...
	'mat');
C = R.C;
if ~issparse(C)
	fail('CC = C^H * C is too slow if not sparse')
end
CC = C' * C;

% apply mask to data if necessary
if length(w) ~= sum(arg.maskR(:))
	w = w(arg.maskR,:);
	y = y(arg.maskR,:,:);
end

[np,nc,n] = size(y);
time = zeros(arg.niter+1,1);
cost = zeros(arg.niter+1,1);

% check the data size matches the echo times
if n ~= size(delta,2), fail 'need delta to be [1 n]', end

%% calculate the magnitude and angles used in the data-fit curvatures
abss = abs(smap);
sjtotal = sum(abss.^2, 2); %[np,1]
angy = angle(y);
angs = angle(smap);
if arg.df
    Gamma = phiInv(arg.relamp, arg.df, delta); %[L,L]
end
set = 1; 
nset = cumsum(1:n-1);nset = nset(end);
wj_mag = zeros(np,nset,nc,nc);
d2 = zeros(1,nset);
ang2 = zeros(np,nset,nc,nc);

for j=1:n % for each pair of scans
    for i=1:n
        if i<j % only need one pair of differences
            d2(set) = delta(i) - delta(j);
            for c = 1:nc
                for d = 1:nc
                    wj_mag(:,set,c,d) = smap(:,c) .* conj(smap(:,d)) .*...
                        conj(y(:,c,i)) .* y(:,d,j); 
                    wj_mag(:,set,c,d) = abs(wj_mag(:,set,c,d));
                    % difference of the echo times and angles
                    ang2(:,set,c,d) = angs(:,c) - angs(:,d) + ...
                         angy(:,d,j) - angy(:,c,i);
                    if arg.df
                        wj_mag(:,set,c,d) = wj_mag(:,set,c,d)*abs(Gamma(i,j));
                        ang2(:,set,c,d) = ang2(:,set,c,d) + angle(Gamma(i,j));
                    end 
                end
            end
            set = set+1;
        end
    end
end
% compute |s_c s_d' y_dj' y_ci| /L/s * (tj - ti)^2
sjtotal(sjtotal==0) = 1; %avoid outside mask = 0 
wj_mag = wj_mag./sjtotal;
if ~arg.df
    wj_mag = wj_mag/n;
end
wm_deltaD = wj_mag .* d2;
wm_deltaD2 = wj_mag .* (d2.^2);
ang2(isnan(ang2))=0; %avoid atan causing nan 

% prepare outpute variables
out.ws = zeros(length(w(:)), arg.niter+1);
out.ws(:,1) = w;

if strcmp(arg.precon,'diag')
	dCC = diag(CC);
end

% initialize projections and NCG variables
CCw = CC * w;
oldinprod = 0;
warned.dir = 0;
warned.step = 0;

%% begin iterations
% start timing the iterations
tt = tic;
fprintf('\n ********** ite_solve: NCG-MLS **********\n')
for iter=1:arg.niter
	% compute the gradient of the cost function and curvatures
    [hderiv,hcurv,sm] = Adercurv(d2,ang2,wm_deltaD,wm_deltaD2,w);
    
	grad = hderiv + CCw;
	ngrad = -grad;
    cost(iter) = sum(wj_mag.*(1-cos(sm)),'all') + norm(C*w,'fro');

    fprintf(' ite: %d , cost: %f3\n', iter-1, cost(iter)) 

	% apply preconditioner
	switch arg.precon
	case 'diag'
		H = hcurv + dCC;
		npregrad = ngrad ./ H;
	case 'chol'
		%spparms('spumoni',2) % this will let you see if sparse Cholesky is used
		H = spdiag(hcurv) + CC;
        L = chol(H, 'lower');
		npregrad = L' \ (L \ ngrad);
	case 'ichol'
		H = spdiag(hcurv) + CC;
		alpha = max(max(sum(abs(H),2) ./ diag(H)) - 2,0); % cl max with 0
        if arg.tol
            L = ichol(H,struct('type','ict','droptol',arg.tol*max(H,[],'all'),'diagcomp',alpha));
        else
            L = ichol(H);
        end
        npregrad = L' \ (L \ ngrad);
    otherwise
		npregrad = ngrad;
    end
	% compute CG direction
	newinprod = ngrad' * npregrad;
	newinprod = real(newinprod); %should be real, but just in case

	if oldinprod == 0 || mod(iter, arg.reset) == 0
		ddir = npregrad;
		gamma = 0;
		ngradO = ngrad;
	else
		if strcmp(arg.gammaType,'FR') % Fletcher-Reeves
			gamma = newinprod / oldinprod;
			ddir = npregrad + gamma * ddir;

		elseif strcmp(arg.gammaType,'PR') % Polack-Ribeir
			gamma = real((ngrad - ngradO)' * npregrad) / oldinprod;
			ngradO = ngrad;

			if (gamma < 0)
				printm('RESETTING GAMMA, iter=%d', iter)
				gamma = 0;
			end

			ddir = npregrad + gamma * ddir;

		end
	end
	oldinprod = newinprod;

	% check if correct descent direction
	if ddir' * grad > 0
		if ~warned.dir
			warned.dir = 1;
			warn 'wrong direction so resetting'
			printm('<ddir,grad>=%g, |ddir|=%g, |grad|=%g', ...
				ddir' * grad, norm(ddir), norm(grad))
		end
		% reset direction if not descending
		ddir = npregrad;
		oldinprod = 0;
	end

	% step size in search direction
	Cdir = C * ddir; % caution: can be a big array for 3D problems

	% compute the monotonic line search using quadratic surrogates
	CdCd = Cdir'*Cdir;
	CdCw = ddir'*CCw;
	step = 0;
    for is=1:arg.stepper{2}
        
        % compute the curvature and derivative for subsequent steps
        if step ~= 0
            [hderiv,hcurv] = Adercurv(d2,ang2,wm_deltaD,wm_deltaD2,w + step * ddir);
        end

        % compute numer and denom of the Huber's algorithm based line search
        denom = (ddir.^2)' * hcurv + CdCd;
        numer = ddir' * hderiv + (CdCw + step * CdCd); 

        if denom == 0
            warn 'found exact solution??? step=0 now!?'
            step = 0;
        else
            % update line search
            step = step - numer / denom;
        end

    end

	% update the estimate and the finite differences of the estimate
	CCw = CCw + step * C' * Cdir;
	w = w + step * ddir;

	% save any iterations that are required (with times)
	out.ws(:,iter+1) = w;
    time(iter+1) = toc(tt);
	% display counter
	if (mod(iter,500) == 1)
		printm([num2str(iter) ' of ' num2str(arg.niter)])
	end
end

sm = w * d2 + ang2;
cost(iter+1) = sum(wj_mag.*(1-cos(sm)),'all') + norm(C*w,'fro');

fprintf(' ite: %d , cost: %f3\n', iter, cost(iter+1)) 

%output water & fat images
if arg.df
    x = decomp(w,arg.relamp,arg.df,delta,smap,y);
    out.xw = x(1,:).';
    out.xf = x(2,:).';
end

end

function [hderiv, hcurv, sm] = Adercurv(d2,ang2,wm_deltaD,wm_deltaD2,w)
% compute the data-fit derivatives and curvatures as in Funai paper
sm = w * d2 + ang2;
hderiv = 2 * sum(wm_deltaD .* sin(sm), [2:4]);

srm = mod(sm + pi,2*pi) - pi;
hcurv = 2 * sum(wm_deltaD2 .* ir_sinc_nopi(srm), [2:4]);
end

function Gamma = phiInv(relamp,df,delta)
n = length(delta);
phi = [ones(n,1) sum(relamp.*exp(1i*delta(:)*df),2)]; %[n,2]
Gamma = phi*inv(phi'*phi)*phi';
end

function x = decomp(w,relamp,df,delta,smap,y)
[np,nc,n] = size(y);
phi = [ones(n,1) sum(relamp.*exp(1i*delta(:)*df),2)]; %[n,2]
x = zeros(2,np);
for ip = 1:np
    B = phi.*col(exp(1i*w(ip)*delta(:))); %[n,2]
    B = kron(col(smap(ip,:)),B); %[n*nc,2]
    yc = permute(y(ip,:,:),[1,3,2]);
    x(:,ip) = B\yc(:);
end
end

function wmap = fmap_est_pcg_ls_test(type, varargin)
% test example
printm 'simulate noisy multicoil 3d data'
etime = [0 2 10] * 1e-3; % echo times
SNR = 30; % dB
ne = length(etime);
wtrue = 2*pi * ir_get_data(fullfile('mri','2001-phase-data','fieldmap128.fld')); % "true" fieldmap
mag = ir_get_data(fullfile('mri','2001-phase-data','mag128.fld')); % "true" magnitude
[nx ny] = size(mag); % 128
nz = 2; nc = 4;
wtrue = repmat(wtrue,[1,1,nz]);
mag = repmat(mag,[1,1,nz]);
mask = mag > 0.05 * max(mag(:));
for iz = 1:nz
    mask(:,:,iz) = bwconvhull(mask(:,:,iz), 'union');
end
smap = double(ir_mri_sensemap_sim('nx', nx, 'ny', ny, 'nz', nz, 'ncoil', nc));

im plc 2 3
im(1, mag, 'true mag'), cbar
im(2, mask,'mask'), cbar
im(3, smap, 'sense map'), cbar
clim = [-40, 128];
im(4, wtrue/(2*pi), 'true field map', clim), cbar('Hz')

image_power = 10*log10(sum(mag.^2,'all')/(nx*ny*nz));
noise_power = image_power - SNR;
noise_std = sqrt(10^(noise_power/10));
noise_std = noise_std / 2; % because complex
yik = zeros(nx,ny,nz,nc,ne);
for kk=1:ne
    yik(:,:,:,:,kk) = mag ...
        .* exp(1i * wtrue * (etime(kk) - etime(1))) ...
        .* smap; 
end
rng(0)
yik = yik + noise_std * (randn(size(yik)) + 1i * randn(size(yik)));

yik_c = reshape(yik,[nx*ny*nz,nc,ne]);
smap_c = reshape(smap,[nx*ny*nz,nc]);
yik_sos = reshape(sum(yik.*reshape(conj(smap),[nx,ny,nz,nc]),4),[nx,ny,nz,ne]);%coil combine
winit = angle(stackpick(yik_sos,2) .* conj(stackpick(yik_sos,1))) ...
    / (etime(2) - etime(1));

printm 'estimate field map'
[out,cost,time] = fmap_est_pcg_ls(winit(mask),yik_c(mask,:,:),etime, ...
    smap_c(mask,:),'maskR', mask,'l2b',-3,'niter',20,'order',1);
wmap = embed(out.ws(:,end),mask);

	finit = winit(mask) / (2*pi);
	ftrue = wtrue(mask) / (2*pi);
	fpcg = wmap(mask) / (2*pi);
	rmse_init = sqrt(sum((finit - ftrue).^2) / sum(mask(:)));
	rmse_pcg = sqrt(sum((fpcg - ftrue).^2) / sum(mask(:)));
	im(5, embed(finit,mask), 'initial field map', clim), cbar('Hz')
	titlef('initial field map, RMSE %3.1f Hz', rmse_init)
	im(6, embed(fpcg,mask), 'regularized field map', clim), cbar('Hz')
	titlef('regularized field map, RMSE %3.1f Hz', rmse_pcg)
end
