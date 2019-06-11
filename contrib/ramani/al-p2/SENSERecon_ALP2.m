% by Sathish Ramani
function [xnew, CAL, TAL, EAL, ERR] = ...
	SENSERecon_ALP2(img, SP3, Data, recon_F, d, z, xnew, params)

%% Compute SER at current estimate
dispitr = params.dispitr;
dispitrnum = params.dispitrnum;

dispfig = params.dispfig;
Npix = params.Npix;
mn = params.mn;
mx = params.mx;
rs = params.rs;
cs = params.cs;
figno = params.figno;

maxitr = params.maxitr;
maxitr_in = params.maxitr_in;
dxtol = params.dxtol;
dcosttol = params.dcosttol;

xinf = params.xinf; %% Solution obtained after infinite number of iterations
xinfnorm = params.xinfnorm;

recon_sos = params.recon_sos;

smap = params.smap;

mu = params.AL.mu;
nu1 = params.AL.nu1;
nu2 = params.AL.nu2;
nuratio = nu2/nu1;

iPpmu = params.AL.iPpmu;
iSpnu2 = params.AL.iSpnu2;
iRpnu2nu1 = params.AL.iRpnu2nu1;

sv = params.AL.sALP2;
ev = params.AL.eALP2;

%% Compute Error in estimate
err = xnew - img;
errnorm = sqrt(sum(abs(err(:)).^2)/Npix)/mx;

%% Main while loop for varying sigma
tE1 = zeros(maxitr_in, 8);
tE = zeros(1, 2);
time_elapse = 0;
itr = 1;

%% Prepare cost and thresholds
cost_new = compute_Cost_SENSE(SP3, Data, xnew, params); % Reinitialize cost_old for current sigma
cost_old = 2*cost_new;
diffcost = cost_new - cost_old;

normx = sqrt(sum(abs(xnew(:)).^2));
diffx = normx;

CAL(itr) = cost_new;
TAL(itr) = time_elapse;
EAL(itr) = sqrt(sum(abs(xnew(:) - xinf(:)).^2))/xinfnorm;
ERR(itr) = errnorm;

%% Show images
if (dispfig)
	im(params.subplot_img, abs(xnew))
	titlef('iter = %g', itr)
	xlabelf('nrme = %g dB', 20*log10(ERR(itr)))
	im(params.subplot_err, abs(xnew - img))
end

if 0 && dispfig
	figure(figno); clf;
	subplot(2,2,1); imagesc(abs(img), [mn mx]); colormap(gray); colorbar; axis off; axis equal; title('Unknown original = x');
	subplot(2,2,2); imagesc(abs(err)); colormap(gray); colorbar; axis off; axis equal; title(cat(2,'x - x', int2str(itr), '; ', num2str(ERR(itr))));
	subplot(2,2,3); imagesc(abs(recon_sos), [mn mx]); colormap(gray); colorbar; axis off; axis equal;
	subplot(2,2,4); imagesc(abs(xnew), [mn mx]); colormap(gray); colorbar; axis off; axis equal; title(cat(2, 'Current estimate x', int2str(itr)));
end

%% Display cost and other parameters
if(dispitr)
    disp('--------------------------------------------------------------------');
    disp(['ItrALP2 = ' int2str(itr-1), '; C = ' num2str(cost_new) '; DC = ' num2str(diffcost) '; l2D = ' num2str(EAL(itr)) ...
          '; Dx/||x|| = ' num2str(diffx/normx)  '; NRMSE =' num2str(errnorm) '; Time = ' num2str(time_elapse)]);
end

%% Inner while loop for minimizing cost corresponding to a given sigma
Sx = smap .* repmat(xnew, [1 1 ev(1)-sv(1)+1]);
z(:, :, sv(2):ev(2)) = xnew;
RW = 0;

while((itr <= maxitr) && (diffx/normx > dxtol) && (diffcost < dcosttol))
    itr_in = 1; % Inner iterations for minimization over z's
    cost_old = cost_new;
    xold = xnew;
    
    while(itr_in <= maxitr_in)
        % Compute u0 (Refer revised manuscript)
        tS = tic;
        rhsz = recon_F + mu*(Sx + d(:, :, sv(1):ev(1)));
        z(:, :, sv(1):ev(1)) = ifft2(iPpmu.*fft2(rhsz));
        tE1(itr_in, 1) = toc(tS);
        
	if 0 && (~mod(itr-1, dispitrnum))
		figure(10); im( z( :, :, sv(1):ev(1))); pause(1);
	end
        
        % Compute u1 (Refer revised manuscript)
        [RW, tE1(itr_in, 2)] = doR(z(:, :, sv(2):ev(2)), params); % First do Ru2
        
        tS = tic;
        Rw = RW + d(:, :, sv(3):ev(end));
        tE1(itr_in, 3) = toc(tS);
        
        [z(:, :, sv(3):ev(end)), tE1(itr_in, 4)] = thresholdRw(Rw, params, nu1*mu); % Perform thresholding

        % Compute u2 (Refer revised manuscript)
        tS = tic;
        vmd = z(:, :, sv(3):ev(end)) - d(:, :, sv(3):ev(end));  % Construct v-dv
        tE1(itr_in, 5) = toc(tS);
        
        [rhsv, tE1(itr_in, 6)] = doRAdj(vmd, params);                             % doRAdj(v-dv) first
        
        tS = tic;
        rhsw = nuratio*(xnew + d(:, :, sv(2):ev(2))) + rhsv;              % Construct the rhs for the system
        z(:, :, sv(2):ev(2)) = ifft2(iRpnu2nu1.*fft2(rhsw));
        tE1(itr_in, 7) = toc(tS);

        % Compute x (Refer notes)
        tS = tic;
        rhsx = nu2*(z(:, :, sv(2):ev(2)) - d(:, :, sv(2):ev(2))) + sum(conj(smap) .* (z(:, :, sv(1):ev(1)) - d(:, :, sv(1):ev(1))), 3);
        xnew = iSpnu2.*rhsx;
        
        % Book Keeping
        Sx = smap .* repmat(xnew, [1 1 ev(1)-sv(1)+1]);
        RW = doR(z(:, :, sv(2):ev(2)), params);
        tE1(itr_in, 8) = toc(tS);

        itr_in = itr_in + 1;
    end
    tE(1) = sum(tE1(:));
    
    %% Update d
    tS = tic;
    d(:, :, sv(1):ev(1)) = d(:, :, sv(1):ev(1)) - (z(:, :, sv(1):ev(1)) - Sx);
    d(:, :, sv(2):ev(2)) = d(:, :, sv(2):ev(2)) - (z(:, :, sv(2):ev(2)) - xnew);
    d(:, :, sv(3):ev(end)) = d(:, :, sv(3):ev(end)) - (z(:, :, sv(3):ev(end)) - RW);
    tE(2) = toc(tS);
    
    %% End Timer
    time_elapse = sum(tE(:));
    itr = itr + 1;

    %% Compute Cost at current estimate
    cost_new = compute_Cost_SENSE(SP3, Data, xnew, params); % Cost at current iterate
    diffcost = cost_new - cost_old;

    %% Compute error
    err = xnew - img;
    errnorm = sqrt(sum(abs(err(:)).^2)/Npix)/mx;

    CAL(itr) = cost_new;
    TAL(itr) = time_elapse;
    EAL(itr) = sqrt(sum(abs(xnew(:) - xinf(:)).^2))/xinfnorm;
    ERR(itr) = errnorm;

    %% Compute the norm difference between successive iterates
    diffx = sqrt(sum(abs(xnew(:)-xold(:)).^2));

	%% Show updates of cost, estimate, etc

	%% Show images
	if dispfig && ~mod(itr-1, dispitrnum)
		im(params.subplot_img, abs(xnew))
		titlef('iter = %g', itr)
		xlabelf('nrme = %g dB', 20*log10(ERR(itr)))
		im(params.subplot_err, abs(xnew - img))
	end

	if 0 && dispfig && mod(itr-1, dispitrnum)
            figure(figno); clf;
            subplot(2,2,1); imagesc(abs(img), [mn mx]); colormap(gray); colorbar; axis off; axis equal; title('Unknown original = x');
            subplot(2,2,2); imagesc(abs(err)); colormap(gray); colorbar; axis off; axis equal; title(cat(2,'x - x', int2str(itr), '; ', num2str(ERR(itr))));
            subplot(2,2,3); imagesc(abs(recon_sos), [mn mx]); colormap(gray); colorbar; axis off; axis equal; title('Initial estimate x0');
            subplot(2,2,4); imagesc(abs(xnew), [mn mx]); colormap(gray); colorbar; axis off; axis equal; title(cat(2, 'Current estimate x', int2str(itr)));
            pause(0.1);
	end
    
	%% Display cost and other parameters
	if dispitr && ~mod(itr-1, dispitrnum)
		disp('--------------------------------------------------------------------');
		disp(['ItrALP2 = ' int2str(itr-1), '; C = ' num2str(cost_new) '; DC = ' num2str(diffcost) '; l2D = ' num2str(EAL(itr)) ...
                '; Dx/||x|| = ' num2str(diffx/normx)  '; NRMSE =' num2str(errnorm) '; Time = ' num2str(time_elapse)]);
	end
end
