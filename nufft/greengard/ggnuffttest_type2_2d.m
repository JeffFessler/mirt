% 1-d example of leslie greengard's gaussian kernel-based nufft
% (type 2)

R = 4; % oversampling ratio
Msp = 12; % spreading parameter - the number of neighbors to spread
          % to in each direction
M = 32; % dimension of uniform grid on which we evaluate nufft
tau = 1/M/M*pi/R/(R-0.5)*Msp; % width of gaussian
Mr = R*M; % size of oversampled grid

% precomputations
l = [-Msp+1:Msp]';

E3 = exp(-((pi*l/Mr).^2)/tau); % part of gaussian interpolation
                               % kernel
k = [-M/2:M/2-1]';
[kdc1,kdc2] = ndgrid(k);
E4 = exp(tau*(kdc1.^2+kdc2.^2));

% make up irregular grid points
t = [-M/2:M/2-1]'/M*2;
[xj,yj] = ndgrid(pi*sign(t).*t.^2);

% make up data points
Fk = rand(M,M);

% find the nearest grid points
m1 = floor(xj/2/pi*Mr);
m2 = floor(yj/2/pi*Mr);
xi1 = 2*pi*m1/Mr;
xi2 = 2*pi*m2/Mr;

% calculate irregular-point-dependent terms
E1 = exp(-((xj-xi1).^2 + (yj-xi2).^2)/4/tau);
E2x = exp(pi*(xj-xi1)./Mr./tau);
E2y = exp(pi*(yj-xi2)./Mr./tau);

disp 'in-script nufft'

% first deconvolve
Fmintau = pi/tau*E4.*Fk;

% then take iDFT to get fmintau
padval = (Mr-M)/2;
fmintau = fftn_fast(fftshift(pad(Fmintau,padval,padval,padval,padval,0)))/Mr/Mr;

fj = zeros(M,M);
for j = 1:M*M
  for l2 = -Msp+1:Msp
    for l1 = -Msp+1:Msp
      fj(j) = fj(j) + fmintau(mod(m1(j)+l1,Mr)+1,mod(m2(j)+l2,Mr)+1)*(E2y(j)^l2)*(E2x(j)^l1)*E3(l1+Msp)*E3(l2+Msp);
    end
  end
end
fj = fj.*E1;

% now do object, to test
disp 'creating object'
G = gg_Gnufft([xj(:),yj(:)],M,Msp,R);
disp 'object nufft'
fjob = G*Fk(:);fjob = reshape(fjob,M,M);

% now do slow dft, for comparison
disp 'slow dft'
[K1,K2] = ndgrid(k);
dftmtx = exp(-1i*(xj(:)*K1(:)'+yj(:)*K2(:)'));
fex = reshape(dftmtx*Fk(:),M,M);

% report error
disp(sprintf('Type 2 2D NRMSE, script: %0.2d',nrmse(fex,fj,ones(size(fex)))))
disp(sprintf('Type 2 2D MAD, script: %0.2d',max(abs(fex(:)-fj(:)))));
disp(sprintf('Type 2 2D NRMSE, object: %0.2d',nrmse(fex,fjob,ones(size(fex)))))
disp(sprintf('Type 2 2D MAD, object: %0.2d',max(abs(fex(:)-fjob(:)))));
