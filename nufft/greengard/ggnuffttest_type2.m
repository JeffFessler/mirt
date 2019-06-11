% 1-d example of leslie greengard's gaussian kernel-based nufft
% (type 2)

R = 4; % oversampling ratio
Msp = 12; % spreading parameter - the number of neighbors to spread
          % to in each direction
M = 512; % dimension of uniform grid on which we evaluate nufft
tau = 1/M/M*pi/R/(R-0.5)*Msp;
Mr = R*M;

% precomputations
l = [-Msp+1:Msp]';

E3 = exp(-((pi*l/Mr).^2)/tau); % part of gaussian interpolation
                               % kernel
k = [-M/2:M/2-1]';
E4 = exp(tau*k.^2);

% make up irregular grid points
t = [0:M-1]'/M;
xj = 2*pi*t.^2;

% make up data points
Fk = rand(M,1);

% find the nearest grid points
m1 = floor(xj/2/pi*Mr);
xi1 = 2*pi*m1/Mr;

% calculate irregular-point-dependent terms
E1 = exp(-((xj-xi1).^2)/4/tau);
E2x = exp(pi*(xj-xi1)./Mr./tau);

% first deconvolve
Fmintau = sqrt(pi/tau)*E4.*Fk;

% then take iDFT to get fmintau
fmintau = ifft(ifftshift([zeros((Mr-M)/2,1); Fmintau; zeros((Mr-M)/2,1)]));

fj = zeros(M,1);
for j = 1:M
  for l1 = -Msp+1:Msp
    fj(j) = fj(j) + fmintau(mod(m1(j)+l1,Mr)+1)*(E2x(j)^l1)*E3(l1+Msp);
  end
end
fj = fj.*E1;

% now do slow dft, for comparison
dftmtx = exp(1i*xj(:)*k(:)');
fex = dftmtx*Fk;

% report error
disp(sprintf('Type 2 1D NRMSE: %0.2d',nrmse(fex,fj,ones(size(fex)))))
disp(sprintf('Type 2 1D MAD: %0.2d',max(abs(fex-fj))));
