% 1-d example of leslie greengard's gaussian kernel-based nufft,
% type 1.

R = 4; % oversampling ratio
Msp = 12; % spreading parameter - the number of neighbors to spread
          % to in each direction
M = 512; % dimension of uniform grid on which we evaluate nufft
tau = 1/M/M*pi/R/(R-0.5)*Msp; % width of gridding kernel (see greengard)
Mr = R*M; % points in oversampled grid

% precomputations
l = [-Msp+1:Msp]';

E3 = exp(-((pi*l/Mr).^2)/tau); % part of gaussian interpolation
                               % kernel
k = [-M/2:M/2-1]';
kdeconv = [-Mr/2:Mr/2-1]';
E4 = exp(tau*kdeconv.^2);

% make up irregular grid points
t = [0:M-1]'/M;
xj = 2*pi*t.^2;

% make up data points
fj = rand(M,1);

% find the nearest grid points
m1 = floor(xj/2/pi*Mr);
xi1 = 2*pi*m1/Mr;

% calculate irregular-point-dependent terms
E1 = exp(-((xj-xi1).^2)/4/tau);
E2x = exp(pi*(xj-xi1)./Mr./tau);

V0 = fj.*E1;
ftau = zeros(Mr,1);
for j = 1:M
  for l1 = -Msp+1:Msp
    ftau(mod(m1(j)+l1,Mr)+1) = ftau(mod(m1(j)+l1,Mr)+1) + V0(j)* ...
        (E2x(j)^l1)*E3(l1+Msp);
  end
end

% take fft of ftau
Ftau = fftshift(fft(ftau));

% deconvolve to get F
tmp = Ftau.*E4;
F = sqrt(pi/tau)*tmp(Mr/2-M/2+1:Mr/2+M/2)/Mr;

% now do slow dft, for comparison
dftmtx = exp(-1i*k(:)*xj(:)');
Fex = dftmtx*fj;

% report error
disp(sprintf('Type 1 1D NRMSE: %0.2d',nrmse(Fex,F,ones(size(Fex)))))
disp(sprintf('Type 1 1D MAD: %0.2d',max(abs(Fex-F))));

