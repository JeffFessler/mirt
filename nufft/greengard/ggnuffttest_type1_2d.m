% 2-d example of leslie greengard's gaussian kernel-based nufft

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
%k = [0:M-1]';
[kdc1,kdc2] = ndgrid(k);
E4 = exp(tau*(kdc1.^2+kdc2.^2));

% make up irregular grid points
t = [-M/2:M/2-1]'/M*2;
[xj,yj] = ndgrid(pi*sign(t).*t.^2);

% make up data points
fj = rand(M,M);

% find the nearest grid points
m1 = floor(xj/2/pi*Mr);
m2 = floor(yj/2/pi*Mr);
xi1 = 2*pi*m1/Mr;
xi2 = 2*pi*m2/Mr;

% calculate irregular-point-dependent terms
E1 = exp(-((xj-xi1).^2+(yj-xi2).^2)/4/tau);
E2x = exp(pi/Mr/tau*(xj-xi1));
E2y = exp(pi/Mr/tau*(yj-xi2));

disp 'script nufft'
V0 = fj.*E1;
ftau = zeros(Mr,Mr);
for j = 1:M*M
  for l2 = -Msp+1:Msp
    Vy = V0(j)*(E2y(j)^l2)*E3(l2+Msp);
    for l1 = -Msp+1:Msp
      ftau(mod(m1(j)+l1,Mr)+1,mod(m2(j)+l2,Mr)+1) = ftau(mod(m1(j)+l1,Mr)+1,mod(m2(j)+l2,Mr)+1) + ...
          Vy*(E2x(j)^l1)*E3(l1+Msp);
    end
  end
end

% now do fast gridding
ftau2 = back_grid(complexify(V0),Msp,m1,m2,Mr,E2x,E2y,circshift(E3,13),1024);

% take fft of ftau
Ftau = ifftshift(ifftn_fast(ftau));

% deconvolve to get F
F = pi/tau*E4.*Ftau(Mr/2-M/2+1:Mr/2+M/2,Mr/2-M/2+1:Mr/2+M/2);

% now do object, to test
disp 'creating object'
G = gg_Gnufft([xj(:),yj(:)],M,Msp,R);
disp 'object nufft'
Fob = G'*fj(:);

% now do slow dft, for comparison
disp 'slow dft'
k = [-M/2:M/2-1]';
[K1,K2] = ndgrid(k);
dftmtx = exp(1i*(K1(:)*xj(:)'+K2(:)*yj(:)'));
Fex = reshape(dftmtx*fj(:),M,M);

% report error
disp(sprintf('Type 1 2D NRMSE, script: %0.2d',nrmse(Fex,F,ones(size(Fex)))))
disp(sprintf('Type 1 2D MAD, script: %0.2d',max(abs(Fex(:)-F(:)))));
disp(sprintf('Type 1 2D NRMSE, object: %0.2d',nrmse(Fex,Fob,ones(size(Fex)))))
disp(sprintf('Type 1 2D MAD, object: %0.2d',max(abs(Fex(:)-Fob(:)))));
