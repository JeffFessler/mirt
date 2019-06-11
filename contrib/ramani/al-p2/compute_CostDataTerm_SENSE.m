function [DT, tE] = compute_CostDataTerm_SENSE(Data, SP3, x, params)
%% Compute the data part of the cost function
% 1/2 * (y-Bx)^H * W * (y-Bx); B = S * F

smap = params.smap;
ncoils = params.ncoils;

%% Compute the Data Term
tS = tic;
Sx = smap .* repmat(x, [1, 1, ncoils]);
FSx = fft2(Sx) .* SP3;
dif = Data - FSx;

DT = sum(abs(dif(:)).^2) / 2;
tE = toc(tS);
