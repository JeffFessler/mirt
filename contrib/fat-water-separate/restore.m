function [w f] = restore(data,map,F_fw,te)

% data = source image
% map = estimated field map
% F_fw = water fat shift (omega)
% te = time

n = size(data,1);
nsets = size(data,3);

A = zeros(nsets,2);
A(:,1) = ones(nsets,1);
A(:,2) = exp(j*F_fw*te');

M = inv(A'*A)*A';
w = zeros(n,n);
f = zeros(n,n);
E_data = zeros(n,n,nsets);

for k = 1:nsets
    E_data(:,:,k) = data(:,:,k).*exp(-j*map*te(k));
    w = w + M(1,k)*reshape(E_data(:,:,k),n,n);
    f = f + M(2,k)*reshape(E_data(:,:,k),n,n);
end