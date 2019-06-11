function trv = trfdv(v)
% Compute "time"-reversed vertical finite differences on the 
% current estimate x 
%
% Use periodic boundaries

%% Application of "time"-reversed vertical finite difference
trv = [v(end,:); v];
trv = imfilter(trv, [1 -1]', 'circular');
trv = trv(1:end-1, :);