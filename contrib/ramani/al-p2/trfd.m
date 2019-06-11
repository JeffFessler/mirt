function [trh, trv, tE] = trfd(h, v)
% Compute "time"-reversed horizontal and vertical finite differences on the 
% current estimate x 
%
% Use periodic boundaries

%% Application of "time"-reversed horizontal finite difference
tS = tic;

trh = [h(:,end) h];
trh = imfilter(trh, [1 -1], 'circular');
trh = trh(:, 1:end-1);

%% Application of "time"-reversed vertical finite difference
trv = [v(end,:); v];
trv = imfilter(trv, [1 -1]', 'circular');
trv = trv(1:end-1, :);

tE = toc(tS);
