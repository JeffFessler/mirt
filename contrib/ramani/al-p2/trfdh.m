function trh = trfdh(h)
% Compute "time"-reversed horizontal finite differences on the 
% current estimate x 
%
% Use periodic boundaries

%% Application of "time"-reversed horizontal finite difference
trh = [h(:,end) h];
trh = imfilter(trh, [1 -1], 'circular');
trh = trh(:, 1:end-1);
