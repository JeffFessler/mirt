function [m] = region_sub(map,A,afterest,arr,col)
% this function is estimate initial guess for pixel independet method
% at field map(arr,col) by region growing method
% in 
%   map         [n,n]   initial guess field map
%   A           [3*n,3] matrix of [1 arr col]
%   afterest    [n,n]   afterest(arr,col)=1 means we estimate
%                       field map(arr,col) by region growing method
%   arr,col             the position of the pixel now we try to estimate

%% initialize
window = zeros(512,512);
% make 41*41 window whose center is (arr,col)th pixel for region growing 
% method
sa = arr-20;    % start array
ea = arr+20;    % end array
sc = col-20;    % start col
ec = col+20;    % end col

if sa<1
    sa = 1;
    ea = 41;
elseif ea>512
    ea = 512;
    sa = 472;
end
if sc<1
    sc = 1;
    ec = 41;
elseif ec>512
    ec = 512;
    sc = 472;
end

window(sa:ea,sc:ec) = 1;
maskgre = afterest.*window;
% we use only pixels which is already estimated and near(41*41) from
% (arr,col)th pixel

choose = find(maskgre==1);
map2 = map(choose);
w= map2.^2; % wieght

B = zeros(size(choose,1),3);

for i = 1:3
    B(:,i) = A(choose,i);
    C(:,i) = B(:,i).*w;
end



theta = inv(C'*B)*C'*map2; % the result of MSE 

m = theta(1) + theta(2)*arr+theta(3)*col; 
% m is new initla guess value for fieldmap(arr,col)






 




