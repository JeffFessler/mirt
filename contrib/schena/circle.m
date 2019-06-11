function [circle_mask,nnz_crcl] = circle(N);

% N is the circle diameter;
centro=ceil(N/2);

xax = (1:N);
x = repmat(xax,N,1) ; 
% x coordinates, the y coordinates are rot90(x)
y = rot90(x);

raggio = centro;

circle_mask=( ((x-centro).^2+(y-centro).^2)<=raggio^2 );

figure, imshow(circle_mask,[ ]); title(' simple mask')
nnz_crcl=nnz(circle_mask);
