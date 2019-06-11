function [DISx, DISy] = calc_displacement(A, W, refx, refy, gridrefx, gridrefy)
% calculate displacemnets in x and y
[size_n, size_m] = size(gridrefx);
n = length(refx);
U = zeros(n,1);
gridholx = zeros(size_n, size_m);
gridholy = zeros(size_n, size_m);
DISx=0;DISy=0;
for i=1:size_n,
    for j=1:size_m,
        x = gridrefx(i,j); y = gridrefy(i,j);
        [gridholx(i,j) gridholy(i,j)] = map_points(A, W, refx, refy, x, y);
        DISx = DISx + abs(gridholx(i,j) - x);
        DISy = DISy + abs(gridholy(i,j) - y);
    end;
end;
