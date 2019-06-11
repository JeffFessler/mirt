 function x = back_rot_view0(y, G0, na)
%function x = back_rot_view0(y, G0, na)
%	fast (?) backprojector based on rotation
%	and view 0 sparse matrix

Gt = G0';
nb = size(G0, 1);
mag = sqrt(size(G0,2)) / size(x,1);

%	magnify image
xs = single(imresize(x, mag, 'nearest'));

y = zeros(nb, na);
for ia=1:na
	angle = (ia-1)/na * 180;
	tmp = rotmex(xs, angle, 1);
	tmp = double(tmp(:))';
%	tmp = tmp(:);
%	tmp = tmp';
	tmp = tmp * Gt;
	tmp = tmp';
	y(:,ia) = tmp;
end
