 function y = proj_rot_view0(x, G0, na, filter)
%function y = proj_rot_view0(x, G0, na, filter)
% fast (?) projector based on rotation
% and view 0 sparse matrix

Gt = G0';
nb = size(G0, 1);
mag = sqrt(size(G0,2)) / size(x,1);

% magnify image
if mag ~= 1
	x = imresize(x, mag, 'nearest');
end

y = zeros(nb, na);
for ia=1:na
	angle = (ia-1)/na * 180;
	tmp = rotmex(single(x), angle, int32(filter));
	tmp = double(tmp(:))';
%	tmp = tmp(:);
%	tmp = tmp';
	tmp = tmp * Gt;
	tmp = tmp';
	y(:,ia) = tmp;
end
