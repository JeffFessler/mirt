 function y = proj_rotmex_view0(x, G0, na, rotate_key)
%function y = proj_rotmex_view0(x, G0, na, rotate_key)
% fast (?) projector based on rotation
% and view 0 sparse matrix
% rotate_key is from previous call to rotmex('init', ...)
% orbit, orbit_start, filter)

Gt = G0';
nb = size(G0, 1);
mag = sqrt(size(G0,2)) / size(x,1);

% magnify image
if mag ~= 1
	x = imresize(x, mag, 'nearest');
end

if 0
	orbit = 180;
	orbit_start = 0;
	filter = 3;
	chat = 0;
rotate_key = rotmex('init', int32(size(x,1)), int32(size(x,2)), int32(na), ...
	orbit, orbit_start, int32(filter), int32(chat))
end

y = zeros(nb, na);
ticker reset
for ia=1:na
	ticker(mfilename, ia, na)
	tmp = rotmex('spin', single(x), int32(ia-1), rotate_key);
	tmp = double(tmp(:))';
%	tmp = tmp(:);
%	tmp = tmp';
	tmp = tmp * Gt;
	tmp = tmp';
	y(:,ia) = tmp;
end

% rotmex('free', rotate_key)
