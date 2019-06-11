% test_pattern1
% test pattern for monitors

nx = 1024;
ny = nx;
pat = zeros(nx, ny, 'uint8');

npat = 2;
for jx=1:npat
	for jy=1:npat
		mx = nx/npat;
		my = ny/npat;
		ix = 1:mx;
		iy = 1:my;
		tx = round(ix / jx);
		ty = round(iy / jy);
		tmp = ((-1).^tx(:)) * ((-1).^ty); 
		tmp = uint8(255 * (1 + tmp)/2); % 0,255
		pat(ix+mx*(jx-1),iy+my*(jy-1)) = tmp;
	end
end
im(pat(1:40,1:40))
im(pat)

%imwrite(pat, 'tmp.gif')
