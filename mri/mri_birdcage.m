% mri_birdcage
% examine response of an idealized birdcage coil
% using biot-savart law for an infinitely long current carrying wire
% caution: written by a signal processor so the emag could be all wrong...

ig = image_geom('nx', 20, 'fov', 20);

nc = 4; % # of legs
if 0 % image sides
    rc = 14;
    xlist = rc * cos(2*pi*(0:nc-1)/nc); % position of legs
    ylist = rc * sin(2*pi*(0:nc-1)/nc);
else % image corners
    rc = ig.fov / sqrt(2) + 2*ig.dx;
    xlist = rc * cos(2*pi*(0:nc-1)/nc + pi/4); % position of legs
    ylist = rc * sin(2*pi*(0:nc-1)/nc + pi/4);
end
plist = -(0:nc-1)/nc * 2*pi; % quadrature phase for each leg

bcoil = zeros(ig.nx, ig.ny, nc);
for ic=1:nc
	x1 = xlist(ic);
	y1 = ylist(ic);
	xo = ig.xg - x1;
	yo = ig.yg - y1;
	rr = sqrt(xo.^2 + yo.^2); % distance of leg to grid points
	phi = atan2(xo, -yo) + plist(ic); % phase
%	phi = atan2(yo, xo) + plist(ic); % phase - which is correct?

	bcoil(:,:,ic) = (1 ./ rr) .* exp(1i * phi); % 1/r falloff and phase
end

bsum = sum(bcoil, 3);
bsum = bsum / max(abs(bsum(:))); % normalize to unit magnitude for display

im plc 2 2
im('row', 2, abs(bcoil), 'Coil magnitudes'), cbar
im(abs(bsum), [0 1], 'Magnitude of sum'), cbar
im('row', 2, angle(bcoil), 'Coil phases'), cbar
colormap(gca, hsv)

im subplot 4
quiver(ig.xg, ig.yg, real(bsum), imag(bsum))
titlef('Phase of sum')
axis([-1 1 -1 1]*12)
axis square

if 0
	t = linspace(0, 2*pi, 25);
	r = 4;
	x = r * cos(t);
	y = r * sin(t);

%	phi = pi/2 + atan2(y,x);
	phi = atan2(x,-y);
	u = cos(phi);
	v = sin(phi);

	quiver(x,y,u,v)
	axis square
end
