% test_2z_block.m
% todo: replace dd with wtmex or Gtomo2_table, or several of them

if ~isvar('G2'), printm 'G2'
	f.down = 8;
	ig = image_geom('nx', 512, 'fov', 500, 'down', f.down);
	ig.mask = ig.circ > 0;
	sg = sino_geom('ge1', 'down', f.down);

	G2 = Gtomo_dd(sg, ig);
end

if ~isvar('b3')
	x2 = ig.mask;
	y2 = G2 * x2; % single 2d slice baseline
	b2 = G2' * y2;

	nz = 4;
	x3 = repmat(x2, [1 1 nz]);
	mask3 = repmat(ig.mask, [1 1 nz]);

	y3 = repmat(y2, [1 1 nz]);
	b3 = repmat(b2, [1 1 nz]);
	im(b2)
prompt
end

% block_fatrix
if ~isvar('bz'), printm 'bz'
	Gz = block_fatrix({G2}, 'type', 'kron', 'Mkron', nz);

	yz = Gz * x3(mask3);
	yz = reshape(yz, [sg.dim nz]); % ideal forw
	jf_equal(yz, y3);

	bz = Gz' * y3(:); % ideal back
	bz = reshape(bz, [ig.np nz]);
	bz = ig.embed(bz);
	jf_equal(bz, b3);
prompt
end

% Gblock
if ~isvar('bb'), printm 'bb'
	nblock = 3;
	Gb = Gblock(Gz, nblock);

	yb = Gb * x3(mask3);
	yb = reshape(yb, [sg.dim nz]); % ideal
	jf_equal(yb, y3)

	bb = Gb' * y3(:);
	bb = reshape(bb, [ig.np nz]);
	bb = ig.embed(bb);
	jf_equal(bb, b3);
prompt
end

if 1
	tmp = Gb{1};
	y1 = tmp * x3(mask3);
	y1 = reshape(y1, [sg.nb sg.na/nblock nz]); % ideal
	ia = 1:nblock:sg.na;
	jf_equal(y1, yz(:,ia,:))

	y2s = zeros(size(y2));
	y2s(:,ia) = y2(:,ia);
	b2s = G2' * y2s; % reference backproj of subset
end

if 1 % todo: check backproj
%	y3s = zeros(size(y3));
%	y3s(:,ia,:) = y3(:,ia,:);
	% todo: not implemented properly yet!
	y3s = y3(:,ia,:);
	b3s = tmp' * y3s(:);
	b3s = reshape(b3s, [ig.np nz]);
	b3s = ig.embed(b3s);
printm 'ok'
return
end
