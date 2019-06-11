% Robject_test.m
% compare new Robject to old Rbuild

ig = image_geom('nx', 512, 'ny', 500, 'fov', 500);
ig.mask = ig.circ(200) > 0;
nx = ig.nx;
ny = ig.ny;
mask = ig.mask;
if 0
	nx = 40;
	ny = 38;
	mask = ellipse_im(ig, [0 0 14 17 0 1], 'oversample', 3) > 0;
end

f.l2b = 3;
f.nbrs = 8;
f.delta = 0.5;
f.pot = 'hyper3';
f.pot = 'huber';

if 0 % check denom
	R2 = Robject(mask, 'order', 2, 'beta', 1, 'type_denom', 'matlab', ...
		'potential', 'huber', 'delta', 1e4);
	den = R2.denom(R2, ig.zeros(ig.mask));
	den = embed(den, mask);
	pred = 4 * 4 * (2*1 + 2*1/sqrt(2));
%	im(pred.*mask - den)
	max_percent_diff pred max(den(:))
return
end

if 0 % compare to C2sparse
	rj = 10;
	Rq = Robject(mask, 'edge_type', 'leak', 'beta', 2^f.l2b, ...
		'offsets', [1 nx nx+1 nx-1], 'user_wt', rj, ...
		'type_denom', 'aspire');
	[Cs wjk] = C2sparse('leak', mask, 8, 0, 1);
	Cs = Cs(:,mask);
	Cs = sqrt(2^f.l2b) * spdiag(sqrt(rj .* wjk), 'nowarn') * Cs;
	Cq = Rq.C;

	x = mask;
	ys = Cs * x(mask);
	ys = reshape(ys, [nx ny 4]);
	yq = Cq * x(mask);
	yq = reshape(yq, [nx ny 4]);
	im clf
	im(221, ys, 'ys'), cbar
	im(222, yq, 'yq'), cbar
	im(223, yq-ys, 'yq-ys'), cbar
return
end

if ~isvar('R1')
	R1 = Robject(mask, 'edge_type', 'leak', 'beta', 2^f.l2b, ...
		'potential', f.pot, 'delta', f.delta, 'type_denom', 'aspire');
end

if 0 % test: see how much time it saves to have a concrete value in there
	R1.pot.potk = @(dum,t) 0.5^2 * (sqrt(1 + abs(t / 0.5).^2) - 1);
	R1.pot.wpot = @(dum,t) 1 ./ sqrt(1 + abs(t / 0.5).^2);
	R1.pot.dpot = @(dum,t) t ./ sqrt(1 + abs(t / 0.5).^2);
end

if 0 % test: time with implicits
	R1.pot.potk = @(dum,t) t.^2/2;
	R1.pot.wpot = @(dum,t) ones(size(t));
	R1.pot.dpot = @(dum,t) t;
end

if 0 % timing
	x = single(ig.mask);
%	[tmp1 tmp2] = R1.dercurv(R1, R1.C1*x);
	Cx = R1.C1 * x(ig.mask);

	t = Cx;
	w = Cx;
	n = 20;

	pot_time_test(t);

%	profile on
	cpu etic
	for ii=1:n
		[tmp1 tmp2] = feval(R1.dercurv, R1, Cx);
	end
	cpu etoc 'R1.dercurv'
%	profile report

return
end

if 0 && ~isvar('R0') % obsolete
%	R0 = Rbuild('leak', mask, f.nbrs, 2^f.l2b, f.pot, f.delta, 1);

	x = mask;
	rng(0)
	x = rand(size(mask));
	x = x(mask);
	max_percent_diff(R0.penal(R0, x), R1.penal(R1, x))
	max_percent_diff(R0.diag(R0), R1.diag(R1))
	max_percent_diff(R0.cgrad(R0, x), R1.cgrad(R1, x))
	max_percent_diff(R0.denom(R0, x), R1.denom(R1, x))
	if 0
		d0 = R0.denom(R0, x);
		d1 = R1.denom(R1, x);
		d0 = embed(d0, mask);
		d1 = embed(d1, mask);
		im plc 2 2
		im(1, d0)
		im(2, d1)
		im(3, d1-d0)
	end
end

if 1 % test 3D
	nx = 9; ny = 7; nz = 11;
	mask = ones(nx, ny, nz);
	t = [1 nx nx+1 nx-1 nx*ny nx*ny-1 nx*ny+1 nx*ny+nx nx*ny-nx];
	t = '3d:26';
	R = Robject(mask, 'offsets', t);
	ej = 0 * mask; ej((end+1)/2, (end+1)/2, (end+1)/2) = 1;
	t = R.C1' * (R.C1 * ej);
	t = fftshift(fftn(ifftshift(t)));
	t = reale(t);
	im(t), cbar
return
end
