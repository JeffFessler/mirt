% dp_test_yzz
% test distance_power bug (?) for yzz

if ~isvar('R2')
	ig = image_geom('nx', 22, 'ny', 20, 'nz', 8, 'dx', 2, 'dz', 1);
	ej = ig.unitv;
	l2b = 0;
	offsets = ''; % default is x,y,z
	offsets = '3d:26'; % full 3d neighborhood
	offsets = [ig.nx+1]; % test a diagonal
	arg = {ig.mask, 'edge_type', 'tight', 'beta', 2^l2b, ...
		'type_denom', 'matlab', ...
		'offsets', offsets, 'distance_power'};
%	Rgen = @(dp) Robject(arg{:}, dp);
	Rgen = @(dp) Reg1(arg{:}, dp);
	R0 = Rgen(0);
	R1 = Rgen(1);
	R2 = Rgen(2);
end

if 1
	pr R0.penal(R0, ej(ig.mask))
	pr R1.penal(R1, ej(ig.mask))
	pr R2.penal(R2, ej(ig.mask))
end

if 0
	w0 = reshape(R0.wt, [ig.dim length(R0.offsets)]);
	w1 = reshape(R1.wt, [ig.dim length(R0.offsets)]);
	w2 = reshape(R2.wt, [ig.dim length(R0.offsets)]);
	im(w2), cbar
return
end

%tmp = R0.C * ej(ig.mask);
%return

psf0 = ig.embed(R0.C' * (R0.C * ej(ig.mask)));
psf1 = ig.embed(R1.C' * (R1.C * ej(ig.mask)));
psf2 = ig.embed(R2.C' * (R2.C * ej(ig.mask)));
minmax(psf0)
minmax(psf1)
minmax(psf2)

try0 = ig.embed(R0.cgrad(R0, 1e-2 * ej(R0.mask(:))) / 1e-2);
try1 = ig.embed(R1.cgrad(R1, 1e-2 * ej(R1.mask(:))) / 1e-2);
try2 = ig.embed(R2.cgrad(R2, 1e-2 * ej(R2.mask(:))) / 1e-2);

im pl 2 3
im(1, psf0)
im(2, psf1)
im(3, psf2)
im(4, try0)
im(5, try1)
im(6, try2)
