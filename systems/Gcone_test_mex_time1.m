% Gcone_test_mex_time1.m
% test timing of subset versions of Gcone
% much slower than expected due to matlab mex overhead?

if ~isvar('A')
	cg = ct_geom('ge1', 'na', 642, 'nt', 64);
	ig = image_geom('nx', 128, 'fov', 500, 'zfov', 64*0.625, 'nz', 64);
	A = Gcone(cg, ig, 'type', 'sf2', 'is_ns_nt', 1);

end
	x = single(ig.circ);
	x = x(ig.mask);

if 1
%	profile on
	cpu etic
	y1 = A * x;
	cpu etoc A
%	profile report
end
if 0
	cpu etic
		y1 = A * x;
	cpu etoc A2
end

%nblock = 6;
%nblock = 10;
%nblock = 107;
nblock = floor(642/8 / 2)
B = Gblock(A, nblock);

%profile on
cpu etic
for ib=1:nblock
	y2 = B{ib} * x;
end
cpu etoc B
%profile report
