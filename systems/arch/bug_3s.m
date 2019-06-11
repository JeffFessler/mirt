% bug_3s.m
% explore an imperfect adjoint that used to exist in 3s

if 0 % rotmex bug, fixed by AY on 2006-3-8:
	im(rotmex(ones(200, 'single'), 40.97561, int32(3)))
end

if ~isvar('G')
	nx = 12;
	ny = 12;
	nz = 8;
	na = 6;

%	mask = []; % for 3s object!
	f.sys_type = sprintf('3s@1,1,1,360,0,%d%s@%s@%s@-%d,%d,%d', ...
			1, ',none', '-', '-', nx, nz, na)

	G = Gtomo3(f.sys_type, [], nx, ny, nz);

	ig = image_geom('nx', nx, 'ny', ny, 'nz', nz, 'dx', 1);
	ig.mask = G.arg.mask;
prompt
end

if ~isvar('A')
	[A Aa] = test_adjoint(G);
	im(131, A), im(132, Aa'), im(133, Aa'-A)
end

jj = find(sum(abs(Aa'-A)) > 0.1);
if isempty(jj), printm 'no bug now!', return, end
jj = jj(4)

a1 = reshape(A(:,jj), nx, nz, na, []);
a2 = reshape(Aa(jj,:)', nx, nz, na, []);

im pl 3 1
im(1, a1, 'A ej'), cbar
im(2, a2, 'A'''' ej'), cbar
im(3, a2-a1, 'diff'), cbar
