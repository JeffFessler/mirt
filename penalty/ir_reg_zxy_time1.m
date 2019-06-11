% ir_reg_zxy_time1.m
% timing of new zxy regularization

if ~has_mex_jf, warn 'no mex', return, end

if ~isvar('times')

f.down = 1;
nthread = jf('ncore');
%nthread_list = nthread:-1:(nthread-2);
%nthread_list = [nthread nthread/2 1];
%nthread_list = [2*nthread 1.5*nthread nthread:-1:1];
%nthread_list = [(1.5*nthread):-1:1];
nthread_list = [24 16 12 8 4 2 1];
nthread_list(nthread_list > 1.5*nthread) = [];
nlist = numel(nthread_list);

pot_list = {'qgg2', 'gf1', 'quad'}; % 'hyper3';
npot = numel(pot_list);

times = zeros(nlist, npot, 2);

zxy = @(x) permute(x, [3 1 2]);
xyz = @(x) permute(x, [2 3 1]);

f.offsets = '3d:26'; % works for both xyz and zxy orderings
f.l2b = 3;
f.order = 1;
f.delta = 1;

ig = image_geom('nx', 512+16, 'ny', 512-16, 'nz', 768, 'fov', 500, 'down', f.down);

%tmp = ig.mask;
tmp = ig.circ > 0; % circ mask
tmp([1 end],:,:) = 0; tmp(:, [1 end],:) = 0; % zero x,y border for 3d
ig.mask = tmp; clear tmp
if 0 % save mask for hct2 timing
	tmp = sprintf('mask-%d-%d.fld', ig.nx, ig.ny)
	fld_write(tmp, ig.mask(:,:,1), 'type', 'uint8')
	im(ig.mask(:,:,1))
return
end

kappa = single(ig.mask);

for ip=1:npot
	f.pot = pot_list{ip};
	switch f.pot
	case 'quad'
		f.pot_arg = {f.pot};
	case 'qgg2'
		f.pot_arg = {f.pot, f.delta, 1.2};
	case 'hyper3'
		f.pot_arg = {f.pot, f.delta};
	case 'gf1'
		f.pot_arg = {f.pot, [f.delta, 2, 3]}; % fake
	otherwise
		fail 'bug'
	end

 for it=1:nlist
	nthread = nthread_list(it);

	printm('Rx,Rz nthread=%d', nthread)

	f.arg = {'beta', 2^f.l2b, 'pot_arg', f.pot_arg, ...
		'edge_type', 'tight', 'order', f.order, ...
		'offsets', f.offsets, 'nthread', nthread};
	Rx = Reg1(kappa, f.arg{:}, 'type_penal', 'mex', 'control', 2);
	Rz = Reg1(zxy(kappa), f.arg{:}, 'type_penal', 'zxy');

	rng(0)
	x = rand(ig.dim) .* ig.mask;
%	x = ig.unitv;

	printm 'time cgrad / denom'
	tmp = zxy(x);
	cpu etic
	[g5 d5] = feval(Rz.cgrad_denom, Rz, tmp);
	times(it,ip,1) = cpu('etoc');
	g5 = xyz(g5);
	d5 = xyz(d5);
	printm('Rz cgrad/denom %d time: %g', nthread, times(it,ip,1))

	if 1
		cpu etic
		g3 = Rx.cgrad(Rx, x);
		times(it,ip,2) = cpu('etoc');
		equivs(g3, g5)
		printm('Rx cgrad %d time: %g', nthread, times(it,ip,2))
	end

	if 0
		cpu etic
		d3 = Rx.denom(Rx, x);
		cpu etoc 'Rx denom time'

		equivs(d3, d5)
	end
 end % it
end % ip

pr times

end % if

clf, pl=@(i) subplot(220+i);
for ip=1:npot
	pl(ip)
	it = nthread_list;
	i1 = find(it == 1); time1 = times(i1,ip,1);
	plot(	it, it, ':', ...
		it, time1 ./ times(:,ip,1), '-o', ...
		it, time1 ./ times(:,ip,2), '-x')
	xlabel 'nthread', ylabel 'speedup '
	legend('ideal', 'zxy', 'xyz', 2)
	axis([1 max(it) 1 1+jf('ncore')]), grid
	titlef('speedup for (%d,%d,%d) %s', ig.dim, pot_list{ip})
end
pl(4)
ic = find(it == jf('ncore'))
bar(squeeze(times(ic,:,:)))
tmp = strvcat(pot_list{:});
set(gca, 'xticklabel', tmp)
titlef('wall time for nthread=%d', nthread_list(ic))
legend('zxy', 'xyz')
