% ir_shrink1_test

delta = 10; tmax = 4 * delta; reg = 1.5 * delta;
t = tmax * linspace(-1, 1, 2001)';
zs = 3 * max(delta + reg, delta * (1+reg)) * linspace(-1, 1, 301)';
dt = delta / 20;
tmp = potential_fun('qgg2', delta, 1.2);
table = tmp.dpot([1:1e4] * dt);
param = [dt, table];
pot = potential_fun('table1', delta, param);

xs = pot.shrink(zs, reg);
% plot(zs, zs, '-', zs, xs, ':')

thresh = 0; % not used
nthread = jf('ncore');
nthread = 1
chat = 100;
niter = 10;
s = @(x) single(x);
regs = reg * ones(size(zs));
xc = ir_shrink1_mex(s(zs), s(regs), s(table), s(dt), ...
	int32(niter), s(thresh), int32(nthread), int32(chat));

equivs(xc, xs, 'fail', 0)

plot(zs, zs, '-', zs, xs, ':', zs, xc, '.--')
legend('identity', 'mat', 'mex', 2)
axis([-1 1 -1 1] * max(zs))
axis square
