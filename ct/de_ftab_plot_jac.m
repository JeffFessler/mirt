  function dummy = de_ftab_plot_jac(ftab)
%|function dummy = de_ftab_plot_jac(ftab)
%| show Jacobian transformation

dummy = [];

sl = ftab.sls.sl;
sll = ndgrid_jf('mat', sl{:});
ftmp = ftab.fit.fmfun(sll);

zz = de_ftab_xform(ftab, ftmp);
z1 = zz(:,:,1);
z2 = zz(:,:,2);

if 1 % array of dots showing effect of transformation
	clf, plot(z1, z2, '.')
%	axis([0 max(sl{1}) 0 max(sl{2})])
	axis tight
	xlabel 's_1 [cm^2/g]', ylabel 's_2 [cm^2/g]'
	xlabel 'f^*_1(s_1,s_2)', ylabel 'f^*_2(s_1,s_2)'
	title 'Effect of Linearizing Transformation'
%	ir_savefig('c', fdir, 'fig_f_linearized')
prompt
end

% nice figure showing inverse and departure from linearity thereof
if 1 && usejava('jvm')
	i1 = 1:4:length(sl{1});
	i2 = 1:4:length(sl{2});
	z1 = z1(i1, i2);
	z2 = z2(i1, i2);
	s1 = sll(i1, i2, 1);
	s2 = sll(i1, i2, 2);

	m1 = ftab.sls.max(1);
	m2 = ftab.sls.max(2);

	clf, pl = @(n) subplot(220 + n);
%	plot3(z1, z2, s1-0*z1, '.')
	pl(1), mesh(z1, z2, s1-0*z1), title 's_1'
	axis([0 m1 0 m2 0 m1])
	xtick([]), ytick([]), zwhite
%	xlabel f^*_1, ylabel f^*_2

	pl(2), mesh(z1, z2, s2-0*z2), title 's_2'
	axis([0 m1 0 m2 0 m2])
	xtick([]), ytick([]), ztick, zwhite
	% xlabel f^*_1, ylabel f^*_2

	pl(3), mesh(z1, z2, s1-1*z1), title 's_1 - f^*_1'
	axis([0 m1 0 m2 -3 12])
	xtick, ytick, ztick, zwhite
	% xlabel f^*_1, ylabel f^*_2

	pl(4), mesh(z1, z2, s2-1*z2), title 's_2 - f^*_2'
	axis([0 m1 0 m2 -3 12])
	xtick, ytick, ztick, zwhite, xlabel f^*_1, ylabel f^*_2

%	ir_savefig('c', fdir, 'fig_fs_inverse')
prompt
end
