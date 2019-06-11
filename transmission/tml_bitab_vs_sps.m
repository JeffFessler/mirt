% tml_bitab_vs_sps.m
% compare TML-BITAB and T-ML-OS-SPS
% Copyright Apr 2000, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('yi')
	f.randpercent = 0; % because BITAB "requires" ri=0
	tr_test_setup
	Gt = G';
prompt
end

f.niter = 8;

xinit = 0.01+0*xfbp;
f.pixmax = 0.012;

%
% loop over 1 subset case (monotone) and multi-subset case (nonmonotone)
%
im plc 4 2, pl=420;
subset_list = [1 4];
for is=1:length(subset_list)
	f.nblock = subset_list(is);

	Gb = Gblock(G, f.nblock);

	%
	% os-sps
	%
	if 1
		if f.nblock == 1, f.curv = 'oc'; else f.curv = 'pc'; end

		xosps = tpl_os_sps(xinit(ig.mask), Gb, yi, bi, ri, ...
			[], f.niter, f.pixmax, f.curv);
		xosps = ig.embed(xosps);

		im(is, xosps, sprintf('T-ML-OSPS-%d', f.nblock))

		t = tpl_obj(xosps, G, yi(:), bi(:), ri(:), [], ig.mask);
		if im
			subplot(pl+2+is), plot(0:f.niter-1, t, '-o')
		end
	end

	%
	% bitab
	%
	if 1
		xbitab = tml_bitab(xinit(ig.mask), Gb, yi, bi, ri, ...
			f.niter, f.pixmax, 1e5);
		xbitab = ig.embed(xbitab);

		im(4+is, xbitab, sprintf('T-ML-BITAB-%d', f.nblock))

		t = tpl_obj(xbitab, G, yi(:), bi(:), ri(:), [], ig.mask);
		if im
			subplot(pl+6+is), plot(0:f.niter-1, t, '-o')
		end
	end
end
