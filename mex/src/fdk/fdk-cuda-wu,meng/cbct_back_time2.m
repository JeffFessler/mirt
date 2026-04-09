 % cbct_back_time
% compare run times of various versions of fdk cone-beam CT back-projector

% do down=8 first as a "warm up" so that subsequent timings are accurate
dlist = [8 4 2 1]; % use this later
dlist = [8]; % for now just quick tests
for down = dlist % various levels of down-sampling, fastest first

	cg = ct_geom('ge1', 'nt', 800, 'down', down);
	ig = image_geom('nx', 512, 'ny', 496, 'nz', 480, 'fov', 500, ...
		'mask', 'all-but-edge', 'down', down);

	pr [ig.nx ig.ny ig.nz ig.dx ig.dy ig.dz]
	pr [cg.ns cg.nt cg.na cg.ds cg.dt]
	printm('image Mbytes %d', ig.nx*ig.ny*ig.nz*4 / 1024^2)
	printm('proj Mbytes %d', cg.ns*cg.nt*cg.na*4 / 1024^2)

	printm 'ellipse proj' % somewhat realistic phantom object
	ell = [ ...
		[30 10 10	150 150 280	0 0 1000]; % 30cm diam
		[80 10 10	50 50 30	0 0 300]; % bone-like
		[-10 -40 75	40 40 40	0 0 600];
		[-10 80 -20	30 30 30	0 0 600];
	];
%	xtrue = ellipsoid_im(ig, ell); im(xtrue); return

	proj = ellipsoid_proj(cg, ell);
	proj = fdk_filter(proj, 'ramp', cg.dsd, cg.dfs, cg.ds);
	if 0 % zero outer edges of projection for testing
		proj([1 cg.ns], :, :) = 0;
		proj(:, [1 cg.nt], :) = 0;
	end
%	im(proj, 'true projections'), cbar, return
if exist('fdk_mex') == 3 %, printm 'found fdk_mex'
		back2 = cbct_back(args{:}, 'use_mex', 1, 'back_call', @fdk_mex);
end
	args = {proj, cg, ig, 'ia_skip', 1}; % increase 1 for faster debugging
	if 0 && down >= 4 % test nthread scaling
		nthreads = 2 .^ fliplr(0:log2(jf('ncore')));
		for in = 1:length(nthreads)
			nthread = nthreads(in);
			pr nthread
			cpu etic
			tmp{in} = cbct_back(args{:}, 'use_mex', 1, ...
				'back_call', @jf_mex, 'nthread', nthread);
			times(in) = cpu('etoc');
			printm('time %d = %g', nthread, times(in))
		end
		pr times
		pr times(end) ./ times % speedup
	continue
	return
    end
    ttest = 1;
    for mengtest = 1:ttest
	cpu etic
	back1 = cbct_back(args{:}, 'use_mex', 1, 'back_call', @jf_mex);
	time1(mengtest) = cpu('etoc');

	if exist('fdk_mex') == 3 %, printm 'found fdk_mex'
		cpu etic
		back2 = cbct_back(args{:}, 'use_mex', 1, 'back_call', @fdk_mex);
		time2(mengtest) = cpu('etoc');

		printm('times jf_mex=%g fdk_mex=%g', time1(mengtest), time2(mengtest))
     %   max_percent_diff(back1, back2);
    end
    end
    max_percent_diff(back1, back2);
    printm('downsampling rate=%g',dlist );
    printm('cpu_time=%g',sum(time1)/ttest);
    printm('gpu_time=%g',sum(time2)/ttest);
   printm('ratio=%g', sum(time1)/sum(time2));
    

	if 0
		im plc 2 2
		im(1, 'mid3', back{1})
		im(2, 'mid3', back{2})
		err = back{2} - back{1};
		im(3, 'mid3', err), ylabel 'err mid'
		im(4, 'mip3', abs(err)), ylabel '|err| mip'
	keyboard
	end

end % down
