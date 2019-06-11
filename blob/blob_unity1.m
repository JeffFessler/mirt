% blob_unity1.m
% find blob parameters that make the best "partition of unity"
% jeff fessler

% 1d
% J m alpha	alpha/J	(max/min-1)%
% 2 2 2.502	1.251	0.0192		1 min
% 3 2 2.472	0.824	0.0021		2 min
% 3 2 7.464	2.488	0.0587776%	"
% 4 2 2.461	0.61525 0.000315658%	3 min
% 4 2 8.687	2.17175	0.0102284%	"
% 4 2 11.17	2.7925	0.00177428%	"

% 2d
% J m alpha	alpha/J	(max/min-1)%
% 4 2 7.888	1.972	0.0411645
% 4 2 10.83	2.7075	0.0180456%

% 3d
% J m alpha	alpha/J	(max/min-1)%
% 4 2 10.4	2.6	0.0290578	and mae and rmse < 0.01%.
% for CT this is a fraction of a HU!

if ~isvar('J')
	J = 4; % a = J/2

	alf_list = linspace(0.5*J, 3.0*J, 26+0*101);
	alf_list = linspace(7, 12, 5*10+1); % 2d
	alf_list = linspace(10.5, 11, 5*10+1); % 2d
	alf_list = linspace(2, 12, 10001); % 1d
	alf_list = linspace(7.0, 8.0, 1001); % 1d fine, J=3
alf_list = linspace(0.1, 20, 1001); % 1d coarse
	alf_list = linspace(8.0, 9.0, 1001); % 1d fine, J=4
%	alf_list = linspace(10.3, 10.5, 3); % 3d
%	alf_list = 10.4; % 3d
%	m_list = linspace(1.9, 2.1, 21);
	m_list = 2.0;
	%alf_list = 2.34 * J;
	%m_list = 2;
	[aa mm] = ndgrid(alf_list, m_list);

	x1 = linspace(0, J/2, 101)';
	x2 = linspace(0, J/2, 101)';
	x3 = linspace(0, J/2, 101)';
	x2 = 0; % 1d
	x3 = 0; % 2d
	[xx1 xx2 xx3] = ndgrid(x1, x2, x3);

	% get all j that affect x in [0,J/2]
	j1max = ceil(J-1);
	j1min = floor(-J/2+1);
	ndim = 1 + (length(x2) > 1);
	nj1 = j1max-j1min+1;

	if ndim >= 3
		nj3 = nj1;
		j3min = j1min;
		j3max = j1max;
	else
		nj3 = 1;
		j3min = 0;
		j3max = 0;
	end

	if ndim >= 2
		nj2 = nj1;
		j2min = j1min;
		j2max = j1max;
	else
		nj2 = 1;
		j2min = 0;
		j2max = 0;
	end
end


%
% precompute r samples for each j offset
%
if ~isvar('rr'), disp 'do rr'
	rr = zeros([length(x1) length(x2) length(x3) nj1 nj2 nj3]);
	for j3=j3min:j3max
		for j2=j2min:j2max
			for j1=j1min:j1max
				rr(:,:,:,j1-j1min+1,j2-j2min+1,j3-j3min+1) = ...
				sqrt((xx1-j1).^2 + (xx2-j2).^2 + (xx3-j3).^2);
			end
		end
	end
	%im(x1, x2, rr)
end

if ~isvar('bad'), disp 'do bad'
	bad = zeros(size(mm));
	bb = zeros(length(x1), length(x2), length(x3), nj1*nj2*nj3);

	for ii=1:numel(aa)
		ticker(mfilename, ii, numel(aa))
		kb_a = aa(ii);
		kb_m = mm(ii);

		bb = kaiser_bessel(rr, J, kb_a, kb_m);
		bsum = sum(sum(sum(bb, 3), 4), 5);
		bad(ii) = max(bsum(:)) / min(bsum(:));
	end
	bad = abs(bad-1);
end

if 1
	im clf, im(121, alf_list/J, m_list, bad), cbar
	axis normal
	hold on
	ibest = imin(bad, 2);
	plot(alf_list(ibest(1))/J, m_list(ibest(2)), '*')
	hold off
	xlabel '\alpha/J', ylabel 'm'
end

if ~isvar('bmean')
	kb_a = alf_list(ibest(1));
	kb_m = m_list(ibest(2));

	bb = kaiser_bessel(rr, J, kb_a, kb_m);
	bsum = sum(sum(sum(bb, 3), 4), 5);
	bmean = mean(bsum(:));
end

printf('J m alpha alpha/J ratio-1')
printf('%g %g %g %g %g%%', ...
	J, kb_m, kb_a, kb_a / J, (max(bsum(:)) / min(bsum(:))-1)*100)
printf('mea=%g%%', mean(abs(bsum(:)-bmean))/bmean * 100)
printf('rmse=%g%%', sqrt(mean(abs(bsum(:)-bmean).^2))/bmean * 100)

if 1
	subplot(122)
	if ndim == 1
		plot(x1, reshape(bb, length(x1), []), '--', x1, bsum, '-')
		clf, semilogy(alf_list, bad, '.-'), xlabel '\alpha'
	else
		im(x1, x2, bsum/bmean, 'bsum/bmean'), cbar
	end
end

if 0
	clf
	for j1=1:nj1
		for j2=1:nj2
			subplot(nj1, nj2, j1 + (j2-1)*nj1)
			if ndim == 1
				plot(x1, bb(:,1,j1,j2))
			else
				im(x1, x2, bb(:,:,j1,j2))
			end
		end
	end
end
