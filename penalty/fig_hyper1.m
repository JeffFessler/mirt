% fig_hyper1
% show hyperbola for various delta values

dlist = [1000 10 4 1];
tmax = 4;
ptype = 'hyper3';
plist = {'quad', 'huber2', 'hyper3', 'lange1', 'lange3', ...
	'cauchy', 'qgg2'};
t = tmax * linspace(-1, 1, 401)';
for ii=1:length(dlist)
%	ptype = plist{ii};
	delta = dlist(ii);
	leg{ii} = [ptype ' \delta = ' num2str(delta)];
	param = [];

	pot = potential_fun(ptype, delta, param);
	pp(:,ii) = pot.potk(t);
	pw(:,ii) = pot.wpot(t);
	pd(:,ii) = pot.dpot(t);
end

if im
	clf
	plot(t, pp), title 'potk'
	legend(leg)
	xlabel 't'
	ylabel 'p(t)'
end
