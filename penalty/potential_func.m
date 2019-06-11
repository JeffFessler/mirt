 function pot = potential_func(type, delta, param)
%function pot = potential_func(type, delta, param)
%|
%| Define roughness penalty potential functions (using function handles).
%|
%| The penalty will have the form
%|	R(x) = sum_k w_k * potential_k([Cx]_k, delta_k)
%| where w_k is provided elsewhere, not here!
%|
%| in
%|	type		quad, huber, hyper2, hyper3, cauchy, lange1, lange3, ...
%|			recommended: 'hyper3'
%|	delta		scalar or image-sized array;
%|			"cutoff" parameter for edge-preserving regularization
%|	param		optional additional parameter(s) for some choices
%| out
%|	pot.delta	local copy
%|	anonymous functions:
%|	pot.potk(pot, C*x)	potential function value
%|	pot.wpot(pot, C*x)	potential 'weights' (aka half-quad. curvatures)
%|	pot.dpot(pot, C*x)	potential derivative
%|
% Copyright 2004-5-18, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(type, 'test') potential_func_test, return, end

pot.delta = [];
pot.param = [];
if isvar('delta')
	pot.delta = delta(:);
	clear delta
end
if isvar('param')
	pot.param = param(:);
	clear param
end

% trick: huber2 is just a huber with delta / 2
% so that weighting function drops to 1/2 at delta, like hyper3 etc.
if streq(type, 'huber2')
	type = 'huber';
	pot.delta = pot.delta / 2;
end

% trick: hyper3 is just a hyperbola with delta scaled by sqrt(3)
% to approximately "match" the properties of 'cauchy' (old erroneous 'hyper')
if streq(type, 'hyper3')
	type = 'hyper2';
	pot.delta = pot.delta / sqrt(3);
end

switch type

%
% quadratic potential function
%
case 'quad'
	potk = '(abs(t).^2)/2';
	wpot = 'ones(size(t))';
	dpot = 't';

%
% huber potential function
%
case 'huber'
	potk = 'huber_pot(t, pot.delta)';
	wpot = 'huber_wpot(t, pot.delta)';
	dpot = 'huber_dpot(t, pot.delta)';

%
% cauchy penalty: d^2 / 2 * log(1 + (t/d)^2)
% Not convex!
%
case 'cauchy'
	potk = 'pot.delta.^2 / 2 .* log(1 + abs(t ./ pot.delta).^2)';
	wpot = '1 ./ (1 + abs(t ./ pot.delta).^2)';
	dpot = 't ./ (1 + abs(t ./ pot.delta).^2)';

%
% Geman&McClure penalty: d^2 / 2 * (t/d)^2 / (1 + (t/d)^2)
% Not convex!
%
case 'geman&mcclure'
	potk = 'pot.delta.^2 / 2 .* (t/pot.delta).^2 ./ (1 + abs(t ./ pot.delta).^2)';
	wpot = '1 ./ (1 + abs(t ./ pot.delta).^2).^2';
	dpot = 't ./ (1 + abs(t ./ pot.delta).^2).^2';

%
% hyperbola penalty: d^2 * [ sqrt(1 + (t/d)^2) - 1 ]
%
case 'hyper2'
	potk = 'pot.delta.^2 .* (sqrt(1 + abs(t ./ pot.delta).^2) - 1)';
	wpot = '1 ./ sqrt(1 + abs(t ./ pot.delta).^2)';
	dpot = 't ./ sqrt(1 + abs(t ./ pot.delta).^2)';

case 'hyper'
	error 'use "cauchy" or "hyper3" not "hyper" now'

%
% Lange1 penalty
%
case 'lange1'
	potk = 't.^2 / 2 ./ (1+abs(t./pot.delta))';
	wpot = '(1 + abs(t ./ pot.delta) / 2) ./ (1 + abs(t ./ pot.delta)).^2';
	dpot = ['t .* (' wpot ')'];

%
% Lange3 penalty
%
case 'lange3'
	potk = 'pot.delta.^2 .* (abs(t./pot.delta) - log(1+abs(t./pot.delta)))';
	wpot = '1 ./ (1 + abs(t ./ pot.delta))';
	dpot = 't ./ (1 + abs(t ./ pot.delta))';

%
% li98cfs
%
case 'li98cfs'
	% f = @(x) atan(x) / x - 0.5; fsolve(f, 2.3)
	pot.delta = pot.delta / 2.3311;
	potk = 'pot.delta.^2 .* ((t ./ pot.delta) .* atan(t ./ pot.delta) - 0.5 * log(1 + (t ./ pot.delta).^2))';
	wpot = 'ir_li98cfs_wpot(t, pot.delta)';
	dpot = ['t .* ' wpot];

%
% qgg2: q-generalized gaussian for p=2, due to Thibault, Sauer, Bouman
% q = "param", same as lange1 when q=1
%
case 'qgg2'
	potk = 't.^2 / 2 ./ (1+abs(t./pot.delta).^(2-pot.param))';
	wpot = ['(1 + abs(t ./ pot.delta).^(2-pot.param) * pot.param / 2) ' ...
		' ./ (1 + abs(t ./ pot.delta).^(2-pot.param)).^2'];
	dpot = ['t .* (' wpot ')'];

%
% genhub : "generalized" Huber, same as huber when p=2 and q=1
% p = param(1), q = param(2)
%
case 'genhub'
	potk = ['0.5 * abs(t) .^ p .* (abs(t) <= d) + ' ...
		'0.5 * (p ./ q .* d .^ (p-q) .* abs(t) .^ q + (1 - p ./ q) .* d .^ p) .* (abs(t) > d)'];
	reps = {'p', 'zp', 'q', 'zq', 'd', 'zd', ...
		'zp', 'pot.param(1)', 'zq', 'pot.param(2)', 'zd', 'pot.delta'};
	potk = strreps(potk, reps{:});
	wpot = ['p / 2 .* (abs(t) .^ (p-2)) .* (abs(t) <= d) + ' ...
		'p / 2 .* (d .^ (p-q)) .* ((abs(t) + (abs(t)==0)) .^ (q-2)) .* (abs(t) > d)']; % trick to avoid 0^-:
	wpot = strreps(wpot, reps{:});
	dpot = ['t .* (' wpot ')'];

%
% stevenson:94:dpr
% p = param(1), q = param(2), same as huber when p=2 and q=1 ???
%
case 'stevenson94dpr'
	potk = ['0.5 * abs(t) .^ p .* (abs(t) <= d) + ' ...
		'0.5 * ( (p .* d .^ (p-1) .* abs(t) - p .* d .^ p + (1 ./ q) .^ (1 ./ (q-1)) ) .^ q + d .^ p - (1 ./ q) .^ (q ./ (q-1)) ) .* (abs(t) > d)'];
	reps = {'p', 'zp', 'q', 'zq', 'd', 'zd', ...
		'zp', 'pot.param(1)', 'zq', 'pot.param(2)', 'zd', 'pot.delta'};
	potk = strreps(potk, reps{:});
	wpot = 'ones(size(t))'; % fix: fake for now
	dpot = ['t .* (' wpot ')'];

otherwise
	try % trick: try new strum version
		p = potential_fun(type, pot.delta, pot.param);
		pot.potk = @(dum, t) p.potk(t);
		pot.wpot = @(dum, t) p.wpot(t);
		pot.dpot = @(dum, t) p.dpot(t);
	catch
		fail('Unknown potential "%s"', type)
	end
end

if ~isfield(pot, 'potk')
	pot.potk = ir_str2func(['@(pot,t) ' potk]);
	pot.wpot = ir_str2func(['@(pot,t) ' wpot]);
	pot.dpot = ir_str2func(['@(pot,t) ' dpot]);
end


% test routine
% examine potential functions after rescaling.
function potential_func_test

delta = 10; tmax = 4;
plist = {'quad', 'huber2', 'hyper3', 'lange1', 'lange3', ...
	'cauchy', 'qgg2', 'gf1'};
%plist = {'quad', 'li98cfs', 'hyper3', 'huber2'}; % show li98cfs roughly hyper3
%plist = {'huber', 'genhub', 'quad'};%, 'stevenson94dpr'};
%plist = {'hyper3', 'qgg2'}; delta = 20; tmax = 10;
%plist = {'huber', 'hyper3', 'qgg2'}; delta = 10; tmax = 8;
%plist = {'qgg2', 'gf1'};
t = tmax * linspace(-delta, delta, 401)';
for ii=1:length(plist)
	type = plist{ii};
	leg{ii} = [plist{ii} ' \delta = ' num2str(delta)];
	if streq(type, 'qgg2')
		param = 1.2;
		leg{ii} = [leg{ii} ' q = ' num2str(param)];
	elseif streq(type, 'gf1')
		param = [0.0558 1.6395];
	elseif streq(type, 'genhub')
		param = [1.9 1.1];
		leg{ii} = [leg{ii} sprintf('p=%g q=%g', param(1), param(2))];
	elseif streq(type, 'stevenson94dpr')
		param = [2 2.01];
		leg{ii} = [leg{ii} sprintf('p=%g q=%g', param(1), param(2))];
	else
		param = [];
	end
	pot = potential_func(type, delta, param);
	pp(:,ii) = pot.potk(pot, t);
	pw(:,ii) = pot.wpot(pot, t);
	pd(:,ii) = pot.dpot(pot, t);
end

if im
	clf
	subplot(411), plot(t, pp), title 'potk'
	axis tight
	axisy([-0.0 2.5] * delta^2)
	legend(leg, 'location', 'north')
	subplot(412), plot(t, pw), title 'wpot'
	axisy(0, 1.1)
	subplot(413), plot(t, pd), title 'dpot'
	axisy([-1 1] * 1.1 * delta)
	% check derivatives
	subplot(414)
	plot(t, pd)
	title 'dpot check'
	hold on
	d = diffc(pp) / (t(2)-t(1));
	plot(t(1:end-1), d(1:end-1,:), '--')
	hold off
	axisy([-1 1] * 1.1 * delta)
end
