 function pot_time_test(t)
%function pot_time_test(t)
% compare timing: code < function < implicit

if nargin < 1
	t = linspace(0, 100, 512*512*4)';
end

w = t;
w = 1;
nrep = 20;

wpot = @(t, d) 1 ./ sqrt(1 + abs(t / d).^2);
dpot = @(t, d) t ./ sqrt(1 + abs(t / d).^2);

d = 0.5;
% warm up
if 1
	wpot(0,d);
	dpot(0,d);
	hyper_dpot(0, d);
	hyper_wpot(0, d);
end

cpu etic
for ii=1:nrep
	a = w .* (t ./ sqrt(1 + abs(t / d).^2));
	b = w .* (1 ./ sqrt(1 + abs(t / d).^2));
end
cpu etoc 'code'

cpu etic
for ii=1:nrep
	a = w .* hyper_dpot(t, d);
	b = w .* hyper_wpot(t, d);
end
cpu etoc 'func'

cpu etic
for ii=1:nrep
	a = w .* dpot(t, d);
	b = w .* wpot(t, d);
end
cpu etoc 'implicit'



 function wpot = hyper_wpot(t, d)
%function wpot = hyper_wpot(t, d)
wpot = 1 ./ sqrt(1 + abs(t / d).^2);


 function dpot = hyper_dpot(t, d)
%function dpot = hyper_dpot(t, d)
dpot = t ./ sqrt(1 + abs(t / d).^2);
