function h = huber_pot(t, d)
%| function h = huber_pot(t, d)
%| huber potential function

if nargin < 1, help(mfilename), error(mfilename), end
if nargin < 2, huber_pot_test, return, end

h = t.^2 / 2;
ii = abs(t) > d;
h(ii) = d * abs(t(ii)) - d.^2/2;

function huber_pot_test
t = linspace(-9,9,101);
delta = 4;
plot(t, huber_pot(t, delta), '-', delta, huber_pot(delta, delta), 'o')
