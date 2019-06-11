 function g = huber_dpot(t, d)
%function g = huber_dpot(t, d)
% huber potential derivative function
g = t;
ii = abs(t) > d;
g(ii) = d * sign(t(ii));
