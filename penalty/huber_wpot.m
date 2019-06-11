 function w = huber_wpot(t, d)
%function w = huber_wpot(t, d)
% huber potential 'weighting' (i.e. curvature) function
w = ones(size(t));
ii = abs(t) > d;
w(ii) = d ./ abs(t(ii));
