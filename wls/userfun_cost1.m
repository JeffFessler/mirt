  function cost = userfun_cost1(x, varargin)
%|function cost = userfun_cost1(x, varargin)
%|
%| cost(x) = (y-Ax)' W (y-Ax)/2 + R(x)
%|
%| Copyright 2011-4, Jeff Fessler, University of Michigan

A = evalin('caller', 'A');
W = evalin('caller', 'W');
R = evalin('caller', 'R');
yi = evalin('caller', 'yi');
%x = x(:,end);
cost = pwls_cost(x, A, W, yi, R);
