  function out = wls_info_step
%|function out = wls_info_step
%| a 'userfun' for the pcg type of routines that also prints
%| Copyright 2008-09-27, Jeff Fessler, University of Michigan

gamma = evalin('caller', 'gamma');
step = evalin('caller', 'step');
ddir = evalin('caller', 'ddir');
change = step * ddir(:);
x = evalin('caller', 'x');
out = [gamma step cpu('toc')];

printf('change: %g %g, %g%%, x: %g %g', ...
	min(change), max(change), ...
	100 * max(abs(change)) / max(x(:)), ...
	min(x(:)), max(x(:)))
