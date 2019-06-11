  function out = ir_iter_fld_write(x, iter, varargin)
%|function out = ir_iter_fld_write(x, iter, varargin)
%|
%| fld_write image arg.x_fun(x) to file arg.file_fun(iter)
%|
%| option
%|	'file_fun'	default @(iter) '%04d.fld'
%|	'x_fun'		default @(x) x;
%|
%| Jeff Fessler

arg.file_fun = @(iter) '%04d.fld';
arg.x_fun = @(x) x;

arg = vararg_pair(arg, varargin);

out = 0;
file = arg.file_fun(iter);
x = arg.x_fun(x);

fld_write(file, x)
file = sprintf('%s/%04d', iter);
