 function pn = jf_protected_names
%function pn = jf_protected_names
%|
%| A serious drawback of the original matlab language was its lack
%| of a protected or local namespace.  Every m-file that is in the path
%| is available to all functions (except those in "private" subdirectories).
%| Users who have their own m-files that happen to have the same names as
%| any of the routines in a toolbox like this one will have problems.
%|
%| To try to overcome this limitation, I created this function in late 2009
%| to serve as a repository of simple functions.
%| To use any of these functions, one types something like
%|	pn = jf_protected_names;
%| and then one can call the functions using
%|	out = pn.fun(arg1, arg2, ...);
%|
%| Mathworks eventually added the "+package/import"  capability...
%|
%| Copyright 2009-11-21, Jeff Fessler, University of Michigan


pn = strum(struct, { ...
	'struct_recurse', @jf_struct_recurse, '()';
	'color_order', @jf_color_order, '()';
	'diary', @jf_diary, '(file)';
	'has_hct2', @jf_has_hct2, '()';
	'hct_arg', @jf_hct_arg, '(cg, ig)';
	'ind2sub', @jf_ind2sub, '(siz, ind)';
	'mid3', @jf_mid3, '(im_3d, [dim])';
	'normcdf', @jf_normcdf, '(x, mu, sigma)';
	'prctile', @jf_prctile, '(x, p, [dim])';
	'cases', @jf_cases, '(x, v0, v1, ...)';
	'test', @jf_test, '()';
	});

end % jf_protected_names()


% jf_cases()
% its purpose is to act like the '?' operator in C
% jf_cases(x, 'value if x is 0', 'value if x is 1', ...)
function out = jf_cases(st, x, varargin)
if x+1 > numel(varargin)
	fail('x=%d but num=%d', x, numel(varargin))
end
out = varargin{x+1};
end % jf_cases()


% jf_diary()
% this version prompts if file exists!
function dummy = jf_diary(st, file)
if streq(file, 'off')
	diary('off');
return
end

if exist(file, 'file')
	fail 'file exists'
else
	printm('starting diary for "%s"', file)
	diary(file);
end
dummy = [];
end % jf_diary()


% jf_ind2sub()
% version with a single matrix output, one dimension per column
function subs = jf_ind2sub(st, Nd, ind)
ind = ind(:);
subs = zeros(numel(ind), numel(Nd));
switch numel(Nd)
case 2
	[subs(:,1) subs(:,2)] = ind2sub(Nd, ind);
case 3
	[subs(:,1) subs(:,2) subs(:,3)] = ind2sub(Nd, ind);
case 4
	[subs(:,1) subs(:,2) subs(:,3) subs(:,4)] = ind2sub(Nd, ind);
otherwise
        fail 'not done'
end
end % jf_ind2sub()


% jf_test()
function yn = jf_test(st, varargin)
yn = jf_equal(st.cases(1, 20, 30), 30)
end % jf_test()
