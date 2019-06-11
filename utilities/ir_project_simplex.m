 function x = ir_project_simplex(y)
%function x = ir_project_simplex(y)
%|
%| Project each n-dim column in y onto the simplex Dn,
%| Dn = { x : x n-dim, 0 <= x <= 1, sum(x) = 1}.
%|
%| in
%|	y	[n m]	each column is a vector in R^n
%|
%| out
%|	x	[n m]	each column is a vector in R^n
%|
%| Based on algorithm by Xiaojing Ye from paper with Yunmei Chen (chen:11:poa)
%| http://arxiv.org/abs/1101.6081
%|
%| See also duchi:08:epo doi: 10.1145/1390156.1390191
%|
%| 2014-01-09, Jeff Fessler, University of Michigan

% The original version of this was based on code by
% (c) Xiaojing Ye
% xyex19@gmail.com
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1101.6081
% or
% http://ufdc.ufl.edu/IR00000353/
%
% Jan. 14, 2011.
%
% 2012-06-08, modified by JF to include built-in help and example.
% 2014-01-09, rewritten by JF to vectorize.

if nargin < 1, ir_usage, end
if streq(y, 'test'), ir_project_simplex_test, return, end

x = ir_project_simplex_jf(y);


% ir_project_simplex_jf()
% vectorized version, applied to each column
function x = ir_project_simplex_jf(y)

s = sort(y, 1, 'descend'); % [n m] each column sorted
tmp = cumsum(s);

[n m] = size(y);
ii = [1:n]';
tmax = (tmp - 1) ./ repmat(ii, [1 m]);
tmp = [tmax(1:end-1,:); inf(1,m)] >= [s(2:end,:); zeros(1,m)];
tmp = cumsum(tmp, 1) >= 1;
tmp = n + 1 - sum(tmp); % [n 1] first index in each col where tmax >= s
tmax = tmax(sub2ind([n m], tmp, 1:m)); % [n 1] corresponding value of tmax
x = max(y - repmat(tmax, [n 1]), 0);


% ir_project_simplex_orig()
% original version by Xiaojing Ye, uses loops and handles only one column
function x = ir_project_simplex_orig(y)
m = length(y); bget = false;

s = sort(y, 'descend');

tmpsum = 0;

for ii = 1:m-1
	tmpsum = tmpsum + s(ii);
	tmax = (tmpsum - 1) / ii;
	if tmax >= s(ii+1)
		bget = true;
	break
	end
end
 
if ~bget, tmax = (tmpsum + s(m) - 1) / m; end;

x = max(y-tmax,0);


% ir_project_simplex_test2()
function ir_project_simplex_test2

mm = 2^20; % # of trials to compute
m = 30; % # of trials to show
rng(0)
yy = 4*(rand(2,mm)-0.5);
args = {[0 1], [1 0], '-'};
cpu etic
x1 = ir_project_simplex(yy);
cpu etoc 'old'
if has_mex_jf % test new mex version based on condat:16:fpo
	cpu etic
	x2 = jf_mex('simplex', single(yy), int32(1));
	cpu etoc 'new'
	equivs(x1, x2)
end
for ii=1:m
	yi = yy(:,ii);
	xi = ir_project_simplex_orig(yi);
	jf_equal(x1(:,ii), xi)
	args = {args{:}, [xi(1) yi(1)], [xi(2) yi(2)], ':'};
end
if im
	clf
	plot(args{:})
	axis([-1 1 -1 1]*2), xtick([0 1]), ytick([0 1]), grid
	axis equal
end


% ir_project_simplex_test3()
function ir_project_simplex_test3

m = 80; % # of trials
rng(0)
zz = 4*(rand(3,m)-0.5);
args = {[1 0 0 1], [0 1 0 0], [0 0 1 0], '-'};
tmp = ir_project_simplex(zz);
for ii=1:m
	z = zz(:,ii);
	w = ir_project_simplex_orig(z);
	jf_equal(tmp(:,ii), w)
	args = {args{:}, [w(1) z(1)], [w(2) z(2)], [w(3) z(3)], ':'};
end
if im
	clf
	plot3(args{:})
	axis([-1 1 -1 1 -1 1]*2), xtick([0 1]), ytick([0 1]), ztick([0 1]), grid
	axis equal
end


% ir_project_simplex_test()
function ir_project_simplex_test
ir_project_simplex_test2
prompt
ir_project_simplex_test3
