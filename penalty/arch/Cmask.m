 function [out] = Cmask(type, mask)
%function [out] = Cmask(type, mask)
% select appropriate portions of penalty_mex output based on support mask.
% Cmask() .* penalty_mex() should be same as C * x.
% in
%	type	string	
%	mask	[nx,ny[,nz]]	support
% out
%	out	[nx,ny,4]	binary
%
% Copyright 2004-5-16, Jeff Fessler, The University of Michigan

% default is to show arguments and run test case(s)
if nargin < 1, help(mfilename), Cmask_test, return, end

[nx ny nz] = size(mask);
if nz > 1, error '3d not done', end

if streq(type, 'tight,2d,hvd')
	out = zeros([size(mask) 4]);

	out(2:end,:,1) = mask(2:end,:) & mask(1:end-1,:);
	out(:,2:end,2) = mask(:,2:end) & mask(:,1:end-1);
	out(2:end,2:end,3) = mask(2:end,2:end) & mask(1:end-1,1:end-1);
	out(1:end-1,2:end,4) = mask(1:end-1,2:end) & mask(2:end,1:end-1);

else
	help Cmask
	error([type ' unknown'])
end



% run test
function Cmask_test
nx = 16; ny = 14;
nx = 512; ny = 500;
mask = conv2(ellipse_im(nx, ny, [0 0 [nx ny]/2-5 0 1]) > 0, ones(2), 'same') > 0;
beta = Cmask('tight,2d,hvd', mask);
im clf
im(121, mask, 'mask')
im(122, beta, 'beta')

%prompt

tic
C = C2sparse('tight', mask, 8);
printf('make C time %g', toc)
rng(0)
x = rand(nx, ny);
x = dsingle(x);
x = x .* mask;

tic
d1 = C * x(:); 
d1 = dsingle(d1);
printf('Cx time %g', toc)
d1 = reshape(d1, [nx ny 4]); 

tic
beta = Cmask('tight,2d,hvd', mask);
printf('make beta time %g', toc)

tic
d2 = penalty_mex('forw,2d,hvd', single(x));
d2 = double(d2);
d2 = d2 .* beta;
printf('mex time %g', toc)

if 1
	im clf
	im(131, d1, 'C * x')
	im(132, d2, 'penalty--mex')
	im(133, d2-d1, 'err')
	printf('Cx vs penalty_mex: %g%%', max_percent_diff(d1, d2))
end

d = d1;
%d = zeros(nx,ny,4);
%d(end/2,end/2,1) = 1;
tic
x1 = C' * d(:);
printf('C''d time %g', toc)
x1 = dsingle(x1);
x1 = reshape(x1, nx, ny);

tic
x2 = penalty_mex('back,2d,hvd', single(d));
x2 = double(x2);
x2 = x2 .* mask;
printf('mex time %g', toc)

printf('C''x vs penalty_mex: %g%%', max_percent_diff(x1, x2))
if 1
	im clf
	im(131, x1, 'C''x'), cbar
	im(132, x2, 'penalty--mex'), cbar
	im(133, x2-x1, 'err'), cbar
end
