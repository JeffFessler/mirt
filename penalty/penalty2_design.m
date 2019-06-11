 function [C, rj] = penalty2_design(type, varargin)
%|function [C, rj] = penalty2_design(type, ['leak'|'tight'], wang, ang, mask)
%|
%| Design the penalty matrix "C" (and penalty coefficients "rj")
%| for a quadratic penalty R(x) = 1/2 x' C * C * x with:
%| 1st-order differences and a 2nd-order neighborhood (8 neighbors).
%| For 2D parallel-beam tomography with shift-invariant blur.
%| Design based on Fourier method in fessler:03:aat (IEEE MIC, 2003).
%|
%| in
%|	type			'test' or 'quad,d1,n2' or 'quad,d1,n1'
%|				(n2 for usual 2nd order neighborhood)
%|	wang	[na nx ny]	angular weights for each pixel
%|	ang	[na]		angles
%|	mask	[nx ny]		(logical) reconstruction support
%| out:
%|	C	[4*nx*ny nx*ny] "modified" C containing sqrt{r_j} factors
%|	rj	[nx ny 4]	horiz,vert,diag1,diag2 coefficients
%|
%| Copyright 2003-5-23, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end

% run a self-test to compare my analytical solution
% to the numerical solution computed using NNLS.
if nargin == 1 && streq(type, 'test')
	penalty2_design_test
return
end

% analytical design of 1st-order difference, 2nd-order neighborhood
if streq(type, 'quad,d1', 7)
	if length(varargin) < 1 || length(varargin) > 3
		help(mfilename), error(mfilename)
	end

	if ischar(varargin{1})
		Ctype = varargin{1};
		varargin = varargin{2:end};
	else
		Ctype = 'leak';
	end

	wang = varargin{1};
	[na nx ny] = size(wang);
	np = nx*ny;
	wang = reshape(wang, [na np]);

	if length(varargin) >= 2
		ang = col(varargin{2});
		if na ~= length(ang), error 'angle dim', end
	else
		ang = [0:(na-1)]'/na * pi; % [0,pi)
	end

	if length(varargin) == 3
		mask = varargin{3};
		if any(size(mask) ~= [nx ny]), error 'mask dims', end
	else
		mask = true(nx,ny);
	end

	% dot product of each basis with wj
	% the "mean" takes care of normalization
	basis = [ones(size(ang)) cos(2*ang) sin(2*ang)];
	dot_products = zeros(ncol(basis),np);
	for ib=1:ncol(basis)
		b = basis(:,ib);
		b = repmat(basis(:,ib), 1, np);
		dot_products(ib,:) = mean(b .* wang, 1);
	end, clear ib b

	sptmp = @(tmp) spdiag(tmp, 'nowarn'); % for now

	if streq(type, 'quad,d1,n1') % 1st-order neighborhood
		rj = penalty2_design_d1_n1(dot_products); % [2 np]
		rj = permute(reshape(rj, [2 nx ny]), [2 3 1]); % [nx ny 2]
		C = C2sparse('leak', mask, 4, 0, 0);
		i1 = 1:(nx*ny);
		C = [	sptmp(sqrt(col(rj(:,:,1)))) * C(0*nx*ny+i1,:);
			sptmp(sqrt(col(rj(:,:,2)))) * C(1*nx*ny+i1,:)];

	elseif streq(type, 'quad,d1,n2') % 2nd-order neighborhood
		rj = penalty2_design_d1_n2(dot_products); % [4 np]
		rj = permute(reshape(rj, [4 nx ny]), [2 3 1]); % [nx ny 4]
		% fix: this must depend on flip_y !!
		printf('Warn: penalty2_design.m may fail if flip_y = -1')
		rj = rj(:,:,[1 2 4 3]); % re-order 45,135 (empirical)

		C = C2sparse('leak', mask, 8, 0, 0);
		i1 = 1:(nx*ny);
		% analysis assumes 1/sqrt(2) for diagonal 1st differences!
		C = [	sptmp(sqrt(col(rj(:,:,1))))	* C(0*nx*ny+i1,:);
			sptmp(sqrt(col(rj(:,:,2))))	* C(1*nx*ny+i1,:);
			sptmp(sqrt(col(rj(:,:,3))/2))	* C(2*nx*ny+i1,:);
			sptmp(sqrt(col(rj(:,:,4))/2))	* C(3*nx*ny+i1,:)];
	else
		error 'bad type'
	end

else
	error(['unknown type: ' type])
end


%
% in:	iprod	[3 np]		inner products of angular weights with
%					[1 cos(2a) sin(2a)]
% out:	rj	[2 np]		horiz,vert
%
function rj = penalty2_design_d1_n1(iprod)

d1 = iprod(1,:);
d2 = iprod(2,:);

if any(d1 < 0), error 'bad inner products: d1<0', end
if any(abs(d2) > d1), error 'bad inner products: d2>d1', end

rj = [d1+2*d2; d1-2*d2];

ii = d2 > d1/2;
r1 = 4/3 * (d1 + d2);
rj(1,ii) = r1(ii);
rj(2,ii) = 0;

ii = d2 < -d1/2;
r2 = 4/3 * (d1 - d2);
rj(1,ii) = 0;
rj(2,ii) = r2(ii);



%
% in:
%	iprod	[3 np]		inner products of angular weights with
%					[1 cos(2a) sin(2a)]
%					typically np = # of pixels
% out:
%	rj	[4 np]		horiz,vert,diag1,diag2 coefficients
%
% caution: diagonals coefficients are designed to be used with
% diagonal differences that are normalized by 1/\sqrt{2}.
%
function rj = penalty2_design_d1_n2(iprod)

d1 = iprod(1,:);
dp2 = iprod(2,:);
dp3 = iprod(3,:);
d2 = max(abs(dp2), abs(dp3)); % d2 >= d3
d3 = min(abs(dp2), abs(dp3));
if any(d1 < 0), error 'bad inner products: d1<0', end
if any(d2 > d1), error 'bad inner products: d2>d1', end

rj = zeros(4,length(d1));

ii = (d2 >= d1/2) & (d3 <= 2/3*d2 - d1/3); % case 1
r1 = 4/3 * (d1 + d2);
rj(1,ii) = r1(ii);

ii = (d3 >= 2/3*d2 - d1/3) & (d3 + d2 >= d1/2); % case 2
r1 = 8/5 * (d1/2 + 3/2*d2 - d3);
r3 = 12/5 * (d3 - 2/3*d2 + 1/3*d1);
rj(1,ii) = r1(ii);
rj(3,ii) = r3(ii);

ii = (d3 + d2 <= d1/2) & (d2 >= d1/4); % case 3
rj(1,ii) = 4*d2(ii);
r3 = d1 - 2*d2 + 2*d3;
r4 = d1 - 2*d2 - 2*d3;
rj(3,ii) = r3(ii);
rj(4,ii) = r4(ii);

ii = (d2 <= d1/4); % case 4
r1 = d1/2 + 2*d2;
r2 = d1/2 - 2*d2;
r3 = d1/2 + 2*d3;
r4 = d1/2 - 2*d3;
rj(1,ii) = r1(ii);
rj(2,ii) = r2(ii);
rj(3,ii) = r3(ii);
rj(4,ii) = r4(ii);

% symmetries
ii = abs(dp3) > abs(dp2);
rj(:,ii) = rj([3 4 1 2],ii);

ii = dp3 < 0;
rj([3 4],ii) = rj([4 3],ii);

ii = dp2 < 0;
rj([1 2],ii) = rj([2 1],ii);


%
% brute-force NNLS approach
% dot	[3 np]
% rj	[4 np]
%
function rj = penalty2_design_nnls(Amat, dot)

ww = [dot(1,:); sqrt(2)*dot([2 3],:)];
opt = optimset('lsqnonneg');
opt = optimset(opt, 'tolx', 10*eps);

rj = zeros(4, ncol(ww));
for jj=1:ncol(ww)
	rj(:,jj) = lsqnonneg(Amat, ww(:,jj));
	if ~rem(jj,100), printf('%d of %d', jj, ncol(ww)), end
end



% penalty2_design_test()
function penalty2_design_test
Amat = 0.5 * [1 1 1 1;
		1/sqrt(2) -1/sqrt(2) 0 0 ;
		0 0 1/sqrt(2) -1/sqrt(2)];
n2 = 41;
n3 = 43;
d2 = linspace(-1,1,n2);
d3 = linspace(-1,1,n3);
dot = zeros(3,n2*n3);
[t2, t3] = ndgrid(d2, d3);
dot(1,:) = 1;
dot(2,:) = t2(:)';
dot(3,:) = t3(:)';

rj_anal = penalty2_design_d1_n2(dot);
rj_nnls = penalty2_design_nnls(Amat, dot); % [4 n2*n3]

penalty2_design_test_figure(Amat, dot, d2, d3, rj_anal, rj_nnls)


function penalty2_design_test_figure(Amat, dot, d2, d3, ra, rn)

n2 = length(d2);
n3 = length(d3);

ww = [dot(1,:); sqrt(2)*dot([2 3],:)];
erra = permute(reshape(Amat * ra - ww, [3 n2 n3]), [2 3 1]); % [n2 n3 3]
errn = permute(reshape(Amat * rn - ww, [3 n2 n3]), [2 3 1]);

ra = permute(reshape(ra, [4 n2 n3]), [2 3 1]); % [n2 n3 4]
rn = permute(reshape(rn, [4 n2 n3]), [2 3 1]);

if im
	clf, pl=440;
	for ib=1:4
		im(pl+ib+0, d2, d3, ra(:,:,ib)), axis xy, cbar
		title(sprintf('Anal. ib=%d', ib))
		im(pl+ib+4, d2, d3, rn(:,:,ib)), axis xy, cbar
		title NNLS
	end

	subplot(4,4,9)
	im(d2, d3, sum(ra > 0,3), 'Anal. sum(r > 0)'), axis xy, cbar
	subplot(4,4,10)
	im(d2, d3, sum(rn > 0,3), 'NNLS sum(r > 0)'), axis xy, cbar

	% hold on, plot(d2, 2/3*(d2 - 1/sqrt(2)), 'c-'), hold off
	% hold on, plot(d2, sqrt(2)/3 - 2/3*d2, 'r-'), hold off

	tmp = sqrt(sum(erra.^2, 3));
	subplot(4, 4, 11)
	im(d2, d3, tmp, 'Anal. |err|'), axis xy, cbar
	subplot(4, 4, 15)
	im(d2, d3, tmp < 10*eps, 'Anal. |err|=0'), axis xy, cbar

	tmp = sqrt(sum(errn.^2, 3));
	subplot(4, 4, 12)
	im(d2, d3, tmp, 'NNLS |err|'), axis xy, cbar
	subplot(4, 4, 16)
	im(d2, d3, tmp < 10*eps, 'NNLS |err|=0'), axis xy, cbar

	subplot(4, 4, 13)
	%im(d2, d3, sum(erra.^2, 3) < sum(errn.^2, 3)+10*eps, 'Anal<NNLS')
	im(d2, d3, sum(errn.^2, 3) - sum(erra.^2, 3), 'NNLS-Anal')
	axis xy, cbar
end
