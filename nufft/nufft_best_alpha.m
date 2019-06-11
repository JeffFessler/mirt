 function [alpha, beta, ok] = nufft_best_alpha(J, L, K_N)
%function [alpha, beta, ok] = nufft_best_alpha(J, L, K_N)
%
%	return previously numerically optimized alpha and beta
%	use L=0 to return best available choice of any L
%
%	Copyright 2001-12-17	Jeff Fessler	The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin < 2, L=2; end
if nargin < 3, K_N=2; end

Jlist1 = 6;		% list of which J's
alpha1 = [...		% last colum is best beta
	1 -0.46		0.19;	% J=6
];

Jlist2 = 2:10;		% list of which J's
alpha2 = [...		% last colum is best beta
	1 -0.200 -0.04	0.34;	% J=2
	1 -0.485 0.090	0.48;
	1 -0.470 0.085	0.56;	% J=4
	1 -0.4825 0.12	0.495;
	1 -0.57 0.14	0.43;	% J=6
	1 -0.465 0.07	0.65;
	1 -0.540 0.16	0.47;	% J=8
	1 -0.625 0.14	0.325;	% J=9
	1 -0.57 0.185	0.43;	% J=10 5.9707e-07
];

Jlist3 = [4 6];		% list of which J's
alpha3 = [...		% last colum is best beta
	1 -0.5319 0.1522 -0.0199	0.6339;	% J=4 2.5953e-04
	1 -0.6903 0.2138 -0.0191	0.2254;	% J=6 1.0097e-04
];


if K_N == 2

	if L==0
		if any(J == Jlist3)
			L = 3;
		else
			L = 2;
		end	% current best
	end

	if L==1
		alpha = alpha1;
		Jlist = Jlist1;

	elseif L==2
		alpha = alpha2;
		Jlist = Jlist2;

	elseif L==3
		alpha = alpha3;
		Jlist = Jlist3;

	else
		warning 'L not done'
		alpha = nan;
		beta = nan;
		ok = 0;
		return
	end

else
	warning 'K_N not done'
	alpha = 1;
	beta = 0.5;
	ok = true;
	return
end

if any(J == Jlist)
	j = find(J == Jlist);
	beta = alpha(j,end);
	alpha = alpha(j,1:(end-1));
	ok = true;
else
	ok = false;
	alpha = nan;
	beta = nan;
end
