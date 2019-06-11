 function [yik, scalefactor] = ir_mri_field_map_reg_scale(yik, etime, varargin)
%function [yik, scalefactor] = ir_mri_field_map_reg_scale(yik, etime, varargin)
%|
%| Scale images to account for R2* effects and differences in absolute value
%| using median(ri), where ri = sum_j sum_k |yik_j yik_k|^2 (t_k - t_j)^2
%|
%| in
%|	yik	[N nset]	scan images
%|	etime	[1 nset]	echo times (units of sec if fieldmap is in Hz)
%|
%| option
%|	fmax			threshold for absolute yik value (default 0.1)
%|	dmax			threshold for absolute rj value (default 0.1)
%|	show	0|1		1 to show result (default 0)
%|
%| out
%|	yik	[N nset]	scaled scan images
%|	scalefactor		sqrt(median(rj))
%|
%| MJ Allison
%| 2015-06-27 JF cosmetic changes

if nargin < 4, help(mfilename), error(mfilename), end

arg.fmax = 0.1;
arg.dmax = 0.1;
arg.show = 0;
arg = vararg_pair(arg, varargin);

dim_yik = size(yik);
yik = reshapee(yik, [], dim_yik(end)); % [*N nset]

[nn nset] = size(yik);

% Scale by median of first set of data to get rid of large mag_j
% effects (not actually needed, but left in for consistency)
if arg.fmax > 0
	y1 = abs(yik(:,1));
	scalefactor = median(y1(y1(:) > arg.fmax * max(y1(:))));
	if scalefactor == 0
		fail 'median is zero?'
	end
	yik = yik / scalefactor;
else
	scalefactor = 1;
end


% Try to compensate for R2 effects on effective regularization.

d = zeros(nn,1);
for j=1:nset
	for k=1:nset
		tmp = abs(yik(:,j) .* yik(:,k)).^2 * (etime(k)-etime(j)).^2;
		d = d + tmp;
	end
end

% divide by numerator of wj^mn -> sum(abs(y)^2)
d = div0(d, sum(abs(yik).^2,2));

if arg.show
	im(reshape(d, [dim_yik(1:2)]))
end

% compute dtyp
dtyp = median(d(d(:) > arg.dmax * max(d(:))));

% now uniformly scale by the square root of dtyp
for j = 1:nset
	yik(:,j) = div0(yik(:,j),sqrt(dtyp));
end

scalefactor = sqrt(dtyp);
yik = reshape(yik, dim_yik); % [(N) nset]
