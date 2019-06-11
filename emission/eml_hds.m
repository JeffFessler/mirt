 function gam = eml_hds(Gb, ci, ri, hds)
%function gam = eml_hds(Gb, ci, ri, hds)
% precompute hidden data space factor, based on 1995 PML SAGE paper

if nargin < 4, help(mfilename), error(mfilename), end

if iscell(ci)
	gam = inf;
	nblock = block_op(Gb, 'n');
	for ii=1:nblock
		gam_new = eml_hds_1(Gb{ii}, ci{ii}, ri{ii}, hds);
		gam = min(gam, gam_new);
	end
else
	gam = eml_hds_1(Gb, ci, ri, hds);
end


%
% eml_hds_1()
%
function gam = eml_hds_1(Gb, ci, ri, hds)

if hds == 1
	gam = 0;

elseif hds == 3
	ai = ci(:) .* col(Gb * ones(size(Gb,2),1)); % sum_j \aij
	gam = min(ri(ai~=0) ./ ai(ai~=0));

else
	error 'bad hds'
end
