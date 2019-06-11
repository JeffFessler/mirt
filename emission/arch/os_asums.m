 function Asums = os_asums(Gt, ci, nb, nsubset)
%function Asums = os_asums(Gt, ci, nb, nsubset)
%
%	for ordered subsets (block iterative, incremental gradient) algorithms,
%	we need to precompute the following backprojections:
%		\sum_{i \in subset_m} a_ij = \sum_{i \in subset_m} c_i g_ij

if nargin ~= 4, help(mfilename), error(mfilename), end

na = ncol(Gt) / nb;
if na * nb ~= ncol(Gt)
	error('block size not a divisor of projection data size')
end

[starts, nsubset] = subset_start(nsubset);

if ~isvar('ci') || isempty(ci)
	ci = ones(nb,na);
end

for iset=1:nsubset
        ia = starts(iset):nsubset:na;
        ii = outer_sum(1:nb,(ia-1)*nb);         % relevant row indices
        Asums(:,iset) = Gt(:,ii) * col(ci(:,ia));
end
