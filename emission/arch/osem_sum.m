 function Asums = osem_sum(G, yi, ci, ri, block_order)
%	precompute some sums needed for OSEM
%	NOT DONE!!!!!!!!!!!
%	block_order is a vector of length nblock
%	that specifies the starting angle for each block
%	simplest case is just block_order = 1:nblock,
%	but that choice is suboptimal
%	if block_order is an empty matrix, then 1 block is used
%	if block_order is a scalar (nblock) power of 2 != 1,
%	then the "bit-reversal ordering" is used
%	(nblock can be 1 to na)
%
%	Copyright Apr 1999, Jeff Fessler

if nargin < 3, help osem_sum, error args, end

[nb, na] = size(yi);

if (nargin < 4 || isempty(ci))
	ci = ones(size(yi));
end
if (nargin < 5 || isempty(ri))
	ri = zeros(size(yi));
end
if (nargin < 6 || isempty(block_order))
	block_order = 1;
end
nblock = length(block_order);
if length(block_order) == 1
	nblock = block_order;
	block_order = 1:nblock;
end
if (nargin < 7 || isempty(Asums))
	% compute Asum

	not done!!!!!!!1
	for ii=1:nblock
		ib = block_order(ii);
	end
	Asum = zeros([size(x) nblock]);
end
if (nargin < 7 || isempty(nblock))
	nblock = 1;
end

Gt = G';	% easier to work with transpose

for ib=1:nblock
	ia = ib:nblock:na;
	Asum = zeros(size(x));
	eterm = zeros(size(x));

	%	Asum = G' * ci
	for ia=ib:nblock:na
		Asum = Asum + Gt(:,[1:nb]+(ia-1)*na) * ci(ia,:)';
	end
	return
	yp = ci .* (G * x) + ri;	% predicted measurements
	eterm = G' * (ci .* (yi ./ yp));
	x = x .* eterm ./ Asum;
end
