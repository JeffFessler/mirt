  function sp = aspire_buff_parse(buff)
%|function sp = aspire_buff_parse(buff)
%|
%| parse 'buffer' read from aspire .wtf
%| primarily for testing
%|
%| in
%|	buff	byte	read using wtfmex()
%|
%| out
%|	sp	struct	header information etc.
%|
%| Copyright 2011-3-4, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end

[buff sp.head] = aspire_buff_skip_ff(buff);

% methods to cast buffer bytes into appropriate types
geti = @(i, b) typecast(b((4*(min(i)-1)+1):4*(max(i)-1)+4), 'uint32');
getf = @(i, b) typecast(b((4*(min(i)-1)+1):4*(max(i)-1)+4), 'single');

% first 128 bytes is sparse header

sp.group_by = geti(1, buff);
sp.index_by = geti(2, buff);
sp.value_by = geti(3, buff);
sp.nrow_used = geti(4, buff);
sp.ncol_used = geti(5, buff);
sp.nwt = geti(6, buff);
sp.max_wt = getf(7, buff);
sp.min_wt = getf(8, buff);
sp.total = getf(9, buff);
sp.max_index = geti(10, buff);
sp.min_index = geti(11, buff);
sp.nx = geti(12, buff);
sp.ny = geti(13, buff);
sp.nb = geti(14, buff);
sp.na = geti(15, buff);
sp.fill = geti(16:32, buff);

buff = buff(129:end);

% mask
sp.mask = reshape(buff(1:(sp.nx*sp.ny)), sp.nx, sp.ny);
buff = buff((sp.nx*sp.ny+1):end);

switch sp.group_by
case 0 % by_row
	sp.ngroup = sp.nb * sp.na;
case 1 % by_col
	sp.ngroup = sp.nx * sp.ny;
otherwise
	fail('unknown group_by %d', sp.group_by)
end

% length & offset
sp.length = geti(1:sp.ngroup, buff);
buff = buff((4*sp.ngroup+1):end);
sp.offset = geti(1:sp.ngroup, buff);
buff = buff((4*sp.ngroup+1):end);

switch sp.index_by
%case 0 % by_uint2
case 1 % by_uint4
	sp.index = geti(1:sp.nwt, buff);
	buff = buff((4*sp.nwt+1):end);
otherwise
	fail('unknown index_by %d', sp.index_by)
end

switch sp.value_by
case 0 % by_float4
	sp.value = getf(1:sp.nwt, buff); % final "sp.nwt" entries are matrix values
	buff = buff((4*sp.nwt+1):end);
otherwise
	fail('unknown value_by %d', sp.value_by)
end

if length(buff)
	fail 'buffer not empty at end'
end


% aspire_buff_skip_ff()
function [buff head] = aspire_buff_skip_ff(buff)
tmp = find(buff == 12); % \f
f1 = min(tmp);
if buff(f1+1) ~= 12, error 'two form feeds', end
head = char(buff(1:(f1-1))'); % ascii header
buff = buff((f1+2):end);
