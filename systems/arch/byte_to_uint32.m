 function val = byte_to_uint32(buff)
%function val = byte_to_uint32(buff)
% convert 4-byte buffer into uint32

try
	val = typecast(buff, 'uint32')
catch
	u = @(x) uint32(x);
	buff = u(buff);
	s = u(256);
	val = buff(4) + s * (buff(3) + s * (buff(2) + s * buff(1)));
end
