 function y = mtimes(ob, x)
%function y = mtimes(ob, x)
% y = G * x	or x = G' * y

% full operation
if ~ob.i_block || ob.nblock == 1
	y = ob.base * x;

% block operation
% caution: i_block is really "i_start" (index of starting block in 1...nblock)
else
	y = mtimes_block(ob.base, x, ob.i_block, ob.nblock);
end
