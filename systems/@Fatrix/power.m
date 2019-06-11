 function ob = power(ob, sup)
%function ob = power(ob, sup)
% array "power" method (G.^2)
% matrix power (mpower) would be for G^2

if isempty(ob.handle_power)
	error 'no power method for this object'
end

ob = ob.handle_power(ob, sup);
