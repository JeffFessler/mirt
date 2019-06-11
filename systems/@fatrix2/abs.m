 function ob = abs(ob)
%function ob = abs(ob)

if isempty(ob.handle_abs)
	error 'no abs() method for this object'
end

ob = ob.handle_abs(ob);
