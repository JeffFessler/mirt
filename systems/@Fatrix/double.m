 function ob = double(ob)
%function ob = double(ob)

if isfield(ob, 'handle_double')
	ob = ob.handle_double(ob);
else
	tmp = full(ob);
	if streq(class(tmp), 'Fatrix')
		fail 'bug: still a Fatrix?'
	else
		ob = double(tmp);
	end
end
