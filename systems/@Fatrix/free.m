 function free(ob)
%function free(ob)

if isempty(ob.handle_free)
	warning 'attemp to free Fatrix object with no free() method'
	return
end

ob.handle_free(ob.arg)
