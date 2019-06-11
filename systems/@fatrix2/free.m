 function free(ob)
%function free(ob)

if isempty(ob.handle_free)
	warn('attemp to free %s object with no free() method', class(ob))
return
end

ob.handle_free(ob.arg)

% todo: use onCleanup
