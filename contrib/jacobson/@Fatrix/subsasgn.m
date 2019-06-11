 function ob = subsasgn(ob, S, arg)
%function ob = subsasgn(ob, sub, arg)
% method for "ob.sub = arg" etc.


index_string=indexstr(S);

eval(['ob' index_string  '=arg;']);