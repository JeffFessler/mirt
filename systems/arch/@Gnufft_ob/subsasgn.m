 function ob = subsasgn(ob, sub, arg)
%function ob = subsasgn(ob, sub, arg)

% trick: just convert it to a structure and then back!
c = class(ob);
ob = struct(ob);
ob = subsasgn(ob, sub, arg);
ob = class(ob, c);
