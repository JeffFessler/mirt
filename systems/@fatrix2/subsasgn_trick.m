 function ob = subsasgn_trick(ob, field, arg)
%function ob = subsasgn_trick(ob, sub, arg)
% method for ob.field = arg
% tricky way for dangerously replacing elements of object

ob = struct(ob);
ob.(field) = arg;
ob = fatrix2(ob);
