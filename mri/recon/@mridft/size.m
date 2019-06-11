function dim = size(ob)
%function dim = size(ob)
%       "size" method for "mridft" class

dim = [size(ob.t,1) size(ob.we,1)];

