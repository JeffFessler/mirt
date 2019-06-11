 function out = isfreemat
%function out = isfreemat
% determine if this is freemat or matlab!
tmp = version;
out = streq(tmp, '4.0');
