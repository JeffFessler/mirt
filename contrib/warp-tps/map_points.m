function [newholx, newholy] = map_points(A, W, refx, refy, newrefx, newrefy) 
% map (refx, refy) according to mapping specified by A, W, refx, refy
n=length(refx); %number of control points
U = zeros(n,1);
for k=1:n,
    r_2 = ( newrefx - refx(k) )^2 + (newrefy - refy(k) )^2;
    if r_2 > 0,    U(k) = r_2*log(r_2); end;
    if r_2 == 0,   U(k) = 0; end;
end;
newholx =  A(1,1)+ A(1,2)*newrefx+A(1,3)*newrefy+W(1,:)*U;
newholy =  A(2,1)+ A(2,2)*newrefx+A(2,3)*newrefy+W(2,:)*U;
