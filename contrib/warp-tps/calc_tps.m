function [L, W, A, BE] = calc_tps(refx, refy, holx, holy, calc_BE) %coord in column vector style
%calculate 2D tps warp
n = length(refx);
P = [ ones(n,1) refx refy]; 
V = [holx'; holy'];
K = zeros(n,n);
BE = 0 ;
for i=1:n,
    for j=1:n,
        r_2= (refx(i) - refx(j))^2 + (refy(i) -refy(j))^2 ;
        if r_2 > 0, K(i,j) = r_2 * log(r_2); end;
    end;
end;
for i=1:n, K(i,i)=0; end;
L = [K P;P' zeros(3,3)];
W=zeros(2,n);
A=zeros(2,3);
Y=[V zeros(2,3)]';
temp = ( inv(L)*Y )';
W = temp(1:2,1:n);
A = temp(1:2,n+1:n+3);
if calc_BE == 0, return; end;

temp = zeros(n,n); % server as Ln-1
temp1 = inv(L);
Ln_inv = temp1(1:n,1:n);
[bas_u, s, bas_v] = svd(Ln_inv);%col of bas_v basis vectors
B = zeros(2,n-3); %storage of coefs of decomp of V
for i=1:n-3,
    B(1,i) = V(1,:)*bas_v(:,i);
    B(2,i) = V(2,:)*bas_v(:,i);
end;
BE = 0;
for i=1:n-3,
    BE = BE + s(i,i)*( B(1,i)^2 + B(2,i)^2 );
end;
