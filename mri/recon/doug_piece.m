addpath('/net/waters/home/bpsutton/Reconitmr');

%Need we1 in rad/s
% fov = 22?;
% N = 64?;

          [bin_vals, bin_cens] = hist(we1(:),512);
           we_histo = [bin_cens', bin_vals'];
 %Af = fast_mr(kx,ky,fov,N,K,J,tt,we,flag_swt,L,int_opt,we_histo)
Af = fast_mr(kx,ky,fov,N,2*N,6,t,we1,0,6,2,we_histo);