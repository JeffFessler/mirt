% Demonstrated Grid2Glinear routine
% Demostrates how to go from the traditional back-projection Grid matrix used 
% by iradon and other FBP routines to the G sparse matrix used by the iterative methods
% i.e. - reconciliates FBP and ART methods
% uses iradon and imshow then needs Image Proc toolbox 

clear all; close all;

% read phantom 
P = phantom(64);
theta=[0:2:179]; 
np=length(theta);


% generates projection data
R = radon(P,theta);

%interp='nearest neighbor';
interp='linear';
flt='shepp-logan' ; % The Shepp-Logan filter 
%flt='ram-lak'

tic
% IRADON !
[I,H1] = iradon(R,[theta],interp,flt);
toc
imshow(P); title(' original phantom')
figure; imshow(I,[]); title('iradon')

% FBP STEP BY STEP

% check sinogram size with required image size and pad with 0 if necessary
len=size(R,1); ctrIdx = ceil(len/2)  ;  
N = 2*floor(size(R,1)/(2*sqrt(2))); % default size
imgDiag = 2*ceil(N/sqrt(2))+1 ; % largest distance through image.
if size(R,1) < imgDiag 
    rz = imgDiag - size(R,1) ; % how many rows of zeros
    R = [zeros(ceil(rz/2),size(R,2)); p; zeros(floor(rz/2),size(R,2))];
    ctrIdx = ctrIdx+ceil(rz/2);
    display('pad necessary ! ')
end

len=size(R,1); 

% Design filter for the sinogram 
d=1;   
H =fbkp_filter_def(flt, len, d);


% back projection Grid
Grid=[];Grid_dec=[]; % floor part and decimal part to save memory 
[Grid,Grid_dec] =fbkp_grid([theta],N,interp);

% FBP into the Grid
tic
[Img] =fbkp_image(R,[theta],H,Grid,Grid_dec,N,interp);
toc

figure; imshow(Img,[]); title('ricostructed by Fbkp');

nx=N, ny=N, , 	nb = size(R, 1),
ray_pix=1 , 

% ART
% mask
%mask=logical(ones(nx,ny)); % defaul mask
[mask,nnz_mask] = circle(N); % simple circle mask

% FBP reconstructed image taken as a guess for iterative methods ? ... 
Iguess=[ ]; % NO
%Iguess=Img(:); Iguess=Iguess(mask); % SI

% see memory available to you for number of projection to be used each time !!!!
% with large memory you can use: 
%szblk=np; %
szblk=1, % min memory req. but poor solution

szblk=3; % 

np=length(theta),
iteraz=5;

for block=[1:szblk:np];
   szblk=min(szblk, (np-block+1) );
  % transforms the Grid into the associated sparse matrix G
    G = Grid2GlinearROI(Grid,Grid_dec,nb,ctrIdx,block,interp,szblk,mask(:)',nnz_mask);
    % setting up an under-determined system of equations 
    
    % solving the system in the LS sense
    yi=R(:, [block: min(block+szblk-1,np)] );
    
    % no regularization for this simple jet simplistic example
    I=lsqr(G,yi(:),[],iteraz,[],[],Iguess(:)); % use J.Fessler pcg routines for a faster solution
    % use Fessler C2sparse for regularisation
    Iguess=I;
end
% show results
I = embed(I(:,end), mask); % embed by J. Fessler
figure, imshow(I,[]); title('lsqr & Grid2Glinear');

return
