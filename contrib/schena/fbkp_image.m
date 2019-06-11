function [img] =fbkp_image(p,theta,H,T,T_decimal,N,interp)
%fbkp_image computes inverse Radon transform.
%   I = fbkp_image(p,theta) reconstructs the image I from projection 
%   data in the 2-D array p.  The columns of R are parallel beam 
%   projection data.  
%   fbkp_image assumes that the center of rotation
%   is the center point of the projections, which is defined as
%   ceil(size(p,1)/2).
%
%   THETA describes the angles (in degrees) at which the projections    
%   were taken. THETA is a vector, it must contain angles with 
%   equal spacing between them.  
%   H is the filter
%   T is the geometry
% 
%   N is a scalar that specifies the number of rows and columns in the 
%   reconstructed image.  If N is not specified, the size is determined   
%   from the length of the projections:
%
%       N = 2*floor(size(p,1)/(2*sqrt(2)))

len=size(p,1);   
p(length(H),1)=0;  % Zero pad projections 

% In the code below, I continuously reuse the array p so as to
% save memory.  This makes it harder to read, but the comments
% explain what is going on.

p = fft(p);    % p holds fft of projections

for i = 1:size(p,2)
   p(:,i) = p(:,i).*H; % frequency domain filtering
end

p = real(ifft(p));     % p is the filtered projections
p(len+1:end,:) = [];   % Truncate the filtered projections
img = zeros(N);        % Allocate memory for the image.

% Define the x & y axes for the reconstructed image so that the origin
% (center) is in the spot which RADON would choose.

ctrIdx = ceil(len/2) ;  % index of the center of the projections

% Zero pad the projections to size 1+2*ceil(N/sqrt(2)) if this
% quantity is greater than the length of the projections
imgDiag = 2*ceil(N/sqrt(2))+1 ; % largest distance through image.
if size(p,1) < imgDiag 
   rz = imgDiag - size(p,1) ; % how many rows of zeros
   p = [zeros(ceil(rz/2),size(p,2)); p; zeros(floor(rz/2),size(p,2))];
   ctrIdx = ctrIdx+ceil(rz/2)
end
% Backprojection - vectorized in (x,y), looping over theta
   
if strcmp(interp, 'nearest neighbor')

   for i=1:length(theta) 
      proj = p(:,i);
    %  t = round(x*costheta(i) + y*sintheta(i));
      t=T(:,:,i); t=double(t); % t is originaly INT16 for economy of memory
       img = img + proj(t+ctrIdx);
      
  end
elseif strcmp(interp, 'linear')

      for i=1:length(theta) 
    proj = p(:,i);
      %t = x.*costheta(i) + y.*sintheta(i); 
      a = T(:,:,i); a=double(a); % parte intera  (a = floor(t)in Geometry)
      t= a + double( T_decimal(:,:,i) )/100.;
      img = img + (t-a).*proj(a+1+ctrIdx) + (a+1-t).*proj(a+ctrIdx);
   end
   
end

img = img*pi/(2*length(theta));
