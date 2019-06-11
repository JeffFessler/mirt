function [Geometry,Geom_res] =fbkp_grid(theta,N,interp)

% theta are imput as angles in degrees
% Define the x & y axes for the reconstructed image so that the origin
% (center) is in the spot which RADON would choose.

%   N is a scalar that specifies the number of rows and columns in the 
%   reconstructed image. N is determined from the length of the projections: !!!!
%
%       N = 2*floor(size(P,1)/(2*sqrt(2)))
%
theta = pi*theta/180;

xax = (1:N)-ceil(N/2);
x = repmat(xax,N,1) ; 
% x coordinates, the y coordinates are rot90(x)
y = rot90(x);
costheta = cos(theta);
sintheta = sin(theta);

% Backprojection - vectorized in (x,y), looping over theta
sa=size(x,1); sb=size(x,2); sc=length(theta);


if strcmp(interp, 'nearest neighbor')
        Geometry=repmat(int16(0),[sa,sb,sc]);
        Geom_res=[];
    for i=1:length(theta)
        t=[];
        t = round( x*costheta(i) + y*sintheta(i) );
        Geometry(:,:,i)=t;
    end
    
elseif strcmp(interp, 'linear')
    Geometry=repmat(int16(0),[sa,sb,sc]); % parte intera
    Geom_res=repmat(uint8(0),[sa,sb,sc]); % parte decimale
    
    for i=1:length(theta)  
        t = x.*costheta(i) + y.*sintheta(i); 
        a = floor(t);  
        % img = img + (t-a).*proj(a+1+ctrIdx) + (a+1-t).*proj(a+ctrIdx);
        Geometry(:,:,i)=a ; % floor
        Geom_res(:,:,i)=uint8( abs(t-a)*100) ; % decimal i.e. double - floor

    end
end
