function a = myiswt2(SWC, lo_R, hi_R)
%ISWT2 Inverse discrete stationary wavelet transform 2-D.
%   ISWT2 performs a multilevel 2-D stationary wavelet 
%   reconstruction using either a specific orthogonal wavelet  
%   ('wname', see WFILTERS for more information) or specific
%   reconstruction filters (Lo_R and Hi_R).
%
%   Filters lo_R and hi_R must be normalized by dividing by sqrt(2) (necessary for implmentation of Stationary Discrete WT)
%
%   For X = MYISWT2(SWC,Lo_R,Hi_R)
%   SWC contains H(1:N) V(1:N) D(1:N) A(N)
%   Lo_R is the reconstruction low-pass filter.
%   Hi_R is the reconstruction high-pass filter.
%
%   See also IDWT2, SWT2, WAVEREC2.

%   Sathish Ramani, October 11, 2009

%% Set DWT_Mode to 'per'.
old_modeDWT = dwtmode('status','nodisp');
modeDWT = 'per';
dwtmode(modeDWT,'nodisp');

%% Load coefficients
n = (size(SWC, 3)-1)/3; %% Number of levels
h = SWC(:,:,1:n);
v = SWC(:,:,n+1:2*n);
d = SWC(:,:,2*n+1:3*n);
a = SWC(:,:,3*n+1);

[rx, cx, dump] = size(h);
for k = n:-1:1
    step = 2^(k-1);
    last = step;
    for first1 = 1:last
        iRow = first1:step:rx;
        lR   = length(iRow);
        for first2 = 1:last
            iCol = first2:step:cx;
            lC   = length(iCol);
            sR   = iRow;
            sC   = iCol;
            a(iRow,iCol) = idwt2LOC(a(sR,sC),h(sR,sC,k),v(sR,sC,k),d(sR,sC,k), lo_R,hi_R, [lR lC]);
        end
    end
end

% Restore DWT_Mode.
dwtmode(old_modeDWT,'nodisp');


%===============================================================%
% INTERNAL FUNCTIONS
%===============================================================%
function y = idwt2LOC(a,h,v,d,lo_R,hi_R,sy)

y = upconvLOC(a,lo_R,lo_R,sy)+ ... % Approximation.
    upconvLOC(h,hi_R,lo_R,sy)+ ... % Horizontal Detail.
    upconvLOC(v,lo_R,hi_R,sy)+ ... % Vertical Detail.
    upconvLOC(d,hi_R,hi_R,sy);     % Diagonal Detail.
%---------------------------------------------------------------%
function y = upconvLOC(x,f1,f2,s)

lf = length(f1);
% y  = dyadup(x,'mat',0,1);
y = x;
y  = wextend('2D','per',y,[lf/2,lf/2]);
y  = wconv2('col',y,f1);
y  = wconv2('row',y,f2);
y  = wkeep2(y,s,[lf lf]);
%===============================================================%
