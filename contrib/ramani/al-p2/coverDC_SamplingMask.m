function SP = coverDC_SamplingMask(SP, Ncentx, Ncenty)
%% Ensure a small neighbourhood of size 2N x 2N around origin is sampled

%% Origin at center of image
SP = fftshift(SP);
[rs, cs] = size(SP);

%% Center of image
if(mod(rs,2))
    rsby2 = (rs+1)/2;
else
    rsby2 = rs/2;
end

if(mod(cs,2))
    csby2 = (cs+1)/2;
else
    csby2 = cs/2;
end

%% Small Neighbourhood around origin sampled
SP(rsby2 - Ncentx +1 : rsby2 + Ncentx, csby2 - Ncenty + 1 : csby2 + Ncenty) = ones(2*Ncentx, 2*Ncenty);

%% Origin at (1,1) of the array
SP = fftshift(SP);
SP(1,1) = 1;