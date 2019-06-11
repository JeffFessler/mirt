function fj = gg_Gnufft_forw(st,Fk)

% first deconvolve
Fmintau = pi/st.tau*st.E4.*Fk;

% then take iDFT to get fmintau
padval = (st.Mr-st.M)/2;
fmintau = fftn_fast(fftshift(pad(reshape(Fmintau,st.M,st.M),padval, ...
                              padval,padval,padval,0)))/st.Mr/st.Mr;

% now do fast gridding
fj = forw_grid(fmintau,st.Msp,st.m1,st.m2,st.Mr,st.E2x,st.E2y,st.E3,size(st.om,1));

% now apply final gaussian weights
fj = fj.*st.E1;
