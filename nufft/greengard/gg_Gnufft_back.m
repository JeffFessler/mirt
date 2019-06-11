function Fk = gg_Gnufft_back(st,fj)

% apply j-dependent gaussian weights
V0 = fj.*st.E1;

% now do fast gridding
ftau = back_grid(complexify(V0),st.Msp,st.m1,st.m2,st.Mr,st.E2x,st.E2y,st.E3,size(st.om,1));

% take fft of ftau
Ftau = ifftshift(ifftn_fast(ftau));

% deconvolve to get F
Fk = reshape(pi/st.tau*st.E4.*col(Ftau(st.Mr/2-st.M/2+1:st.Mr/2+st.M/2,st.Mr/2-st.M/2+1:st.Mr/2+st.M/2)),st.M,st.M);

