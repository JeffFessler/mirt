function [WRA, WF, WR] = compute_SWTFiltResp(params)
%% Load parameters
rs = params.rs;
cs = params.cs;
nlev = params.Wavelet.nlev;

lod = params.Wavelet.lod; %% Normalize filter for easier SWT implementation
hid = params.Wavelet.hid; %% Normalize filter for easier SWT implementation

% Get the shift-invariant wavelet filters
indicator = zeros(rs, cs);
indicator(1,1) = 1;
SWC = myswt2(indicator, nlev, lod, hid);

WF = zeros(rs, cs, 3*nlev+1);
WR = WF;
WRA = WR;

for ib = 1:3
    for ilev = 1:nlev
        WF(:,:,(ib-1)*nlev+ilev) = SWC(:,:,(ib-1)*nlev+ilev);
        WR(:,:,(ib-1)*nlev+ilev) = fft2(WF(:,:,(ib-1)*nlev+ilev));
        WRA(:,:,(ib-1)*nlev+ilev) = abs(WR(:,:,(ib-1)*nlev+ilev)).^2;
    end
end
WF(:,:,3*nlev+1) = SWC(:,:,3*nlev+1);
WR(:,:,3*nlev+1) = fft2(WF(:,:,3*nlev+1));
WRA(:,:,3*nlev+1) = abs(WR(:,:,3*nlev+1)).^2;

% DCR = 0;
% for ib = 1:3*nlev
%     DCR = DCR + WRA(:,:,ib);
% end
% dif = DCR-ACR;
% max(abs(dif(:)))
