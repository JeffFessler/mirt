function[RR] = compute_RR(params)
%% Load params
rs = params.rs;
cs = params.cs;
Operator = params.Operator;
redundancy = params.Wavelet.redundancy;
includeApprox = params.Wavelet.includeApprox;
nlev = params.Wavelet.nlev;

%% Get the frequency responses of wavelet filters at various levels
[WRA] = compute_SWTFiltResp(params);
DCR = ones(rs,cs) - WRA(:,:,3*nlev+1);

%% Get the frequency response of FD filters
Rh = fft2([-1 1], rs, cs);
Rv = fft2([-1 1]', rs, cs);
Rabs = abs(Rh).^2 + abs(Rv).^2;

%% Compute W^T*W and R^T*R
switch(Operator)
    case{'FD'}
        RR = Rabs;
                        
    case{'W'}
        switch(redundancy)
            case{'undecimated'}
                if(includeApprox)
                    RR = 1;
                else
                    RR = DCR;
                end
        end
    case{'WFD'}
        switch(redundancy)
            case{'undecimated'}
                if(includeApprox)
                    RR = Rabs + 1;
                else
                    RR = Rabs + DCR;
                end
        end
end
