function [Rz, tE] = doR(z, params)
% Do R*z, where R is a linear operator that signifies a wavelet transform or finite differences

%% Load parameters
nlev = params.Wavelet.nlev;
lod = params.Wavelet.lod;
hid = params.Wavelet.hid;
redundancy = params.Wavelet.redundancy;
includeApprox = params.Wavelet.includeApprox;

Operator = params.Operator;

%% Compute v (Refer notes)
switch(Operator)
    case{'FD'} % Finite difference
        [h, v, tE(1)] = fd(z);
        Rz(:, :, 1) = h;
        Rz(:, :, 2) = v;
        
    case{'W'} % Wavelets
        switch(redundancy)
            case{'undecimated'} % Shift-invariant wavelets
                tS = tic;
                SWC = myswt2(z, nlev, lod, hid);
                tE(1) = toc(tS);
                if(includeApprox)
                    tS = tic;
                    Rz = SWC;
                    tE(2) = toc(tS);
                else
                    tS = tic;
                    Rz = SWC(:,:,1:3*nlev);
                    tE(2) = toc(tS);
                end
            otherwise
                error('Not coded for shift-variant wavelets');
        end
        
    case{'WFD'}
        % And finite difference
        [h, v, tE(1)] = fd(z);
        Rz(:, :, 1) = h;
        Rz(:, :, 2) = v;
        
        switch(redundancy)
            case{'undecimated'} % Shift-invariant wavelets
                tS = tic;
                SWC = myswt2(z, nlev, lod, hid);
                tE(2) = toc(tS);
                if(includeApprox)
                    tS = tic;
                    Rz(:, :, 3:3*nlev+3) = SWC;
                    tE(3) = toc(tS);
                else
                    tS = tic;
                    Rz(:, :, 3:3*nlev+2) = SWC(:,:,1:3*nlev);
                    tE(3) = toc(tS);
                end
            otherwise
                error('Not coded for shift-variant wavelets');
        end
    otherwise
        error('Only available options are FD, W, and WFD');
end
tE = sum(tE(:));