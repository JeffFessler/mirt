function [Rz, tE] = doRAdj(z, params)
% Do R*z, where R is a linear operator that signifies a wavelet transform or finite differences

%% Load parameters
rs = params.rs;
cs = params.cs;

nlev = params.Wavelet.nlev;
lor = params.Wavelet.lor;
hir = params.Wavelet.hir;
redundancy = params.Wavelet.redundancy;
includeApprox = params.Wavelet.includeApprox;

Operator = params.Operator;

%% Compute v (Refer notes)
switch(Operator)
    case{'FD'} % Finite difference
        [h, v, tE(1)] = trfd(z(:, :, 1), z(:, :, 2));
        Rz = h + v;
        
    case{'W'} % Wavelets
        switch(redundancy)
            case{'undecimated'} % Shift-invariant wavelets
                if(includeApprox)
                    tS = tic;
                    z1 = z;
                    tE(1) = toc(tS);
                else
                    tS = tic;
                    z1 = z;
                    z1(:, :, end+1) = zeros(rs, cs);
                    tE(1) = toc(tS);
                end
                tS = tic;
                Rz = myiswt2(z1, lor, hir);
                tE(2) = toc(tS);
            otherwise
                error('Not coded for shift-variant wavelets');
        end
        
    case{'WFD'}
        [h, v, tE(1)] = trfd(z(:, :, 1), z(:, :, 2));
        Rz = h + v;
        
        switch(redundancy)
            case{'undecimated'} % Shift-invariant wavelets
                if(includeApprox)
                    tS = tic;
                    z1 = z(:, :, 3:3*nlev+3);
                    tE(2) = toc(tS);
                else
                    tS = tic;
                    z1 = z(:, :, 3:3*nlev+2);
                    z1(:, :, end+1) = zeros(rs, cs);
                    tE(2) = toc(tS);
                end
                tS = tic;
                Rz = Rz + myiswt2(z1, lor, hir);
                tE(3) = toc(tS);
            otherwise
                error('Not coded for shift-variant wavelets');
        end
                
    otherwise
        error('Only available options are FD, W, and WFD');
end
tE = sum(tE(:));