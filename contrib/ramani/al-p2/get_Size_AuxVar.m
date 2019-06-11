function[lzALP2 sALP2 eALP2 sr er] = get_Size_AuxVar(params)
% Get the size of auxiliary variables
% sALP2 - stores starting indices of parts of auxiliary variables
% eALP2 - stores ending  indices of parts of auxiliary variables
% sr - stores starting indices of parts of regularization variables
% er - stores starting indices of parts of regularization variables
% lz - total number of auxiliary variables

%% Load parameters
ncoils = params.ncoils;

nlev = params.Wavelet.nlev;
redundancy = params.Wavelet.redundancy;
includeApprox = params.Wavelet.includeApprox;

Operator = params.Operator;

%% Number of variables for dealing with weights, i.e., u0 = Sx in ALP2
sALP2(1) = 1;
eALP2(1) = sALP2(1) + ncoils - 1;

%% Number of variables for the constraint u2 = x in ALP2
sALP2(2) = eALP2(1) + 1;
eALP2(2) = sALP2(2);

%% Number of variables for regularization, i.e., u1 = Rx in ALP1 and ALP2
switch(Operator)
    case{'FD'} % Finite difference
        sALP2(3) = eALP2(2) + 1;
        eALP2(3) = sALP2(3) + 1;
               
    case{'W'} % Wavelets
        sALP2(3) = eALP2(2) + 1;
        
        switch(redundancy)
            case{'undecimated'} % Shift-invariant wavelets
                if(includeApprox)
                    eALP2(3) = sALP2(3) + 3*nlev;
                else
                    eALP2(3) = sALP2(3) + 3*nlev-1;
                end
            otherwise
                error('Not coded for shift-variant wavelets');
        end
        
    case{'WFD'}
        % For finite difference
        sALP2(3) = eALP2(2) + 1;
        eALP2(3) = sALP2(3) + 1;
        
        % For Wavelets
        sALP2(4) = eALP2(3) + 1;
        
        switch(redundancy)
            case{'undecimated'} % Shift-invariant wavelets
                if(includeApprox)
                    eALP2(4) = sALP2(4) + 3*nlev;
                else
                    eALP2(4) = sALP2(4) + 3*nlev-1;
                end
            otherwise
                error('Not coded for shift-variant wavelets');
        end
    otherwise
        error('Only available options are FD, W, and WFD');
end
sALP1 = sALP2(3:end)-sALP2(3)+1;
eALP1 = eALP2(3:end)-sALP2(3)+1;

sADMM = sALP2;
eADMM = eALP2;

sr = sALP1;
er = eALP1;

lzADMM = eADMM(end) - sADMM(1) + 1;
lzALP2 = eALP2(end) - sALP2(1) + 1;
lzALP1 = eALP1(end) - sALP1(1) + 1;