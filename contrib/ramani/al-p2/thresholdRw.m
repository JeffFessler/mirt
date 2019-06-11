function [v, tE] = thresholdRw(Rw, params, mu)
%% Load Parameters
sr = params.AL.sr;
er = params.AL.er;
Operator = params.Operator;
redundancy = params.Wavelet.redundancy;

Prior = params.Prior.PriorType;
lambda = params.lambda;

%% Perform thresholding
switch(Operator)
    case{'FD'} % Finite difference
        tS = tic;
        Dtemp = repmat(sqrt(sum(abs(Rw(:, :, sr(1):er(1))).^2,3)), [1 1 er(1)-sr(1)+1]);
        tE(1) = toc(tS);

    case{'W'} % Wavelets
        switch(redundancy)
            case{'undecimated'} % Shift-invariant wavelets
                tS = tic;
                Dtemp = abs(Rw(:, :, sr(1):er(1)));
                tE(1) = toc(tS);
            otherwise
                error('Not coded for shift-variant wavelets');
        end

    case{'WFD'}
        tS = tic;
        Dtemp(:, :, sr(1):er(1)) = repmat(sqrt(sum(abs(Rw(:, :, sr(1):er(1))).^2,3)), [1 1 er(1)-sr(1)+1]);
        tE(1) = toc(tS);

        switch(redundancy)
            case{'undecimated'} % Shift-invariant wavelets
                tS = tic;
                Dtemp(:, :, sr(2):er(2)) = abs(Rw(:, :, sr(2):er(2)));
                tE(2) = toc(tS);
            otherwise
                error('Not coded for shift-variant wavelets');
        end
        
    otherwise
        error('Only available options are FD, W, and WFD');
end

%% Compute diagonal weighting matrix corresponding to the concave prior 'rho' in the paper "Highly Undersampled MRI Reconstruction" by Tryasko et al, IEEE TMI
tE1 = zeros(1, length(sr));
for iD = 1:length(sr)
    switch (Prior{iD})
        case{'TV'} % l1 norm of the sparsifying transform coefficients
            tS = tic;
            D = Dtemp(:, :, sr(iD):er(iD))-lambda(iD)/mu;
            D = (D + abs(D))/2;
            v(:, :, sr(iD):er(iD)) = Rw(:, :, sr(iD):er(iD))./Dtemp(:, :, sr(iD):er(iD)).*D;
            tE1(iD) = toc(tS);
            
        otherwise
            error('Unknown type of prior');
    end
end

%% Ensure v does not contain NaN obtained from 0/0 operations
v(isnan(v)) = 0;
tE = sum(tE(:)) + sum(tE1(:));
