function [cost, tE] = compute_Cost_SENSE(SP3, Data, x, params)
% Compute the cost function sum_i ||yi - TSi x||^2 + lambda * sum(rho(|grad(x)|))

%% Load parameters
lambda = params.lambda;

redundancy = params.Wavelet.redundancy;

Operator = params.Operator;
Prior = params.Prior.PriorType;

sr = params.AL.sr;
er = params.AL.er;

%% Compute the Data Term
[costDT, tE(1)] = compute_CostDataTerm_SENSE(Data, SP3, x, params);

%% Compute v (Refer notes)
[Rx, tE(2)] = doR(x, params); % First do Rx
switch(Operator)
    case{'FD'} % Finite difference
        tS = tic;
        Dtemp{1} = sqrt(sum(abs(Rx(:, :, sr(1):er(1))).^2,3));
        tE(3) = toc(tS);

    case{'W'} % Wavelets
        switch(redundancy)
            case{'undecimated'} % Shift-invariant wavelets
                tS = tic;
                Dtemp{1} = abs(Rx(:, :, sr(1):er(1)));
                tE(3) = toc(tS);
            otherwise
                error('Not coded for shift-variant wavelets');
        end
        
    case{'WFD'}
        tS = tic;
        Dtemp{1} = sqrt(sum(abs(Rx(:, :, sr(1):er(1))).^2,3));
        tE(3) = toc(tS);
        
        switch(redundancy)
            case{'undecimated'} % Shift-invariant wavelets
                tS = tic;
                Dtemp{2} = abs(Rx(:, :, sr(2):er(2)));
                tE(4) = toc(tS);
            otherwise
                error('Not coded for shift-variant wavelets');
        end
    otherwise
        error('Only available options are FD, W, and WFD');
end

%% Compute regularization corresponding to the concave prior 'rho' in the paper "Highly Undersampled MRI Reconstruction" by Tryasko et al, IEEE TMI
costRT = 0;
D = cell(1, length(Dtemp));
tE1 = zeros(1, length(Dtemp));
tE2 = zeros(1, length(Dtemp));
for iD = 1:length(Dtemp)
    switch (Prior{iD})
        case{'TV'}
            tS = tic;
            D{iD} = Dtemp{iD};
            tE1(iD) = toc(tS);

        case{'Log'}
            tS = tic;
            D{iD} = rho_Log(Dtemp{iD}, params);
            tE1(iD) = toc(tS);

        otherwise
            error('Unknown type of prior');
    end
    tS = tic;
    costRT = costRT + lambda(iD)*sum(D{iD}(:));
    tE2(iD) = toc(tS);
end

%% Total cost
tS = tic;
cost = costDT + costRT;
tE = toc(tS) + sum(tE(:)) + sum(tE1(:)) + sum(tE2(:));