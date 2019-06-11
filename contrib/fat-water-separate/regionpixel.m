function [est_map,cost,iter] = regionpixel(data,arr,col,map,te,c,d,A1,A2,A3,niter,cal)
% data = source image
% arr, col = pixel position
% map = initial field map (=0) at (arr,col)th pixel
% te = time
% c,d = cos(F_fw*te), sin(F_fw*te)
% A1 = according to Appendix A
% A2 = inv(A1'*A1)*A1'
% A3 = [1 exp(j*F_fw*te)]
% max = maximum for delta_field map
% niter = iteration number

%% initialize

iter = 0;
nsets = size(data,3); % number of te  = 3
delta_map = 11;
if cal == 1
    % cost function at (arr,col)th pixel
    cost = zeros(2,niter*2); % cost(1,:) = field map, cost(2,:) = cost value
end
% source images at (arr,col)th pixel
orig_s = reshape(data(arr,col,:),nsets,1);


%% Iteration code

% while (abs(delta_map)>1 && iter<niter && abs(delta_map)<230)
while (abs(delta_map)>1 && iter<niter)
    
    iter = iter+1; % number of iteration
    
    %% Apendix A
    
    s =  orig_s.*(exp(-j*2*pi*map*te).'); % shat (m is freq)
    shat = [real(s);imag(s)]; 
    % shat = [s1_R s2_R s3_R s1_I s2_I s3_I]
    % A2 = inv(A1'*A1)*A1'
    rowhat = A2*shat; % rowhat vector = [water_R water_I fat_R fat_I]
  
   if cal == 1 
        %% cost function of Appendix A
        water = rowhat(1) + j*rowhat(2);
        fat = rowhat(3) + j*rowhat(4);
        cost(1,2*iter-1) = map;
        cost(2,2*iter-1) = norm(orig_s - A3*[water;fat].*(exp(j*2*pi*map*te).'));
        % cost(1,:) = field map values
        % cost(2,:) = cost function valuses
   end
    %% Apendix B
    
    % c_1n = 1, d_1n = 0 for all t_n, c_2n = c, d_2n = d
    f1 = -rowhat(1) -rowhat(3)*c +rowhat(4)*d;
    f2 = -rowhat(2) -rowhat(3)*d -rowhat(4)*c;
    
    shathat = shat + [f1;f2]; % according to  appendix B
    
    g = 2*pi*[te';te'].*[f2;-f1];
    
    B = [g A1];
    y = inv(B'*B)*B'*shathat; % y(1)is delta_fieldmap

    % update the fieldmap = prefieldmap + delta_field map
    delta_map = y(1);
    if abs(delta_map) > 30
        map = map + delta_map/abs(delta_map)*30; 
    else
        map = map + delta_map;
    end
    
    if cal == 1
        % cost function of Appendix B
        cost(1,2*(iter)) = map;
        cost(2,2*(iter)) = norm(orig_s - A3*[water;fat].*(exp(j*2*pi*map*te).'));
    end
end
if cal == 1
    est_map = cost(1,iter*2-1);
else
    est_map = map;
    cost = 0;
end
            
            
            