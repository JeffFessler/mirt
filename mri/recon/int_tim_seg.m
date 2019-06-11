function AA = int_tim_seg(tt,L,we,type, we_histo)
%function AA = int_tim_seg(tt,L,we,type, we_histo)
%  to be used with function fast_mr.m in @fast_mr
%  tt is time vector
%  L is number of segments, i.e. number of interpolation points
%  we is fieldmap, not required for linear and hanning types.
%  type is type of interpolator as follows:
%        1   is for ideal min-max using field map
%        2   is for histogram approximation, using we_histo
%  we_histo: only used when type 2 is selected
%            has two columns:
%                    column 1:  bin centers for histogram
%                    column 2:  histogram values at those bin centers


mint = min(tt(:));
maxt = max(tt(:));
rangt = maxt-mint;
minwe = min(we(:));
maxwe = max(we(:));
ndat = length(tt);
N = numel(we);

tau = (rangt+eps)/(L);

AA = zeros(L+1,ndat);

if L==0
  AA = ones(1,ndat); 
  return
end

tt = tt-mint;

if (type == 1) % Exact LS interpolator
    gg = exp(i*we(:)*tau)*ones(1,L);  %CHECK MINUS SIGN HERE
    lll = ones(size(gg(:,1)))*[1:L];
    gl = gg.^lll;
    G = [ones(size(gg(:,1))),gl];
    glsum = sum(gl);
    GTG = zeros(L+1,L+1);
    GTG = GTG+diag(N*ones(L+1,1),0);
    for kk = 1:L
       GTG = GTG+diag(glsum(kk)*ones(L+1-kk,1),kk);
       GTG = GTG+diag(conj(glsum(kk))*ones(L+1-kk,1),-kk);
    end
       if (rcond(GTG)>10*eps)
           iGTGGT = inv(GTG)*G';
       else
           iGTGGT = pinv(GTG)*G';
           sprintf('used pinv instead')
      end 
 %  for yy = 1:ndat
           cc = exp(i*we(:)*(tt(:))');    
          % cc = exp(i*we(:)*tt(yy));    
          % AA(:,yy) = (iGTGGT*cc)'.';
           AA = (iGTGGT*cc)'.';
 %   end
end   
if (type == 2) % Approx Histogram Interpolator using we_histo
     %find bin size to give an integer value for KK
     rangwe = max(we_histo(:,1))-min(we_histo(:,1)); 
     minwe = min(we_histo(:,1));
     num_bins = length(we_histo(:,1));
     KK = floor(2*pi/((rangwe/num_bins)*tau));
     % bin size
     dwn = 2*pi/(KK*tau);
     bin_centers = we_histo(:,1); 
     N_ap = we_histo(:,2).';
     %keyboard
     ftwe_ap = fft(N_ap.*(exp(i*dwn*[0:(num_bins-1)]*tau*(L))),KK);
     ftwe_ap = exp(-i*(minwe+dwn/2)*tau*([0:(KK-1)]-(L))).*ftwe_ap;
     GTGap_ap = zeros(L+1,L+1);
     for kk = 1:(2*L+1)
        GTGap_ap = GTGap_ap+diag(ftwe_ap(kk)*ones(L+1-abs((L+1)-kk),1),-(L+1)+kk);
     end
       if (rcond(GTGap_ap.')>10*eps)
           iGTGap_ap = inv(GTGap_ap.');
       else
           iGTGap_ap = pinv(GTGap_ap.');
           sprintf('used pinv instead')
      end 
     for yy = 1:ndat
        ftc_ap= fft(N_ap.*(exp(i*[0:(num_bins-1)]*dwn*tt(yy))),KK);
        ftc_ap = exp(i*(minwe+dwn/2)*(tt(yy)-tau*[0:KK-1])).*(ftc_ap);
        GTc_ap = ftc_ap(1:L+1).';
        AA(:,yy) = (iGTGap_ap*GTc_ap)'.';
     end
end
