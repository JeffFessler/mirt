function xinew = zero_pad_pow2(xi);
% Pad the given slices using mirror boundary conditions to size them up to a multiple of power of 2
[rs cs ncoils] = size(xi);

%% Size of stack and Round rs and cs to nearest integer that is divisible by a power of 2 greater than 1
powcx = floor(log2(rs)) + 1;
powvx = floor(log2(cs)) + 1;

modcx = mod(rs, 2.^[1:powcx]);
modvx = mod(cs, 2.^[1:powvx]);

powcx = min(find(modcx>0))-1;
powvx = min(find(modvx>0))-1;

if(powcx <= 1)
    if(mod(rs, 2))
        rsnew = rs + 2 + modcx(2);
    else
        rsnew = rs + min(modcx(modcx>0));
    end
else
    rsnew = rs;
end

if(powvx <= 1)
    if(mod(cs, 2))
        csnew = cs + 2 + modvx(2);
    else
        csnew = cs + min(modvx((modvx>0)));
    end
else
    csnew = cs;
end

dcleft = floor((rsnew-rs)/2);
dcright = rsnew - rs - dcleft;

dvtop = floor((csnew-cs)/2);
dvbot = csnew - cs - dvtop;

if(max([dcleft dcright dvtop dvbot]) == 0)
    disp('Size is appropriate, no need for zero-padding');
    xinew = xi;
else
    disp(['Zero-padding with ' int2str(dcleft) ' ' int2str(dcright) ' ' int2str(dvtop) ' ' int2str(dvbot)]);
    for icoil = 1:ncoils
        xir = xi(:, :, icoil);
        
        % Pad the images to match new image dimensions
        xirnew = zeros(rsnew, csnew);
        
        xirnew(1:dcleft, 1:dvtop) = 0;           % Left - Top
        xirnew(1:dcleft, dvtop+1:cs+dvtop) = 0;  % Left - Middle
        xirnew(1:dcleft, cs+dvtop+1:csnew) = 0;  % Left - Bottom
        
        xirnew(dcleft+1:rs+dcleft, 1:dvtop) = 0;             % Middle - Top
        xirnew(dcleft+1:rs+dcleft, dvtop+1:cs+dvtop) = xir;  % Middle - Middle
        xirnew(dcleft+1:rs+dcleft, cs+dvtop+1:csnew) = 0;    % Middle - Bottom
        
        xirnew(rs+dcleft+1:rsnew, 1:dvtop) = 0;           % Right - Top
        xirnew(rs+dcleft+1:rsnew, dvtop+1:cs+dvtop) = 0;  % Right - Middle
        xirnew(rs+dcleft+1:rsnew, cs+dvtop+1:csnew) = 0;  % Right- Bottom
        
        xinew(:, :, icoil) = xirnew;
    end
end
