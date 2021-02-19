function [fun] = ir_mex_fun2()
% two input version of ir_mex_fun()

if ismac
    % Code to run on Mac platform
    fun = @(f1, f2) mex('-O', 'CFLAGS="\$CFLAGS -std=c99"', f1, f2);
elseif isunix
    % Code to run on Linux platform
    fun = @(f1, f2) mex('-O', 'CFLAGS="\$CFLAGS -std=c99"', f1, f2);
elseif ispc
    % Code to run on Windows platform
    fun = @(f1, f2) mex('-O', 'COMPFLAGS="\$COMPFLAGS -std=c99"', f1, f2);
else
    disp('Platform not supported')
end

end
