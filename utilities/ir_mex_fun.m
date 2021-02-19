function [fun] = ir_mex_fun()
% create mex files from variable # of c-files, should work for any OS

if ismac
    % Code to run on Mac platform
    fun = @(varargin) mex('-O', 'CFLAGS="\$CFLAGS -std=c99"', varargin{:});
elseif isunix
    % Code to run on Linux platform
    fun = @(varargin) mex('-O', 'CFLAGS="\$CFLAGS -std=c99"', varargin{:});
elseif ispc
    % Code to run on Windows platform
    fun = @(varargin) mex('-O', 'COMPFLAGS="\$COMPFLAGS -std=c99"', varargin{:});
else
    disp('Platform not supported')
end

end
