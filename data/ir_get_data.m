function [data] = ir_get_data(data_filename, varargin)
% ir_get_data either returns the data at the specifiled filename, or grabs
% it from the data/downloads folder. ir_get_data wil download the data from
% https://github.com/JeffFessler/MIRTdata/ using websave if needed.
%
% Optional inputs:
%   pick - optional second input to ir_read_mat
%
% 2021 - Melissa Haskell, mhask@umich.edu

%% Set defaults and parse variable input arguments
arg.pick = [];
arg = vararg_pair(arg, varargin);


%% Breakdown full filename into parts
data_dir = [fileparts(which('ir_get_data.m')),filesep,'downloads',filesep];
full_fn = [data_dir, data_filename];

sub_dir_fn_parts = strsplit(full_fn,filesep);
filename = sub_dir_fn_parts{end};


%% Download data if it doesn't exist on the local machine
if ~exist(full_fn,'file')
    
    disp(['The file ', full_fn])
    disp(['was not found on local machine, downloading from',...
        ' https://github.com/JeffFessler/MIRTdata/'])
    
    % create new directory if needed
    filename_dir_parts = strsplit(full_fn,filename);
    filename_dir = filename_dir_parts{1};
    if ~exist(filename_dir,'dir')
        mkdir(filename_dir)
    end
    
    % download data from web
    url_start = 'https://github.com/JeffFessler/MIRTdata/raw/main/';
    full_url = [url_start, strrep(data_filename,filesep,'/')];
    if ispc, full_url = strrep(full_url,',','%'); end
    websave(full_fn,full_url);
end


%% Get type of file and load data
[~, ~, fExt] = fileparts(filename);
switch lower(fExt)
    case '.mat'
        try
            if arg.pick
                data = ir_read_mat(full_fn, arg.pick);
            else
                data = ir_read_mat(full_fn);
            end

            return
        catch
            warning(['Warning! Could not load ',full_fn])
        end        
    case '.fld'
        try
            data = fld_read(full_fn);
            return
        catch
            warning(['Warning! Could not load ',data_filename])
        end        
    otherwise 
        error('Unexpected file extension: %s', fExt);
end

end
