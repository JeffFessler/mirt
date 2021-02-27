function [data] = ir_get_data(filename_full, varargin)
% ir_get_data either returns the data at the specified filename, or grabs
% if from the data/downloads folder. ir_get_data wil download the data from
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
data_dir = fileparts(which('ir_get_data.m'));
fn_parts = strsplit(filename_full,data_dir);
sub_dir_fn = fn_parts{2};

sub_dir_fn_parts = strsplit(sub_dir_fn,'/');
filename = sub_dir_fn_parts{end};


%% Download data if it doesn't exist on the local machine
if ~exist(filename_full,'file')
    
    disp(['The file ', filename_full])
    disp(['was not found on local machine, downloading from',...
        ' https://github.com/JeffFessler/MIRTdata/'])
    
    % create new directory if needed
    filename_dir_parts = strsplit(filename_full,filename);
    filename_dir = filename_dir_parts{1};
    if ~exist(filename_dir,'dir')
        mkdir(filename_dir)
    end
    
    % download data from web
    url_start = 'https://github.com/JeffFessler/MIRTdata/raw/main';
    full_url = [url_start, sub_dir_fn];
    websave(filename_full,full_url);
end


%% Get type of file and load data
[~, ~, fExt] = fileparts(filename);
switch lower(fExt)
    case '.mat'
        try
            if arg.pick
                data = ir_read_mat(filename_full, arg.pick);
            else
                data = ir_read_mat(filename_full);
            end

            return
        catch
            warning(['Warning! Could not load ',filename_full])
        end        
    case '.fld'
        try
            data = fld_read(filename_full);
            return
        catch
            warning(['Warning! Could not load ',filename_full])
        end        
    otherwise 
        error('Unexpected file extension: %s', fExt);
end



end

