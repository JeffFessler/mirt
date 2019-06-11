 function x = read_zubal_emis(varargin)
%function x = read_zubal_emis(options)
%| read in zubal emission phantom from data directory
%| options
%|	'nx'		desired size
%|	'ny'
%|	'ddir'		data directory

if nargin==1 && streq(varargin{1}, 'test'), read_zubal_emis_test, return, end

if nargout == 0, help(mfilename), error(mfilename), end

arg.ddir = '';
arg.file = 'zubal,emis.raw';
arg.nx = 128;
arg.ny = [];
arg = vararg_pair(arg, varargin);
if isempty(arg.ny), arg.ny = arg.nx; end

x = read_zubal_emis_do(arg.ddir, arg.file, arg.nx, arg.ny);


function x = read_zubal_emis_do(ddir, file, nx, ny)

% guess .../data directory by looking parallel to '.../emission' directory
if ~isvar('ddir') || isempty(ddir)
	t = path_find_dir([filesep 'emission']);
	ddir = strrep(t, 'emission', 'data');
end

file = [ddir filesep file];
if ~exist(file, 'file')
%	os_run(sprintf([ddir filesep 'do,emis,zubal %s'], file)) % here is how to create
	error 'cannot find zubal phantom raw data'
end

fp = fopen(file, 'rb');
if (fp == -1), error 'open file', end
x = fread(fp, [128 128], 'uint8');
if fclose(fp), error 'close file', end

x = x(:,[(end-10):end 1:(end-11)]); % center it nicely

x = ir_phantom_resize(x, nx, ny);


function read_zubal_emis_test
x = read_zubal_emis('nx', 120, 'ny', 90);
im(x, 'zubal phantom: emission'), cbar
