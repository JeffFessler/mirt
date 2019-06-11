 function x = read_zubal_attn(varargin)
%function x = read_zubal_attn(options)
%| read in zubal attenuation phantom from data directory,
%| and assign it attenuation coefficients in inverse mm units.
%| options
%|	'nx'		desired size
%|	'ny'
%|	'ddir'		data directory

if nargin==1 && streq(varargin{1}, 'test'), read_zubal_attn_test, return, end

if nargout == 0, help(mfilename), error(mfilename), end

arg.ddir = '';
arg.file = 'zubal,attn.raw';
arg.nx = 128;
arg.ny = [];
arg = vararg_pair(arg, varargin);
if isempty(arg.ny), arg.ny = arg.nx; end

x = read_zubal_attn_do(arg.ddir, arg.file, arg.nx, arg.ny);


function x = read_zubal_attn_do(ddir, file, nx, ny)

% guess .../data directory by looking parallel to '.../transmission' directory
if ~isvar('ddir') || isempty(ddir)
	t = path_find_dir([filesep 'transmission']);
	ddir = strrep(t, 'transmission', 'data');
end

if ~exist(ddir, 'dir')
	warning(sprintf('cannot find data directory %s', ddir))
	error(sprintf('edit path in %s.m for your installation', mfilename))
end

file = [ddir filesep file];
if ~exist(file, 'file')
%	os_run(sprintf([ddir filesep 'do,attn,zubal %s'], file)) % here is how to create!
	error 'cannot find zubal phantom raw data'
end

fp = fopen(file, 'rb');
if (fp == -1), error 'open file', end
x = fread(fp, [128 128], 'uint8');
if fclose(fp), error 'close file', end

x = x(:,[(end-10):end 1:(end-11)]); % center it nicely

% assign attenuation coefficients in inverse mm units
mulist = [0.002 0.0096 0.0120];
for ii=1:length(mulist)
	x(x == ii) = mulist(ii);
end

x = ir_phantom_resize(x, nx, ny);


function read_zubal_attn_test
x = read_zubal_attn('nx', 120, 'ny', 90);
im(x, 'zubal phantom: transmission'), cbar
