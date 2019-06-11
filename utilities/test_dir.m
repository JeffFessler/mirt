 function tdir = test_dir(tdir)
%function tdir = test_dir(tdir)
%
% Return path to test directory.
% Users can/should customize this file.

persistent pdir
if ~isvar('pdir') || isempty(pdir)
	try
		pdir.name = getenv('USER');
	catch
		pdir.name = 'user';
	end
	pdir.name = ['test,'  pdir.name];
	pdir.top = '/tmp/';
	pdir.path = [pdir.top pdir.name '/']; % trailing '/' is essential!
end

if nargin
	pdir.path = tdir;
end

if ~exist(pdir.path, 'dir')
	printf('Directory "%s" not found, so it will be created', pdir.path)
%	printf('run the following command!:')
%	'!mkdir /tmp/test'
	try
		mkdir(pdir.top, pdir.name) % compat with matlab 6.x
	catch
		try
			mkdir(pdir.path) % matlab 7.x accepts full path
		catch
			error 'give up'
		end
	end
end

tdir = pdir.path;
