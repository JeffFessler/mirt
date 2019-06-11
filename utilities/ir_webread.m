 function st = ir_webread(urlsuff, varargin)
%function st = ir_webread(urlsuff, varargin)
%|
%| read a .mat file from url = [arg.urlroot urlsuff] and return as struct
%|
%| 2015-08-27, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if streq(urlsuff, 'test'), ir_webread_test, return, end

arg.urlroot = 'http://web.eecs.umich.edu/~fessler/irt/data/';
arg.url = '';
arg.tmpfile = 'ir_webread_tmp1.mat';
arg = vararg_pair(arg, varargin);

if isempty(arg.url) % default is to form url from root and suffix
	arg.url = [arg.urlroot urlsuff];
end

if exist(arg.tmpfile, 'file')
	fail('%s exists', arg.tmpfile)
end

if ir_is_octave
	st = ir_webread_oct(arg.url, arg.tmpfile);
else
	st = ir_webread_mat(arg.url, arg.tmpfile);
end


% ir_webread_mat()
function st = ir_webread_mat(url, tmpfile)

if streq(url(end-3:end), '.mat')
	options = weboptions('ContentType', 'binary');
	tmp = webread(url, options); % reads as a long uint8 array!
	fid = fopen(tmpfile, 'w');
	if fid == -1, fail('open %s', tmpfile), end
	count = fwrite(fid, tmp, 'char*1');
	if count ~= numel(tmp), fail('count'), end
	if fclose(fid) == -1, fail('fclose'), end
	st = load(tmpfile);
	delete(tmpfile)
else
	fail('not done for "%s"', url)
end


% ir_webread_oct()
function st = ir_webread_oct(url, tmpfile)

if streq(url(end-3:end), '.mat')
	tmp = urlread(url); % reads it as a long char string!
	fid = fopen(tmpfile, 'w');
	if fid == -1, fail('open %s', tmpfile), end
	count = fwrite(fid, tmp, 'char*1');
	if count ~= numel(tmp), fail('count'), end
	if fclose(fid) == -1, fail('fclose'), end
	st = load(tmpfile);
	delete(tmpfile)
else
	fail('not done for "%s"', url)
end


% ir_webread_test
function ir_webread_test

st = ir_webread('test/test15.mat');
x = st.x;
jf_equal(x, uint32(1:5))
