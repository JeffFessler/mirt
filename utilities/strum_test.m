  function strum_test
%|function strum_test
%| test strum object
%| Copyright 2006-1-19, Jeff Fessler, University of Michigan

dat.a = 'x';
dat.b = 'y';

imp10 = @(st) st.b;
imp11 = @(st, arg) [st.a arg];
imp12 = @(st, arg1, arg2) sum(arg1) + sum(arg2);

list = { ...
	'fun10', @strum_test_fun10, '()'; ...
	'fun11', @strum_test_fun11, '(arg)'; ...
	'fun12', @strum_test_fun12, '(arg1, arg2)'; ...
	'imp11', imp11, '(arg)'; ...
	'imp12', imp12, '(arg1, arg2)'; ...
};

% no longer attempting anything but exactly 1 output so comment out these:
%	'fun00', @strum_test_fun00, '()'; ...
%	'fun22', @strum_test_fun22, '(arg1, arg2)'; ...

for ii=1:2
	if ii==1
		st = strum(dat, list); % with comment
	else
		tmp = col(list(:,1:2)');
		st = strum(dat, tmp); % without comment
	end
	st.a;
	st.b;
	jf_equal(st.fun10(), dat.a)
	jf_equal(st.fun11('d'), [dat.a 'd'])
	jf_equal(st.fun12('d', 'e'), [dat.a 'd' 'e'])
	jf_equal(st.imp11('d'), imp11(dat, 'd'))
	jf_equal(st.imp12(1:2, 3:4), imp12(dat, 1:2, 3:4))

	% augment base strum with more data and methods
	new.c = {10, 20};
	new.d = 'z';
	if ii==1
		s2 = strum(new, {'imp10', imp10, 'imp10()'}, 'base', st);
	else
		s2 = strum(new, {'imp10', imp10}, 'base', st);
	end
	jf_equal(s2.d, 'z')
	jf_equal(s2.imp10(), st.b)

	% cell arguments
	new.c{1};
	s2.c{1};
	jf_equal(s2.c{1}, new.c{1})
	jf_equal(s2.c{2}, new.c{2})

%	s2.c{1:2} % fails due to matlab limitation
%	s2.c{:} % ""
	tmp = s2.c; jf_equal(new.c, tmp) % this is the workaround
end


if 0 % old attempt to see if two outputs would work (it does not)
	st = strum(dat, {'fun22', @strum_test_fun22});
	jf_equal(st.fun22('c', 'd'), 'xc')
%	[a b] = st.fun22('c', 'd'); % "Too many output arguments" says matlab

	try
		[a b] = st.fun22(3, 4)
		printm 'multiple outputs worked!'
	catch
		warn 'darn, matlab cannot handle multiple outputs, as matt said'
	end
end


function strum_test_fun00(ob) % 0 out, 0 in 
printm 'ok'

function a = strum_test_fun10(ob) % 1 out, 0 in
a = ob.a;

function out = strum_test_fun11(ob, arg) % 1 out, 1 in
out = [ob.a arg];

function out = strum_test_fun12(ob, arg1, arg2) % 1 out, 2 in
out = [ob.a arg1 arg2];

function [a, b] = strum_test_fun22(ob, arg1, arg2) % 2 out, 2 in
a = [ob.a arg1];
b = [ob.b arg2];
