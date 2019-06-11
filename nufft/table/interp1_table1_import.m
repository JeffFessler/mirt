function interp1_table1_import

libname = 'interp1_table1.a';

type11r = ['double[K1] r_ck, double[K1] i_ck, int32 K1, ' ...
	'double[J1*L1+1] r_h1, int32 J1, int32 L1, ' ...
	'double[M] p_tm, int32 M, double[M] &r_fm, double[M] &i_fm'];

myimport(libname, 'interp1_table1_real_per', 'void', type11r);

function myimport(libname, name, type_return, type_call)
import(libname, name, name, type_return, type_call);
