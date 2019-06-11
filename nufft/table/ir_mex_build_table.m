% ir_mex_build_table
% run matlab's "mex" command to "compile" the table interpolation code
% into mex files.
% only users on unsupported systems, e.g., PCs, will need to do this

dir_current = pwd;
dir_nufft = path_find_dir('nufft');
dir_table = [dir_nufft filesep 'table'];
cd(dir_table)

% mex compile command with optimization and c99 flag:
fun = @(f1, f2) mex('-O', 'CFLAGS="\$CFLAGS -std=c99"', f1, f2)

fun('interp1_table_adj_mex.c',	'interp1_table1_adj.c')
fun('interp1_table_mex.c',	'interp1_table1_for.c')
fun('interp2_table_adj_mex.c',	'interp2_table1_adj.c')
fun('interp2_table_mex.c',	'interp2_table1_for.c')
fun('interp3_table_adj_mex.c',	'interp3_table1_adj.c')
fun('interp3_table_mex.c',	'interp3_table1_for.c')

%{ old way
mex interp1_table_adj_mex.c	interp1_table1_adj.c
mex interp1_table_mex.c	interp1_table1_for.c
mex interp2_table_adj_mex.c	interp2_table1_adj.c
mex interp2_table_mex.c	interp2_table1_for.c
mex interp3_table_adj_mex.c	interp3_table1_adj.c
mex interp3_table_mex.c	interp3_table1_for.c
%}

cd(dir_current)
