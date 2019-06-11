# Makefile
# one can reduce time by compiling the initialization routine

suf_i686 = mexglx
suf_ppc = mexppc
suf_ = suf_i686; # default is linux
suf = $(suf_($Arch))

goal:	\
	interp1_table_mex.$(suf) \
	interp2_table_mex.$(suf) \
	interp3_table_mex.$(suf) \
	interp1_table_adj_mex.$(suf) \
	interp2_table_adj_mex.$(suf) \
	interp3_table_adj_mex.$(suf)

CFLAGS = -ansi -Wall -Wshadow -Wpointer-arith -W -Werror
mex = mex CFLAGS='$(CFLAGS)'

interp1_table_mex.$(suf):	interp1_table_mex.c
	$(mex) $?

interp2_table_mex.$(suf):	interp2_table_mex.c
	$(mex) $?

interp3_table_mex.$(suf):	interp3_table_mex.c
	$(mex) $?

interp1_table_adj_mex.$(suf):	interp1_table_adj_mex.c
	$(mex) $?

interp2_table_adj_mex.$(suf):	interp2_table_adj_mex.c
	$(mex) $?

interp3_table_adj_mex.$(suf):	interp3_table_adj_mex.c
	$(mex) $?

mat = /users/apps/matlab5.3/toolbox/matlab
mfile = \
	$(mat)/datatypes/ndgrid.m \
	$(mat)/elfun/mod.m \
	$(mat)/ops/kron.m

goal = nufft2_init.mexaxp 
main = nufft2_init.m

$(goal):	$(main) $(mfile)
		mcc -x -h $(main) $(mfile)

reset:
	rm -i $(goal)


#plist = nufft2.m nufft2_init.m nufft_diric.m nufft_R.m
plist = *.m
en7=print,code
comm=$(en7) `linesplit $? | sort`; lpq; touch $@
.print:	$(plist)
	$(comm)
