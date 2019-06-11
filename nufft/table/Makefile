# Makefile
# Build the .mex* files needed for fast NUFFT interpolation.

# assumes that Mathworks scripts "mex" and "mexext" are in your shell path
# determine mex file suffix and location of "mex" compiling script
suf := $(shell mexext)
mw_bin := $(shell dirname `which mexext`)
mexcom := $(mw_bin)/mex

# where to install the compiled mex files:
dir = ../../mex/v7/

mexo = $(mexcom) -output
mexa = chmod 755

# architecture specific compile flags:
aflag_x86_64 = -fPIC
aflag_mi64 = -arch x86_64
aflag = $(aflag_$(Arch))

# this -D_GNU_SOURCE flag came from looking at "mex -v" output...
CFLAGS = -D_GNU_SOURCE \
	-std=c99 -Wall -Wshadow -Wpointer-arith -W -Werror $(aflag)
mex = $(mexcom) CFLAGS='$(CFLAGS)'

def = def,table.h

goal:	\
	$(dir)interp1_table_mex.$(suf) \
	$(dir)interp2_table_mex.$(suf) \
	$(dir)interp3_table_mex.$(suf) \
	$(dir)interp4_table_mex.$(suf) \
	$(dir)interp1_table_adj_mex.$(suf) \
	$(dir)interp2_table_adj_mex.$(suf) \
	$(dir)interp3_table_adj_mex.$(suf) \
	$(dir)interp4_table_adj_mex.$(suf)
	@echo working on $(suf)


i1c = interp1_table1_for.c interp1_table_mex.c
$(dir)interp1_table_mex.$(suf):	$(i1c) $(def) def,table1.h
	$(mexo) $@ $(i1c)
	$(mexa) $@

i2c = interp2_table1_for.c interp2_table_mex.c
$(dir)interp2_table_mex.$(suf):	$(i2c) $(def) def,table2.h
	$(mexo) $@ $(i2c)
	$(mexa) $@

i3c = interp3_table1_for.c interp3_table_mex.c
$(dir)interp3_table_mex.$(suf):	$(i3c) $(def) def,table3.h
	$(mexo) $@ $(i3c)
	$(mexa) $@

i4c = interp4_table1_for.c interp4_table_mex.c
$(dir)interp4_table_mex.$(suf):	$(i4c) $(def) def,table4.h
	$(mexo) $@ $(i4c)
	$(mexa) $@


a1c = interp1_table1_adj.c interp1_table_adj_mex.c
$(dir)interp1_table_adj_mex.$(suf):	$(a1c) $(def) def,table1.h
	$(mexo) $@ $(a1c)
	$(mexa) $@

a2c = interp2_table1_adj.c interp2_table_adj_mex.c
$(dir)interp2_table_adj_mex.$(suf):	$(a2c) $(def) def,table2.h
	$(mexo) $@ $(a2c)
	$(mexa) $@

a3c = interp3_table1_adj.c interp3_table_adj_mex.c
$(dir)interp3_table_adj_mex.$(suf):	$(a3c) $(def) def,table3.h
	$(mexo) $@ $(a3c)
	$(mexa) $@

a4c = interp4_table1_adj.c interp4_table_adj_mex.c
$(dir)interp4_table_adj_mex.$(suf):	$(a4c) $(def) def,table4.h
	$(mexo) $@ $(a4c)
	$(mexa) $@


# how to compile for freemat import?
todo:
	cc -dynamic -fPIC -c interp1_table1.c
	libtool -o interp1_table1.a interp1_table1.o
	#libtool -dynamic -noall_load -o interp1_table1.so interp1_table1.o

reset:
	/bin/rm -i $(dir)interp*table*mex.$(suf)
