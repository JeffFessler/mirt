# makefile for compiling bspline tools 

# assumes that Mathworks scripts "mex" and "mexext" are in your shell path
# determine mex file suffix and location of "mex" compiling script
suf := $(shell mexext)
mw_bin := $(shell dirname `which mexext`)
mexcom := $(mw_bin)/mex

# most likely you will need to edit this path
# MPATH = /usr/local/matlab/bin/
# MPATH = /Volumes/s2/fessler//l/matlab/bin/

# architecture specific compile flags:
aflag_x86_64 = -fPIC
aflag_mi64 = -arch x86_64
aflag = $(aflag_$(Arch))

CFLAGS = -O3 -std=c99 -Wall -Wshadow -Wpointer-arith -W -Werror $(aflag)

mxc = $(mexcom) CFLAGS='$(CFLAGS)' -output

#suf = $(shell $(MPATH)mexext)
dir = ../mex

all:	\
	$(dir)/BsplExpand.$(suf) \
	$(dir)/BsplReduce.$(suf) \
	$(dir)/BsplVal2CoMirr.$(suf) \
	$(dir)/BsplCo2ValMirr.$(suf) \
	$(dir)/BsplCo2GdXMirr.$(suf) \
	$(dir)/BsplCo2GdYMirr.$(suf) \
	$(dir)/BsplCo2GdZMirr.$(suf) \
	$(dir)/BsplCo2ValTranMirr.$(suf) \
	$(dir)/BsplCo2GdXTranMirr.$(suf) \
	$(dir)/BsplCo2GdYTranMirr.$(suf) \
	$(dir)/BsplCo2GdZTranMirr.$(suf) \
	$(dir)/BsplVal2CoZero.$(suf) \
	$(dir)/BsplCo2ValZero.$(suf) \
	$(dir)/BsplCo2GdXZero.$(suf) \
	$(dir)/BsplCo2GdYZero.$(suf) \
	$(dir)/BsplCo2GdZZero.$(suf) \
	$(dir)/BsplCo2ValTranZero.$(suf) \
	$(dir)/BsplCo2GdXTranZero.$(suf) \
	$(dir)/BsplCo2GdYTranZero.$(suf) \
	$(dir)/BsplCo2GdZTranZero.$(suf) \
	$(dir)/BsplCo2ValZeroFilt.$(suf) \
	$(dir)/BsplCo2ValTranZeroFilt.$(suf)

PYRAM = BIG/pyramids/messagedisplay.c BIG/pyramids/pyramidtools.c BIG/pyramids/pyramidfilters.c BIG/pyramidtools3D.c

$(dir)/BsplExpand.$(suf): BsplExpand.c $(PYRAM)
	$(mxc) $@ BsplExpand.c $(PYRAM)

$(dir)/BsplReduce.$(suf): BsplReduce.c $(PYRAM)
	$(mxc) $@ BsplReduce.c $(PYRAM)

INTERP = batchinterpol2D.c batchinterpol3D.c interpol2D.c interpol3D.c

$(dir)/BsplVal2CoMirr.$(suf): BsplVal2CoMirr.c $(INTERP)
	$(mxc) $@ BsplVal2CoMirr.c $(INTERP)

$(dir)/BsplCo2ValMirr.$(suf): BsplCo2ValMirr.c $(INTERP)
	$(mxc) $@ BsplCo2ValMirr.c $(INTERP)

$(dir)/BsplCo2GdXMirr.$(suf): BsplCo2GdXMirr.c $(INTERP)
	$(mxc) $@ BsplCo2GdXMirr.c $(INTERP)

$(dir)/BsplCo2GdYMirr.$(suf): BsplCo2GdYMirr.c $(INTERP)
	$(mxc) $@ BsplCo2GdYMirr.c $(INTERP)

$(dir)/BsplCo2GdZMirr.$(suf): BsplCo2GdZMirr.c $(INTERP)
	$(mxc) $@ BsplCo2GdZMirr.c $(INTERP)

TINTERP =  batchinterpol2Dtran.c batchinterpol3Dtran.c interpol2Dtran.c interpol3Dtran.c

$(dir)/BsplCo2ValTranMirr.$(suf): BsplCo2ValTranMirr.c $(TINTERP)
	$(mxc) $@ BsplCo2ValTranMirr.c $(TINTERP)

$(dir)/BsplCo2GdXTranMirr.$(suf): BsplCo2GdXTranMirr.c $(TINTERP)
	$(mxc) $@ BsplCo2GdXTranMirr.c $(TINTERP)

$(dir)/BsplCo2GdYTranMirr.$(suf): BsplCo2GdYTranMirr.c $(TINTERP)
	$(mxc) $@ BsplCo2GdYTranMirr.c $(TINTERP)

$(dir)/BsplCo2GdZTranMirr.$(suf): BsplCo2GdZTranMirr.c $(TINTERP)
	$(mxc) $@ BsplCo2GdZTranMirr.c $(TINTERP)


$(dir)/BsplVal2CoZero.$(suf): BsplVal2CoZero.c
	$(mxc) $@ BsplVal2CoZero.c

$(dir)/BsplCo2ValZero.$(suf): BsplCo2ValZero.c $(INTERP)
	$(mxc) $@ BsplCo2ValZero.c $(INTERP)

$(dir)/BsplCo2GdXZero.$(suf): BsplCo2GdXZero.c $(INTERP)
	$(mxc) $@ BsplCo2GdXZero.c $(INTERP)

$(dir)/BsplCo2GdYZero.$(suf): BsplCo2GdXZero.c $(INTERP)
	$(mxc) $@ BsplCo2GdYZero.c $(INTERP)

$(dir)/BsplCo2GdZZero.$(suf): BsplCo2GdZZero.c $(INTERP)
	$(mxc) $@ BsplCo2GdZZero.c $(INTERP)

$(dir)/BsplCo2ValTranZero.$(suf): BsplCo2ValTranZero.c $(TINTERP)
	$(mxc) $@ BsplCo2ValTranZero.c $(TINTERP)

$(dir)/BsplCo2GdXTranZero.$(suf): BsplCo2GdXTranZero.c $(TINTERP)
	$(mxc) $@ BsplCo2GdXTranZero.c $(TINTERP)

$(dir)/BsplCo2GdYTranZero.$(suf): BsplCo2GdYTranZero.c $(TINTERP)
	$(mxc) $@ BsplCo2GdYTranZero.c $(TINTERP)

$(dir)/BsplCo2GdZTranZero.$(suf): BsplCo2GdZTranZero.c $(TINTERP)
	$(mxc) $@ BsplCo2GdZTranZero.c $(TINTERP)

FINTERP = BIG/BsplnTrfKer.c BIG/reg/convolve.c BIG/reg/getPut.c BIG/phil.c

$(dir)/BsplCo2ValZeroFilt.$(suf): BsplCo2ValZeroFilt.c filtinterpol2D.c filtinterpol3D.c $(FINTERP)
	$(mxc) $@ BsplCo2ValZeroFilt.c filtinterpol2D.c filtinterpol3D.c $(FINTERP)

$(dir)/BsplCo2ValTranZeroFilt.$(suf): BsplCo2ValTranZeroFilt.c filtinterpol2Dtran.c filtinterpol3Dtran.c $(FINTERP)
	$(mxc) $@ BsplCo2ValTranZeroFilt.c filtinterpol2Dtran.c filtinterpol3Dtran.c $(FINTERP)

reset:
	rm -i ../mex/*.$(suf)
