# mri/Makefile

# compile MRI-related .mex* files
# assumes that Mathworks scripts "mex" and "mexext" are in your shell path

# determine mex file suffix and location of "mex" compiling script
suf := $(shell mexext)
mw_bin := $(shell dirname `which mexext`)
mexcom := $(mw_bin)/mex

dir = ../mex/v7

goal1 = $(dir)/mri_exp_mult_mex.$(suf)
goal2 = $(dir)/exp_xform_mex.$(suf)

all:
	@echo "mw_bin = $(mw_bin)"
	@echo "mexcom = $(mexcom)"
	$(MAKE) $(goal1) # date=-DDate='\"date='`date +%Y-%m-%d,%T`'\"'
	$(MAKE) $(goal2) # date=-DDate='\"date='`date +%Y-%m-%d,%T`'\"'

aflag_x86_64 = -fPIC
aflag_mi64 = -arch x86_64
aflag = $(aflag_$(Arch))

CFLAGS = -D_GNU_SOURCE $(date) \
	-std=c99 -Wall -Wshadow -Wpointer-arith -W -Werror $(aflag)

com = $(mexcom) CFLAGS='$(CFLAGS) -fexceptions'

$(goal1):	mri_exp_mult_mex.c mri_exp_mult_template.c
	echo $(date)
	$(com) mri_exp_mult_mex.c
	mv mri_exp_mult_mex.$(suf) $@
	chmod 755 $@

$(goal2):	exp_xform_mex.c exp_xform_template.c
	$(com) exp_xform_mex.c
	mv exp_xform_mex.$(suf) $@
	chmod 755 $@

reset:
	@echo suf = $(suf)
	rm -i $(goal1) $(goal2)
