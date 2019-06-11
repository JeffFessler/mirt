# bloch/Makefile
# compile Bloch simulator

# We need a unique name for the binary for each OS type.
# We base it on matlab's "mexext" command.

ending = $(shell mexext)

goal = blochmulti_spsp-$(ending)

all:	$(goal)
	@echo "done $(goal) for ending=$(ending)"

$(goal):	blochmulti_spsp_2008_12_03.c
		gcc -o $@ $? -lm
