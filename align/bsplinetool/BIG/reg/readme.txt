This distribution is dated April 8, 2002.

List of files:
-  main.c <-- descriptive only, should _not_ be used as is
-  readme.txt <-- this file
- BsplnTrf.c
- BsplnTrf.h
- BsplnWgt.c
- BsplnWgt.h
- convolve.c
- convolve.h
- getPut.c
- getPut.h
- phil.h
- pyrFilt.c
- pyrGetSz.c
- quant.c
- quant.h
- regFlt3d.c <-- procedure to call
- regFlt3d.h
- register.h
- reg0.c
- reg0.h
- reg1.c
- reg1.h
- reg2.c
- reg2.h
- reg3.c
- reg3.h
- svdcmp.c
- svdcmp.h

Important: the file 'main.c' is supposed to get you started in the
process of integrating this set of registration routines into your own
private system. It gives a description of the variables and of their
use. What it _doesn't_ do is to provide you with a working code for
calling the registration procedure. It gives just the outline; you'll
have to adapt the values to your own case. Link error: the procedure
'message()' is left voluntarily unimplemented, since IO's are strongly
system-dependent.

You'll be free to use this software for research purposes, but you
should not redistribute it without our consent. In addition, we expect
you to include a citation or acknowlegment whenever you present or
publish results that are based on it.

---
P. Thevenaz and M. Unser, "A Pyramid Approach to Subpixel Registration
Based on Intensity," IEEE Transactions on Image Processing, vol. 7, no.
1, pp. 27-41, January 1998.

P. Thevenaz, U.E. Ruttimann and M. Unser, "Iterative Multi-Scale
Registration without Landmarks," Proc. IEEE International Conference on
Image Processing, Washington, DC, USA, October 23-26, 1995, vol. III,
pp. 228-231.

M. Unser, A. Aldroubi and C.R. Gerfen, "A Multiresolution Image
Registration Procedure Using Spline Pyramids," in Proc. SPIE, San Diego,
CA, USA, July 15-16, 1993, vol. 2034, Mathematical Imaging: Wavelet
Applications in Signal and Image Processing, pp. 160-170.
---

Philippe Thevenaz
Swiss Federal Institute of Technology, Lausanne
EPFL/STI/IOA/LIB
Bldg. BM-Ecublens 4.137
Station 17
CH-1015 Lausanne VD
Switzerland

+41(21)693.51.61 (phone)
+41(21)693.37.01 (fax)
philippe.thevenaz@epfl.ch
http://bigwww.epfl.ch/
