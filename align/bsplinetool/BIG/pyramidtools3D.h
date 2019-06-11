/*******************************************************************************
        Trivial extension of pyramidtools.h in spline pyramids package at
        http://bigwww.epfl.ch/algorithms.html to 3D reduce and expand

        Coded by Se Young Chun, May 25, 2007, The Univ. of Michigan

*******************************************************************************/
					
extern int 	Reduce_3D (	
			float *In, long NxIn, long NyIn, long NzIn,
			float *Out,
			double w[], long nw,
			short FlagCentered
		);
				
extern int 	Expand_3D (	
			float *In, long NxIn, long NyIn, long NzIn,
			float *Out,
			double w[], long nw,
			short FlagCentered
		);
