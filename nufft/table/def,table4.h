// def,table4.h
// declarations
// Extended 2013-08-02 to 4D by David Johnson, The Ohio State University Wexner Medical Center

#include <string.h> // for memset

#define mymod(k,K) ((k) - (K) * floor((k) / (double) (K)))
#define iround(x) floor(x + 0.5)

// interp4_table1_for.c

extern void interp4_table0_complex_per(
const double *r_ck, /* [K1,K2,K3,K4] in */
const double *i_ck,
const int K1,
const int K2,
const int K3,
const int K4,
const double *r_h1, /* [J1*L1+1,1] in */
const double *i_h1,
const double *r_h2, /* [J2*L2+1,1] in */
const double *i_h2,
const double *r_h3, /* [J3*L3+1,1] in */
const double *i_h3,
const double *r_h4, /* [J4*L4+1,1] in */
const double *i_h4,
const int J1,
const int J2,
const int J3,
const int J4,
const int L1,
const int L2,
const int L3,
const int L4,
const double *p_tm, /* [M,4] in */
const int M,
double *r_fm, /* [M,1] out */
double *i_fm);

typedef void interp4_table_real_per_t(
const double *r_ck, /* [K1,K2,K3,K4] in */
const double *i_ck,
const int K1,
const int K2,
const int K3,
const int K4,
const double *r_h1, /* [J1*L1+1,1] in */
const double *r_h2, /* [J2*L2+1,1] in */
const double *r_h3, /* [J3*L3+1,1] in */
const double *r_h4, /* [J4*L4+1,1] in */
#ifdef Provide_flip
const int flip1, /* sign flips every K? */
const int flip2,
const int flip3,
const int flip4,
#endif
const int J1,
const int J2,
const int J3,
const int J4,
const int L1,
const int L2,
const int L3,
const int L4,
const double *p_tm, /* [M,4] in */
const int M,
double *r_fm, /* [M,1] out */
double *i_fm);

extern interp4_table_real_per_t interp4_table0_real_per;
extern interp4_table_real_per_t interp4_table1_real_per;


// interp4_table1_adj.c

extern void interp4_table0_complex_per_adj(
double *r_ck,		/* [K1,K2,K3,K4] out */
double *i_ck,
const int K1,
const int K2,
const int K3,
const int K4,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *i_h1,
const double *r_h2,	/* [J2*L2+1,1] in */
const double *i_h2,
const double *r_h3,	/* [J3*L3+1,1] in */
const double *i_h3,
const double *r_h4,	/* [J4*L4+1,1] in */
const double *i_h4,
const int J1,
const int J2,
const int J3,
const int J4,
const int L1,
const int L2,
const int L3,
const int L4,
const double *p_tm,	/* [M,4] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm);

typedef void interp4_table_real_per_adj_t(
double *r_ck,		/* [K1,K2,K3,K4] out */
double *i_ck,
const int K1,
const int K2,
const int K3,
const int K4,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *r_h2,	/* [J2*L2+1,1] in */
const double *r_h3,	/* [J3*L3+1,1] in */
const double *r_h4,	/* [J4*L4+1,1] in */
#ifdef Provide_flip
const int flip1,	/* sign flips every K? */
const int flip2,
const int flip3,
const int flip4,
#endif
const int J1,
const int J2,
const int J3,
const int J4,
const int L1,
const int L2,
const int L3,
const int L4,
const double *p_tm,	/* [M,4] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm);

extern interp4_table_real_per_adj_t interp4_table0_real_per_adj;
extern interp4_table_real_per_adj_t interp4_table1_real_per_adj;
