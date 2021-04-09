
#ifndef MWSIZE_MAX
    #define  mwIndex        int
    #define  mwSignedIndex  int
    #define  mwSize         int
#endif

extern int dgemv_(char *trans, int *m, int *n, double *
	alpha, double *a, int *lda, double *x, int *incx, 
	double *beta, double *y, int *incy);

#if !defined(_WIN32)
#define dgemv dgemv_
#endif

#define mwIndex ptrdiff_t
#define mwSize  ptrdiff_t