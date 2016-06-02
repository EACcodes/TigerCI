#ifndef LAPACK_DEFS_H
#define	LAPACK_DEFS_H


// We use a typedef to handle different possible integer sizes in our LAPACK
// library. Note that in principle, we need 64bit integers our largest calculations
#ifdef TIGER_LAPACK_INT4
typedef int32_t int_lapack;
#else
typedef long long int_lapack;
#endif
extern "C" void dpstf2_(char *uplo, int_lapack *n, double* a, int_lapack *lda,
        int_lapack* piv, int_lapack *rank, double *tol,
        double *work, int_lapack *info);

extern "C" void dgemv_(char *trans, int_lapack *m, int_lapack *n, double* alpha,
        double *A, int_lapack *lda, double* x, int_lapack* incx,
        double *beta, double *y, int_lapack *incy);

extern "C" void dswap_(int_lapack* n, double* x, int_lapack* incx, double * y, int_lapack* incy);


#endif	/* LAPACK_DEFS_H */

