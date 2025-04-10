#ifndef MATRIX_IDZ_H
#define MATRIX_IDZ_H

#include <stddef.h> /* size_t */
#include <R_ext/Complex.h> /* Rcomplex */

void iswap2(size_t,      int *, size_t,      int *, size_t);
void dswap2(size_t,   double *, size_t,   double *, size_t);
void zswap2(size_t, Rcomplex *, size_t, Rcomplex *, size_t);

void iswap1(size_t,      int *, size_t, size_t, int,      int *, size_t, size_t, int);
void dswap1(size_t,   double *, size_t, size_t, int,   double *, size_t, size_t, int);
void zswap1(size_t, Rcomplex *, size_t, size_t, int, Rcomplex *, size_t, size_t, int);

void icopy2(size_t,      int *, size_t, const      int *, size_t);
void dcopy2(size_t,   double *, size_t, const   double *, size_t);
void zcopy2(size_t, Rcomplex *, size_t, const Rcomplex *, size_t);

void icopy1(size_t,      int *, size_t, size_t, int, const      int *, size_t, size_t, int);
void dcopy1(size_t,   double *, size_t, size_t, int, const   double *, size_t, size_t, int);
void zcopy1(size_t, Rcomplex *, size_t, size_t, int, const Rcomplex *, size_t, size_t, int);

void irowperm2(     int *, const      int *, int, int, const int *, int, int);
void drowperm2(  double *, const   double *, int, int, const int *, int, int);
void zrowperm2(Rcomplex *, const Rcomplex *, int, int, const int *, int, int);

void isymperm2(     int *, const      int *, int, char, const int *, int, int);
void dsymperm2(  double *, const   double *, int, char, const int *, int, int);
void zsymperm2(Rcomplex *, const Rcomplex *, int, char, const int *, int, int);

void isymperm1(     int *, const      int *, int, char, const int *, int, int);
void dsymperm1(  double *, const   double *, int, char, const int *, int, int);
void zsymperm1(Rcomplex *, const Rcomplex *, int, char, const int *, int, int);

void ipack2(     int *, const      int *, size_t, char, char, char);
void dpack2(  double *, const   double *, size_t, char, char, char);
void zpack2(Rcomplex *, const Rcomplex *, size_t, char, char, char);

void ipack1(     int *, const      int *, size_t, char, char, char);
void dpack1(  double *, const   double *, size_t, char, char, char);
void zpack1(Rcomplex *, const Rcomplex *, size_t, char, char, char);

void iforce2(     int *, const      int *, size_t, char, char, char);
void dforce2(  double *, const   double *, size_t, char, char, char);
void zforce2(Rcomplex *, const Rcomplex *, size_t, char, char, char);

void iforce1(     int *, const      int *, size_t, char, char, char);
void dforce1(  double *, const   double *, size_t, char, char, char);
void zforce1(Rcomplex *, const Rcomplex *, size_t, char, char, char);

void itrans2(     int *, const      int *, size_t, size_t, char);
void dtrans2(  double *, const   double *, size_t, size_t, char);
void ztrans2(Rcomplex *, const Rcomplex *, size_t, size_t, char);

void itrans1(     int *, const      int *, size_t, char, char);
void dtrans1(  double *, const   double *, size_t, char, char);
void ztrans1(Rcomplex *, const Rcomplex *, size_t, char, char);

void iband2(     int *, const      int *, size_t, size_t, size_t, size_t);
void dband2(  double *, const   double *, size_t, size_t, size_t, size_t);
void zband2(Rcomplex *, const Rcomplex *, size_t, size_t, size_t, size_t);

void iband1(     int *, const      int *, size_t, char, size_t, size_t);
void dband1(  double *, const   double *, size_t, char, size_t, size_t);
void zband1(Rcomplex *, const Rcomplex *, size_t, char, size_t, size_t);

int ntest2(const      int *, size_t, char, char, char);
int ltest2(const      int *, size_t, char, char, char);
int itest2(const      int *, size_t, char, char, char);
int dtest2(const   double *, size_t, char, char, char);
int ztest2(const Rcomplex *, size_t, char, char, char);

int ntest1(const      int *, size_t, char, char, char);
int ltest1(const      int *, size_t, char, char, char);
int itest1(const      int *, size_t, char, char, char);
int dtest1(const   double *, size_t, char, char, char);
int ztest1(const Rcomplex *, size_t, char, char, char);

void ntspaggr(int *, int *,      int *, const int *, const int *, const      int *, int, int, int *, int *,      int *);
void ltspaggr(int *, int *,      int *, const int *, const int *, const      int *, int, int, int *, int *,      int *);
void itspaggr(int *, int *,      int *, const int *, const int *, const      int *, int, int, int *, int *,      int *);
void dtspaggr(int *, int *,   double *, const int *, const int *, const   double *, int, int, int *, int *,   double *);
void ztspaggr(int *, int *, Rcomplex *, const int *, const int *, const Rcomplex *, int, int, int *, int *, Rcomplex *);

void ntspsort(int *, int *,      int *, const int *, const int *, const      int *, int, int, int *, int *,      int *);
void ltspsort(int *, int *,      int *, const int *, const int *, const      int *, int, int, int *, int *,      int *);
void itspsort(int *, int *,      int *, const int *, const int *, const      int *, int, int, int *, int *,      int *);
void dtspsort(int *, int *,   double *, const int *, const int *, const   double *, int, int, int *, int *,   double *);
void ztspsort(int *, int *, Rcomplex *, const int *, const int *, const Rcomplex *, int, int, int *, int *, Rcomplex *);

void ncspsort(int *, int *,      int *, int, int, int *,      int *);
void lcspsort(int *, int *,      int *, int, int, int *,      int *);
void icspsort(int *, int *,      int *, int, int, int *,      int *);
void dcspsort(int *, int *,   double *, int, int, int *,   double *);
void zcspsort(int *, int *, Rcomplex *, int, int, int *, Rcomplex *);

void ncsptrans(int *, int *,      int *, const int *, const int *, const      int *, int, int, char, int *);
void lcsptrans(int *, int *,      int *, const int *, const int *, const      int *, int, int, char, int *);
void icsptrans(int *, int *,      int *, const int *, const int *, const      int *, int, int, char, int *);
void dcsptrans(int *, int *,   double *, const int *, const int *, const   double *, int, int, char, int *);
void zcsptrans(int *, int *, Rcomplex *, const int *, const int *, const Rcomplex *, int, int, char, int *);

void zvreal(Rcomplex *, const Rcomplex *, size_t);
void zvimag(Rcomplex *, const Rcomplex *, size_t);
void zvconj(Rcomplex *, const Rcomplex *, size_t);

#endif /* MATRIX_IDZ_H */
