#ifndef MATRIX_IDZ_H
#define MATRIX_IDZ_H

#include <Rinternals.h>

void iswap2(size_t,      int *, ptrdiff_t,      int *, ptrdiff_t);
void dswap2(size_t,   double *, ptrdiff_t,   double *, ptrdiff_t);
void zswap2(size_t, Rcomplex *, ptrdiff_t, Rcomplex *, ptrdiff_t);

void iswap1(size_t,      int *, ptrdiff_t, ptrdiff_t,      int *, ptrdiff_t, ptrdiff_t);
void dswap1(size_t,   double *, ptrdiff_t, ptrdiff_t,   double *, ptrdiff_t, ptrdiff_t);
void zswap1(size_t, Rcomplex *, ptrdiff_t, ptrdiff_t, Rcomplex *, ptrdiff_t, ptrdiff_t);

void icopy2(size_t,      int *, ptrdiff_t, const      int *, ptrdiff_t);
void dcopy2(size_t,   double *, ptrdiff_t, const   double *, ptrdiff_t);
void zcopy2(size_t, Rcomplex *, ptrdiff_t, const Rcomplex *, ptrdiff_t);

void icopy1(size_t,      int *, ptrdiff_t, ptrdiff_t, const      int *, ptrdiff_t, ptrdiff_t);
void dcopy1(size_t,   double *, ptrdiff_t, ptrdiff_t, const   double *, ptrdiff_t, ptrdiff_t);
void zcopy1(size_t, Rcomplex *, ptrdiff_t, ptrdiff_t, const Rcomplex *, ptrdiff_t, ptrdiff_t);

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

void iband2(     int *, const      int *, size_t, size_t, ptrdiff_t, ptrdiff_t);
void dband2(  double *, const   double *, size_t, size_t, ptrdiff_t, ptrdiff_t);
void zband2(Rcomplex *, const Rcomplex *, size_t, size_t, ptrdiff_t, ptrdiff_t);

void iband1(     int *, const      int *, size_t, char, ptrdiff_t, ptrdiff_t);
void dband1(  double *, const   double *, size_t, char, ptrdiff_t, ptrdiff_t);
void zband1(Rcomplex *, const Rcomplex *, size_t, char, ptrdiff_t, ptrdiff_t);

void zvreal(Rcomplex *, size_t);
void zvimag(Rcomplex *, size_t);
void zvconj(Rcomplex *, size_t);

#endif /* MATRIX_IDZ_H */
