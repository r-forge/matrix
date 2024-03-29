#ifndef MATRIX_IDZ_H
#define MATRIX_IDZ_H

#include <Rinternals.h>

#define IDZ \
TEMPLATE(i,      int); \
TEMPLATE(d,   double); \
TEMPLATE(z, Rcomplex);

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
rowperm2(_CTYPE_ *, int, int, int *, int, int)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
symperm2(_CTYPE_ *, int, char, int *, int, int)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
pack2(_CTYPE_ *, const _CTYPE_ *, int, char, char)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
unpack1(_CTYPE_ *, const _CTYPE_ *, int, char, char)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
trans2(_CTYPE_ *, const _CTYPE_ *, int, int, char)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
trans1(_CTYPE_ *, const _CTYPE_ *, int, char, char)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
syforce2(_CTYPE_ *, int, char, char)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
trforce2(_CTYPE_ *, int, int, char, char)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
band2(_CTYPE_ *, int, int, int, int)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
band1(_CTYPE_ *, int, int, int, char)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
dcopy2(_CTYPE_ *, const _CTYPE_ *, int, R_xlen_t, char, char)
IDZ
#undef TEMPLATE

#define TEMPLATE(_PREFIX_, _CTYPE_) \
void _PREFIX_ ## \
dcopy1(_CTYPE_ *, const _CTYPE_ *, int, R_xlen_t, char, char, char)
IDZ
#undef TEMPLATE

void zdreal2(Rcomplex *, int);
void zdreal1(Rcomplex *, int, char);

void zvreal(Rcomplex *, R_xlen_t);
void zvimag(Rcomplex *, R_xlen_t);
void zvconj(Rcomplex *, R_xlen_t);

#undef IDZ

#endif /* MATRIX_IDZ_H */
