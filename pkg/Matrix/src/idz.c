#include "Mdefines.h"
#include "M5.h"
#include "idz.h"

#define TEMPLATE(c) \
void \
c##swap2(size_t n, \
         c##TYPE *x, size_t dx, \
         c##TYPE *y, size_t dy) \
{ \
	c##TYPE tmp; \
	while (n--) { \
		tmp = *x; \
		*x = *y; \
		*y = tmp; \
		x += dx; \
		y += dy; \
	} \
	return; \
} \
 \
void \
c##swap1(size_t n, \
         c##TYPE *x, size_t dx, size_t addx, int nddx, \
         c##TYPE *y, size_t dy, size_t addy, int nddy) \
{ \
	c##TYPE tmp; \
	while (n--) { \
		tmp = *x; \
		*x = *y; \
		*y = tmp; \
		x += dx; \
		y += dy; \
		if (nddx) dx -= addx; else dx += addx; \
		if (nddy) dy -= addy; else dy += addy; \
	} \
	return; \
} \
 \
void \
c##copy2(size_t n, \
               c##TYPE *x, size_t dx, \
         const c##TYPE *y, size_t dy) \
{ \
	while (n--) { \
		*x = *y; \
		x += dx; \
		y += dy; \
	} \
	return; \
} \
 \
void \
c##copy1(size_t n, \
               c##TYPE *x, size_t dx, size_t addx, int nddx, \
         const c##TYPE *y, size_t dy, size_t addy, int nddy) \
{ \
	while (n--) { \
		*x = *y; \
		x += dx; \
		y += dy; \
		if (nddx) dx -= addx; else dx += addx; \
		if (nddy) dy -= addy; else dy += addy; \
	} \
	return; \
} \
 \
static void \
c##symswapr2(c##TYPE *x, \
             size_t n, char uplo, size_t i, size_t j) \
{ \
	if (n == 0 || i == j) \
		return; \
	if (i > j) { \
		size_t k = i; \
		i = j; \
		j = k; \
	} \
	c##TYPE *xi = x + i * n, *xj = x + j * n, tmp; \
	tmp = xi[i]; \
	xi[i] = xj[j]; \
	xj[j] = tmp; \
	if (uplo == 'U') { \
		c##swap2(i, xi, 1, xj, 1); \
		c##swap2(j - i - 1, xi + i + n, n, xj + i + 1, 1); \
		c##swap2(n - j - 1, xj + i + n, n, xj + j + n, n); \
	} else { \
		c##swap2(i, x + i, n, x + j, n); \
		c##swap2(j - i - 1, xi + i + 1, 1, xi + j + n, n); \
		c##swap2(n - j - 1, xi + j + 1, 1, xj + j + 1, 1); \
	} \
	return; \
} \
 \
static void \
c##symswapr1(c##TYPE *x, \
             size_t n, char uplo, size_t i, size_t j) \
{ \
	if (n == 0 || i == j) \
		return; \
	if (i >= j) { \
		size_t k = i; \
		i = j; \
		j = k; \
	} \
	c##TYPE *xi, *xj, tmp; \
	if (uplo == 'U') { \
		xi = x + PACKED_AR21_UP(0, i); \
		xj = x + PACKED_AR21_UP(0, j); \
		tmp = xi[i]; \
		xi[i] = xj[j]; \
		xj[j] = tmp; \
		c##swap1(i, xi, 1, 0, 0, xj, 1, 0, 0); \
		c##swap1(j - i - 1, xi + i + (i + 1), i + 2, 1, 0, xj + i +       1,     1, 0, 0); \
		c##swap1(n - j - 1, xj + i + (j + 1), j + 2, 1, 0, xj + j + (j + 1), j + 2, 1, 0); \
	} else { \
		xi = x + PACKED_AR21_LO(i, i, n); \
		xj = x + PACKED_AR21_LO(j, j, n); \
		tmp = xi[0]; \
		xi[0] = xj[0]; \
		xj[0] = tmp; \
		c##swap1(i, x + i, n - 1, 1, 1, x + j, n - 1, 1, 1); \
		c##swap1(j - i - 1, xi + (i - i) + 1, 1, 0, 0, xi + (j - i) + (n - i - 1), n - i - 2, 1, 1); \
		c##swap1(n - j - 1, xi + (j - i) + 1, 1, 0, 0, xj + (j - j) +           1,         1, 0, 0); \
	} \
	return; \
} \
 \
static void \
c##symcopyr2(c##TYPE *x, const c##TYPE *y, \
             size_t n, char uplo, size_t i, size_t j) \
{ \
	if (n == 0) \
		return; \
	      c##TYPE *xi = x + i * n, *xj = x + j * n; \
	const c##TYPE *yi = y + i * n, *yj = y + j * n; \
	xi[i] = yj[j]; \
	if (uplo == 'U') { \
		if (i <= j) { \
			xj[i] = yj[i]; \
			c##copy2(i, xi, 1, yj, 1); \
			if (i < j) \
			c##copy2(j - i - 1, xi + i + n, n, yj + i + 1, 1); \
			c##copy2(n - j - 1, xj + i + n, n, yj + j + n, n); \
		} else { \
			xi[j] = yi[j]; \
			c##copy2(j, xi, 1, yj, 1); \
			c##copy2(i - j - 1, xi + j + 1, 1, yj + j + n, n); \
			c##copy2(n - i - 1, xi + i + n, n, yi + j + n, n); \
		} \
	} else { \
		if (i <= j) { \
			xi[j] = yi[j]; \
			c##copy2(i, x + i, n, y + j, n); \
			if (i < j) \
			c##copy2(j - i - 1, xi + i + 1, 1, yi + j + n, n); \
			c##copy2(n - j - 1, xi + j + 1, 1, yj + j + 1, 1); \
		} else { \
			xj[i] = yj[i]; \
			c##copy2(j, x + i, n, y + j, n); \
			c##copy2(i - j - 1, xj + i + n, n, yj + j + 1, 1); \
			c##copy2(n - i - 1, xi + i + 1, 1, yj + i + 1, 1); \
		} \
	} \
	return; \
} \
 \
static void \
c##symcopyr1(c##TYPE *x, const c##TYPE *y, \
             size_t n, char uplo, size_t i, size_t j) \
{ \
	if (n == 0) \
		return; \
	      c##TYPE *xi, *xj; \
	const c##TYPE *yi, *yj; \
	if (uplo == 'U') { \
		xi = x + PACKED_AR21_UP(0, i); \
		xj = x + PACKED_AR21_UP(0, j); \
		yi = y + PACKED_AR21_UP(0, i); \
		yj = y + PACKED_AR21_UP(0, j); \
		xi[i] = yj[j]; \
		if (i <= j) { \
			xj[i] = yj[i]; \
			c##copy1(i, xi, 1, 0, 0, yj, 1, 0, 0); \
			if (i < j) \
			c##copy1(j - i - 1, xi + i + (i + 1), i + 2, 1, 0, yj + i +       1,     1, 0, 0); \
			c##copy1(n - j - 1, xj + i + (j + 1), j + 2, 1, 0, yj + j + (j + 1), j + 2, 1, 0); \
		} else { \
			xi[j] = yi[j]; \
			c##copy1(j, xi, 1, 0, 0, yj, 1, 0, 0); \
			c##copy1(i - j - 1, xi + j + 1, 1, 0, 0,           yj + j + (j + 1), j + 2, 1, 0); \
			c##copy1(n - i - 1, xi + i + (i + 1), 0, i + 2, 1, yi + j + (i + 1), i + 2, 1, 0); \
		} \
	} else { \
		xi = x + PACKED_AR21_LO(i, i, n); \
		xj = x + PACKED_AR21_LO(j, j, n); \
		yi = y + PACKED_AR21_LO(i, i, n); \
		yj = y + PACKED_AR21_LO(j, j, n); \
		xi[0] = yj[0]; \
		if (i <= j) { \
			xi[j] = yi[j]; \
			c##copy1(i, x + i, n - 1, 1, 1, y + j, n - 1, 1, 1); \
			if (i < j) \
			c##copy1(j - i - 1, xi + (i - i) + 1, 1, 0, 0, yi + (j - i) + (n - i - 1), n - i - 2, 1, 1); \
			c##copy1(n - j - 1, xi + (j - i) + 1, 1, 0, 0, yj + (j - j) +           1,         1, 0, 0); \
		} else { \
			xj[i] = yj[i]; \
			c##copy1(j, x + i, n - 1, 1, 1, y + j, n - 1, 1, 1); \
			c##copy1(i - j - 1, xj + (i - j) + (n - j - 1), n - j - 2, 1, 1, yj + (j - j) + 1, 1, 0, 0); \
			c##copy1(n - i - 1, xi + (i - i) +           1,         1, 0, 0, yj + (i - j) + 1, 1, 0, 0); \
		} \
	} \
	return; \
} \
 \
void \
c##rowperm2(c##TYPE *x, const c##TYPE *y, \
            int m, int n, const int *p, int off, int invert) \
{ \
	if (m <= 0 || n <= 0) \
		return; \
	size_t m_ = (size_t) m, n_ = (size_t) n; \
	if (!p) { \
		if (y) \
			memcpy(x, y, sizeof(c##TYPE) * m_ * n_); \
	} else { \
		if (y) { \
			int i, j; \
			if (!invert) \
			for (j = 0; j < n; ++j, x += m, y += m) \
				for (i = 0; i < m; ++i) \
					x[i] = y[p[i] - off]; \
			else \
			for (j = 0; j < n; ++j, x += m, y += m) \
				for (i = 0; i < m; ++i) \
					x[p[i] - off] = y[i]; \
		} else { \
			int i, k0, k1, *q; \
			Matrix_Calloc(q, m, int); \
			if (!invert) \
			for (i = 0; i < m; ++i) \
				q[i] = -(p[i] - off + 1); \
			else \
			for (i = 0; i < m; ++i) \
				q[p[i] - off] = -(i + 1); \
			for (i = 0; i < m; ++i) { \
				if (q[i] > 0) \
					continue; \
				k0 = i; \
				q[k0] = -q[k0]; \
				k1 = q[k0] - 1; \
				while (q[k1] < 0) { \
					c##swap2(n_, x + k0, m_, x + k1, m_); \
					k0 = k1; \
					q[k0] = -q[k0]; \
					k1 = q[k0] - 1; \
				} \
			} \
			Matrix_Free(q, m); \
		} \
	} \
	return; \
} \
 \
void \
c##symperm2(c##TYPE *x, const c##TYPE *y, \
            int n, char uplo, const int *p, int off, int invert) \
{ \
	if (n <= 0) \
		return; \
	size_t n_ = (size_t) n; \
	if (!p) { \
		if (y) \
			memcpy(x, y, sizeof(c##TYPE) * n_ * n_); \
	} else { \
		if (y) { \
			int j; \
			if (!invert) \
			for (j = 0; j < n; ++j) \
				c##symcopyr2(x, y, n_, uplo, (size_t) j, (size_t) (p[j] - off)); \
			else \
			for (j = 0; j < n; ++j) \
				c##symcopyr2(x, y, n_, uplo, (size_t) (p[j] - off), (size_t) j); \
		} else { \
			int j, k0, k1, *q; \
			Matrix_Calloc(q, n, int); \
			if (!invert) \
			for (j = 0; j < n; ++j) \
				q[j] = -(p[j] - off + 1); \
			else \
			for (j = 0; j < n; ++j) \
				q[p[j] - off] = -(j + 1); \
			for (j = 0; j < n; ++j) { \
				if (q[j] > 0) \
					continue; \
				k0 = j; \
				q[k0] = -q[k0]; \
				k1 = q[k0] - 1; \
				while (q[k1] < 0) { \
					c##symswapr2(x, n_, uplo, (size_t) k0, (size_t) k1); \
					k0 = k1; \
					q[k0] = -q[k0]; \
					k1 = q[k0] - 1; \
				} \
			} \
			Matrix_Free(q, n); \
		} \
	} \
	return; \
} \
 \
void \
c##symperm1(c##TYPE *x, const c##TYPE *y, \
            int n, char uplo, const int *p, int off, int invert) \
{ \
	if (n <= 0) \
		return; \
	size_t n_ = (size_t) n; \
	if (!p) { \
		if (y) \
			memcpy(x, y, sizeof(c##TYPE) * PACKED_LENGTH(n_)); \
	} else { \
		if (y) { \
			int j; \
			if (!invert) \
			for (j = 0; j < n; ++j) \
				c##symcopyr1(x, y, n_, uplo, (size_t) j, (size_t) (p[j] - off)); \
			else \
			for (j = 0; j < n; ++j) \
				c##symcopyr1(x, y, n_, uplo, (size_t) (p[j] - off), (size_t) j); \
		} else { \
			int j, k0, k1, *q; \
			Matrix_Calloc(q, n, int); \
			if (!invert) \
			for (j = 0; j < n; ++j) \
				q[j] = -(p[j] - off + 1); \
			else \
			for (j = 0; j < n; ++j) \
				q[p[j] - off] = -(j + 1); \
			for (j = 0; j < n; ++j) { \
				if (q[j] > 0) \
					continue; \
				k0 = j; \
				q[k0] = -q[k0]; \
				k1 = q[k0] - 1; \
				while (q[k1] < 0) { \
					c##symswapr1(x, n_, uplo, (size_t) k0, (size_t) k1); \
					k0 = k1; \
					q[k0] = -q[k0]; \
					k1 = q[k0] - 1; \
				} \
			} \
			Matrix_Free(q, n); \
		} \
	} \
	return; \
} \
 \
void \
c##pack2(c##TYPE *x, const c##TYPE *y, \
         size_t n, char uplo, char trans, char diag) \
{ \
	size_t i, j; \
	c##TYPE *tmp = x; \
	if (uplo == 'U') { \
	for (j = 0; j < n; y += n - (++j)) \
		for (i = 0; i <= j; ++i) \
			*(x++) = *(y++); \
	if (diag != '\0') { \
		if (diag != 'N') \
		for (j = 0, x = tmp; j < n; x += (++j) + 1) \
			c##SET_UNIT(*x); \
	} else if (trans == 'C') \
		for (j = 0, x = tmp; j < n; x += (++j) + 1) \
			c##SET_PROJ_REAL(*x); \
	} else { \
	for (j = 0; j < n; y += (++j)) \
		for (i = j; i < n; ++i) \
			*(x++) = *(y++); \
	if (diag != '\0') { \
		if (diag != 'N') \
		for (j = 0, x = tmp; j < n; x += n - (j++)) \
			c##SET_UNIT(*x); \
	} else if (trans == 'C') \
		for (j = 0, x = tmp; j < n; x += n - (j++)) \
			c##SET_PROJ_REAL(*x); \
	} \
	return; \
} \
 \
void \
c##pack1(c##TYPE *x, const c##TYPE *y, \
         size_t n, char uplo, char trans, char diag) \
{ \
	size_t i, j; \
	if (diag != '\0') { \
		c##TYPE *tmp = x; \
		if (uplo == 'U') { \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i <= j; ++i) \
				*(x++) = *(y++); \
			for (i = j + 1; i < n; ++i) \
				*(x++) = c##ZERO; \
		} \
		} else { \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i < j; ++i) \
				*(x++) = c##ZERO; \
			for (i = j; i < n; ++i) \
				*(x++) = *(y++); \
		} \
		} \
		if (diag != 'N') { \
		size_t dx = n + 1; \
		for (j = 0, x = tmp; j < n; ++j, x += dx) \
			c##SET_UNIT(*x); \
		} \
	} else if (trans == 'C') { \
		c##TYPE *u = x, *l = x; \
		if (uplo == 'U') { \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i < j; ++i) { \
				c##ASSIGN_IDEN(*u, *y); \
				c##ASSIGN_CONJ(*l, *y); \
				u += 1; y += 1; l += n; \
			} \
			c##ASSIGN_PROJ_REAL(*u, *y); \
			u += 1; y += 1; \
			u += n - j - 1; l = x + j + 1; \
		} \
		} else { \
		for (j = 0; j < n; ++j) { \
			l += j; u = l + n; \
			c##ASSIGN_PROJ_REAL(*l, *y); \
			l += 1; y += 1; \
			for (i = j + 1; i < n; ++i) { \
				c##ASSIGN_IDEN(*l, *y); \
				c##ASSIGN_CONJ(*u, *y); \
				l += 1; y += 1; u += n; \
			} \
		} \
		} \
	} else { \
		c##TYPE *u = x, *l = x; \
		if (uplo == 'U') { \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i < j; ++i) { \
				c##ASSIGN_IDEN(*u, *y); \
				c##ASSIGN_IDEN(*l, *y); \
				u += 1; y += 1; l += n; \
			} \
			c##ASSIGN_IDEN(*u, *y); \
			u += 1; y += 1; \
			u += n - j - 1; l = x + j + 1; \
		} \
		} else { \
		for (j = 0; j < n; ++j) { \
			l += j; u = l + n; \
			c##ASSIGN_IDEN(*l, *y); \
			l += 1; y += 1; \
			for (i = j + 1; i < n; ++i) { \
				c##ASSIGN_IDEN(*l, *y); \
				c##ASSIGN_IDEN(*u, *y); \
				l += 1; y += 1; u += n; \
			} \
		} \
		} \
	} \
	return; \
} \
 \
void \
c##force2(c##TYPE *x, const c##TYPE *y, \
          size_t n, char uplo, char trans, char diag) \
{ \
	size_t i, j; \
	if (diag < '\0') { \
		if (diag != -'N') { \
			size_t dx = n + 1; \
			memset(x, 0, sizeof(c##TYPE) * n * n); \
			for (j = 0; j < n; ++j, x += dx) \
				c##SET_UNIT(*x); \
		} else if (y) { \
			size_t dx = n + 1; \
			size_t dy = (trans) ? 1 : dx; \
			memset(x, 0, sizeof(c##TYPE) * n * n); \
			for (j = 0; j < n; ++j, x += dx, y += dy) \
				*x = *y; \
		} else { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) \
					*(x++) = c##ZERO; \
				x++; \
				for (i = j + 1; i < n; ++i) \
					*(x++) = c##ZERO; \
			} \
		} \
	} else if (diag > '\0') { \
		c##TYPE *tmp = x; \
		if (uplo == 'U') { \
		for (j = 0; j < n; ++j) { \
			if (!y) \
				x += j + 1; \
			else { \
				for (i = 0; i <= j; ++i) \
					*(x++) = *(y++); \
				y += n - j - 1; \
			} \
			for (i = j + 1; i < n; ++i) \
				*(x++) = c##ZERO; \
		} \
		} else { \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i < j; ++i) \
				*(x++) = c##ZERO; \
			if (!y) \
				x += n - j; \
			else { \
				y += j; \
				for (i = j; i < n; ++i) \
					*(x++) = *(y++); \
			} \
		} \
		} \
		if (diag != 'N') { \
		size_t dx = n + 1; \
		for (j = 0, x = tmp; j < n; ++j, x += dx) \
			c##SET_UNIT(*x); \
		} \
	} else if (trans == 'C') { \
		c##TYPE *u = x, *l = x; \
		if (uplo == 'U') { \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i < j; ++i) { \
				if (!y) \
					c##ASSIGN_CONJ(*l, *u); \
				else { \
					c##ASSIGN_IDEN(*u, *y); \
					c##ASSIGN_CONJ(*l, *y); \
					y += 1; \
				} \
				u += 1; l += n; \
			} \
			if (!y) \
				c##SET_PROJ_REAL(*u); \
			else { \
				c##ASSIGN_PROJ_REAL(*u, *y); \
				y += 1; \
				y += n - j - 1; \
			} \
			u += 1; \
			u += n - j - 1; l = x + j + 1; \
		} \
		} else { \
		for (j = 0; j < n; ++j) { \
			l += j; u = l + n; \
			if (!y) \
				c##SET_PROJ_REAL(*l); \
			else { \
				y += j; \
				c##ASSIGN_PROJ_REAL(*l, *y); \
				y += 1; \
			} \
			l += 1; \
			for (i = j + 1; i < n; ++i) { \
				if (!y) \
					c##ASSIGN_CONJ(*u, *l); \
				else { \
					c##ASSIGN_IDEN(*l, *y); \
					c##ASSIGN_CONJ(*u, *y); \
					y += 1; \
				} \
				l += 1; u += n; \
			} \
		} \
		} \
	} else { \
		c##TYPE *u = x, *l = x; \
		if (uplo == 'U') { \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i < j; ++i) { \
				if (!y) \
					c##ASSIGN_IDEN(*l, *u); \
				else { \
					c##ASSIGN_IDEN(*u, *y); \
					c##ASSIGN_IDEN(*l, *y); \
					y += 1; \
				} \
				u += 1; l += n; \
			} \
			if (y) { \
				c##ASSIGN_IDEN(*u, *y); \
				y += 1; \
				y += n - j - 1; \
			} \
			u += 1; \
			u += n - j - 1; l = x + j + 1; \
		} \
		} else { \
		for (j = 0; j < n; ++j) { \
			l += j; u = l + n; \
			if (y) { \
				y += j; \
				c##ASSIGN_IDEN(*l, *y); \
				y += 1; \
			} \
			l += 1; \
			for (i = j + 1; i < n; ++i) { \
				if (!y) \
					c##ASSIGN_IDEN(*u, *l); \
				else { \
					c##ASSIGN_IDEN(*l, *y); \
					c##ASSIGN_IDEN(*u, *y); \
					y += 1; \
				} \
				l += 1; u += n; \
			} \
		} \
		} \
	} \
	return; \
} \
 \
void \
c##force1(c##TYPE *x, const c##TYPE *y, \
          size_t n, char uplo, char trans, char diag) \
{ \
	size_t i, j; \
	if (uplo == 'U') { \
	if (diag < '\0') { \
		if (diag != -'N') { \
			memset(x, 0, sizeof(c##TYPE) * PACKED_LENGTH(n)); \
			for (j = 0; j < n; x += (++j) + 1) \
				*x = c##UNIT; \
		} else if (y) { \
			memset(x, 0, sizeof(c##TYPE) * PACKED_LENGTH(n)); \
			for (j = 0; j < n; ++j, x += j + 1, y += (trans) ? 1 : j + 1) \
				*x = *y; \
		} else { \
			for (j = 0; j < n; ++j) { \
				for (i = 0; i < j; ++i) \
					*(x++) = c##ZERO; \
				x++; \
			} \
		} \
	} else if (diag > '\0') { \
		if (y) \
			memcpy(x, y, sizeof(c##TYPE) * PACKED_LENGTH(n)); \
		if (diag != 'N') \
		for (j = 0; j < n; x += (++j) + 1) \
			c##SET_UNIT(*x); \
	} else if (trans == 'C') { \
		if (y) \
			memcpy(x, y, sizeof(c##TYPE) * PACKED_LENGTH(n)); \
		for (j = 0; j < n; x += (++j) + 1) \
			c##SET_PROJ_REAL(*x); \
	} \
	} else { \
	if (diag < '\0') { \
		if (diag != -'N') { \
			memset(x, 0, sizeof(c##TYPE) * PACKED_LENGTH(n)); \
			for (j = 0; j < n; x += n - (j++)) \
				*x = c##UNIT; \
		} else if (y) { \
			memset(x, 0, sizeof(c##TYPE) * PACKED_LENGTH(n)); \
			for (j = 0; j < n; x += n - j, y += (trans) ? 1 : n - j, j++) \
				*x = *y; \
		} else { \
			for (j = 0; j < n; ++j) { \
				x++; \
				for (i = j + 1; i < n; ++i) \
					*(x++) = c##ZERO; \
			} \
		} \
	} else if (diag > '\0') { \
		if (y) \
			memcpy(x, y, sizeof(c##TYPE) * PACKED_LENGTH(n)); \
		if (diag != 'N') \
		for (j = 0; j < n; x += n - (j++)) \
			c##SET_UNIT(*x); \
	} else if (trans == 'C') { \
		if (y) \
			memcpy(x, y, sizeof(c##TYPE) * PACKED_LENGTH(n)); \
		for (j = 0; j < n; x += n - (j++)) \
			c##SET_PROJ_REAL(*x); \
	} \
	} \
	return; \
} \
 \
void \
c##trans2(c##TYPE *x, const c##TYPE *y, \
          size_t m, size_t n, char trans) \
{ \
	size_t i, j, mn = m * n, mn1s = mn - 1; \
	if (trans == 'N') \
		memcpy(x, y, sizeof(c##TYPE) * mn); \
	else if (trans == 'C') \
		for (j = 0; j < m; ++j, y -= mn1s) \
			for (i = 0; i < n; ++i, x += 1, y += m) \
				c##ASSIGN_CONJ(*x, *y); \
	else \
		for (j = 0; j < m; ++j, y -= mn1s) \
			for (i = 0; i < n; ++i, x += 1, y += m) \
				c##ASSIGN_IDEN(*x, *y); \
	return; \
} \
 \
void \
c##trans1(c##TYPE *x, const c##TYPE *y, \
          size_t n, char uplo, char trans) \
{ \
	c##TYPE tmp; \
	size_t i, j; \
	if (trans == 'N') { \
		memcpy(x, y, sizeof(c##TYPE) * PACKED_LENGTH(n)); \
	} else if (trans == 'C') { \
		if (uplo == 'U') \
		for (j = 0; j < n; ++j) { \
			for (i = j; i < n; ++i) { \
				tmp = *(y + PACKED_AR21_UP(j, i)); \
				c##ASSIGN_CONJ(*x, tmp); \
				x += 1; \
			} \
		} \
		else \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i <= j; ++i) { \
				tmp = *(y + PACKED_AR21_LO(j, i, n)); \
				c##ASSIGN_CONJ(*x, tmp); \
				x += 1; \
			} \
		} \
	} else { \
		if (uplo == 'U') \
		for (j = 0; j < n; ++j) { \
			for (i = j; i < n; ++i) { \
				tmp = *(y + PACKED_AR21_UP(j, i)); \
				c##ASSIGN_IDEN(*x, tmp); \
				x += 1; \
			} \
		} \
		else \
		for (j = 0; j < n; ++j) { \
			for (i = 0; i <= j; ++i) { \
				tmp = *(y + PACKED_AR21_LO(j, i, n)); \
				c##ASSIGN_IDEN(*x, tmp); \
				x += 1; \
			} \
		} \
	} \
	return; \
} \
 \
void \
c##band2(c##TYPE *x, const c##TYPE *y, \
         size_t m, size_t n, size_t a, size_t b) \
{ \
	if (m == 0 || n == 0) \
		return; \
	size_t d; \
	if (a > b || a >= m + n || b == 0) { \
		d = m * n; \
		memset(x, 0, sizeof(c##TYPE) * d); \
		return; \
	} \
	a = (a == 0) ? 1 : a; \
	b = (b >= m + n) ? m + n - 1 : b; \
	size_t i, j, i0, i1, \
		j0 = (a < m) ? 0 : a - m, \
		j1 = (b < n) ? b : n; \
	if (j0 > 0) { \
		d = m * j0; \
		memset(x, 0, sizeof(c##TYPE) * d); \
		x += d; if (y) y += d; \
	} \
	for (j = j0; j < j1; ++j) { \
		i0 = (b <= m + j) ? m + j - b     : 0; \
		i1 = (a >= j + 1) ? m + j - a + 1 : m; \
		for (i =  0; i < i0; ++i) \
			x[i] = c##ZERO; \
		if (y) { \
		for (i = i0; i < i1; ++i) \
			x[i] = y[i]; \
		y += m; \
		} \
		for (i = i1; i <  m; ++i) \
			x[i] = c##ZERO; \
		x += m; \
	} \
	if (j1 < n) { \
		d = m * (n - j1); \
		memset(x, 0, sizeof(c##TYPE) * d); \
		x += d; if (y) y += d; \
	} \
	return; \
} \
 \
void \
c##band1(c##TYPE *x, const c##TYPE *y, \
         size_t n, char uplo, size_t a, size_t b) \
{ \
	if (n == 0) \
		return; \
	if (uplo == 'U') \
		a = (a < n) ? n : a; \
	else \
		b = (b > n) ? n : b; \
	size_t d; \
	if (a > b || a >= n + n || b == 0) { \
		d = PACKED_LENGTH(n); \
		memset(x, 0, sizeof(c##TYPE) * d); \
		return; \
	} \
	if (uplo == 'U') \
		b = (b >= n + n) ? n + n - 1 : b; \
	else \
		a = (a == 0) ? 1 : a; \
	size_t i, j, i0, i1, \
		j0 = (a < n) ? 0 : a - n, \
		j1 = (b < n) ? b : n; \
	if (uplo == 'U') { \
		if (j0 > 0) { \
			d = PACKED_LENGTH(j0); \
			memset(x, 0, sizeof(c##TYPE) * d); \
			x += d; if (y) y += d; \
		} \
		for (j = j0; j < j1; ++j) { \
			i0 = (b <= n + j) ? n + j - b     :     0; \
			i1 = (a >= n    ) ? n + j - a + 1 : j + 1; \
			for (i =  0; i < i0; ++i) \
				x[i] = c##ZERO; \
			if (y) { \
			for (i = i0; i < i1; ++i) \
				x[i] = y[i]; \
			y += j + 1; \
			} \
			for (i = i1; i <= j; ++i) \
				x[i] = c##ZERO; \
			x += j + 1; \
		} \
		if (j1 < n) { \
			d = PACKED_LENGTH(n) - PACKED_LENGTH(j1); \
			memset(x, 0, sizeof(c##TYPE) * d); \
			x += d; if (y) y += d; \
		} \
	} else { \
		if (j0 > 0) { \
			d = PACKED_LENGTH(n) - PACKED_LENGTH(j0); \
			memset(x, 0, sizeof(c##TYPE) * d); \
			x += d; if (y) y += d; \
		} \
		for (j = j0; j < j1; ++j) { \
			i0 = (b <= n    ) ? n + j - b     : j; \
			i1 = (a >= j + 1) ? n + j - a + 1 : n; \
			for (i =  j; i < i0; ++i) \
				x[i - j] = c##ZERO; \
			if (y) { \
			for (i = i0; i < i1; ++i) \
				x[i - j] = y[i - j]; \
			y += n - j; \
			} \
			for (i = i1; i <  n; ++i) \
				x[i - j] = c##ZERO; \
			x += n - j; \
		} \
		if (j1 < n) { \
			d = PACKED_LENGTH(n - j1); \
			memset(x, 0, sizeof(c##TYPE) * d); \
			x += d; if (y) y += d; \
		} \
	} \
	return; \
} \
 \
void \
c##sptrans(int *p0, int *i0, c##TYPE *x0, \
           int *p1, int *i1, c##TYPE *x1, \
           int m, int n, char trans) \
{ \
	int i, j, k, kend, nnz = p0[n], *p_ = NULL; \
	Matrix_Calloc(p_, m, int); \
	p0++; *(p1++) = 0; \
	memset(pp1, 0, sizeof(int) * (size_t) m); \
	for (k = 0; k < nnz; ++k) \
		++p1[i0[k]]; \
	for (i = 0; i < m; ++i) \
		p1[i] += (p_[i] = p1[i - 1]); \
	if (!x1) \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = p0[j]; \
			while (k < kend) { \
				i = i0[k]; \
				i1[p_[i]++] = j; \
				++k; \
			} \
		} \
	else if (trans == 'C') \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = p0[j]; \
			while (k < kend) { \
				i = i0[k]; \
				c##ASSIGN_CONJ(x1[p_[i]], x0[k]); \
				i1[p_[i]++] = j; \
				++k; \
			} \
		} \
	else \
		for (j = 0, k = 0; j < n; ++j) { \
			kend = p0[j]; \
			while (k < kend) { \
				i = i0[k]; \
				c##ASSIGN_IDEN(x1[p_[i]], x0[k]); \
				i1[p_[i]++] = j; \
				++k; \
			} \
		} \
	Matrix_Free(p_, m); \
	return; \
} \


TEMPLATE(i)
TEMPLATE(d)
TEMPLATE(z)

void zvreal(Rcomplex *x, const Rcomplex *y, size_t n)
{
	if (y)
	while (n--) {
		(* x   ).r = (*(y++)).r;
		(*(x++)).i = 0.0;
	}
	else
	while (n--)
		(*(x++)).i = 0.0;
	return;
}

void zvimag(Rcomplex *x, const Rcomplex *y, size_t n)
{
	if (y)
	while (n--) {
		(* x   ).r = 0.0;
		(*(x++)).i = (*(y++)).i;
	}
	else
	while (n--)
		(*(x++)).r = 0.0;
	return;
}

void zvconj(Rcomplex *x, const Rcomplex *y, size_t n)
{
	if (y)
	while (n--) {
		(* x   ).r =  (* y   ).r;
		(*(x++)).i = -(*(y++)).i;
	}
	else
	while (n--) {
		(*x).i = -(*x).i;
		x++;
	}
	return;
}
