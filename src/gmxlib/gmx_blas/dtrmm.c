/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include <math.h>

#include <types/simple.h>

#include "gmx_blas.h"

void 
F77_FUNC(dtrmm,DTRMM)(const char *side, 
       const char *uplo, 
       const char *transa, 
       const char *diag, 
       int *m__, 
       int *n__, 
       double *alpha__, 
       double *a, 
       int *lda__, 
       double *b, 
       int *ldb__)
{
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    int m = *m__;
    int n = *n__;
    int lda = *lda__;
    int ldb = *ldb__;
    double alpha = *alpha__;
    
    /* Local variables */
    int i__, j, k, info;
    double temp;
    int lside;
    int nrowa;
    int upper;
    int nounit;
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    lside = (*side=='L' || *side=='l');
    if (lside) {
	nrowa = m;
    } else {
	nrowa = n;
    }
    nounit = (*diag=='N' || *diag=='n');
    upper = (*uplo=='U' || *uplo=='u');

    info = 0;

    if (n == 0) {
	return;
    }
    if (fabs(alpha)<GMX_DOUBLE_MIN) {
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = 0.;
	    }
	}
	return;
    }
    if (lside) {
	if (*transa=='N' || *transa=='n') {
	    if (upper) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = m;
		    for (k = 1; k <= i__2; ++k) {
			if (fabs(b[k + j * b_dim1])>GMX_DOUBLE_MIN) {
			    temp = alpha * b[k + j * b_dim1];
			    i__3 = k - 1;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * a_dim1];
			    }
			    if (nounit) {
				temp *= a[k + k * a_dim1];
			    }
			    b[k + j * b_dim1] = temp;
			}
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    for (k = m; k >= 1; --k) {
			if (fabs(b[k + j * b_dim1])>GMX_DOUBLE_MIN) {
			    temp = alpha * b[k + j * b_dim1];
			    b[k + j * b_dim1] = temp;
			    if (nounit) {
				b[k + j * b_dim1] *= a[k + k * a_dim1];
			    }
			    i__2 = m;
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
			    }
			}
		    }
		}
	    }
	} else {

	    if (upper) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = m; i__ >= 1; --i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__2 = i__ - 1;
			for (k = 1; k <= i__2; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
			}
			b[i__ + j * b_dim1] = alpha * temp;
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__3 = m;
			for (k = i__ + 1; k <= i__3; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
			}
			b[i__ + j * b_dim1] = alpha * temp;
		    }
		}
	    }
	}
    } else {
	if (*transa=='N' || *transa=='n') {

	    if (upper) {
		for (j = n; j >= 1; --j) {
		    temp = alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__1 = m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
		    }
		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			if (fabs(a[k + j * a_dim1])>GMX_DOUBLE_MIN) {
			    temp = alpha * a[k + j * a_dim1];
			    i__2 = m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    temp = alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
		    }
		    i__2 = n;
		    for (k = j + 1; k <= i__2; ++k) {
			if (fabs(a[k + j * a_dim1])>GMX_DOUBLE_MIN) {
			    temp = alpha * a[k + j * a_dim1];
			    i__3 = m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		}
	    }
	} else {

	    if (upper) {
		i__1 = n;
		for (k = 1; k <= i__1; ++k) {
		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if (fabs(a[j + k * a_dim1])>GMX_DOUBLE_MIN) {
			    temp = alpha * a[j + k * a_dim1];
			    i__3 = m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if (fabs(temp-1.0)>GMX_DOUBLE_EPS) {
			i__2 = m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
			}
		    }
		}
	    } else {
		for (k = n; k >= 1; --k) {
		    i__1 = n;
		    for (j = k + 1; j <= i__1; ++j) {
			if (fabs(a[j + k * a_dim1])>GMX_DOUBLE_MIN) {
			    temp = alpha * a[j + k * a_dim1];
			    i__2 = m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if (fabs(temp-1.0)>GMX_DOUBLE_EPS) {
			i__1 = m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
			}
		    }
		}
	    }
	}
    }

    return;

}


