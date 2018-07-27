/* --- --- ---
 * Copyright (C) 2008--2010 Idiap Research Institute (.....@idiap.ch)
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* Eigen decomposition code for symmetric 3x3 matrices, copied from the public
   domain Java Matrix library JAMA. */

#include <math.h>

#include <ctime>						// clock
#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <fstream>						// file I/O
#include <cmath>						// math includes
#include <iostream>						// I/O streams


#ifdef MAX
#undef MAX
#endif

#define MAX(a, b) ((a)>(b)?(a):(b))

static double hypot2(double x, double y) {
	return sqrt(x*x+y*y);
}

// Symmetric Householder reduction to tridiagonal form.
static void tred2(double V[2][2], double d[2], double e[2], int n) {

	//  This is derived from the Algol procedures tred2 by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	int i, j, k;

	for ( j = 0; j < n; j++) {
		d[j] = V[n-1][j];
	}

	// Householder reduction to tridiagonal form.
	for ( i = n-1; i > 0; i--) {

		// Scale to avoid under/overflow.

		double scale = 0.0;
		double h = 0.0;
		for ( k = 0; k < i; k++) {
			scale = scale + fabs(d[k]);
		}
		if (scale == 0.0) {
			e[i] = d[i-1];
			for ( j = 0; j < i; j++) {
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		} 
		else {

			// Generate Householder vector.

			for ( k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i-1];
			double g = sqrt(h);
			if (f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i-1] = f - g;
			for ( j = 0; j < i; j++) {
				e[j] = 0.0;
			}

			// Apply similarity transformation to remaining columns.

			for ( j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (int k = j+1; k <= i-1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for ( j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for ( j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for ( j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i-1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.
	for ( i = 0; i < n-1; i++) {
		V[n-1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i+1];
		if (h != 0.0) {
			for ( k = 0; k <= i; k++) {
				d[k] = V[k][i+1] / h;
			}
			for ( j = 0; j <= i; j++) {
				double g = 0.0;
				for ( k = 0; k <= i; k++) {
					g += V[k][i+1] * V[k][j];
				}
				for ( k = 0; k <= i; k++) {
					  V[k][j] -= g * d[k];
				}
			}
		}
		for ( k = 0; k <= i; k++) {
			V[k][i+1] = 0.0;
		}
	}

	for ( j = 0; j < n; j++) {
		d[j] = V[n-1][j];
		V[n-1][j] = 0.0;
	}
	
	V[n-1][n-1] = 1.0;
	e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.

static void tql2(double V[2][2], double d[2], double e[2], int n) {

	//  This is derived from the Algol procedures tql2, by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	int i, k, m, l;

	for ( i = 1; i < n; i++) {
		e[i-1] = e[i];
	}
	e[n-1] = 0.0;

	double f = 0.0;
	double tst1 = 0.0;
	//double eps = pow(2.0,-35.0);
	double eps = 1.0e-30;

	for ( l = 0; l < n; l++) {

		// Find small subdiagonal element

		tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
		m = l;
		while (m < n) {
			if (fabs(e[m]) <= eps*tst1) {
				break;
			}
			m++;
		}

		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.

		if (m > l) {
			int iter = 0;
			do {
				iter = iter + 1;  // (Could check iteration count here.)

				// Compute implicit shift

				double g = d[l];
				double p = (d[l+1] - g) / (2.0 * e[l]);
				double r = hypot2(p,1.0);
				if (p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l+1] = e[l] * (p + r);
				double dl1 = d[l+1];
				double h = g - d[l];
				for ( i = l+2; i < n; i++) {
					d[i] -= h;
				}
				f = f + h;

				// Implicit QL transformation.

				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l+1];
				double s = 0.0;
				double s2 = 0.0;
				for ( i = m-1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot2(p,e[i]);
					e[i+1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i+1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.

					for ( k = 0; k < n; k++) {
						h = V[k][i+1];
						V[k][i+1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;

				// Check for convergence.

			} while (fabs(e[l]) > eps*tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}

	// Sort eigenvalues and corresponding vectors.
	
	/*
	int j;
	for ( i = 0; i < n-1; i++) {
		k = i;
		double p = d[i];
		for ( j = i+1; j < n; j++) {
			if (d[j] > p) {
				k = j;
				p = d[j];
			}
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for ( j = 0; j < n; j++) {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
	*/
}

// Symmetric Householder reduction to tridiagonal form.
static void tred3(double V[3][3], double d[3], double e[3], int n) {

	//  This is derived from the Algol procedures tred2 by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	int i, j, k;

	for ( j = 0; j < n; j++) {
		d[j] = V[n-1][j];
	}

	// Householder reduction to tridiagonal form.
	for ( i = n-1; i > 0; i--) {

		// Scale to avoid under/overflow.

		double scale = 0.0;
		double h = 0.0;
		for ( k = 0; k < i; k++) {
			scale = scale + fabs(d[k]);
		}
		if (scale == 0.0) {
			e[i] = d[i-1];
			for ( j = 0; j < i; j++) {
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		} 
		else {

			// Generate Householder vector.

			for ( k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i-1];
			double g = sqrt(h);
			if (f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i-1] = f - g;
			for ( j = 0; j < i; j++) {
				e[j] = 0.0;
			}

			// Apply similarity transformation to remaining columns.

			for ( j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (int k = j+1; k <= i-1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for ( j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for ( j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for ( j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i-1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.
	for ( i = 0; i < n-1; i++) {
		V[n-1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i+1];
		if (h != 0.0) {
			for ( k = 0; k <= i; k++) {
				d[k] = V[k][i+1] / h;
			}
			for ( j = 0; j <= i; j++) {
				double g = 0.0;
				for ( k = 0; k <= i; k++) {
					g += V[k][i+1] * V[k][j];
				}
				for ( k = 0; k <= i; k++) {
					  V[k][j] -= g * d[k];
				}
			}
		}
		for ( k = 0; k <= i; k++) {
			V[k][i+1] = 0.0;
		}
	}

	for ( j = 0; j < n; j++) {
		d[j] = V[n-1][j];
		V[n-1][j] = 0.0;
	}
	
	V[n-1][n-1] = 1.0;
	e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.

static void tql3(double V[3][3], double d[3], double e[3], int n) {

	//  This is derived from the Algol procedures tql2, by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	int i, k, m, l;

	for ( i = 1; i < n; i++) {
		e[i-1] = e[i];
	}
	e[n-1] = 0.0;

	double f = 0.0;
	double tst1 = 0.0;
	//double eps = pow(2.0,-35.0);
	double eps = 1.0e-30;

	for ( l = 0; l < n; l++) {

		// Find small subdiagonal element

		tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
		m = l;
		while (m < n) {
			if (fabs(e[m]) <= eps*tst1) {
				break;
			}
			m++;
		}

		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.

		if (m > l) {
			int iter = 0;
			do {
				iter = iter + 1;  // (Could check iteration count here.)

				// Compute implicit shift

				double g = d[l];
				double p = (d[l+1] - g) / (2.0 * e[l]);
				double r = hypot2(p,1.0);
				if (p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l+1] = e[l] * (p + r);
				double dl1 = d[l+1];
				double h = g - d[l];
				for ( i = l+2; i < n; i++) {
					d[i] -= h;
				}
				f = f + h;

				// Implicit QL transformation.

				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l+1];
				double s = 0.0;
				double s2 = 0.0;
				for ( i = m-1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot2(p,e[i]);
					e[i+1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i+1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.

					for ( k = 0; k < n; k++) {
						h = V[k][i+1];
						V[k][i+1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;

				// Check for convergence.

			} while (fabs(e[l]) > eps*tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}

	// Sort eigenvalues and corresponding vectors.
	
	/*
	int j;
	for ( i = 0; i < n-1; i++) {
		k = i;
		double p = d[i];
		for ( j = i+1; j < n; j++) {
			if (d[j] > p) {
				k = j;
				p = d[j];
			}
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for ( j = 0; j < n; j++) {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
	*/
}

// Symmetric Householder reduction to tridiagonal form.
static void tred4(double V[4][4], double d[4], double e[4], int n) {

	//  This is derived from the Algol procedures tred2 by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	int i, j, k;

	for ( j = 0; j < n; j++) {
		d[j] = V[n-1][j];
	}

	// Householder reduction to tridiagonal form.
	for ( i = n-1; i > 0; i--) {

		// Scale to avoid under/overflow.

		double scale = 0.0;
		double h = 0.0;
		for ( k = 0; k < i; k++) {
			scale = scale + fabs(d[k]);
		}
		if (scale == 0.0) {
			e[i] = d[i-1];
			for ( j = 0; j < i; j++) {
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		} 
		else {

			// Generate Householder vector.

			for ( k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i-1];
			double g = sqrt(h);
			if (f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i-1] = f - g;
			for ( j = 0; j < i; j++) {
				e[j] = 0.0;
			}

			// Apply similarity transformation to remaining columns.

			for ( j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (int k = j+1; k <= i-1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for ( j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for ( j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for ( j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i-1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.
	for ( i = 0; i < n-1; i++) {
		V[n-1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i+1];
		if (h != 0.0) {
			for ( k = 0; k <= i; k++) {
				d[k] = V[k][i+1] / h;
			}
			for ( j = 0; j <= i; j++) {
				double g = 0.0;
				for ( k = 0; k <= i; k++) {
					g += V[k][i+1] * V[k][j];
				}
				for ( k = 0; k <= i; k++) {
					  V[k][j] -= g * d[k];
				}
			}
		}
		for ( k = 0; k <= i; k++) {
			V[k][i+1] = 0.0;
		}
	}

	for ( j = 0; j < n; j++) {
		d[j] = V[n-1][j];
		V[n-1][j] = 0.0;
	}
	
	V[n-1][n-1] = 1.0;
	e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.

static void tql4(double V[4][4], double d[4], double e[4], int n) {

	//  This is derived from the Algol procedures tql2, by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	int i, k, m, l;

	for ( i = 1; i < n; i++) {
		e[i-1] = e[i];
	}
	e[n-1] = 0.0;

	double f = 0.0;
	double tst1 = 0.0;
	//double eps = pow(2.0,-35.0);
	double eps = 1.0e-30;

	for ( l = 0; l < n; l++) {

		// Find small subdiagonal element

		tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
		m = l;
		while (m < n) {
			if (fabs(e[m]) <= eps*tst1) {
				break;
			}
			m++;
		}

		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.

		if (m > l) {
			int iter = 0;
			do {
				iter = iter + 1;  // (Could check iteration count here.)

				// Compute implicit shift

				double g = d[l];
				double p = (d[l+1] - g) / (2.0 * e[l]);
				double r = hypot2(p,1.0);
				if (p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l+1] = e[l] * (p + r);
				double dl1 = d[l+1];
				double h = g - d[l];
				for ( i = l+2; i < n; i++) {
					d[i] -= h;
				}
				f = f + h;

				// Implicit QL transformation.

				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l+1];
				double s = 0.0;
				double s2 = 0.0;
				for ( i = m-1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot2(p,e[i]);
					e[i+1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i+1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.

					for ( k = 0; k < n; k++) {
						h = V[k][i+1];
						V[k][i+1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;

				// Check for convergence.

			} while (fabs(e[l]) > eps*tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}

	// Sort eigenvalues and corresponding vectors.
	
	/*
	int j;
	for ( i = 0; i < n-1; i++) {
		k = i;
		double p = d[i];
		for ( j = i+1; j < n; j++) {
			if (d[j] > p) {
				k = j;
				p = d[j];
			}
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for ( j = 0; j < n; j++) {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
	*/
}


// Symmetric Householder reduction to tridiagonal form.
static void tred_12(double V[12][12], double d[12], double e[12], int n) {

	//  This is derived from the Algol procedures tred2 by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	int i, j, k;

	for ( j = 0; j < n; j++) {
		d[j] = V[n-1][j];
	}

	// Householder reduction to tridiagonal form.
	for ( i = n-1; i > 0; i--) {

		// Scale to avoid under/overflow.

		double scale = 0.0;
		double h = 0.0;
		for ( k = 0; k < i; k++) {
			scale = scale + fabs(d[k]);
		}
		if (scale == 0.0) {
			e[i] = d[i-1];
			for ( j = 0; j < i; j++) {
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		} 
		else {

			// Generate Householder vector.

			for ( k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i-1];
			double g = sqrt(h);
			if (f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i-1] = f - g;
			for ( j = 0; j < i; j++) {
				e[j] = 0.0;
			}

			// Apply similarity transformation to remaining columns.

			for ( j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (int k = j+1; k <= i-1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for ( j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for ( j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for ( j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i-1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.
	for ( i = 0; i < n-1; i++) {
		V[n-1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i+1];
		if (h != 0.0) {
			for ( k = 0; k <= i; k++) {
				d[k] = V[k][i+1] / h;
			}
			for ( j = 0; j <= i; j++) {
				double g = 0.0;
				for ( k = 0; k <= i; k++) {
					g += V[k][i+1] * V[k][j];
				}
				for ( k = 0; k <= i; k++) {
					  V[k][j] -= g * d[k];
				}
			}
		}
		for ( k = 0; k <= i; k++) {
			V[k][i+1] = 0.0;
		}
	}

	for ( j = 0; j < n; j++) {
		d[j] = V[n-1][j];
		V[n-1][j] = 0.0;
	}
	
	V[n-1][n-1] = 1.0;
	e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.

static void tql_12(double V[12][12], double d[12], double e[12], int n) {

	//  This is derived from the Algol procedures tql2, by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	int i, k, m, l;

	for ( i = 1; i < n; i++) {
		e[i-1] = e[i];
	}
	e[n-1] = 0.0;

	double f = 0.0;
	double tst1 = 0.0;
	//double eps = pow(2.0,-35.0);
	double eps = 1.0e-30;

	for ( l = 0; l < n; l++) {

		// Find small subdiagonal element

		tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
		m = l;
		while (m < n) {
			if (fabs(e[m]) <= eps*tst1) {
				break;
			}
			m++;
		}

		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.

		if (m > l) {
			int iter = 0;
			do {
				iter = iter + 1;  // (Could check iteration count here.)

				// Compute implicit shift

				double g = d[l];
				double p = (d[l+1] - g) / (2.0 * e[l]);
				double r = hypot2(p,1.0);
				if (p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l+1] = e[l] * (p + r);
				double dl1 = d[l+1];
				double h = g - d[l];
				for ( i = l+2; i < n; i++) {
					d[i] -= h;
				}
				f = f + h;

				// Implicit QL transformation.

				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l+1];
				double s = 0.0;
				double s2 = 0.0;
				for ( i = m-1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot2(p,e[i]);
					e[i+1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i+1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.

					for ( k = 0; k < n; k++) {
						h = V[k][i+1];
						V[k][i+1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;

				// Check for convergence.

			} while (fabs(e[l]) > eps*tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}

	// Sort eigenvalues and corresponding vectors.
	
	
	int j;
	for ( i = 0; i < n-1; i++) {
		k = i;
		double p = d[i];
		for ( j = i+1; j < n; j++) {
			if (d[j] > p) {
				k = j;
				p = d[j];
			}
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for ( j = 0; j < n; j++) {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
	
}

void eigen2_decomposition(double *A, double *V, double *d) {
	/*	
	double xx = A[0], xy = A[1], yy = A[2];
	double b = -xx-yy;
	double xy2 = xy*xy;
	double c = xx*yy - xy2;

	double b4ac = sqrt(b*b/4.0-c);
	double b2 = -b/2.0;
	
	d[0] = b2 - b4ac;
	d[1] = b2 + b4ac;

	double xx_r;

	xx_r = d[0] - xx;

	double v_norm = sqrt(xx_r*xx_r + xy2);
	
	V[3] = xy/v_norm;
	V[0] = -V[3];
	V[1] = V[2] = -xx_r/v_norm;

	return;
	*/

	/*

	xx_r = yy - d[1];

	v_norm = sqrt(xx_r*xx_r + xy2);

	V[1] = xx_r/v_norm;
	V[3] = -xy/v_norm;
	*/
	

	
	double e[2];
	double _V[2][2];
	
	_V[0][0] = A[0];
	_V[0][1] = A[1];
	_V[1][0] = A[1];
	_V[1][1] = A[2];
	
	tred2(_V, d, e, 2);
	tql2(_V, d, e, 2);

	V[0] = _V[0][0];
	V[1] = _V[1][0];
	V[2] = _V[0][1];
	V[3] = _V[1][1];
		
}

void eigen3_decomposition(double *A, double *V, double *d) {

	double e[3];
	double _V[3][3];
	
	_V[0][0] = A[0];
	_V[0][1] = A[1];
	_V[0][2] = A[2];
	_V[1][0] = A[1];
	_V[1][1] = A[3];
	_V[1][2] = A[4];
	_V[2][0] = A[2];
	_V[2][1] = A[4];
	_V[2][2] = A[5];

	tred3(_V, d, e, 3);
	tql3(_V, d, e, 3);

	V[0] = _V[0][0];
	V[1] = _V[1][0];
	V[2] = _V[2][0];
	V[3] = _V[0][1];
	V[4] = _V[1][1];
	V[5] = _V[2][1];
	V[6] = _V[0][2];
	V[7] = _V[1][2];
	V[8] = _V[2][2];
}

void eigen4_decomposition(double *A, double *V, double *d) {

	double e[4];
	double _V[4][4];
	
	_V[0][0] = A[0];
	_V[0][1] = A[1];
	_V[0][2] = A[2];
	_V[0][3] = A[3];
	_V[1][0] = A[1];
	_V[1][1] = A[4];
	_V[1][2] = A[5];
	_V[1][3] = A[6];
	_V[2][0] = A[2];
	_V[2][1] = A[5];
	_V[2][2] = A[7];
	_V[2][3] = A[8];
	_V[3][0] = A[3];
	_V[3][1] = A[6];
	_V[3][2] = A[8];
	_V[3][3] = A[9];

	tred4(_V, d, e, 4);
	tql4(_V, d, e, 4);

	V[0] = _V[0][0];
	V[1] = _V[1][0];
	V[2] = _V[2][0];
	V[3] = _V[3][0];
	V[4] = _V[0][1];
	V[5] = _V[1][1];
	V[6] = _V[2][1];
	V[7] = _V[3][1];
	V[8] = _V[0][2];
	V[9] = _V[1][2];
	V[10] = _V[2][2];
	V[11] = _V[3][2];
	V[12] = _V[0][3];
	V[13] = _V[1][3];
	V[14] = _V[2][3];
	V[15] = _V[3][3];
}


void eigen_n_decomposition(double *A, double *V, double *d, int n) {
	if ( n > 12 ) {
		printf("The dimension of matrix to be decomposied must be less than 12!\n");
		exit(1);
	}

	double e[12];
	double _V[12][12];

	double *tmp_A = A;
	int i, j;

	for ( i = 0 ; i < n ; i++ )
		for ( j = 0 ; j < n ; j++ )
			_V[i][j] = *tmp_A++;

	tred_12(_V, d, e, n);
	tql_12(_V, d, e, n);

	double *tmp_V = V;
	for ( i = 0 ; i < n ; i++ )
		for ( j = 0 ; j < n ; j++ )
			*tmp_V++ = _V[j][i];
}


