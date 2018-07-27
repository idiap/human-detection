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
// SymmetricMatrix.h: interface for the CCovarianceMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_COVARIANCE_MATRIX_H_)
#define _COVARIANCE_MATRIX_H_

#include <cstdio>
#include "stdlib.h"
#include <memory>
#include <cmath>

#include "sym_eig234.h"

#include <ctime>						// clock
#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <fstream>						// file I/O
#include <cmath>						// math includes
#include <iostream>						// I/O streams

#define MIN_DOUBLE_NUM	1.0E-10
#define MAX_WEIGHTED_MEAN_ITR_NUM	8

#define MIN_IDENTITY_VALUE  1.0E-4
#define PARENT_MIN_IDENTITY_VALUE	1.0
//#define MIN_IDENTITY_VALUE  1.0

using namespace std;

#include "OpencvDataConversion.h"

#ifndef SQRT2
#define	SQRT2	1.4142135623730949
#endif

#ifndef sgn
#define	sgn(x)	( (x) >= 0 ? 1 : -1 )
#endif

#include "cv.h"

class CCovarianceMatrix
{
public:
    void Print(const char* file_name=NULL, bool stdprint=true);
	void GeneratedRand(CvRNG rng_state, double range=1.0);
	int* GetRowStartPos();
	void Resize(int mat_size, bool used_avg=true);
	void GetSubMatrix(int *rows, int rows_num, CCovarianceMatrix *sub_cov_mat, int *row_start_pos=NULL);
	void GetNormDiagonalElems(double *norm_diag_elems);

	void Import(const char* file_name);

	void Export(const char* log_msg_fn = NULL, bool upper=false);

	void SetIdentityMatrix(double identity_value = 1.0);

	void SetZero();

	void NormalizeIllumination(double *parent_win_diag_elems, double *parent_avg);

	void GetDiagonalElems(double *diag_elems);

	/* remove local illumination variant w.r.t. the parent window */
	void NormalizeIllumination(CCovarianceMatrix *parent_win_mat);

	/* vec_X ( log_X ( Y ) ) */
	void OrthogonalVectorLog(CCovarianceMatrix *Y_mat, CCovarianceMatrix *vec_log = NULL, CCovarianceMatrix *this_mat_sqrt=NULL, CCovarianceMatrix *this_mat_inv_sqrt=NULL);

	/* scaling a matrix, X = s * X */
	void Scale(double scale);

	/* sum of two covariance matrices, Z = X + Y */
	void Add(CCovarianceMatrix *Y_mat, CCovarianceMatrix *sum_mat = NULL);

	/* finding weighted mean of a set of points on Riemannian manifolds using Gradient Descent algorithm */
	void WeightedMean(CCovarianceMatrix **Y_mats, int mat_num, double *weights=NULL, int *sel_pts=NULL, int itr_num=10, double *dist_sums=NULL, double itr_stop_value=1.0);

	/* computing sum of geodesic distances of a set of covariance matrices */
	double GeodesicDistanceSum(CCovarianceMatrix **Y_mats, int mat_num, double *weights = NULL, int *sel_pts=NULL, CCovarianceMatrix *this_mat_inv_sqrt=NULL);

	/* computing region covariance matrix from a set of feature image with rectangle of interested and mask image */
	void ComputeRegionCovFeature(IplImage **feature_imgs, CvRect roi, IplImage *mask_img = NULL);

	/* normalize columns and rows corresponding to position features */
	void NormalizePositionFeatures(CvRect roi, int x_feature_idx, int y_feature_idx);

	/* do cholesky factorization, X = F * F' */
	void CholeskyFactorization(CvMat *F);

	/* do cholesky factorization, X = F * F' */
	void CholeskyFactorization_2(CvMat *F);

	/* computing generalized eigenvalues */
	double* GeneralizedEigenvalues(CCovarianceMatrix *Y_mat);

	/* dissimilarity between two covariance matrices */
	double CovDistance(CCovarianceMatrix *Y_mat);

	/* update reference covariance matrix using current covariance matrix with a update rate */
	void Update(CCovarianceMatrix *Y_mat, double update_rate, CCovarianceMatrix *updated_mat = NULL);

	/* orthogonal coordinates of a vector y on the tangent space , vec_X ( Y ) */
	void OrthogonalVector(CCovarianceMatrix *y_tangent_vec, CCovarianceMatrix *orth_vec = NULL, CCovarianceMatrix *this_mat_inv_sqrt=NULL);

	/* matrix multiplication, A = A * D, D is a diagonal matrix */
	void Mul_Diag(CvMat *A, CvMat *D);

	/* do SVD decomposition on a covariance matrix,  X = U * D * U',
	   D - diagonal eigenvalue matrix, U - left orthogonal matrix */
	void SVD(CvMat *D, CvMat *U);

	/* matrix multiplicatoin, Z = X * Y * X */
	void Mul2(CCovarianceMatrix *Y_mat, CCovarianceMatrix  *mul2_mat = NULL);

	/* get matrix data */
	double* GetData(int one_upper_zero_full=1);

	/* set the matrix data */
	void SetData(double *data, int length, double *avg_data=NULL);

	/* geodesic distance between two points on Riemannian manifolds */
	double GeodesicDistance(CCovarianceMatrix *Y_mat, CCovarianceMatrix *this_mat_inv_sqrt=NULL);

	/* copy current matrix data into a new matrix */
	void Copy(CCovarianceMatrix *new_mat);

	/* exponential of a matrix, log(X) = U * exp(D) * U', [U D U'] = SVD(X) */
	void Exp(CCovarianceMatrix *exp_mat = NULL);

	/* logarithm of a matrix, log(X) = U * exp(D) * U', [U D U'] = SVD(X) */
	void Log(CCovarianceMatrix *log_mat = NULL);

	/* exp_X(Y) = X^0.5 * exp ( X^-0.5 * Y * X^-0.5 ) * X^0.5 */
	void Exp_Y(CCovarianceMatrix *Y_mat, CCovarianceMatrix *exp_Y_mat = NULL, CCovarianceMatrix *this_mat_sqrt=NULL, CCovarianceMatrix *this_mat_inv_sqrt=NULL);

	/* log_X(Y) = X^0.5 * log ( X^-0.5 * Y * X^-0.5 ) * X^0.5 */
	void Log_Y(CCovarianceMatrix *Y_mat, CCovarianceMatrix *log_Y_mat = NULL, CCovarianceMatrix *this_mat_sqrt=NULL, CCovarianceMatrix *this_mat_inv_sqrt=NULL);

	/* square root of a matrix, sqrt(X) = U * sqrt(D) * U', [U D U'] = SVD(X)  */
	void Sqrt(CCovarianceMatrix *sqrt_mat = NULL);

	/* square of a matrix  */
	void Square(CCovarianceMatrix *square_mat = NULL);

	/* inverse of a matrix  */
	void Inverse(CCovarianceMatrix *inv_mat = NULL);

	/* matrix multiplication using fast winograd's algorithm implemented by Jian Yao */
	void Mul(CCovarianceMatrix *mat1, CCovarianceMatrix *mat2 = NULL,  CCovarianceMatrix *mul_mat = NULL);

	/* The function cvTrace returns sum of diagonal elements of the matrix */
	double Trace();

	/* release memories for the matrix */
	void Delete();

	double* SM_ptr;		/* the symmetric positive definite covariance matrix, here,
				   it is represented by an upper triangle part */

	double* AVG_ptr;	/* the mean values */

	int SM_length;		/* the length of the upper triangular part, SM_length = (M_size*(M_size+1))/2 */

	int M_size;		/* the size of covariance matrix */
	int M_size2;		/* M_size2 = M_size * M_size */

	CCovarianceMatrix(int mat_size, bool used_avg=true);

	virtual ~CCovarianceMatrix();

private:
	/* logarithm of a diagonal eigenvalue matrix */
	void Log_Diag(CvMat *D);

	/* exponential of a diagonal eigenvalue matrix */
	void Exp_Diag(CvMat *D);

	/* square root of a diagonal eigenvalue matrix */
	void Sqrt_Diag(CvMat*D);

	/* logarithm of a diagonal eigenvalue matrix */
	void Exp_Diag();

	/* exponential of a diagonal eigenvalue matrix */
	void Log_Diag();

	/* square root of a diagonal eigenvalue matrix */
	void Sqrt_Diag();

};

#endif // !defined(_COVARIANCE_MATRIX_H_)

