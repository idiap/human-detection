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
// SymmetricMatrix.cpp: implementation of the CCovarianceMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "CovarianceMatrix.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCovarianceMatrix::CCovarianceMatrix(int mat_size, bool used_avg)
{
	M_size = mat_size;
	SM_length = (mat_size*(mat_size+1))/2;
	M_size2 = M_size*M_size;

	SM_ptr = new double[SM_length];
	AVG_ptr = NULL;
	if ( used_avg )
        AVG_ptr = new double[M_size];

	memset(SM_ptr, 0, sizeof(double)*SM_length);
}

CCovarianceMatrix::~CCovarianceMatrix()
{
	Delete();
}

void CCovarianceMatrix::Delete()
{
	if ( SM_ptr ) {
		delete [] SM_ptr;
		SM_ptr = NULL;
	}
	if ( AVG_ptr ) {
		delete [] AVG_ptr;
		AVG_ptr = NULL;
	}
}

double CCovarianceMatrix::Trace()
{
	double trace = 0;
	double *tmp_ptr = SM_ptr;

	int dy;
	for ( dy = M_size ; dy > 0 ; dy-- ) {
		trace += *tmp_ptr;
		tmp_ptr += dy;
	}

	return trace;
}


void CCovarianceMatrix::Copy(CCovarianceMatrix *new_mat)
{
	memcpy(new_mat->SM_ptr, this->SM_ptr, sizeof(double)*SM_length);
	if ( new_mat->AVG_ptr && this->AVG_ptr )
        memcpy(new_mat->AVG_ptr, this->AVG_ptr, sizeof(double)*M_size);
}

void CCovarianceMatrix::Inverse(CCovarianceMatrix *inv_mat)
{
	if ( this->M_size == 1 ) {
		if ( inv_mat )
			inv_mat->SM_ptr[0] = 1.0/SM_ptr[0];
		else
			SM_ptr[0] = 1.0/SM_ptr[0];
	}
	else if ( this->M_size == 2 ) {
		double a0, a1, a2;

		double *tmp_SM_ptr = this->SM_ptr;

		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;

		if ( inv_mat )
			tmp_SM_ptr = inv_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		double D = a0*a2 - a1*a1;

		*tmp_SM_ptr++ = a2/D;
		*tmp_SM_ptr++ = -a1/D;
		*tmp_SM_ptr++ = a0/D;
	}
	else if ( this->M_size == 3 ) {
		double a0, a1, a2, a3, a4, a5;

		double *tmp_SM_ptr = this->SM_ptr;

		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;

		if ( inv_mat )
			tmp_SM_ptr = inv_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		double t2 = a4*a4;
		double t4 = a0*a3;
		double t7 = a1*a2;
		double t10 = a1*a1;
		double t12 = a2*a2;
		double t15 = 1.0/(t4*a5-a0*t2+2.0*t7*a4-t10*a5-t12*a3);

		*tmp_SM_ptr++ = (a3*a5-t2)*t15;
		*tmp_SM_ptr++ = -(-a2*a4+a1*a5)*t15;
		*tmp_SM_ptr++ = (a1*a4-a2*a3)*t15;
		*tmp_SM_ptr++ = (a0*a5-t12)*t15;
		*tmp_SM_ptr++ = -(a0*a4-t7)*t15;
		*tmp_SM_ptr++ = (t4-t10)*t15;
	}
	else if ( this->M_size == 4 ) {
		double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;
		a6 = *tmp_SM_ptr++;
		a7 = *tmp_SM_ptr++;
		a8 = *tmp_SM_ptr++;
		a9 = *tmp_SM_ptr++;

		if ( inv_mat )
			tmp_SM_ptr = inv_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		double t1 = a8*a8;
		double t5 = a6*a6;
		double t7 = a5*a8;
		double t10 = a5*a5;
		double t13 = a0*t1;
		double t15 = a0*a9;
		double t18 = a0*t5;
		double t20 = a0*a8;
		double t21 = a5*a6;
		double t25 = a2*a8;
		double t26 = a3*a4;
		double t29 = a2*a2;
		double t30 = a9*t29;
		double t32 = a3*a3;
		double t33 = t32*a4;
		double t35 = a1*a6;
		double t40 = a1*a3;
		double t43 = a1*a1;
		double t45 = a3*a6;
		double t46 = a2*a5;
		double t51 = a1*a2;
		double t56 = a1*a8;
		double t57 = a2*a6;
		double t60 = t13*a4-t15*a4*a7+t18*a7-2.0*t20*t21+t15*t10-2.0*t25*t26+t30*a4+t33*
			a7-2.0*t35*a3*a7-t29*t5+2.0*t7*t40-t1*t43+2.0*t45*t46-t10*t32-2.0*a9*a5*t51+a9*
			a7*t43+2.0*t56*t57;
		double t61 = 1.0/t60;
		double t71 = (a1*t1-t7*a3+t46*a9-a1*a9*a7+t45*a7-t57*a8)*t61;
		double t72 = a2*a4;
		double t78 = a1*a5;
		double t82 = (-t72*a9+a8*a3*a4+a2*t5-t56*a6+t78*a9-t45*a5)*t61;
		double t90 = (-t72*a8+t26*a7+t56*a5+t21*a2-a3*t10-t35*a7)*t61;
		double t104 = (-t20*a6+t15*a5+t56*a3-a5*t32-t51*a9+t57*a3)*t61;
		double t105 = a0*a5;
		double t115 = (t105*a8-a0*a6*a7+t29*a6-a3*a5*a2+t40*a7-t56*a2)*t61;
		double t130 = (t20*a4-a8*t43-a2*a3*a4-t105*a6+t78*a3+t51*a6)*t61;

		*tmp_SM_ptr++ = (t1*a4-a9*a4*a7+t5*a7-2.0*t7*a6+a9*t10)*t61;
		*tmp_SM_ptr++ = -t71;
		*tmp_SM_ptr++ = -t82;
		*tmp_SM_ptr++ = t90;
		*tmp_SM_ptr++ = (t13-t15*a7+t32*a7+t30-2.0*t25*a3)*t61;
		*tmp_SM_ptr++ = t104;
		*tmp_SM_ptr++ = -t115;
		*tmp_SM_ptr++ = -(t15*a4-t18-t33+2.0*t35*a3-a9*t43)*t61;
		*tmp_SM_ptr++ = t130;
		*tmp_SM_ptr++ = -t61*(a7*a0*a4-a7*t43-t29*a4-a0*t10+2.0*t78*a2);
	}
	else {
		if ( inv_mat ) {
			CvMat* A = cvCreateMat(M_size, M_size, CV_64F);
			CvMat* invA = cvCreateMat(M_size, M_size, CV_64F);

			double *M_ptr = this->GetData(0);

			memcpy(A->data.ptr, M_ptr, sizeof(double)*this->M_size2);

			delete [] M_ptr;

			cvInv(A, invA);

			inv_mat->SetData((double*)invA->data.ptr, M_size2);


			cvReleaseMat(&A);
			cvReleaseMat(&invA);
		}
		else {
			CvMat* A = cvCreateMat(M_size, M_size, CV_64F);
			CvMat* invA = cvCreateMat(M_size, M_size, CV_64F);

			double *M_ptr = this->GetData(0);

			memcpy(A->data.ptr, M_ptr, sizeof(double)*this->M_size2);

			delete [] M_ptr;

			cvInv(A, invA);

			this->SetData((double*)invA->data.ptr, M_size2);

			cvReleaseMat(&A);
			cvReleaseMat(&invA);
		}
	}
}

void CCovarianceMatrix::SVD(CvMat *D, CvMat *U)
{
	if ( this->M_size == 1 ) {
		((double*)D->data.ptr)[0] = SM_ptr[0];
		((double*)U->data.ptr)[0] = 1.0;
	}
	else if ( this->M_size == 2 ) {
		eigen2_decomposition(this->SM_ptr, (double*)U->data.ptr, (double*)D->data.ptr);
	}
	else if ( this->M_size == 3 ) {
		eigen3_decomposition(this->SM_ptr, (double*)U->data.ptr, (double*)D->data.ptr);
	}
	else if ( this->M_size == 4 ) {
		eigen4_decomposition(this->SM_ptr, (double*)U->data.ptr, (double*)D->data.ptr);
	}
	else if ( this->M_size <= 12 ) {

		double *M_ptr = this->GetData(0);
		eigen_n_decomposition(M_ptr, (double*)U->data.ptr, (double*)D->data.ptr, this->M_size);
		delete [] M_ptr;
	}
	else {
		CvMat* A = cvCreateMat(M_size, M_size, CV_64F);

		double *M_ptr = this->GetData(0);

		memcpy(A->data.ptr, M_ptr, sizeof(double)*this->M_size2);

		delete [] M_ptr;

		double eps = 1.0E-15;

		cvEigenVV(A, U, D, eps);

		cvReleaseMat(&A);
	}
}

void CCovarianceMatrix::Square(CCovarianceMatrix *square_mat)
{
	if ( this->M_size == 1 ) {
		if ( square_mat ) {
			square_mat->SM_ptr[0] = SM_ptr[0]*SM_ptr[0];
		}
		else
			SM_ptr[0] *= SM_ptr[0];
	}
	else if ( this->M_size == 2 ) {
		double a0, a1, a2;

		double *tmp_SM_ptr = this->SM_ptr;

		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;

		if ( square_mat )
			tmp_SM_ptr = square_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		double t2 = a1*a1;

		*tmp_SM_ptr++ = a0*a0+t2;
		*tmp_SM_ptr++ = a0*a1+a1*a2;
		*tmp_SM_ptr++ = t2+a2*a2;
	}
	else if ( this->M_size == 3 ) {
		double a0, a1, a2, a3, a4, a5;

		double *tmp_SM_ptr = this->SM_ptr;

		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;

		if ( square_mat )
			tmp_SM_ptr = square_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		double t2 = a1*a1;
		double t3 = a2*a2;
		double t14 = a4*a4;

		*tmp_SM_ptr++ = a0*a0+t2+t3;
		*tmp_SM_ptr++ = a0*a1+a1*a3+a2*a4;
		*tmp_SM_ptr++ = a0*a2+a1*a4+a2*a5;
		*tmp_SM_ptr++ = t2+a3*a3+t14;
		*tmp_SM_ptr++ = a1*a2+a3*a4+a4*a5;
		*tmp_SM_ptr++ = t3+t14+a5*a5;
	}
	else if ( this->M_size == 4 ) {
		double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;
		a6 = *tmp_SM_ptr++;
		a7 = *tmp_SM_ptr++;
		a8 = *tmp_SM_ptr++;
		a9 = *tmp_SM_ptr++;

		if ( square_mat )
			tmp_SM_ptr = square_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		double t2 = a1*a1;
		double t3 = a2*a2;
		double t4 = a3*a3;
		double t22 = a5*a5;
		double t23 = a6*a6;
		double t36 = a8*a8;

		*tmp_SM_ptr++ = a0*a0+t2+t3+t4;
		*tmp_SM_ptr++ = a0*a1+a1*a4+a2*a5+a3*a6;
		*tmp_SM_ptr++ = a0*a2+a1*a5+a2*a7+a3*a8;
		*tmp_SM_ptr++ = a0*a3+a1*a6+a2*a8+a3*a9;
		*tmp_SM_ptr++ = t2+a4*a4+t22+t23;
		*tmp_SM_ptr++ = a1*a2+a4*a5+a5*a7+a6*a8;
		*tmp_SM_ptr++ = a1*a3+a4*a6+a5*a8+a6*a9;
		*tmp_SM_ptr++ = t3+t22+a7*a7+t36;
		*tmp_SM_ptr++ = a2*a3+a5*a6+a7*a8+a8*a9;
		*tmp_SM_ptr++ = t4+t23+t36+a9*a9;
	}
	else if ( this->M_size == 5 ) {
		double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;
		a6 = *tmp_SM_ptr++;
		a7 = *tmp_SM_ptr++;
		a8 = *tmp_SM_ptr++;
		a9 = *tmp_SM_ptr++;
		a10 = *tmp_SM_ptr++;
		a11 = *tmp_SM_ptr++;
		a12 = *tmp_SM_ptr++;
		a13 = *tmp_SM_ptr++;
		a14 = *tmp_SM_ptr++;

		if ( square_mat )
			tmp_SM_ptr = square_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		double t2 = a1*a1;
		double t3 = a2*a2;
		double t4 = a3*a3;
		double t5 = a4*a4;
		double t32 = a6*a6;
		double t33 = a7*a7;
		double t34 = a8*a8;
		double t55 = a10*a10;
		double t56 = a11*a11;
		double t71 = a13*a13;

		*tmp_SM_ptr++ = a0*a0+t2+t3+t4+t5;
		*tmp_SM_ptr++ = a0*a1+a1*a5+a2*a6+a3*a7+a4*a8;
		*tmp_SM_ptr++ = a0*a2+a1*a6+a2*a9+a3*a10+a4*a11;
		*tmp_SM_ptr++ = a0*a3+a1*a7+a2*a10+a3*a12+a4*a13;
		*tmp_SM_ptr++ = a0*a4+a1*a8+a2*a11+a3*a13+a4*a14;
		*tmp_SM_ptr++ = t2+a5*a5+t32+t33+t34;
		*tmp_SM_ptr++ = a1*a2+a5*a6+a6*a9+a7*a10+a8*a11;
		*tmp_SM_ptr++ = a1*a3+a5*a7+a6*a10+a7*a12+a8*a13;
		*tmp_SM_ptr++ = a1*a4+a5*a8+a6*a11+a7*a13+a8*a14;
		*tmp_SM_ptr++ = t3+t32+a9*a9+t55+t56;
		*tmp_SM_ptr++ = a2*a3+a6*a7+a9*a10+a10*a12+a11*a13;
		*tmp_SM_ptr++ = a2*a4+a6*a8+a9*a11+a10*a13+a11*a14;
		*tmp_SM_ptr++ = t4+t33+t55+a12*a12+t71;
		*tmp_SM_ptr++ = a3*a4+a7*a8+a10*a11+a12*a13+a13*a14;
		*tmp_SM_ptr++ = t5+t34+t56+t71+a14*a14;
	}
	else {
		double *M_ptr = this->GetData(0);
		CvMat* A = cvCreateMat(M_size, M_size, CV_64F);
		CvMat* AA = cvCreateMat(M_size, M_size, CV_64F);

		memcpy(A->data.ptr, M_ptr, sizeof(double)*M_size2);
		cvMatMul(A, A, AA);

		delete [] M_ptr;

		if ( square_mat )
			square_mat->SetData((double*)AA->data.ptr, M_size2);
		else
			this->SetData((double*)AA->data.ptr, M_size2);

		cvReleaseMat(&A);
		cvReleaseMat(&AA);
	}
}

void CCovarianceMatrix::Mul2(CCovarianceMatrix *Y_mat, CCovarianceMatrix *mul2_mat)
{
	if ( this->M_size == 1 ) {
		if ( mul2_mat )
			mul2_mat->SM_ptr[0] = SM_ptr[0] * Y_mat->SM_ptr[0];
		else
			SM_ptr[0] *= Y_mat->SM_ptr[0] * SM_ptr[0];
	}
	else if ( this->M_size == 2 ) {
		double a0, a1, a2, b0, b1, b2;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;

		tmp_SM_ptr = Y_mat->SM_ptr;
		b0 = *tmp_SM_ptr++;
		b1 = *tmp_SM_ptr++;
		b2 = *tmp_SM_ptr++;

		double t2 = a1*b1;
		double t3 = a0*b0+t2;
		double t7 = a0*b1+a1*b2;
		double t15 = a1*b0+a2*b1;
		double t18 = t2+a2*b2;

		if ( mul2_mat )
			tmp_SM_ptr = mul2_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		*tmp_SM_ptr++ = t3*a0+t7*a1;
		*tmp_SM_ptr++ = t3*a1+t7*a2;
		*tmp_SM_ptr++ = t15*a1+t18*a2;
	}
	else if ( this->M_size == 3 ) {
		double a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;

		tmp_SM_ptr = Y_mat->SM_ptr;
		b0 = *tmp_SM_ptr++;
		b1 = *tmp_SM_ptr++;
		b2 = *tmp_SM_ptr++;
		b3 = *tmp_SM_ptr++;
		b4 = *tmp_SM_ptr++;
		b5 = *tmp_SM_ptr++;

		double t2 = a1*b1;
		double t3 = a2*b2;
		double t4 = a0*b0+t2+t3;
		double t9 = a0*b1+a1*b3+a2*b4;
		double t14 = a0*b2+a1*b4+a2*b5;
		double t28 = a1*b0+a3*b1+a4*b2;
		double t31 = a4*b4;
		double t32 = t2+a3*b3+t31;
		double t37 = a1*b2+a3*b4+a4*b5;
		double t51 = a2*b0+a4*b1+a5*b2;
		double t56 = a2*b1+a4*b3+a5*b4;
		double t59 = t3+t31+a5*b5;

		if ( mul2_mat )
			tmp_SM_ptr = mul2_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		*tmp_SM_ptr++ = t4*a0+t9*a1+t14*a2;
		*tmp_SM_ptr++ = t4*a1+t9*a3+t14*a4;
		*tmp_SM_ptr++ = t4*a2+t9*a4+t14*a5;
		*tmp_SM_ptr++ = t28*a1+t32*a3+t37*a4;
		*tmp_SM_ptr++ = t28*a2+t32*a4+t37*a5;
		*tmp_SM_ptr++ = t51*a2+t56*a4+t59*a5;
	}
	else if ( this->M_size == 4 ) {
		double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;
		a6 = *tmp_SM_ptr++;
		a7 = *tmp_SM_ptr++;
		a8 = *tmp_SM_ptr++;
		a9 = *tmp_SM_ptr++;

		tmp_SM_ptr = Y_mat->SM_ptr;
		b0 = *tmp_SM_ptr++;
		b1 = *tmp_SM_ptr++;
		b2 = *tmp_SM_ptr++;
		b3 = *tmp_SM_ptr++;
		b4 = *tmp_SM_ptr++;
		b5 = *tmp_SM_ptr++;
		b6 = *tmp_SM_ptr++;
		b7 = *tmp_SM_ptr++;
		b8 = *tmp_SM_ptr++;
		b9 = *tmp_SM_ptr++;


		double t2 = a1*b1;
		double t3 = a2*b2;
		double t4 = a3*b3;
		double t5 = a0*b0+t2+t3+t4;
		double t11 = a0*b1+a1*b4+a2*b5+a3*b6;
		double t17 = a0*b2+a1*b5+a2*b7+a3*b8;
		double t23 = a0*b3+a1*b6+a2*b8+a3*b9;
		double t45 = a1*b0+a4*b1+a5*b2+a6*b3;
		double t48 = a5*b5;
		double t49 = a6*b6;
		double t50 = t2+a4*b4+t48+t49;
		double t56 = a1*b2+a4*b5+a5*b7+a6*b8;
		double t62 = a1*b3+a4*b6+a5*b8+a6*b9;
		double t84 = a2*b0+a5*b1+a7*b2+a8*b3;
		double t90 = a2*b1+a5*b4+a7*b5+a8*b6;
		double t93 = a8*b8;
		double t94 = t3+t48+a7*b7+t93;
		double t100 = a2*b3+a5*b6+a7*b8+a8*b9;
		double t122 = a3*b0+a6*b1+a8*b2+a9*b3;
		double t128 = a3*b1+a6*b4+a8*b5+a9*b6;
		double t134 = a3*b2+a6*b5+a8*b7+a9*b8;
		double t137 = t4+t49+t93+a9*b9;

		if ( mul2_mat )
			tmp_SM_ptr = mul2_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		*tmp_SM_ptr++ = t5*a0+t11*a1+t17*a2+t23*a3;
		*tmp_SM_ptr++ = t5*a1+t11*a4+t17*a5+t23*a6;
		*tmp_SM_ptr++ = t5*a2+t11*a5+t17*a7+t23*a8;
		*tmp_SM_ptr++ = t5*a3+t11*a6+t17*a8+t23*a9;
		*tmp_SM_ptr++ = t45*a1+t50*a4+t56*a5+t62*a6;
		*tmp_SM_ptr++ = t45*a2+t50*a5+t56*a7+t62*a8;
		*tmp_SM_ptr++ = t45*a3+t50*a6+t56*a8+t62*a9;
		*tmp_SM_ptr++ = t84*a2+t90*a5+t94*a7+t100*a8;
		*tmp_SM_ptr++ = t84*a3+t90*a6+t94*a8+t100*a9;
		*tmp_SM_ptr++ = t122*a3+t128*a6+t134*a8+t137*a9;
	}
	else if ( this->M_size == 5 ) {
		double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;
		a6 = *tmp_SM_ptr++;
		a7 = *tmp_SM_ptr++;
		a8 = *tmp_SM_ptr++;
		a9 = *tmp_SM_ptr++;
		a10 = *tmp_SM_ptr++;
		a11 = *tmp_SM_ptr++;
		a12 = *tmp_SM_ptr++;
		a13 = *tmp_SM_ptr++;
		a14 = *tmp_SM_ptr++;

		tmp_SM_ptr = Y_mat->SM_ptr;
		b0 = *tmp_SM_ptr++;
		b1 = *tmp_SM_ptr++;
		b2 = *tmp_SM_ptr++;
		b3 = *tmp_SM_ptr++;
		b4 = *tmp_SM_ptr++;
		b5 = *tmp_SM_ptr++;
		b6 = *tmp_SM_ptr++;
		b7 = *tmp_SM_ptr++;
		b8 = *tmp_SM_ptr++;
		b9 = *tmp_SM_ptr++;
		b10 = *tmp_SM_ptr++;
		b11 = *tmp_SM_ptr++;
		b12 = *tmp_SM_ptr++;
		b13 = *tmp_SM_ptr++;
		b14 = *tmp_SM_ptr++;

		double t2 = a1*b1;
		double t3 = a2*b2;
		double t4 = a3*b3;
		double t5 = a4*b4;
		double t6 = a0*b0+t2+t3+t4+t5;
		double t13 = a0*b1+a1*b5+a2*b6+a3*b7+a4*b8;
		double t20 = a0*b2+a1*b6+a2*b9+a3*b10+a4*b11;
		double t27 = a0*b3+a1*b7+a2*b10+a3*b12+a4*b13;
		double t34 = a0*b4+a1*b8+a2*b11+a3*b13+a4*b14;
		double t66 = a1*b0+a5*b1+a6*b2+a7*b3+a8*b4;
		double t69 = a6*b6;
		double t70 = a7*b7;
		double t71 = a8*b8;
		double t72 = t2+a5*b5+t69+t70+t71;
		double t79 = a1*b2+a5*b6+a6*b9+a7*b10+a8*b11;
		double t86 = a1*b3+a5*b7+a6*b10+a7*b12+a8*b13;
		double t93 = a1*b4+a5*b8+a6*b11+a7*b13+a8*b14;
		double t125 = a2*b0+a6*b1+a9*b2+a10*b3+a11*b4;
		double t132 = a2*b1+a6*b5+a9*b6+a10*b7+a11*b8;
		double t135 = a10*b10;
		double t136 = a11*b11;
		double t137 = t3+t69+a9*b9+t135+t136;
		double t144 = a2*b3+a6*b7+a9*b10+a10*b12+a11*b13;
		double t151 = a2*b4+a6*b8+a9*b11+a10*b13+a11*b14;
		double t183 = a3*b0+a7*b1+a10*b2+a12*b3+a13*b4;
		double t190 = a3*b1+a7*b5+a10*b6+a12*b7+a13*b8;
		double t197 = a3*b2+a7*b6+a10*b9+a12*b10+a13*b11;
		double t200 = a13*b13;
		double t201 = t4+t70+t135+a12*b12+t200;
		double t208 = a3*b4+a7*b8+a10*b11+a12*b13+a13*b14;
		double t240 = a4*b0+a8*b1+a11*b2+a13*b3+a14*b4;
		double t247 = a4*b1+a8*b5+a11*b6+a13*b7+a14*b8;
		double t254 = a4*b2+a8*b6+a11*b9+a13*b10+a14*b11;
		double t261 = a4*b3+a8*b7+a11*b10+a13*b12+a14*b13;
		double t264 = t5+t71+t136+t200+a14*b14;

		if ( mul2_mat )
			tmp_SM_ptr = mul2_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		*tmp_SM_ptr++ = t6*a0+t13*a1+t20*a2+t27*a3+t34*a4;
		*tmp_SM_ptr++ = t6*a1+t13*a5+t20*a6+t27*a7+t34*a8;
		*tmp_SM_ptr++ = t6*a2+t13*a6+t20*a9+t27*a10+t34*a11;
		*tmp_SM_ptr++ = t6*a3+t13*a7+t20*a10+t27*a12+t34*a13;
		*tmp_SM_ptr++ = t6*a4+t13*a8+t20*a11+t27*a13+t34*a14;
		*tmp_SM_ptr++ = t66*a1+t72*a5+t79*a6+t86*a7+t93*a8;
		*tmp_SM_ptr++ = t66*a2+t72*a6+t79*a9+t86*a10+t93*a11;
		*tmp_SM_ptr++ = t66*a3+t72*a7+t79*a10+t86*a12+t93*a13;
		*tmp_SM_ptr++ = t66*a4+t72*a8+t79*a11+t86*a13+t93*a14;
		*tmp_SM_ptr++ = t125*a2+t132*a6+t137*a9+t144*a10+t151*a11;
		*tmp_SM_ptr++ = t125*a3+t132*a7+t137*a10+t144*a12+t151*a13;
		*tmp_SM_ptr++ = t125*a4+t132*a8+t137*a11+t144*a13+t151*a14;
		*tmp_SM_ptr++ = t183*a3+t190*a7+t197*a10+t201*a12+t208*a13;
		*tmp_SM_ptr++ = t183*a4+t190*a8+t197*a11+t201*a13+t208*a14;
		*tmp_SM_ptr++ = t240*a4+t247*a8+t254*a11+t261*a13+t264*a14;
	}
	else if ( this->M_size == 6 ) {
		double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20;
		double b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;
		a6 = *tmp_SM_ptr++;
		a7 = *tmp_SM_ptr++;
		a8 = *tmp_SM_ptr++;
		a9 = *tmp_SM_ptr++;
		a10 = *tmp_SM_ptr++;
		a11 = *tmp_SM_ptr++;
		a12 = *tmp_SM_ptr++;
		a13 = *tmp_SM_ptr++;
		a14 = *tmp_SM_ptr++;
		a15 = *tmp_SM_ptr++;
		a16 = *tmp_SM_ptr++;
		a17 = *tmp_SM_ptr++;
		a18 = *tmp_SM_ptr++;
		a19 = *tmp_SM_ptr++;
		a20 = *tmp_SM_ptr++;

		tmp_SM_ptr = Y_mat->SM_ptr;
		b0 = *tmp_SM_ptr++;
		b1 = *tmp_SM_ptr++;
		b2 = *tmp_SM_ptr++;
		b3 = *tmp_SM_ptr++;
		b4 = *tmp_SM_ptr++;
		b5 = *tmp_SM_ptr++;
		b6 = *tmp_SM_ptr++;
		b7 = *tmp_SM_ptr++;
		b8 = *tmp_SM_ptr++;
		b9 = *tmp_SM_ptr++;
		b10 = *tmp_SM_ptr++;
		b11 = *tmp_SM_ptr++;
		b12 = *tmp_SM_ptr++;
		b13 = *tmp_SM_ptr++;
		b14 = *tmp_SM_ptr++;
		b15 = *tmp_SM_ptr++;
		b16 = *tmp_SM_ptr++;
		b17 = *tmp_SM_ptr++;
		b18 = *tmp_SM_ptr++;
		b19 = *tmp_SM_ptr++;
		b20 = *tmp_SM_ptr++;

		double t2 = a1*b1;
		double t3 = a2*b2;
		double t4 = a3*b3;
		double t5 = a4*b4;
		double t6 = a5*b5;
		double t7 = a0*b0+t2+t3+t4+t5+t6;
		double t15 = a0*b1+a1*b6+a2*b7+a3*b8+a4*b9+a5*b10;
		double t23 = a0*b2+a1*b7+a2*b11+a3*b12+a4*b13+a5*b14;
		double t31 = a0*b3+a1*b8+a2*b12+a3*b15+a4*b16+a5*b17;
		double t39 = a0*b4+a1*b9+a2*b13+a3*b16+a4*b18+a5*b19;
		double t47 = a0*b5+a1*b10+a2*b14+a3*b17+a4*b19+a5*b20;
		double t91 = a1*b0+a6*b1+a7*b2+a8*b3+a9*b4+a10*b5;
		double t94 = a7*b7;
		double t95 = a8*b8;
		double t96 = a9*b9;
		double t97 = a10*b10;
		double t98 = t2+a6*b6+t94+t95+t96+t97;
		double t106 = a1*b2+a6*b7+a7*b11+a8*b12+a9*b13+a10*b14;
		double t114 = a1*b3+a6*b8+a7*b12+a8*b15+a9*b16+a10*b17;
		double t122 = a1*b4+a6*b9+a7*b13+a8*b16+a9*b18+a10*b19;
		double t130 = a1*b5+a6*b10+a7*b14+a8*b17+a9*b19+a10*b20;
		double t174 = a2*b0+a7*b1+a11*b2+a12*b3+a13*b4+a14*b5;
		double t182 = a2*b1+a7*b6+a11*b7+a12*b8+a13*b9+a14*b10;
		double t185 = a12*b12;
		double t186 = a13*b13;
		double t187 = a14*b14;
		double t188 = t3+t94+a11*b11+t185+t186+t187;
		double t196 = a2*b3+a7*b8+a11*b12+a12*b15+a13*b16+a14*b17;
		double t204 = a2*b4+a7*b9+a11*b13+a12*b16+a13*b18+a14*b19;
		double t212 = a2*b5+a7*b10+a11*b14+a12*b17+a13*b19+a14*b20;
		double t256 = a3*b0+a8*b1+a12*b2+a15*b3+a16*b4+a17*b5;
		double t264 = a3*b1+a8*b6+a12*b7+a15*b8+a16*b9+a17*b10;
		double t272 = a3*b2+a8*b7+a12*b11+a15*b12+a16*b13+a17*b14;
		double t275 = a16*b16;
		double t276 = a17*b17;
		double t277 = t4+t95+t185+a15*b15+t275+t276;
		double t285 = a3*b4+a8*b9+a12*b13+a15*b16+a16*b18+a17*b19;
		double t293 = a3*b5+a8*b10+a12*b14+a15*b17+a16*b19+a17*b20;
		double t337 = a4*b0+a9*b1+a13*b2+a16*b3+a18*b4+a19*b5;
		double t345 = a4*b1+a9*b6+a13*b7+a16*b8+a18*b9+a19*b10;
		double t353 = a4*b2+a9*b7+a13*b11+a16*b12+a18*b13+a19*b14;
		double t361 = a4*b3+a9*b8+a13*b12+a16*b15+a18*b16+a19*b17;
		double t364 = a19*b19;
		double t365 = t5+t96+t186+t275+a18*b18+t364;
		double t373 = a4*b5+a9*b10+a13*b14+a16*b17+a18*b19+a19*b20;
		double t417 = a5*b0+a10*b1+a14*b2+a17*b3+a19*b4+a20*b5;
		double t425 = a5*b1+a10*b6+a14*b7+a17*b8+a19*b9+a20*b10;
		double t433 = a5*b2+a10*b7+a14*b11+a17*b12+a19*b13+a20*b14;
		double t441 = a5*b3+a10*b8+a14*b12+a17*b15+a19*b16+a20*b17;
		double t449 = a5*b4+a10*b9+a14*b13+a17*b16+a19*b18+a20*b19;
		double t452 = t6+t97+t187+t276+t364+a20*b20;

		if ( mul2_mat )
			tmp_SM_ptr = mul2_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		*tmp_SM_ptr++ = t7*a0+t15*a1+t23*a2+t31*a3+t39*a4+t47*a5;
		*tmp_SM_ptr++ = t7*a1+t15*a6+t23*a7+t31*a8+t39*a9+t47*a10;
		*tmp_SM_ptr++ = t7*a2+t15*a7+t23*a11+t31*a12+t39*a13+t47*a14;
		*tmp_SM_ptr++ = t7*a3+t15*a8+t23*a12+t31*a15+t39*a16+t47*a17;
		*tmp_SM_ptr++ = t7*a4+t15*a9+t23*a13+t31*a16+t39*a18+t47*a19;
		*tmp_SM_ptr++ = t7*a5+t15*a10+t23*a14+t31*a17+t39*a19+t47*a20;
		*tmp_SM_ptr++ = t91*a1+t98*a6+t106*a7+t114*a8+t122*a9+t130*a10;
		*tmp_SM_ptr++ = t91*a2+t98*a7+t106*a11+t114*a12+t122*a13+t130*a14;
		*tmp_SM_ptr++ = t91*a3+t98*a8+t106*a12+t114*a15+t122*a16+t130*a17;
		*tmp_SM_ptr++ = t91*a4+t98*a9+t106*a13+t114*a16+t122*a18+t130*a19;
		*tmp_SM_ptr++ = t91*a5+t98*a10+t106*a14+t114*a17+t122*a19+t130*a20;
		*tmp_SM_ptr++ = t174*a2+t182*a7+t188*a11+t196*a12+t204*a13+t212*a14;
		*tmp_SM_ptr++ = t174*a3+t182*a8+t188*a12+t196*a15+t204*a16+t212*a17;
		*tmp_SM_ptr++ = t174*a4+t182*a9+t188*a13+t196*a16+t204*a18+t212*a19;
		*tmp_SM_ptr++ = t174*a5+t182*a10+t188*a14+t196*a17+t204*a19+t212*a20;
		*tmp_SM_ptr++ = t256*a3+t264*a8+t272*a12+t277*a15+t285*a16+t293*a17;
		*tmp_SM_ptr++ = t256*a4+t264*a9+t272*a13+t277*a16+t285*a18+t293*a19;
		*tmp_SM_ptr++ = t256*a5+t264*a10+t272*a14+t277*a17+t285*a19+t293*a20;
		*tmp_SM_ptr++ = t337*a4+t345*a9+t353*a13+t361*a16+t365*a18+t373*a19;
		*tmp_SM_ptr++ = t337*a5+t345*a10+t353*a14+t361*a17+t365*a19+t373*a20;
		*tmp_SM_ptr++ = t417*a5+t425*a10+t433*a14+t441*a17+t449*a19+t452*a20;
	}
	else if ( this->M_size == 7 ) {
		double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27;
		double b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26, b27;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;
		a6 = *tmp_SM_ptr++;
		a7 = *tmp_SM_ptr++;
		a8 = *tmp_SM_ptr++;
		a9 = *tmp_SM_ptr++;
		a10 = *tmp_SM_ptr++;
		a11 = *tmp_SM_ptr++;
		a12 = *tmp_SM_ptr++;
		a13 = *tmp_SM_ptr++;
		a14 = *tmp_SM_ptr++;
		a15 = *tmp_SM_ptr++;
		a16 = *tmp_SM_ptr++;
		a17 = *tmp_SM_ptr++;
		a18 = *tmp_SM_ptr++;
		a19 = *tmp_SM_ptr++;
		a20 = *tmp_SM_ptr++;
		a21 = *tmp_SM_ptr++;
		a22 = *tmp_SM_ptr++;
		a23 = *tmp_SM_ptr++;
		a24 = *tmp_SM_ptr++;
		a25 = *tmp_SM_ptr++;
		a26 = *tmp_SM_ptr++;
		a27 = *tmp_SM_ptr++;

		tmp_SM_ptr = Y_mat->SM_ptr;
		b0 = *tmp_SM_ptr++;
		b1 = *tmp_SM_ptr++;
		b2 = *tmp_SM_ptr++;
		b3 = *tmp_SM_ptr++;
		b4 = *tmp_SM_ptr++;
		b5 = *tmp_SM_ptr++;
		b6 = *tmp_SM_ptr++;
		b7 = *tmp_SM_ptr++;
		b8 = *tmp_SM_ptr++;
		b9 = *tmp_SM_ptr++;
		b10 = *tmp_SM_ptr++;
		b11 = *tmp_SM_ptr++;
		b12 = *tmp_SM_ptr++;
		b13 = *tmp_SM_ptr++;
		b14 = *tmp_SM_ptr++;
		b15 = *tmp_SM_ptr++;
		b16 = *tmp_SM_ptr++;
		b17 = *tmp_SM_ptr++;
		b18 = *tmp_SM_ptr++;
		b19 = *tmp_SM_ptr++;
		b20 = *tmp_SM_ptr++;
		b21 = *tmp_SM_ptr++;
		b22 = *tmp_SM_ptr++;
		b23 = *tmp_SM_ptr++;
		b24 = *tmp_SM_ptr++;
		b25 = *tmp_SM_ptr++;
		b26 = *tmp_SM_ptr++;
		b27 = *tmp_SM_ptr++;

		double t2 = a1*b1;
		double t3 = a2*b2;
		double t4 = a3*b3;
		double t5 = a4*b4;
		double t6 = a5*b5;
		double t7 = a6*b6;
		double t8 = a0*b0+t2+t3+t4+t5+t6+t7;
		double t17 = a0*b1+a1*b7+a2*b8+a3*b9+a4*b10+a5*b11+a6*b12;
		double t20 = a1*b8;
		double t26 = a0*b2+t20+a2*b13+a3*b14+a4*b15+a5*b16+a6*b17;
		double t34 = a0*b3+t20+a2*b14+a3*b18+a4*b19+a5*b20+a6*b21;
		double t43 = a0*b4+a1*b10+a2*b15+a3*b19+a4*b22+a5*b23+a6*b24;
		double t52 = a0*b5+a1*b11+a2*b16+a3*b20+a4*b23+a5*b25+a6*b26;
		double t61 = a0*b6+a1*b12+a2*b17+a3*b21+a4*b24+a5*b26+a6*b27;
		double t73 = t17*a8;
		double t118 = a1*b0+a7*b1+a8*b2+a8*b3+a10*b4+a11*b5+a12*b6;
		double t121 = a8*b8;
		double t123 = a10*b10;
		double t124 = a11*b11;
		double t125 = a12*b12;
		double t126 = t2+a7*b7+t121+a8*b9+t123+t124+t125;
		double t129 = a7*b8;
		double t131 = a8*b14;
		double t135 = a1*b2+t129+a8*b13+t131+a10*b15+a11*b16+a12*b17;
		double t142 = a1*b3+t129+t131+a8*b18+a10*b19+a11*b20+a12*b21;
		double t151 = a1*b4+a7*b10+a8*b15+a8*b19+a10*b22+a11*b23+a12*b24;
		double t160 = a1*b5+a7*b11+a8*b16+a8*b20+a10*b23+a11*b25+a12*b26;
		double t169 = a1*b6+a7*b12+a8*b17+a8*b21+a10*b24+a11*b26+a12*b27;
		double t181 = t126*a8;
		double t226 = a2*b0+a8*b1+a13*b2+a14*b3+a15*b4+a16*b5+a17*b6;
		double t235 = a2*b1+a8*b7+a13*b8+a14*b9+a15*b10+a16*b11+a17*b12;
		double t238 = a14*b14;
		double t239 = a15*b15;
		double t240 = a16*b16;
		double t241 = a17*b17;
		double t242 = t3+t121+a13*b13+t238+t239+t240+t241;
		double t250 = a2*b3+t121+a13*b14+a14*b18+a15*b19+a16*b20+a17*b21;
		double t259 = a2*b4+a8*b10+a13*b15+a14*b19+a15*b22+a16*b23+a17*b24;
		double t268 = a2*b5+a8*b11+a13*b16+a14*b20+a15*b23+a16*b25+a17*b26;
		double t277 = a2*b6+a8*b12+a13*b17+a14*b21+a15*b24+a16*b26+a17*b27;
		double t289 = t235*a8;
		double t334 = a3*b0+a9*b1+a14*b2+a18*b3+a19*b4+a20*b5+a21*b6;
		double t343 = a3*b1+a9*b7+a14*b8+a18*b9+a19*b10+a20*b11+a21*b12;
		double t346 = a9*b8;
		double t352 = a3*b2+t346+a14*b13+a18*b14+a19*b15+a20*b16+a21*b17;
		double t355 = a19*b19;
		double t356 = a20*b20;
		double t357 = a21*b21;
		double t358 = t4+t346+t238+a18*b18+t355+t356+t357;
		double t367 = a3*b4+a9*b10+a14*b15+a18*b19+a19*b22+a20*b23+a21*b24;
		double t376 = a3*b5+a9*b11+a14*b16+a18*b20+a19*b23+a20*b25+a21*b26;
		double t385 = a3*b6+a9*b12+a14*b17+a18*b21+a19*b24+a20*b26+a21*b27;
		double t397 = t343*a8;
		double t442 = a4*b0+a10*b1+a15*b2+a19*b3+a22*b4+a23*b5+a24*b6;
		double t451 = a4*b1+a10*b7+a15*b8+a19*b9+a22*b10+a23*b11+a24*b12;
		double t454 = a10*b8;
		double t460 = a4*b2+t454+a15*b13+a19*b14+a22*b15+a23*b16+a24*b17;
		double t468 = a4*b3+t454+a15*b14+a19*b18+a22*b19+a23*b20+a24*b21;
		double t471 = a23*b23;
		double t472 = a24*b24;
		double t473 = t5+t123+t239+t355+a22*b22+t471+t472;
		double t482 = a4*b5+a10*b11+a15*b16+a19*b20+a22*b23+a23*b25+a24*b26;
		double t491 = a4*b6+a10*b12+a15*b17+a19*b21+a22*b24+a23*b26+a24*b27;
		double t548 = a5*b0+a11*b1+a16*b2+a20*b3+a23*b4+a25*b5+a26*b6;
		double t557 = a5*b1+a11*b7+a16*b8+a20*b9+a23*b10+a25*b11+a26*b12;
		double t560 = a11*b8;
		double t566 = a5*b2+t560+a16*b13+a20*b14+a23*b15+a25*b16+a26*b17;
		double t574 = a5*b3+t560+a16*b14+a20*b18+a23*b19+a25*b20+a26*b21;
		double t583 = a5*b4+a11*b10+a16*b15+a20*b19+a23*b22+a25*b23+a26*b24;
		double t586 = a26*b26;
		double t587 = t6+t124+t240+t356+t471+a25*b25+t586;
		double t596 = a5*b6+a11*b12+a16*b17+a20*b21+a23*b24+a25*b26+a26*b27;
		double t653 = a6*b0+a12*b1+a17*b2+a21*b3+a24*b4+a26*b5+a27*b6;
		double t662 = a6*b1+a12*b7+a17*b8+a21*b9+a24*b10+a26*b11+a27*b12;
		double t665 = a12*b8;
		double t671 = a6*b2+t665+a17*b13+a21*b14+a24*b15+a26*b16+a27*b17;
		double t679 = a6*b3+t665+a17*b14+a21*b18+a24*b19+a26*b20+a27*b21;
		double t688 = a6*b4+a12*b10+a17*b15+a21*b19+a24*b22+a26*b23+a27*b24;
		double t697 = a6*b5+a12*b11+a17*b16+a21*b20+a24*b23+a26*b25+a27*b26;
		double t700 = t7+t125+t241+t357+t472+t586+a27*b27;

		if ( mul2_mat )
			tmp_SM_ptr = mul2_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		*tmp_SM_ptr++ = t8*a0+t17*a1+t26*a2+t34*a3+t43*a4+t52*a5+t61*a6;
		*tmp_SM_ptr++ = t8*a1+t17*a7+t26*a8+t34*a9+t43*a10+t52*a11+t61*a12;
		*tmp_SM_ptr++ = t8*a2+t73+t26*a13+t34*a14+t43*a15+t52*a16+t61*a17;
		*tmp_SM_ptr++ = t8*a3+t73+t26*a14+t34*a18+t43*a19+t52*a20+t61*a21;
		*tmp_SM_ptr++ = t8*a4+t17*a10+t26*a15+t34*a19+t43*a22+t52*a23+t61*a24;
		*tmp_SM_ptr++ = t8*a5+t17*a11+t26*a16+t34*a20+t43*a23+t52*a25+t61*a26;
		*tmp_SM_ptr++ = t8*a6+t17*a12+t26*a17+t34*a21+t43*a24+t52*a26+t61*a27;
		*tmp_SM_ptr++ = t118*a1+t126*a7+t135*a8+t142*a9+t151*a10+t160*a11+t169*a12;
		*tmp_SM_ptr++ = t118*a2+t181+t135*a13+t142*a14+t151*a15+t160*a16+t169*a17;
		*tmp_SM_ptr++ = t118*a3+t181+t135*a14+t142*a18+t151*a19+t160*a20+t169*a21;
		*tmp_SM_ptr++ = t118*a4+t126*a10+t135*a15+t142*a19+t151*a22+t160*a23+t169*a24;
		*tmp_SM_ptr++ = t118*a5+t126*a11+t135*a16+t142*a20+t151*a23+t160*a25+t169*a26;
		*tmp_SM_ptr++ = t118*a6+t126*a12+t135*a17+t142*a21+t151*a24+t160*a26+t169*a27;
		*tmp_SM_ptr++ = t226*a2+t289+t242*a13+t250*a14+t259*a15+t268*a16+t277*a17;
		*tmp_SM_ptr++ = t226*a3+t289+t242*a14+t250*a18+t259*a19+t268*a20+t277*a21;
		*tmp_SM_ptr++ = t226*a4+t235*a10+t242*a15+t250*a19+t259*a22+t268*a23+t277*a24;
		*tmp_SM_ptr++ = t226*a5+t235*a11+t242*a16+t250*a20+t259*a23+t268*a25+t277*a26;
		*tmp_SM_ptr++ = t226*a6+t235*a12+t242*a17+t250*a21+t259*a24+t268*a26+t277*a27;
		*tmp_SM_ptr++ = t334*a3+t397+t352*a14+t358*a18+t367*a19+t376*a20+t385*a21;
		*tmp_SM_ptr++ = t334*a4+t343*a10+t352*a15+t358*a19+t367*a22+t376*a23+t385*a24;
		*tmp_SM_ptr++ = t334*a5+t343*a11+t352*a16+t358*a20+t367*a23+t376*a25+t385*a26;
		*tmp_SM_ptr++ = t334*a6+t343*a12+t352*a17+t358*a21+t367*a24+t376*a26+t385*a27;
		*tmp_SM_ptr++ = t442*a4+t451*a10+t460*a15+t468*a19+t473*a22+t482*a23+t491*a24;
		*tmp_SM_ptr++ = t442*a5+t451*a11+t460*a16+t468*a20+t473*a23+t482*a25+t491*a26;
		*tmp_SM_ptr++ = t442*a6+t451*a12+t460*a17+t468*a21+t473*a24+t482*a26+t491*a27;
		*tmp_SM_ptr++ = t548*a5+t557*a11+t566*a16+t574*a20+t583*a23+t587*a25+t596*a26;
		*tmp_SM_ptr++ = t548*a6+t557*a12+t566*a17+t574*a21+t583*a24+t587*a26+t596*a27;
		*tmp_SM_ptr++ = t653*a6+t662*a12+t671*a17+t679*a21+t688*a24+t697*a26+t700*a27;
	}
	else if ( this->M_size == 8 ) {
		double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35;
		double b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34, b35;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;
		a6 = *tmp_SM_ptr++;
		a7 = *tmp_SM_ptr++;
		a8 = *tmp_SM_ptr++;
		a9 = *tmp_SM_ptr++;
		a10 = *tmp_SM_ptr++;
		a11 = *tmp_SM_ptr++;
		a12 = *tmp_SM_ptr++;
		a13 = *tmp_SM_ptr++;
		a14 = *tmp_SM_ptr++;
		a15 = *tmp_SM_ptr++;
		a16 = *tmp_SM_ptr++;
		a17 = *tmp_SM_ptr++;
		a18 = *tmp_SM_ptr++;
		a19 = *tmp_SM_ptr++;
		a20 = *tmp_SM_ptr++;
		a21 = *tmp_SM_ptr++;
		a22 = *tmp_SM_ptr++;
		a23 = *tmp_SM_ptr++;
		a24 = *tmp_SM_ptr++;
		a25 = *tmp_SM_ptr++;
		a26 = *tmp_SM_ptr++;
		a27 = *tmp_SM_ptr++;
		a28 = *tmp_SM_ptr++;
		a29 = *tmp_SM_ptr++;
		a30 = *tmp_SM_ptr++;
		a31 = *tmp_SM_ptr++;
		a32 = *tmp_SM_ptr++;
		a33 = *tmp_SM_ptr++;
		a34 = *tmp_SM_ptr++;
		a35 = *tmp_SM_ptr++;

		tmp_SM_ptr = Y_mat->SM_ptr;
		b0 = *tmp_SM_ptr++;
		b1 = *tmp_SM_ptr++;
		b2 = *tmp_SM_ptr++;
		b3 = *tmp_SM_ptr++;
		b4 = *tmp_SM_ptr++;
		b5 = *tmp_SM_ptr++;
		b6 = *tmp_SM_ptr++;
		b7 = *tmp_SM_ptr++;
		b8 = *tmp_SM_ptr++;
		b9 = *tmp_SM_ptr++;
		b10 = *tmp_SM_ptr++;
		b11 = *tmp_SM_ptr++;
		b12 = *tmp_SM_ptr++;
		b13 = *tmp_SM_ptr++;
		b14 = *tmp_SM_ptr++;
		b15 = *tmp_SM_ptr++;
		b16 = *tmp_SM_ptr++;
		b17 = *tmp_SM_ptr++;
		b18 = *tmp_SM_ptr++;
		b19 = *tmp_SM_ptr++;
		b20 = *tmp_SM_ptr++;
		b21 = *tmp_SM_ptr++;
		b22 = *tmp_SM_ptr++;
		b23 = *tmp_SM_ptr++;
		b24 = *tmp_SM_ptr++;
		b25 = *tmp_SM_ptr++;
		b26 = *tmp_SM_ptr++;
		b27 = *tmp_SM_ptr++;
		b28 = *tmp_SM_ptr++;
		b29 = *tmp_SM_ptr++;
		b30 = *tmp_SM_ptr++;
		b31 = *tmp_SM_ptr++;
		b32 = *tmp_SM_ptr++;
		b33 = *tmp_SM_ptr++;
		b34 = *tmp_SM_ptr++;
		b35 = *tmp_SM_ptr++;

		double t2 = a1*b1;
		double t3 = a2*b2;
		double t4 = a3*b3;
		double t5 = a4*b4;
		double t6 = a5*b5;
		double t7 = a6*b6;
		double t8 = a7*b7;
		double t9 = a0*b0+t2+t3+t4+t5+t6+t7+t8;
		double t19 = a0*b1+a1*b8+a2*b9+a3*b10+a4*b11+a5*b12+a6*b13+a7*b14;
		double t29 = a0*b2+a1*b9+a2*b15+a3*b16+a4*b17+a5*b18+a6*b19+a7*b20;
		double t39 = a0*b3+a1*b10+a2*b16+a3*b21+a4*b22+a5*b23+a6*b24+a7*b25;
		double t49 = a0*b4+a1*b11+a2*b17+a3*b22+a4*b26+a5*b27+a6*b28+a7*b29;
		double t59 = a0*b5+a1*b12+a2*b18+a3*b23+a4*b27+a5*b30+a6*b31+a7*b32;
		double t69 = a0*b6+a1*b13+a2*b19+a3*b24+a4*b28+a5*b31+a6*b33+a7*b34;
		double t79 = a0*b7+a1*b14+a2*b20+a3*b25+a4*b29+a5*b32+a6*b34+a7*b35;
		double t153 = a1*b0+a8*b1+a9*b2+a10*b3+a11*b4+a12*b5+a13*b6+a14*b7;
		double t156 = a9*b9;
		double t157 = a10*b10;
		double t158 = a11*b11;
		double t159 = a12*b12;
		double t160 = a13*b13;
		double t161 = a14*b14;
		double t162 = t2+a8*b8+t156+t157+t158+t159+t160+t161;
		double t172 = a1*b2+a8*b9+a9*b15+a10*b16+a11*b17+a12*b18+a13*b19+a14*b20;
		double t182 = a1*b3+a8*b10+a9*b16+a10*b21+a11*b22+a12*b23+a13*b24+a14*b25;
		double t192 = a1*b4+a8*b11+a9*b17+a10*b22+a11*b26+a12*b27+a13*b28+a14*b29;
		double t202 = a1*b5+a8*b12+a9*b18+a10*b23+a11*b27+a12*b30+a13*b31+a14*b32;
		double t212 = a1*b6+a8*b13+a9*b19+a10*b24+a11*b28+a12*b31+a13*b33+a14*b34;
		double t222 = a1*b7+a8*b14+a9*b20+a10*b25+a11*b29+a12*b32+a13*b34+a14*b35;
		double t296 = a2*b0+a9*b1+a15*b2+a16*b3+a17*b4+a18*b5+a19*b6+a20*b7;
		double t306 = a2*b1+a9*b8+a15*b9+a16*b10+a17*b11+a18*b12+a19*b13+a20*b14;
		double t309 = a16*b16;
		double t310 = a17*b17;
		double t311 = a18*b18;
		double t312 = a19*b19;
		double t313 = a20*b20;
		double t314 = t3+t156+a15*b15+t309+t310+t311+t312+t313;
		double t324 = a2*b3+a9*b10+a15*b16+a16*b21+a17*b22+a18*b23+a19*b24+a20*b25;
		double t334 = a2*b4+a9*b11+a15*b17+a16*b22+a17*b26+a18*b27+a19*b28+a20*b29;
		double t344 = a2*b5+a9*b12+a15*b18+a16*b23+a17*b27+a18*b30+a19*b31+a20*b32;
		double t354 = a2*b6+a9*b13+a15*b19+a16*b24+a17*b28+a18*b31+a19*b33+a20*b34;
		double t364 = a2*b7+a9*b14+a15*b20+a16*b25+a17*b29+a18*b32+a19*b34+a20*b35;
		double t438 = a3*b0+a10*b1+a16*b2+a21*b3+a22*b4+a23*b5+a24*b6+a25*b7;
		double t448 = a3*b1+a10*b8+a16*b9+a21*b10+a22*b11+a23*b12+a24*b13+a25*b14;
		double t458 = a3*b2+a10*b9+a16*b15+a21*b16+a22*b17+a23*b18+a24*b19+a25*b20;
		double t461 = a22*b22;
		double t462 = a23*b23;
		double t463 = a24*b24;
		double t464 = a25*b25;
		double t465 = t4+t157+t309+a21*b21+t461+t462+t463+t464;
		double t475 = a3*b4+a10*b11+a16*b17+a21*b22+a22*b26+a23*b27+a24*b28+a25*b29;
		double t485 = a3*b5+a10*b12+a16*b18+a21*b23+a22*b27+a23*b30+a24*b31+a25*b32;
		double t495 = a3*b6+a10*b13+a16*b19+a21*b24+a22*b28+a23*b31+a24*b33+a25*b34;
		double t505 = a3*b7+a10*b14+a16*b20+a21*b25+a22*b29+a23*b32+a24*b34+a25*b35;
		double t579 = a4*b0+a11*b1+a17*b2+a22*b3+a26*b4+a27*b5+a28*b6+a29*b7;
		double t589 = a4*b1+a11*b8+a17*b9+a22*b10+a26*b11+a27*b12+a28*b13+a29*b14;
		double t599 = a4*b2+a11*b9+a17*b15+a22*b16+a26*b17+a27*b18+a28*b19+a29*b20;
		double t609 = a4*b3+a11*b10+a17*b16+a22*b21+a26*b22+a27*b23+a28*b24+a29*b25;
		double t612 = a27*b27;
		double t613 = a28*b28;
		double t614 = a29*b29;
		double t615 = t5+t158+t310+t461+a26*b26+t612+t613+t614;
		double t625 = a4*b5+a11*b12+a17*b18+a22*b23+a26*b27+a27*b30+a28*b31+a29*b32;
		double t635 = a4*b6+a11*b13+a17*b19+a22*b24+a26*b28+a27*b31+a28*b33+a29*b34;
		double t645 = a4*b7+a11*b14+a17*b20+a22*b25+a26*b29+a27*b32+a28*b34+a29*b35;
		double t719 = a5*b0+a12*b1+a18*b2+a23*b3+a27*b4+a30*b5+a31*b6+a32*b7;
		double t729 = a5*b1+a12*b8+a18*b9+a23*b10+a27*b11+a30*b12+a31*b13+a32*b14;
		double t739 = a5*b2+a12*b9+a18*b15+a23*b16+a27*b17+a30*b18+a31*b19+a32*b20;
		double t749 = a5*b3+a12*b10+a18*b16+a23*b21+a27*b22+a30*b23+a31*b24+a32*b25;
		double t759 = a5*b4+a12*b11+a18*b17+a23*b22+a27*b26+a30*b27+a31*b28+a32*b29;
		double t762 = a31*b31;
		double t763 = a32*b32;
		double t764 = t6+t159+t311+t462+t612+a30*b30+t762+t763;
		double t774 = a5*b6+a12*b13+a18*b19+a23*b24+a27*b28+a30*b31+a31*b33+a32*b34;
		double t784 = a5*b7+a12*b14+a18*b20+a23*b25+a27*b29+a30*b32+a31*b34+a32*b35;
		double t858 = a6*b0+a13*b1+a19*b2+a24*b3+a28*b4+a31*b5+a33*b6+a34*b7;
		double t868 = a6*b1+a13*b8+a19*b9+a24*b10+a28*b11+a31*b12+a33*b13+a34*b14;
		double t878 = a6*b2+a13*b9+a19*b15+a24*b16+a28*b17+a31*b18+a33*b19+a34*b20;
		double t888 = a6*b3+a13*b10+a19*b16+a24*b21+a28*b22+a31*b23+a33*b24+a34*b25;
		double t898 = a6*b4+a13*b11+a19*b17+a24*b22+a28*b26+a31*b27+a33*b28+a34*b29;
		double t908 = a6*b5+a13*b12+a19*b18+a24*b23+a28*b27+a31*b30+a33*b31+a34*b32;
		double t911 = a34*b34;
		double t912 = t7+t160+t312+t463+t613+t762+a33*b33+t911;
		double t922 = a6*b7+a13*b14+a19*b20+a24*b25+a28*b29+a31*b32+a33*b34+a34*b35;
		double t996 = a7*b0+a14*b1+a20*b2+a25*b3+a29*b4+a32*b5+a34*b6+a35*b7;
		double t1006 = a7*b1+a14*b8+a20*b9+a25*b10+a29*b11+a32*b12+a34*b13+a35*b14;
		double t1016 = a7*b2+a14*b9+a20*b15+a25*b16+a29*b17+a32*b18+a34*b19+a35*b20;
		double t1026 = a7*b3+a14*b10+a20*b16+a25*b21+a29*b22+a32*b23+a34*b24+a35*b25;
		double t1036 = a7*b4+a14*b11+a20*b17+a25*b22+a29*b26+a32*b27+a34*b28+a35*b29;
		double t1046 = a7*b5+a14*b12+a20*b18+a25*b23+a29*b27+a32*b30+a34*b31+a35*b32;
		double t1056 = a7*b6+a14*b13+a20*b19+a25*b24+a29*b28+a32*b31+a34*b33+a35*b34;
		double t1059 = t8+t161+t313+t464+t614+t763+t911+a35*b35;

		if ( mul2_mat )
			tmp_SM_ptr = mul2_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		*tmp_SM_ptr++ = t9*a0+t19*a1+t29*a2+t39*a3+t49*a4+t59*a5+t69*a6+t79*a7;
		*tmp_SM_ptr++ = t9*a1+t19*a8+t29*a9+t39*a10+t49*a11+t59*a12+t69*a13+t79*a14;
		*tmp_SM_ptr++ = t9*a2+t19*a9+t29*a15+t39*a16+t49*a17+t59*a18+t69*a19+t79*a20;
		*tmp_SM_ptr++ = t9*a3+t19*a10+t29*a16+t39*a21+t49*a22+t59*a23+t69*a24+t79*a25;
		*tmp_SM_ptr++ = t9*a4+t19*a11+t29*a17+t39*a22+t49*a26+t59*a27+t69*a28+t79*a29;
		*tmp_SM_ptr++ = t9*a5+t19*a12+t29*a18+t39*a23+t49*a27+t59*a30+t69*a31+t79*a32;
		*tmp_SM_ptr++ = t9*a6+t19*a13+t29*a19+t39*a24+t49*a28+t59*a31+t69*a33+t79*a34;
		*tmp_SM_ptr++ = t9*a7+t19*a14+t29*a20+t39*a25+t49*a29+t59*a32+t69*a34+t79*a35;
		*tmp_SM_ptr++ = t153*a1+t162*a8+t172*a9+t182*a10+t192*a11+t202*a12+t212*a13+t222*a14;
		*tmp_SM_ptr++ = t153*a2+t162*a9+t172*a15+t182*a16+t192*a17+t202*a18+t212*a19+t222*a20;
		*tmp_SM_ptr++ = t153*a3+t162*a10+t172*a16+t182*a21+t192*a22+t202*a23+t212*a24+t222*a25;
		*tmp_SM_ptr++ = t153*a4+t162*a11+t172*a17+t182*a22+t192*a26+t202*a27+t212*a28+t222*a29;
		*tmp_SM_ptr++ = t153*a5+t162*a12+t172*a18+t182*a23+t192*a27+t202*a30+t212*a31+t222*a32;
		*tmp_SM_ptr++ = t153*a6+t162*a13+t172*a19+t182*a24+t192*a28+t202*a31+t212*a33+t222*a34;
		*tmp_SM_ptr++ = t153*a7+t162*a14+t172*a20+t182*a25+t192*a29+t202*a32+t212*a34+t222*a35;
		*tmp_SM_ptr++ = t296*a2+t306*a9+t314*a15+t324*a16+t334*a17+t344*a18+t354*a19+t364*a20;
		*tmp_SM_ptr++ = t296*a3+t306*a10+t314*a16+t324*a21+t334*a22+t344*a23+t354*a24+t364*a25;
		*tmp_SM_ptr++ = t296*a4+t306*a11+t314*a17+t324*a22+t334*a26+t344*a27+t354*a28+t364*a29;
		*tmp_SM_ptr++ = t296*a5+t306*a12+t314*a18+t324*a23+t334*a27+t344*a30+t354*a31+t364*a32;
		*tmp_SM_ptr++ = t296*a6+t306*a13+t314*a19+t324*a24+t334*a28+t344*a31+t354*a33+t364*a34;
		*tmp_SM_ptr++ = t296*a7+t306*a14+t314*a20+t324*a25+t334*a29+t344*a32+t354*a34+t364*a35;
		*tmp_SM_ptr++ = t438*a3+t448*a10+t458*a16+t465*a21+t475*a22+t485*a23+t495*a24+t505*a25;
		*tmp_SM_ptr++ = t438*a4+t448*a11+t458*a17+t465*a22+t475*a26+t485*a27+t495*a28+t505*a29;
		*tmp_SM_ptr++ = t438*a5+t448*a12+t458*a18+t465*a23+t475*a27+t485*a30+t495*a31+t505*a32;
		*tmp_SM_ptr++ = t438*a6+t448*a13+t458*a19+t465*a24+t475*a28+t485*a31+t495*a33+t505*a34;
		*tmp_SM_ptr++ = t438*a7+t448*a14+t458*a20+t465*a25+t475*a29+t485*a32+t495*a34+t505*a35;
		*tmp_SM_ptr++ = t579*a4+t589*a11+t599*a17+t609*a22+t615*a26+t625*a27+t635*a28+t645*a29;
		*tmp_SM_ptr++ = t579*a5+t589*a12+t599*a18+t609*a23+t615*a27+t625*a30+t635*a31+t645*a32;
		*tmp_SM_ptr++ = t579*a6+t589*a13+t599*a19+t609*a24+t615*a28+t625*a31+t635*a33+t645*a34;
		*tmp_SM_ptr++ = t579*a7+t589*a14+t599*a20+t609*a25+t615*a29+t625*a32+t635*a34+t645*a35;
		*tmp_SM_ptr++ = t719*a5+t729*a12+t739*a18+t749*a23+t759*a27+t764*a30+t774*a31+t784*a32;
		*tmp_SM_ptr++ = t719*a6+t729*a13+t739*a19+t749*a24+t759*a28+t764*a31+t774*a33+t784*a34;
		*tmp_SM_ptr++ = t719*a7+t729*a14+t739*a20+t749*a25+t759*a29+t764*a32+t774*a34+t784*a35;
		*tmp_SM_ptr++ = t858*a6+t868*a13+t878*a19+t888*a24+t898*a28+t908*a31+t912*a33+t922*a34;
		*tmp_SM_ptr++ = t858*a7+t868*a14+t878*a20+t888*a25+t898*a29+t908*a32+t912*a34+t922*a35;
		*tmp_SM_ptr++ = t996*a7+t1006*a14+t1016*a20+t1026*a25+t1036*a29+t1046*a32+t1056*a34+t1059*a35;
	}
	else if ( this->M_size == 9 ) {
		double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35;
		double a36,a37,a38,a39,a40,a41,a42,a43,a44;
		double b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34, b35;
		double b36,b37,b38,b39,b40,b41,b42,b43,b44;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;
		a6 = *tmp_SM_ptr++;
		a7 = *tmp_SM_ptr++;
		a8 = *tmp_SM_ptr++;
		a9 = *tmp_SM_ptr++;
		a10 = *tmp_SM_ptr++;
		a11 = *tmp_SM_ptr++;
		a12 = *tmp_SM_ptr++;
		a13 = *tmp_SM_ptr++;
		a14 = *tmp_SM_ptr++;
		a15 = *tmp_SM_ptr++;
		a16 = *tmp_SM_ptr++;
		a17 = *tmp_SM_ptr++;
		a18 = *tmp_SM_ptr++;
		a19 = *tmp_SM_ptr++;
		a20 = *tmp_SM_ptr++;
		a21 = *tmp_SM_ptr++;
		a22 = *tmp_SM_ptr++;
		a23 = *tmp_SM_ptr++;
		a24 = *tmp_SM_ptr++;
		a25 = *tmp_SM_ptr++;
		a26 = *tmp_SM_ptr++;
		a27 = *tmp_SM_ptr++;
		a28 = *tmp_SM_ptr++;
		a29 = *tmp_SM_ptr++;
		a30 = *tmp_SM_ptr++;
		a31 = *tmp_SM_ptr++;
		a32 = *tmp_SM_ptr++;
		a33 = *tmp_SM_ptr++;
		a34 = *tmp_SM_ptr++;
		a35 = *tmp_SM_ptr++;
		a36 = *tmp_SM_ptr++;
		a37 = *tmp_SM_ptr++;
		a38 = *tmp_SM_ptr++;
		a39 = *tmp_SM_ptr++;
		a40 = *tmp_SM_ptr++;
		a41 = *tmp_SM_ptr++;
		a42 = *tmp_SM_ptr++;
		a43 = *tmp_SM_ptr++;
		a44 = *tmp_SM_ptr++;

		tmp_SM_ptr = Y_mat->SM_ptr;

		b0 = *tmp_SM_ptr++;
		b1 = *tmp_SM_ptr++;
		b2 = *tmp_SM_ptr++;
		b3 = *tmp_SM_ptr++;
		b4 = *tmp_SM_ptr++;
		b5 = *tmp_SM_ptr++;
		b6 = *tmp_SM_ptr++;
		b7 = *tmp_SM_ptr++;
		b8 = *tmp_SM_ptr++;
		b9 = *tmp_SM_ptr++;
		b10 = *tmp_SM_ptr++;
		b11 = *tmp_SM_ptr++;
		b12 = *tmp_SM_ptr++;
		b13 = *tmp_SM_ptr++;
		b14 = *tmp_SM_ptr++;
		b15 = *tmp_SM_ptr++;
		b16 = *tmp_SM_ptr++;
		b17 = *tmp_SM_ptr++;
		b18 = *tmp_SM_ptr++;
		b19 = *tmp_SM_ptr++;
		b20 = *tmp_SM_ptr++;
		b21 = *tmp_SM_ptr++;
		b22 = *tmp_SM_ptr++;
		b23 = *tmp_SM_ptr++;
		b24 = *tmp_SM_ptr++;
		b25 = *tmp_SM_ptr++;
		b26 = *tmp_SM_ptr++;
		b27 = *tmp_SM_ptr++;
		b28 = *tmp_SM_ptr++;
		b29 = *tmp_SM_ptr++;
		b30 = *tmp_SM_ptr++;
		b31 = *tmp_SM_ptr++;
		b32 = *tmp_SM_ptr++;
		b33 = *tmp_SM_ptr++;
		b34 = *tmp_SM_ptr++;
		b35 = *tmp_SM_ptr++;
		b36 = *tmp_SM_ptr++;
		b37 = *tmp_SM_ptr++;
		b38 = *tmp_SM_ptr++;
		b39 = *tmp_SM_ptr++;
		b40 = *tmp_SM_ptr++;
		b41 = *tmp_SM_ptr++;
		b42 = *tmp_SM_ptr++;
		b43 = *tmp_SM_ptr++;
		b44 = *tmp_SM_ptr++;

		double t2 = a1*b1;
		double t3 = a2*b2;
		double t4 = a3*b3;
		double t5 = a4*b4;
		double t6 = a5*b5;
		double t7 = a6*b6;
		double t8 = a7*b7;
		double t9 = a8*b8;
		double t10 = a0*b0+t2+t3+t4+t5+t6+t7+t8+t9;
		double t21 = a0*b1+a1*b9+a2*b10+a3*b11+a4*b12+a5*b13+a6*b14+a7*b15+a8*b16;
		double t32 = a0*b2+a1*b10+a2*b17+a3*b18+a4*b19+a5*b20+a6*b21+a7*b22+a8*b23;
		double t43 = a0*b3+a1*b11+a2*b18+a3*b24+a4*b25+a5*b26+a6*b27+a7*b28+a8*b29;
		double t54 = a0*b4+a1*b12+a2*b19+a3*b25+a4*b30+a5*b31+a6*b32+a7*b33+a8*b34;
		double t65 = a0*b5+a1*b13+a2*b20+a3*b26+a4*b31+a5*b35+a6*b36+a7*b37+a8*b38;
		double t76 = a0*b6+a1*b14+a2*b21+a3*b27+a4*b32+a5*b36+a6*b39+a7*b40+a8*b41;
		double t87 = a0*b7+a1*b15+a2*b22+a3*b28+a4*b33+a5*b37+a6*b40+a7*b42+a8*b43;
		double t98 = a0*b8+a1*b16+a2*b23+a3*b29+a4*b34+a5*b38+a6*b41+a7*b43+a8*b44;
		double t190 = a1*b0+a9*b1+a10*b2+a11*b3+a12*b4+a13*b5+a14*b6+a15*b7+a16*b8;
		double t193 = a10*b10;
		double t194 = a11*b11;
		double t195 = a12*b12;
		double t196 = a13*b13;
		double t197 = a14*b14;
		double t198 = a15*b15;
		double t199 = a16*b16;
		double t200 = t2+a9*b9+t193+t194+t195+t196+t197+t198+t199;
		double t211 = a1*b2+a9*b10+a10*b17+a11*b18+a12*b19+a13*b20+a14*b21+a15*b22+a16*b23;
		double t222 = a1*b3+a9*b11+a10*b18+a11*b24+a12*b25+a13*b26+a14*b27+a15*b28+a16*b29;
		double t233 = a1*b4+a9*b12+a10*b19+a11*b25+a12*b30+a13*b31+a14*b32+a15*b33+a16*b34;
		double t244 = a1*b5+a9*b13+a10*b20+a11*b26+a12*b31+a13*b35+a14*b36+a15*b37+a16*b38;
		double t255 = a1*b6+a9*b14+a10*b21+a11*b27+a12*b32+a13*b36+a14*b39+a15*b40+a16*b41;
		double t266 = a1*b7+a9*b15+a10*b22+a11*b28+a12*b33+a13*b37+a14*b40+a15*b42+a16*b43;
		double t277 = a1*b8+a9*b16+a10*b23+a11*b29+a12*b34+a13*b38+a14*b41+a15*b43+a16*b44;
		double t369 = a2*b0+a10*b1+a17*b2+a18*b3+a19*b4+a20*b5+a21*b6+a22*b7+a23*b8;
		double t380 = a2*b1+a10*b9+a17*b10+a18*b11+a19*b12+a20*b13+a21*b14+a22*b15+a23*b16;
		double t383 = a18*b18;
		double t384 = a19*b19;
		double t385 = a20*b20;
		double t386 = a21*b21;
		double t387 = a22*b22;
		double t388 = a23*b23;
		double t389 = t3+t193+a17*b17+t383+t384+t385+t386+t387+t388;
		double t400 = a2*b3+a10*b11+a17*b18+a18*b24+a19*b25+a20*b26+a21*b27+a22*b28+a23*b29;
		double t411 = a2*b4+a10*b12+a17*b19+a18*b25+a19*b30+a20*b31+a21*b32+a22*b33+a23*b34;
		double t422 = a2*b5+a10*b13+a17*b20+a18*b26+a19*b31+a20*b35+a21*b36+a22*b37+a23*b38;
		double t433 = a2*b6+a10*b14+a17*b21+a18*b27+a19*b32+a20*b36+a21*b39+a22*b40+a23*b41;
		double t444 = a2*b7+a10*b15+a17*b22+a18*b28+a19*b33+a20*b37+a21*b40+a22*b42+a23*b43;
		double t455 = a2*b8+a10*b16+a17*b23+a18*b29+a19*b34+a20*b38+a21*b41+a22*b43+a23*b44;
		double t547 = a3*b0+a11*b1+a18*b2+a24*b3+a25*b4+a26*b5+a27*b6+a28*b7+a29*b8;
		double t558 = a3*b1+a11*b9+a18*b10+a24*b11+a25*b12+a26*b13+a27*b14+a28*b15+a29*b16;
		double t569 = a3*b2+a11*b10+a18*b17+a24*b18+a25*b19+a26*b20+a27*b21+a28*b22+a29*b23;
		double t572 = a25*b25;
		double t573 = a26*b26;
		double t574 = a27*b27;
		double t575 = a28*b28;
		double t576 = a29*b29;
		double t577 = t4+t194+t383+a24*b24+t572+t573+t574+t575+t576;
		double t588 = a3*b4+a11*b12+a18*b19+a24*b25+a25*b30+a26*b31+a27*b32+a28*b33+a29*b34;
		double t599 = a3*b5+a11*b13+a18*b20+a24*b26+a25*b31+a26*b35+a27*b36+a28*b37+a29*b38;
		double t610 = a3*b6+a11*b14+a18*b21+a24*b27+a25*b32+a26*b36+a27*b39+a28*b40+a29*b41;
		double t621 = a3*b7+a11*b15+a18*b22+a24*b28+a25*b33+a26*b37+a27*b40+a28*b42+a29*b43;
		double t632 = a3*b8+a11*b16+a18*b23+a24*b29+a25*b34+a26*b38+a27*b41+a28*b43+a29*b44;
		double t724 = a4*b0+a12*b1+a19*b2+a25*b3+a30*b4+a31*b5+a32*b6+a33*b7+a34*b8;
		double t735 = a4*b1+a12*b9+a19*b10+a25*b11+a30*b12+a31*b13+a32*b14+a33*b15+a34*b16;
		double t746 = a4*b2+a12*b10+a19*b17+a25*b18+a30*b19+a31*b20+a32*b21+a33*b22+a34*b23;
		double t757 = a4*b3+a12*b11+a19*b18+a25*b24+a30*b25+a31*b26+a32*b27+a33*b28+a34*b29;
		double t760 = a31*b31;
		double t761 = a32*b32;
		double t762 = a33*b33;
		double t763 = a34*b34;
		double t764 = t5+t195+t384+t572+a30*b30+t760+t761+t762+t763;
		double t775 = a4*b5+a12*b13+a19*b20+a25*b26+a30*b31+a31*b35+a32*b36+a33*b37+a34*b38;
		double t786 = a4*b6+a12*b14+a19*b21+a25*b27+a30*b32+a31*b36+a32*b39+a33*b40+a34*b41;
		double t797 = a4*b7+a12*b15+a19*b22+a25*b28+a30*b33+a31*b37+a32*b40+a33*b42+a34*b43;
		double t808 = a4*b8+a12*b16+a19*b23+a25*b29+a30*b34+a31*b38+a32*b41+a33*b43+a34*b44;
		double t900 = a5*b0+a13*b1+a20*b2+a26*b3+a31*b4+a35*b5+a36*b6+a37*b7+a38*b8;
		double t911 = a5*b1+a13*b9+a20*b10+a26*b11+a31*b12+a35*b13+a36*b14+a37*b15+a38*b16;
		double t922 = a5*b2+a13*b10+a20*b17+a26*b18+a31*b19+a35*b20+a36*b21+a37*b22+a38*b23;
		double t933 = a5*b3+a13*b11+a20*b18+a26*b24+a31*b25+a35*b26+a36*b27+a37*b28+a38*b29;
		double t944 = a5*b4+a13*b12+a20*b19+a26*b25+a31*b30+a35*b31+a36*b32+a37*b33+a38*b34;
		double t947 = a36*b36;
		double t948 = a37*b37;
		double t949 = a38*b38;
		double t950 = t6+t196+t385+t573+t760+a35*b35+t947+t948+t949;
		double t961 = a5*b6+a13*b14+a20*b21+a26*b27+a31*b32+a35*b36+a36*b39+a37*b40+a38*b41;
		double t972 = a5*b7+a13*b15+a20*b22+a26*b28+a31*b33+a35*b37+a36*b40+a37*b42+a38*b43;
		double t983 = a5*b8+a13*b16+a20*b23+a26*b29+a31*b34+a35*b38+a36*b41+a37*b43+a38*b44;
		double t1075 = a6*b0+a14*b1+a21*b2+a27*b3+a32*b4+a36*b5+a39*b6+a40*b7+a41*b8;
		double t1086 = a6*b1+a14*b9+a21*b10+a27*b11+a32*b12+a36*b13+a39*b14+a40*b15+a41*b16;
		double t1097 = a6*b2+a14*b10+a21*b17+a27*b18+a32*b19+a36*b20+a39*b21+a40*b22+a41*b23;
		double t1108 = a6*b3+a14*b11+a21*b18+a27*b24+a32*b25+a36*b26+a39*b27+a40*b28+a41*b29;
		double t1119 = a6*b4+a14*b12+a21*b19+a27*b25+a32*b30+a36*b31+a39*b32+a40*b33+a41*b34;
		double t1130 = a6*b5+a14*b13+a21*b20+a27*b26+a32*b31+a36*b35+a39*b36+a40*b37+a41*b38;
		double t1133 = a40*b40;
		double t1134 = a41*b41;
		double t1135 = t7+t197+t386+t574+t761+t947+a39*b39+t1133+t1134;
		double t1146 = a6*b7+a14*b15+a21*b22+a27*b28+a32*b33+a36*b37+a39*b40+a40*b42+a41*b43;
		double t1157 = a6*b8+a14*b16+a21*b23+a27*b29+a32*b34+a36*b38+a39*b41+a40*b43+a41*b44;
		double t1249 = a7*b0+a15*b1+a22*b2+a28*b3+a33*b4+a37*b5+a40*b6+a42*b7+a43*b8;
		double t1260 = a7*b1+a15*b9+a22*b10+a28*b11+a33*b12+a37*b13+a40*b14+a42*b15+a43*b16;
		double t1271 = a7*b2+a15*b10+a22*b17+a28*b18+a33*b19+a37*b20+a40*b21+a42*b22+a43*b23;
		double t1282 = a7*b3+a15*b11+a22*b18+a28*b24+a33*b25+a37*b26+a40*b27+a42*b28+a43*b29;
		double t1293 = a7*b4+a15*b12+a22*b19+a28*b25+a33*b30+a37*b31+a40*b32+a42*b33+a43*b34;
		double t1304 = a7*b5+a15*b13+a22*b20+a28*b26+a33*b31+a37*b35+a40*b36+a42*b37+a43*b38;
		double t1315 = a7*b6+a15*b14+a22*b21+a28*b27+a33*b32+a37*b36+a40*b39+a42*b40+a43*b41;
		double t1318 = a43*b43;
		double t1319 = t8+t198+t387+t575+t762+t948+t1133+a42*b42+t1318;
		double t1330 = a7*b8+a15*b16+a22*b23+a28*b29+a33*b34+a37*b38+a40*b41+a42*b43+a43*b44;
		double t1422 = a8*b0+a16*b1+a23*b2+a29*b3+a34*b4+a38*b5+a41*b6+a43*b7+a44*b8;
		double t1433 = a8*b1+a16*b9+a23*b10+a29*b11+a34*b12+a38*b13+a41*b14+a43*b15+a44*b16;
		double t1444 = a8*b2+a16*b10+a23*b17+a29*b18+a34*b19+a38*b20+a41*b21+a43*b22+a44*b23;
		double t1455 = a8*b3+a16*b11+a23*b18+a29*b24+a34*b25+a38*b26+a41*b27+a43*b28+a44*b29;
		double t1466 = a8*b4+a16*b12+a23*b19+a29*b25+a34*b30+a38*b31+a41*b32+a43*b33+a44*b34;
		double t1477 = a8*b5+a16*b13+a23*b20+a29*b26+a34*b31+a38*b35+a41*b36+a43*b37+a44*b38;
		double t1488 = a8*b6+a16*b14+a23*b21+a29*b27+a34*b32+a38*b36+a41*b39+a43*b40+a44*b41;
		double t1499 = a8*b7+a16*b15+a23*b22+a29*b28+a34*b33+a38*b37+a41*b40+a43*b42+a44*b43;
		double t1502 = t9+t199+t388+t576+t763+t949+t1134+t1318+a44*b44;


		if ( mul2_mat )
			tmp_SM_ptr = mul2_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		*tmp_SM_ptr++ = t10*a0+t21*a1+t32*a2+t43*a3+t54*a4+t65*a5+t76*a6+t87*a7+t98*a8;
		*tmp_SM_ptr++ = t10*a1+t21*a9+t32*a10+t43*a11+t54*a12+t65*a13+t76*a14+t87*a15+t98*a16;
		*tmp_SM_ptr++ = t10*a2+t21*a10+t32*a17+t43*a18+t54*a19+t65*a20+t76*a21+t87*a22+t98*a23;
		*tmp_SM_ptr++ = t10*a3+t21*a11+t32*a18+t43*a24+t54*a25+t65*a26+t76*a27+t87*a28+t98*a29;
		*tmp_SM_ptr++ = t10*a4+t21*a12+t32*a19+t43*a25+t54*a30+t65*a31+t76*a32+t87*a33+t98*a34;
		*tmp_SM_ptr++ = t10*a5+t21*a13+t32*a20+t43*a26+t54*a31+t65*a35+t76*a36+t87*a37+t98*a38;
		*tmp_SM_ptr++ = t10*a6+t21*a14+t32*a21+t43*a27+t54*a32+t65*a36+t76*a39+t87*a40+t98*a41;
		*tmp_SM_ptr++ = t10*a7+t21*a15+t32*a22+t43*a28+t54*a33+t65*a37+t76*a40+t87*a42+t98*a43;
		*tmp_SM_ptr++ = t10*a8+t21*a16+t32*a23+t43*a29+t54*a34+t65*a38+t76*a41+t87*a43+t98*a44;
		*tmp_SM_ptr++ = t190*a1+t200*a9+t211*a10+t222*a11+t233*a12+t244*a13+t255*a14+t266*a15+t277*a16;
		*tmp_SM_ptr++ = t190*a2+t200*a10+t211*a17+t222*a18+t233*a19+t244*a20+t255*a21+t266*a22+t277*a23;
		*tmp_SM_ptr++ = t190*a3+t200*a11+t211*a18+t222*a24+t233*a25+t244*a26+t255*a27+t266*a28+t277*a29;
		*tmp_SM_ptr++ = t190*a4+t200*a12+t211*a19+t222*a25+t233*a30+t244*a31+t255*a32+t266*a33+t277*a34;
		*tmp_SM_ptr++ = t190*a5+t200*a13+t211*a20+t222*a26+t233*a31+t244*a35+t255*a36+t266*a37+t277*a38;
		*tmp_SM_ptr++ = t190*a6+t200*a14+t211*a21+t222*a27+t233*a32+t244*a36+t255*a39+t266*a40+t277*a41;
		*tmp_SM_ptr++ = t190*a7+t200*a15+t211*a22+t222*a28+t233*a33+t244*a37+t255*a40+t266*a42+t277*a43;
		*tmp_SM_ptr++ = t190*a8+t200*a16+t211*a23+t222*a29+t233*a34+t244*a38+t255*a41+t266*a43+t277*a44;
		*tmp_SM_ptr++ = t369*a2+t380*a10+t389*a17+t400*a18+t411*a19+t422*a20+t433*a21+t444*a22+t455*a23;
		*tmp_SM_ptr++ = t369*a3+t380*a11+t389*a18+t400*a24+t411*a25+t422*a26+t433*a27+t444*a28+t455*a29;
		*tmp_SM_ptr++ = t369*a4+t380*a12+t389*a19+t400*a25+t411*a30+t422*a31+t433*a32+t444*a33+t455*a34;
		*tmp_SM_ptr++ = t369*a5+t380*a13+t389*a20+t400*a26+t411*a31+t422*a35+t433*a36+t444*a37+t455*a38;
		*tmp_SM_ptr++ = t369*a6+t380*a14+t389*a21+t400*a27+t411*a32+t422*a36+t433*a39+t444*a40+t455*a41;
		*tmp_SM_ptr++ = t369*a7+t380*a15+t389*a22+t400*a28+t411*a33+t422*a37+t433*a40+t444*a42+t455*a43;
		*tmp_SM_ptr++ = t369*a8+t380*a16+t389*a23+t400*a29+t411*a34+t422*a38+t433*a41+t444*a43+t455*a44;
		*tmp_SM_ptr++ = t547*a3+t558*a11+t569*a18+t577*a24+t588*a25+t599*a26+t610*a27+t621*a28+t632*a29;
		*tmp_SM_ptr++ = t547*a4+t558*a12+t569*a19+t577*a25+t588*a30+t599*a31+t610*a32+t621*a33+t632*a34;
		*tmp_SM_ptr++ = t547*a5+t558*a13+t569*a20+t577*a26+t588*a31+t599*a35+t610*a36+t621*a37+t632*a38;
		*tmp_SM_ptr++ = t547*a6+t558*a14+t569*a21+t577*a27+t588*a32+t599*a36+t610*a39+t621*a40+t632*a41;
		*tmp_SM_ptr++ = t547*a7+t558*a15+t569*a22+t577*a28+t588*a33+t599*a37+t610*a40+t621*a42+t632*a43;
		*tmp_SM_ptr++ = t547*a8+t558*a16+t569*a23+t577*a29+t588*a34+t599*a38+t610*a41+t621*a43+t632*a44;
		*tmp_SM_ptr++ = t724*a4+t735*a12+t746*a19+t757*a25+t764*a30+t775*a31+t786*a32+t797*a33+t808*a34;
		*tmp_SM_ptr++ = t724*a5+t735*a13+t746*a20+t757*a26+t764*a31+t775*a35+t786*a36+t797*a37+t808*a38;
		*tmp_SM_ptr++ = t724*a6+t735*a14+t746*a21+t757*a27+t764*a32+t775*a36+t786*a39+t797*a40+t808*a41;
		*tmp_SM_ptr++ = t724*a7+t735*a15+t746*a22+t757*a28+t764*a33+t775*a37+t786*a40+t797*a42+t808*a43;
		*tmp_SM_ptr++ = t724*a8+t735*a16+t746*a23+t757*a29+t764*a34+t775*a38+t786*a41+t797*a43+t808*a44;
		*tmp_SM_ptr++ = t900*a5+t911*a13+t922*a20+t933*a26+t944*a31+t950*a35+t961*a36+t972*a37+t983*a38;
		*tmp_SM_ptr++ = t900*a6+t911*a14+t922*a21+t933*a27+t944*a32+t950*a36+t961*a39+t972*a40+t983*a41;
		*tmp_SM_ptr++ = t900*a7+t911*a15+t922*a22+t933*a28+t944*a33+t950*a37+t961*a40+t972*a42+t983*a43;
		*tmp_SM_ptr++ = t900*a8+t911*a16+t922*a23+t933*a29+t944*a34+t950*a38+t961*a41+t972*a43+t983*a44;
		*tmp_SM_ptr++ = t1075*a6+t1086*a14+t1097*a21+t1108*a27+t1119*a32+t1130*a36+t1135*a39+t1146*a40+t1157*a41;
		*tmp_SM_ptr++ = t1075*a7+t1086*a15+t1097*a22+t1108*a28+t1119*a33+t1130*a37+t1135*a40+t1146*a42+t1157*a43;
		*tmp_SM_ptr++ = t1075*a8+t1086*a16+t1097*a23+t1108*a29+t1119*a34+t1130*a38+t1135*a41+t1146*a43+t1157*a44;
		*tmp_SM_ptr++ = t1249*a7+t1260*a15+t1271*a22+t1282*a28+t1293*a33+t1304*a37+t1315*a40+t1319*a42+t1330*a43;
		*tmp_SM_ptr++ = t1249*a8+t1260*a16+t1271*a23+t1282*a29+t1293*a34+t1304*a38+t1315*a41+t1319*a43+t1330*a44;
		*tmp_SM_ptr++ = t1422*a8+t1433*a16+t1444*a23+t1455*a29+t1466*a34+t1477*a38+t1488*a41+t1499*a43+t1502*a44;
	}
	else if ( this->M_size == 10 ) {
		double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35;
		double a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50,a51,a52,a53,a54;
		double b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34, b35;
		double b36,b37,b38,b39,b40,b41,b42,b43,b44,b45,b46,b47,b48,b49,b50,b51,b52,b53,b54;

		double *tmp_SM_ptr = this->SM_ptr;
		a0 = *tmp_SM_ptr++;
		a1 = *tmp_SM_ptr++;
		a2 = *tmp_SM_ptr++;
		a3 = *tmp_SM_ptr++;
		a4 = *tmp_SM_ptr++;
		a5 = *tmp_SM_ptr++;
		a6 = *tmp_SM_ptr++;
		a7 = *tmp_SM_ptr++;
		a8 = *tmp_SM_ptr++;
		a9 = *tmp_SM_ptr++;
		a10 = *tmp_SM_ptr++;
		a11 = *tmp_SM_ptr++;
		a12 = *tmp_SM_ptr++;
		a13 = *tmp_SM_ptr++;
		a14 = *tmp_SM_ptr++;
		a15 = *tmp_SM_ptr++;
		a16 = *tmp_SM_ptr++;
		a17 = *tmp_SM_ptr++;
		a18 = *tmp_SM_ptr++;
		a19 = *tmp_SM_ptr++;
		a20 = *tmp_SM_ptr++;
		a21 = *tmp_SM_ptr++;
		a22 = *tmp_SM_ptr++;
		a23 = *tmp_SM_ptr++;
		a24 = *tmp_SM_ptr++;
		a25 = *tmp_SM_ptr++;
		a26 = *tmp_SM_ptr++;
		a27 = *tmp_SM_ptr++;
		a28 = *tmp_SM_ptr++;
		a29 = *tmp_SM_ptr++;
		a30 = *tmp_SM_ptr++;
		a31 = *tmp_SM_ptr++;
		a32 = *tmp_SM_ptr++;
		a33 = *tmp_SM_ptr++;
		a34 = *tmp_SM_ptr++;
		a35 = *tmp_SM_ptr++;
		a36 = *tmp_SM_ptr++;
		a37 = *tmp_SM_ptr++;
		a38 = *tmp_SM_ptr++;
		a39 = *tmp_SM_ptr++;
		a40 = *tmp_SM_ptr++;
		a41 = *tmp_SM_ptr++;
		a42 = *tmp_SM_ptr++;
		a43 = *tmp_SM_ptr++;
		a44 = *tmp_SM_ptr++;
		a45 = *tmp_SM_ptr++;
		a46 = *tmp_SM_ptr++;
		a47 = *tmp_SM_ptr++;
		a48 = *tmp_SM_ptr++;
		a49 = *tmp_SM_ptr++;
		a50 = *tmp_SM_ptr++;
		a51 = *tmp_SM_ptr++;
		a52 = *tmp_SM_ptr++;
		a53 = *tmp_SM_ptr++;
		a54 = *tmp_SM_ptr++;

		tmp_SM_ptr = Y_mat->SM_ptr;

		b0 = *tmp_SM_ptr++;
		b1 = *tmp_SM_ptr++;
		b2 = *tmp_SM_ptr++;
		b3 = *tmp_SM_ptr++;
		b4 = *tmp_SM_ptr++;
		b5 = *tmp_SM_ptr++;
		b6 = *tmp_SM_ptr++;
		b7 = *tmp_SM_ptr++;
		b8 = *tmp_SM_ptr++;
		b9 = *tmp_SM_ptr++;
		b10 = *tmp_SM_ptr++;
		b11 = *tmp_SM_ptr++;
		b12 = *tmp_SM_ptr++;
		b13 = *tmp_SM_ptr++;
		b14 = *tmp_SM_ptr++;
		b15 = *tmp_SM_ptr++;
		b16 = *tmp_SM_ptr++;
		b17 = *tmp_SM_ptr++;
		b18 = *tmp_SM_ptr++;
		b19 = *tmp_SM_ptr++;
		b20 = *tmp_SM_ptr++;
		b21 = *tmp_SM_ptr++;
		b22 = *tmp_SM_ptr++;
		b23 = *tmp_SM_ptr++;
		b24 = *tmp_SM_ptr++;
		b25 = *tmp_SM_ptr++;
		b26 = *tmp_SM_ptr++;
		b27 = *tmp_SM_ptr++;
		b28 = *tmp_SM_ptr++;
		b29 = *tmp_SM_ptr++;
		b30 = *tmp_SM_ptr++;
		b31 = *tmp_SM_ptr++;
		b32 = *tmp_SM_ptr++;
		b33 = *tmp_SM_ptr++;
		b34 = *tmp_SM_ptr++;
		b35 = *tmp_SM_ptr++;
		b36 = *tmp_SM_ptr++;
		b37 = *tmp_SM_ptr++;
		b38 = *tmp_SM_ptr++;
		b39 = *tmp_SM_ptr++;
		b40 = *tmp_SM_ptr++;
		b41 = *tmp_SM_ptr++;
		b42 = *tmp_SM_ptr++;
		b43 = *tmp_SM_ptr++;
		b44 = *tmp_SM_ptr++;
		b45 = *tmp_SM_ptr++;
		b46 = *tmp_SM_ptr++;
		b47 = *tmp_SM_ptr++;
		b48 = *tmp_SM_ptr++;
		b49 = *tmp_SM_ptr++;
		b50 = *tmp_SM_ptr++;
		b51 = *tmp_SM_ptr++;
		b52 = *tmp_SM_ptr++;
		b53 = *tmp_SM_ptr++;
		b54 = *tmp_SM_ptr++;

		double t2 = a1*b1;
		double t3 = a2*b2;
		double t4 = a3*b3;
		double t5 = a4*b4;
		double t6 = a5*b5;
		double t7 = a6*b6;
		double t8 = a7*b7;
		double t9 = a8*b8;
		double t10 = a9*b9;
		double t11 = a0*b0+t2+t3+t4+t5+t6+t7+t8+t9+t10;
		double t23 = a0*b1+a1*b10+a2*b11+a3*b12+a4*b13+a5*b14+a6*b15+a7*b16+a8*b17+a9*b18;
		double t35 = a0*b2+a1*b11+a2*b19+a3*b20+a4*b21+a5*b22+a6*b23+a7*b24+a8*b25+a9*b26;
		double t47 = a0*b3+a1*b12+a2*b20+a3*b27+a4*b28+a5*b29+a6*b30+a7*b31+a8*b32+a9*b33;
		double t59 = a0*b4+a1*b13+a2*b21+a3*b28+a4*b34+a5*b35+a6*b36+a7*b37+a8*b38+a9*b39;
		double t71 = a0*b5+a1*b14+a2*b22+a3*b29+a4*b35+a5*b40+a6*b41+a7*b42+a8*b43+a9*b44;
		double t83 = a0*b6+a1*b15+a2*b23+a3*b30+a4*b36+a5*b41+a6*b45+a7*b46+a8*b47+a9*b48;
		double t95 = a0*b7+a1*b16+a2*b24+a3*b31+a4*b37+a5*b42+a6*b46+a7*b49+a8*b50+a9*b51;
		double t107 = a0*b8+a1*b17+a2*b25+a3*b32+a4*b38+a5*b43+a6*b47+a7*b50+a8*b52+a9*b53;
		double t119 = a0*b9+a1*b18+a2*b26+a3*b33+a4*b39+a5*b44+a6*b48+a7*b51+a8*b53+a9*b54;
		double t231 = a1*b0+a10*b1+a11*b2+a12*b3+a13*b4+a14*b5+a15*b6+a16*b7+a17*b8+a18*b9;
		double t234 = a11*b11;
		double t235 = a12*b12;
		double t236 = a13*b13;
		double t237 = a14*b14;
		double t238 = a15*b15;
		double t239 = a16*b16;
		double t240 = a17*b17;
		double t241 = a18*b18;
		double t242 = t2+a10*b10+t234+t235+t236+t237+t238+t239+t240+t241;
		double t254 = a1*b2+a10*b11+a11*b19+a12*b20+a13*b21+a14*b22+a15*b23+a16*b24+a17*b25+a18*b26;
		double t266 = a1*b3+a10*b12+a11*b20+a12*b27+a13*b28+a14*b29+a15*b30+a16*b31+a17*b32+a18*b33;
		double t278 = a1*b4+a10*b13+a11*b21+a12*b28+a13*b34+a14*b35+a15*b36+a16*b37+a17*b38+a18*b39;
		double t290 = a1*b5+a10*b14+a11*b22+a12*b29+a13*b35+a14*b40+a15*b41+a16*b42+a17*b43+a18*b44;
		double t302 = a1*b6+a10*b15+a11*b23+a12*b30+a13*b36+a14*b41+a15*b45+a16*b46+a17*b47+a18*b48;
		double t314 = a1*b7+a10*b16+a11*b24+a12*b31+a13*b37+a14*b42+a15*b46+a16*b49+a17*b50+a18*b51;
		double t326 = a1*b8+a10*b17+a11*b25+a12*b32+a13*b38+a14*b43+a15*b47+a16*b50+a17*b52+a18*b53;
		double t338 = a1*b9+a10*b18+a11*b26+a12*b33+a13*b39+a14*b44+a15*b48+a16*b51+a17*b53+a18*b54;
		double t450 = a2*b0+a11*b1+a19*b2+a20*b3+a21*b4+a22*b5+a23*b6+a24*b7+a25*b8+a26*b9;
		double t462 = a2*b1+a11*b10+a19*b11+a20*b12+a21*b13+a22*b14+a23*b15+a24*b16+a25*b17+a26*b18;
		double t465 = a20*b20;
		double t466 = a21*b21;
		double t467 = a22*b22;
		double t468 = a23*b23;
		double t469 = a24*b24;
		double t470 = a25*b25;
		double t471 = a26*b26;
		double t472 = t3+t234+a19*b19+t465+t466+t467+t468+t469+t470+t471;
		double t484 = a2*b3+a11*b12+a19*b20+a20*b27+a21*b28+a22*b29+a23*b30+a24*b31+a25*b32+a26*b33;
		double t496 = a2*b4+a11*b13+a19*b21+a20*b28+a21*b34+a22*b35+a23*b36+a24*b37+a25*b38+a26*b39;
		double t508 = a2*b5+a11*b14+a19*b22+a20*b29+a21*b35+a22*b40+a23*b41+a24*b42+a25*b43+a26*b44;
		double t520 = a2*b6+a11*b15+a19*b23+a20*b30+a21*b36+a22*b41+a23*b45+a24*b46+a25*b47+a26*b48;
		double t532 = a2*b7+a11*b16+a19*b24+a20*b31+a21*b37+a22*b42+a23*b46+a24*b49+a25*b50+a26*b51;
		double t544 = a2*b8+a11*b17+a19*b25+a20*b32+a21*b38+a22*b43+a23*b47+a24*b50+a25*b52+a26*b53;
		double t556 = a2*b9+a11*b18+a19*b26+a20*b33+a21*b39+a22*b44+a23*b48+a24*b51+a25*b53+a26*b54;
		double t668 = a3*b0+a12*b1+a20*b2+a27*b3+a28*b4+a29*b5+a30*b6+a31*b7+a32*b8+a33*b9;
		double t680 = a3*b1+a12*b10+a20*b11+a27*b12+a28*b13+a29*b14+a30*b15+a31*b16+a32*b17+a33*b18;
		double t692 = a3*b2+a12*b11+a20*b19+a27*b20+a28*b21+a29*b22+a30*b23+a31*b24+a32*b25+a33*b26;
		double t695 = a28*b28;
		double t696 = a29*b29;
		double t697 = a30*b30;
		double t698 = a31*b31;
		double t699 = a32*b32;
		double t700 = a33*b33;
		double t701 = t4+t235+t465+a27*b27+t695+t696+t697+t698+t699+t700;
		double t713 = a3*b4+a12*b13+a20*b21+a27*b28+a28*b34+a29*b35+a30*b36+a31*b37+a32*b38+a33*b39;
		double t725 = a3*b5+a12*b14+a20*b22+a27*b29+a28*b35+a29*b40+a30*b41+a31*b42+a32*b43+a33*b44;
		double t737 = a3*b6+a12*b15+a20*b23+a27*b30+a28*b36+a29*b41+a30*b45+a31*b46+a32*b47+a33*b48;
		double t749 = a3*b7+a12*b16+a20*b24+a27*b31+a28*b37+a29*b42+a30*b46+a31*b49+a32*b50+a33*b51;
		double t761 = a3*b8+a12*b17+a20*b25+a27*b32+a28*b38+a29*b43+a30*b47+a31*b50+a32*b52+a33*b53;
		double t773 = a3*b9+a12*b18+a20*b26+a27*b33+a28*b39+a29*b44+a30*b48+a31*b51+a32*b53+a33*b54;
		double t885 = a4*b0+a13*b1+a21*b2+a28*b3+a34*b4+a35*b5+a36*b6+a37*b7+a38*b8+a39*b9;
		double t897 = a4*b1+a13*b10+a21*b11+a28*b12+a34*b13+a35*b14+a36*b15+a37*b16+a38*b17+a39*b18;
		double t909 = a4*b2+a13*b11+a21*b19+a28*b20+a34*b21+a35*b22+a36*b23+a37*b24+a38*b25+a39*b26;
		double t921 = a4*b3+a13*b12+a21*b20+a28*b27+a34*b28+a35*b29+a36*b30+a37*b31+a38*b32+a39*b33;
		double t924 = a35*b35;
		double t925 = a36*b36;
		double t926 = a37*b37;
		double t927 = a38*b38;
		double t928 = a39*b39;
		double t929 = t5+t236+t466+t695+a34*b34+t924+t925+t926+t927+t928;
		double t941 = a4*b5+a13*b14+a21*b22+a28*b29+a34*b35+a35*b40+a36*b41+a37*b42+a38*b43+a39*b44;
		double t953 = a4*b6+a13*b15+a21*b23+a28*b30+a34*b36+a35*b41+a36*b45+a37*b46+a38*b47+a39*b48;
		double t965 = a4*b7+a13*b16+a21*b24+a28*b31+a34*b37+a35*b42+a36*b46+a37*b49+a38*b50+a39*b51;
		double t977 = a4*b8+a13*b17+a21*b25+a28*b32+a34*b38+a35*b43+a36*b47+a37*b50+a38*b52+a39*b53;
		double t989 = a4*b9+a13*b18+a21*b26+a28*b33+a34*b39+a35*b44+a36*b48+a37*b51+a38*b53+a39*b54;
		double t1101 = a5*b0+a14*b1+a22*b2+a29*b3+a35*b4+a40*b5+a41*b6+a42*b7+a43*b8+a44*b9;
		double t1113 = a5*b1+a14*b10+a22*b11+a29*b12+a35*b13+a40*b14+a41*b15+a42*b16+a43*b17+a44*b18;
		double t1125 = a5*b2+a14*b11+a22*b19+a29*b20+a35*b21+a40*b22+a41*b23+a42*b24+a43*b25+a44*b26;
		double t1137 = a5*b3+a14*b12+a22*b20+a29*b27+a35*b28+a40*b29+a41*b30+a42*b31+a43*b32+a44*b33;
		double t1149 = a5*b4+a14*b13+a22*b21+a29*b28+a35*b34+a40*b35+a41*b36+a42*b37+a43*b38+a44*b39;
		double t1152 = a41*b41;
		double t1153 = a42*b42;
		double t1154 = a43*b43;
		double t1155 = a44*b44;
		double t1156 = t6+t237+t467+t696+t924+a40*b40+t1152+t1153+t1154+t1155;
		double t1168 = a5*b6+a14*b15+a22*b23+a29*b30+a35*b36+a40*b41+a41*b45+a42*b46+a43*b47+a44*b48;
		double t1180 = a5*b7+a14*b16+a22*b24+a29*b31+a35*b37+a40*b42+a41*b46+a42*b49+a43*b50+a44*b51;
		double t1192 = a5*b8+a14*b17+a22*b25+a29*b32+a35*b38+a40*b43+a41*b47+a42*b50+a43*b52+a44*b53;
		double t1204 = a5*b9+a14*b18+a22*b26+a29*b33+a35*b39+a40*b44+a41*b48+a42*b51+a43*b53+a44*b54;
		double t1316 = a6*b0+a15*b1+a23*b2+a30*b3+a36*b4+a41*b5+a45*b6+a46*b7+a47*b8+a48*b9;
		double t1328 = a6*b1+a15*b10+a23*b11+a30*b12+a36*b13+a41*b14+a45*b15+a46*b16+a47*b17+a48*b18;
		double t1340 = a6*b2+a15*b11+a23*b19+a30*b20+a36*b21+a41*b22+a45*b23+a46*b24+a47*b25+a48*b26;
		double t1352 = a6*b3+a15*b12+a23*b20+a30*b27+a36*b28+a41*b29+a45*b30+a46*b31+a47*b32+a48*b33;
		double t1364 = a6*b4+a15*b13+a23*b21+a30*b28+a36*b34+a41*b35+a45*b36+a46*b37+a47*b38+a48*b39;
		double t1376 = a6*b5+a15*b14+a23*b22+a30*b29+a36*b35+a41*b40+a45*b41+a46*b42+a47*b43+a48*b44;
		double t1379 = a46*b46;
		double t1380 = a47*b47;
		double t1381 = a48*b48;
		double t1382 = t7+t238+t468+t697+t925+t1152+a45*b45+t1379+t1380+t1381;
		double t1394 = a6*b7+a15*b16+a23*b24+a30*b31+a36*b37+a41*b42+a45*b46+a46*b49+a47*b50+a48*b51;
		double t1406 = a6*b8+a15*b17+a23*b25+a30*b32+a36*b38+a41*b43+a45*b47+a46*b50+a47*b52+a48*b53;
		double t1418 = a6*b9+a15*b18+a23*b26+a30*b33+a36*b39+a41*b44+a45*b48+a46*b51+a47*b53+a48*b54;
		double t1530 = a7*b0+a16*b1+a24*b2+a31*b3+a37*b4+a42*b5+a46*b6+a49*b7+a50*b8+a51*b9;
		double t1542 = a7*b1+a16*b10+a24*b11+a31*b12+a37*b13+a42*b14+a46*b15+a49*b16+a50*b17+a51*b18;
		double t1554 = a7*b2+a16*b11+a24*b19+a31*b20+a37*b21+a42*b22+a46*b23+a49*b24+a50*b25+a51*b26;
		double t1566 = a7*b3+a16*b12+a24*b20+a31*b27+a37*b28+a42*b29+a46*b30+a49*b31+a50*b32+a51*b33;
		double t1578 = a7*b4+a16*b13+a24*b21+a31*b28+a37*b34+a42*b35+a46*b36+a49*b37+a50*b38+a51*b39;
		double t1590 = a7*b5+a16*b14+a24*b22+a31*b29+a37*b35+a42*b40+a46*b41+a49*b42+a50*b43+a51*b44;
		double t1602 = a7*b6+a16*b15+a24*b23+a31*b30+a37*b36+a42*b41+a46*b45+a49*b46+a50*b47+a51*b48;
		double t1605 = a50*b50;
		double t1606 = a51*b51;
		double t1607 = t8+t239+t469+t698+t926+t1153+t1379+a49*b49+t1605+t1606;
		double t1619 = a7*b8+a16*b17+a24*b25+a31*b32+a37*b38+a42*b43+a46*b47+a49*b50+a50*b52+a51*b53;
		double t1631 = a7*b9+a16*b18+a24*b26+a31*b33+a37*b39+a42*b44+a46*b48+a49*b51+a50*b53+a51*b54;
		double t1743 = a8*b0+a17*b1+a25*b2+a32*b3+a38*b4+a43*b5+a47*b6+a50*b7+a52*b8+a53*b9;
		double t1755 = a8*b1+a17*b10+a25*b11+a32*b12+a38*b13+a43*b14+a47*b15+a50*b16+a52*b17+a53*b18;
		double t1767 = a8*b2+a17*b11+a25*b19+a32*b20+a38*b21+a43*b22+a47*b23+a50*b24+a52*b25+a53*b26;
		double t1779 = a8*b3+a17*b12+a25*b20+a32*b27+a38*b28+a43*b29+a47*b30+a50*b31+a52*b32+a53*b33;
		double t1791 = a8*b4+a17*b13+a25*b21+a32*b28+a38*b34+a43*b35+a47*b36+a50*b37+a52*b38+a53*b39;
		double t1803 = a8*b5+a17*b14+a25*b22+a32*b29+a38*b35+a43*b40+a47*b41+a50*b42+a52*b43+a53*b44;
		double t1815 = a8*b6+a17*b15+a25*b23+a32*b30+a38*b36+a43*b41+a47*b45+a50*b46+a52*b47+a53*b48;
		double t1827 = a8*b7+a17*b16+a25*b24+a32*b31+a38*b37+a43*b42+a47*b46+a50*b49+a52*b50+a53*b51;
		double t1830 = a53*b53;
		double t1831 = t9+t240+t470+t699+t927+t1154+t1380+t1605+a52*b52+t1830;
		double t1843 = a8*b9+a17*b18+a25*b26+a32*b33+a38*b39+a43*b44+a47*b48+a50*b51+a52*b53+a53*b54;
		double t1955 = a9*b0+a18*b1+a26*b2+a33*b3+a39*b4+a44*b5+a48*b6+a51*b7+a53*b8+a54*b9;
		double t1967 = a9*b1+a18*b10+a26*b11+a33*b12+a39*b13+a44*b14+a48*b15+a51*b16+a53*b17+a54*b18;
		double t1979 = a9*b2+a18*b11+a26*b19+a33*b20+a39*b21+a44*b22+a48*b23+a51*b24+a53*b25+a54*b26;
		double t1991 = a9*b3+a18*b12+a26*b20+a33*b27+a39*b28+a44*b29+a48*b30+a51*b31+a53*b32+a54*b33;
		double t2003 = a9*b4+a18*b13+a26*b21+a33*b28+a39*b34+a44*b35+a48*b36+a51*b37+a53*b38+a54*b39;
		double t2015 = a9*b5+a18*b14+a26*b22+a33*b29+a39*b35+a44*b40+a48*b41+a51*b42+a53*b43+a54*b44;
		double t2027 = a9*b6+a18*b15+a26*b23+a33*b30+a39*b36+a44*b41+a48*b45+a51*b46+a53*b47+a54*b48;
		double t2039 = a9*b7+a18*b16+a26*b24+a33*b31+a39*b37+a44*b42+a48*b46+a51*b49+a53*b50+a54*b51;
		double t2051 = a9*b8+a18*b17+a26*b25+a33*b32+a39*b38+a44*b43+a48*b47+a51*b50+a53*b52+a54*b53;
		double t2054 = t10+t241+t471+t700+t928+t1155+t1381+t1606+t1830+a54*b54;


		if ( mul2_mat )
			tmp_SM_ptr = mul2_mat->SM_ptr;
		else
			tmp_SM_ptr = this->SM_ptr;

		*tmp_SM_ptr++ = t11*a0+t23*a1+t35*a2+t47*a3+t59*a4+t71*a5+t83*a6+t95*a7+t107*a8+t119*a9;
		*tmp_SM_ptr++ = t11*a1+t23*a10+t35*a11+t47*a12+t59*a13+t71*a14+t83*a15+t95*a16+t107*a17+t119*a18;
		*tmp_SM_ptr++ = t11*a2+t23*a11+t35*a19+t47*a20+t59*a21+t71*a22+t83*a23+t95*a24+t107*a25+t119*a26;
		*tmp_SM_ptr++ = t11*a3+t23*a12+t35*a20+t47*a27+t59*a28+t71*a29+t83*a30+t95*a31+t107*a32+t119*a33;
		*tmp_SM_ptr++ = t11*a4+t23*a13+t35*a21+t47*a28+t59*a34+t71*a35+t83*a36+t95*a37+t107*a38+t119*a39;
		*tmp_SM_ptr++ = t11*a5+t23*a14+t35*a22+t47*a29+t59*a35+t71*a40+t83*a41+t95*a42+t107*a43+t119*a44;
		*tmp_SM_ptr++ = t11*a6+t23*a15+t35*a23+t47*a30+t59*a36+t71*a41+t83*a45+t95*a46+t107*a47+t119*a48;
		*tmp_SM_ptr++ = t11*a7+t23*a16+t35*a24+t47*a31+t59*a37+t71*a42+t83*a46+t95*a49+t107*a50+t119*a51;
		*tmp_SM_ptr++ = t11*a8+t23*a17+t35*a25+t47*a32+t59*a38+t71*a43+t83*a47+t95*a50+t107*a52+t119*a53;
		*tmp_SM_ptr++ = t11*a9+t23*a18+t35*a26+t47*a33+t59*a39+t71*a44+t83*a48+t95*a51+t107*a53+t119*a54;
		*tmp_SM_ptr++ = t231*a1+t242*a10+t254*a11+t266*a12+t278*a13+t290*a14+t302*a15+t314*a16+t326*a17+t338*a18;
		*tmp_SM_ptr++ = t231*a2+t242*a11+t254*a19+t266*a20+t278*a21+t290*a22+t302*a23+t314*a24+t326*a25+t338*a26;
		*tmp_SM_ptr++ = t231*a3+t242*a12+t254*a20+t266*a27+t278*a28+t290*a29+t302*a30+t314*a31+t326*a32+t338*a33;
		*tmp_SM_ptr++ = t231*a4+t242*a13+t254*a21+t266*a28+t278*a34+t290*a35+t302*a36+t314*a37+t326*a38+t338*a39;
		*tmp_SM_ptr++ = t231*a5+t242*a14+t254*a22+t266*a29+t278*a35+t290*a40+t302*a41+t314*a42+t326*a43+t338*a44;
		*tmp_SM_ptr++ = t231*a6+t242*a15+t254*a23+t266*a30+t278*a36+t290*a41+t302*a45+t314*a46+t326*a47+t338*a48;
		*tmp_SM_ptr++ = t231*a7+t242*a16+t254*a24+t266*a31+t278*a37+t290*a42+t302*a46+t314*a49+t326*a50+t338*a51;
		*tmp_SM_ptr++ = t231*a8+t242*a17+t254*a25+t266*a32+t278*a38+t290*a43+t302*a47+t314*a50+t326*a52+t338*a53;
		*tmp_SM_ptr++ = t231*a9+t242*a18+t254*a26+t266*a33+t278*a39+t290*a44+t302*a48+t314*a51+t326*a53+t338*a54;
		*tmp_SM_ptr++ = t450*a2+t462*a11+t472*a19+t484*a20+t496*a21+t508*a22+t520*a23+t532*a24+t544*a25+t556*a26;
		*tmp_SM_ptr++ = t450*a3+t462*a12+t472*a20+t484*a27+t496*a28+t508*a29+t520*a30+t532*a31+t544*a32+t556*a33;
		*tmp_SM_ptr++ = t450*a4+t462*a13+t472*a21+t484*a28+t496*a34+t508*a35+t520*a36+t532*a37+t544*a38+t556*a39;
		*tmp_SM_ptr++ = t450*a5+t462*a14+t472*a22+t484*a29+t496*a35+t508*a40+t520*a41+t532*a42+t544*a43+t556*a44;
		*tmp_SM_ptr++ = t450*a6+t462*a15+t472*a23+t484*a30+t496*a36+t508*a41+t520*a45+t532*a46+t544*a47+t556*a48;
		*tmp_SM_ptr++ = t450*a7+t462*a16+t472*a24+t484*a31+t496*a37+t508*a42+t520*a46+t532*a49+t544*a50+t556*a51;
		*tmp_SM_ptr++ = t450*a8+t462*a17+t472*a25+t484*a32+t496*a38+t508*a43+t520*a47+t532*a50+t544*a52+t556*a53;
		*tmp_SM_ptr++ = t450*a9+t462*a18+t472*a26+t484*a33+t496*a39+t508*a44+t520*a48+t532*a51+t544*a53+t556*a54;
		*tmp_SM_ptr++ = t668*a3+t680*a12+t692*a20+t701*a27+t713*a28+t725*a29+t737*a30+t749*a31+t761*a32+t773*a33;
		*tmp_SM_ptr++ = t668*a4+t680*a13+t692*a21+t701*a28+t713*a34+t725*a35+t737*a36+t749*a37+t761*a38+t773*a39;
		*tmp_SM_ptr++ = t668*a5+t680*a14+t692*a22+t701*a29+t713*a35+t725*a40+t737*a41+t749*a42+t761*a43+t773*a44;
		*tmp_SM_ptr++ = t668*a6+t680*a15+t692*a23+t701*a30+t713*a36+t725*a41+t737*a45+t749*a46+t761*a47+t773*a48;
		*tmp_SM_ptr++ = t668*a7+t680*a16+t692*a24+t701*a31+t713*a37+t725*a42+t737*a46+t749*a49+t761*a50+t773*a51;
		*tmp_SM_ptr++ = t668*a8+t680*a17+t692*a25+t701*a32+t713*a38+t725*a43+t737*a47+t749*a50+t761*a52+t773*a53;
		*tmp_SM_ptr++ = t668*a9+t680*a18+t692*a26+t701*a33+t713*a39+t725*a44+t737*a48+t749*a51+t761*a53+t773*a54;
		*tmp_SM_ptr++ = t885*a4+t897*a13+t909*a21+t921*a28+t929*a34+t941*a35+t953*a36+t965*a37+t977*a38+t989*a39;
		*tmp_SM_ptr++ = t885*a5+t897*a14+t909*a22+t921*a29+t929*a35+t941*a40+t953*a41+t965*a42+t977*a43+t989*a44;
		*tmp_SM_ptr++ = t885*a6+t897*a15+t909*a23+t921*a30+t929*a36+t941*a41+t953*a45+t965*a46+t977*a47+t989*a48;
		*tmp_SM_ptr++ = t885*a7+t897*a16+t909*a24+t921*a31+t929*a37+t941*a42+t953*a46+t965*a49+t977*a50+t989*a51;
		*tmp_SM_ptr++ = t885*a8+t897*a17+t909*a25+t921*a32+t929*a38+t941*a43+t953*a47+t965*a50+t977*a52+t989*a53;
		*tmp_SM_ptr++ = t885*a9+t897*a18+t909*a26+t921*a33+t929*a39+t941*a44+t953*a48+t965*a51+t977*a53+t989*a54;
		*tmp_SM_ptr++ = t1101*a5+t1113*a14+t1125*a22+t1137*a29+t1149*a35+t1156*a40+t1168*a41+t1180*a42+t1192*a43+t1204*a44;
		*tmp_SM_ptr++ = t1101*a6+t1113*a15+t1125*a23+t1137*a30+t1149*a36+t1156*a41+t1168*a45+t1180*a46+t1192*a47+t1204*a48;
		*tmp_SM_ptr++ = t1101*a7+t1113*a16+t1125*a24+t1137*a31+t1149*a37+t1156*a42+t1168*a46+t1180*a49+t1192*a50+t1204*a51;
		*tmp_SM_ptr++ = t1101*a8+t1113*a17+t1125*a25+t1137*a32+t1149*a38+t1156*a43+t1168*a47+t1180*a50+t1192*a52+t1204*a53;
		*tmp_SM_ptr++ = t1101*a9+t1113*a18+t1125*a26+t1137*a33+t1149*a39+t1156*a44+t1168*a48+t1180*a51+t1192*a53+t1204*a54;
		*tmp_SM_ptr++ = t1316*a6+t1328*a15+t1340*a23+t1352*a30+t1364*a36+t1376*a41+t1382*a45+t1394*a46+t1406*a47+t1418*a48;
		*tmp_SM_ptr++ = t1316*a7+t1328*a16+t1340*a24+t1352*a31+t1364*a37+t1376*a42+t1382*a46+t1394*a49+t1406*a50+t1418*a51;
		*tmp_SM_ptr++ = t1316*a8+t1328*a17+t1340*a25+t1352*a32+t1364*a38+t1376*a43+t1382*a47+t1394*a50+t1406*a52+t1418*a53;
		*tmp_SM_ptr++ = t1316*a9+t1328*a18+t1340*a26+t1352*a33+t1364*a39+t1376*a44+t1382*a48+t1394*a51+t1406*a53+t1418*a54;
		*tmp_SM_ptr++ = t1530*a7+t1542*a16+t1554*a24+t1566*a31+t1578*a37+t1590*a42+t1602*a46+t1607*a49+t1619*a50+t1631*a51;
		*tmp_SM_ptr++ = t1530*a8+t1542*a17+t1554*a25+t1566*a32+t1578*a38+t1590*a43+t1602*a47+t1607*a50+t1619*a52+t1631*a53;
		*tmp_SM_ptr++ = t1530*a9+t1542*a18+t1554*a26+t1566*a33+t1578*a39+t1590*a44+t1602*a48+t1607*a51+t1619*a53+t1631*a54;
		*tmp_SM_ptr++ = t1743*a8+t1755*a17+t1767*a25+t1779*a32+t1791*a38+t1803*a43+t1815*a47+t1827*a50+t1831*a52+t1843*a53;
		*tmp_SM_ptr++ = t1743*a9+t1755*a18+t1767*a26+t1779*a33+t1791*a39+t1803*a44+t1815*a48+t1827*a51+t1831*a53+t1843*a54;
		*tmp_SM_ptr++ = t1955*a9+t1967*a18+t1979*a26+t1991*a33+t2003*a39+t2015*a44+t2027*a48+t2039*a51+t2051*a53+t2054*a54;

	}
	else {
		CvMat* A = cvCreateMat(M_size, M_size, CV_64F);
		CvMat* B = cvCreateMat(M_size, M_size, CV_64F);
		CvMat* AB = cvCreateMat(M_size, M_size, CV_64F);
		CvMat* ABA = cvCreateMat(M_size, M_size, CV_64F);

		double *this_M_ptr = this->GetData(0);
		double *Y_mat_M_ptr = Y_mat->GetData(0);

		memcpy(A->data.ptr, this_M_ptr, sizeof(double)*M_size2);
		memcpy(B->data.ptr, Y_mat_M_ptr, sizeof(double)*M_size2);

		delete [] this_M_ptr;
		delete [] Y_mat_M_ptr;

		cvMatMul(A, B, AB);
		cvMatMul(AB, A, ABA);

		if ( mul2_mat ) {
			mul2_mat->SetData((double*)ABA->data.ptr, M_size2);
		}
		else {
			this->SetData((double*)ABA->data.ptr, M_size2);
		}

		cvReleaseMat(&A);
		cvReleaseMat(&B);
		cvReleaseMat(&AB);
		cvReleaseMat(&ABA);
	}
}


void CCovarianceMatrix::Log_Y(CCovarianceMatrix *Y_mat, CCovarianceMatrix *log_Y_mat, CCovarianceMatrix *this_mat_sqrt, CCovarianceMatrix *this_mat_inv_sqrt)
{
	CCovarianceMatrix *log_inv_sqrt_X = new CCovarianceMatrix(M_size);
	CCovarianceMatrix *sqrt_X = this_mat_sqrt;
	CCovarianceMatrix *inv_sqrt_X = this_mat_inv_sqrt;

	if ( !sqrt_X ) {
		sqrt_X = new CCovarianceMatrix(M_size);
		inv_sqrt_X = new CCovarianceMatrix(M_size);
		this->Sqrt(sqrt_X);
		sqrt_X->Inverse(inv_sqrt_X);
	}

	inv_sqrt_X->Mul2(Y_mat, log_inv_sqrt_X);
	log_inv_sqrt_X->Log();

	if ( log_Y_mat ) {
		sqrt_X->Mul2(log_inv_sqrt_X, log_Y_mat);
	}
	else {
		sqrt_X->Mul2(log_inv_sqrt_X, this);
	}

	delete log_inv_sqrt_X;
	if ( !this_mat_sqrt ) {
		delete sqrt_X;
		delete inv_sqrt_X;
	}
}

void CCovarianceMatrix::Exp_Y(CCovarianceMatrix *Y_mat, CCovarianceMatrix *exp_Y_mat, CCovarianceMatrix *this_mat_sqrt, CCovarianceMatrix *this_mat_inv_sqrt)
{
	CCovarianceMatrix *exp_inv_sqrt_X = new CCovarianceMatrix(M_size);
	CCovarianceMatrix *sqrt_X = this_mat_sqrt;
	CCovarianceMatrix *inv_sqrt_X = this_mat_inv_sqrt;

	if ( !sqrt_X ) {
		sqrt_X = new CCovarianceMatrix(M_size);
		inv_sqrt_X = new CCovarianceMatrix(M_size);
		this->Sqrt(sqrt_X);
		sqrt_X->Inverse(inv_sqrt_X);
	}

	inv_sqrt_X->Mul2(Y_mat, exp_inv_sqrt_X);
	exp_inv_sqrt_X->Exp();

	if ( exp_Y_mat ) {
		sqrt_X->Mul2(exp_inv_sqrt_X, exp_Y_mat);
	}
	else {
		sqrt_X->Mul2(exp_inv_sqrt_X, this);
	}

	delete exp_inv_sqrt_X;
	if ( !this_mat_sqrt ) {
		delete sqrt_X;
		delete inv_sqrt_X;
	}
}

void CCovarianceMatrix::Sqrt(CCovarianceMatrix *sqrt_mat)
{
	CCovarianceMatrix *dst_mat;

	if ( sqrt_mat )
		dst_mat = sqrt_mat;
	else
		dst_mat = this;

	if ( this->M_size == 1 ) {
		if ( SM_ptr[0] <= 0 ) {
			printf("SM_ptr[0] = %.20f - sqrt error in sqrt(CCovarianceMatrix *sqrt_mat)!\n", SM_ptr[0]);
			SM_ptr[0] = MIN_DOUBLE_NUM;
		}
		dst_mat->SM_ptr[0] = sqrt(SM_ptr[0]);
	}
	else if ( this->M_size == 2 ) {
		double D[2];
		double U[4];

		eigen2_decomposition(this->SM_ptr, U, D);

		double a0,a1,a2,a3;
		double s0,s1;

		a0 = U[0];
		a1 = U[1];
		a2 = U[2];
		a3 = U[3];

		if ( D[0] <= 0 || D[1] <= 0 ) {
			printf("SM_ptr[0] = %.20f\tSM_ptr[1] = %.20f - sqrt error in sqrt(CCovarianceMatrix *sqrt_mat)!\n", D[0], D[1]);
			D[0] = MAX(D[0], MIN_DOUBLE_NUM);
			D[1] = MAX(D[1], MIN_DOUBLE_NUM);
		}

		s0 = sqrt(D[0]);
		s1 = sqrt(D[1]);

		double t1 = a0*a0;
		double t3 = a2*a2;
		double t10 = a0*s0*a1+a2*s1*a3;
		double t11 = a1*a1;
		double t13 = a3*a3;

		double* tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1;
		*tmp_M_ptr++ = t10;
		*tmp_M_ptr++ = t11*s0+t13*s1;
	}
	else if ( this->M_size == 3 ) {
		double D[3];
		double U[9];

		eigen3_decomposition(this->SM_ptr, U, D);

		double a0,a1,a2,a3,a4,a5,a6,a7,a8;
		double s0,s1,s2;

		a0 = U[0];
		a1 = U[1];
		a2 = U[2];
		a3 = U[3];
		a4 = U[4];
		a5 = U[5];
		a6 = U[6];
		a7 = U[7];
		a8 = U[8];

		if ( D[0] <= 0 || D[1] <= 0 || D[2] <= 0 ) {
			printf("SM_ptr[0] = %.20f\tSM_ptr[1] = %.20f\tSM_ptr[2] = %.20f - sqrt error in sqrt(CCovarianceMatrix *sqrt_mat)!\n", D[0], D[1], D[2]);
			D[0] = MAX(D[0], MIN_DOUBLE_NUM);
			D[1] = MAX(D[1], MIN_DOUBLE_NUM);
			D[2] = MAX(D[2], MIN_DOUBLE_NUM);
		}

		s0 = sqrt(D[0]);
		s1 = sqrt(D[1]);
		s2 = sqrt(D[2]);

		double t1 = a0*a0;
		double t3 = a3*a3;
		double t5 = a6*a6;
		double t8 = a0*s0;
		double t10 = a3*s1;
		double t12 = a6*s2;
		double t14 = t8*a1+t10*a4+t12*a7;
		double t18 = t8*a2+t10*a5+t12*a8;
		double t19 = a1*a1;
		double t21 = a4*a4;
		double t23 = a7*a7;
		double t32 = a1*s0*a2+a4*s1*a5+a7*s2*a8;
		double t33 = a2*a2;
		double t35 = a5*a5;
		double t37 = a8*a8;

		double *tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2;
		*tmp_M_ptr++ = t14;
		*tmp_M_ptr++ = t18;
		*tmp_M_ptr++ = t19*s0+t21*s1+t23*s2;
		*tmp_M_ptr++ = t32;
		*tmp_M_ptr++ = t33*s0+t35*s1+t37*s2;
	}
	else if ( this->M_size == 4 ) {
		double D[4];
		double U[16];

		eigen4_decomposition(this->SM_ptr, U, D);

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
		double s0,s1,s2,s3;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = sqrt(*tmp_M_ptr++);
		s1 = sqrt(*tmp_M_ptr++);
		s2 = sqrt(*tmp_M_ptr++);
		s3 = sqrt(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a4*a4;
		double t5 = a8*a8;
		double t7 = a12*a12;
		double t10 = a0*s0;
		double t12 = a4*s1;
		double t14 = a8*s2;
		double t16 = a12*s3;
		double t18 = t10*a1+t12*a5+t14*a9+t16*a13;
		double t23 = t10*a2+t12*a6+t14*a10+t16*a14;
		double t28 = t10*a3+t12*a7+t14*a11+t16*a15;
		double t29 = a1*a1;
		double t31 = a5*a5;
		double t33 = a9*a9;
		double t35 = a13*a13;
		double t38 = a1*s0;
		double t40 = a5*s1;
		double t42 = a9*s2;
		double t44 = a13*s3;
		double t46 = t38*a2+t40*a6+t42*a10+t44*a14;
		double t51 = t38*a3+t40*a7+t42*a11+t44*a15;
		double t52 = a2*a2;
		double t54 = a6*a6;
		double t56 = a10*a10;
		double t58 = a14*a14;
		double t69 = a2*s0*a3+a6*s1*a7+a10*s2*a11+a14*s3*a15;
		double t70 = a3*a3;
		double t72 = a7*a7;
		double t74 = a11*a11;
		double t76 = a15*a15;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3;
		*tmp_M_ptr++ = t18;
		*tmp_M_ptr++ = t23;
		*tmp_M_ptr++ = t28;
		*tmp_M_ptr++ = t29*s0+t31*s1+t33*s2+t35*s3;
		*tmp_M_ptr++ = t46;
		*tmp_M_ptr++ = t51;
		*tmp_M_ptr++ = t52*s0+t54*s1+t56*s2+t58*s3;
		*tmp_M_ptr++ = t69;
		*tmp_M_ptr++ = t70*s0+t72*s1+t74*s2+t76*s3;
	}
	else if ( this->M_size == 5 ) {
		double D[5];
		double U[25];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 5);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24;
		double s0,s1,s2,s3,s4;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = sqrt(*tmp_M_ptr++);
		s1 = sqrt(*tmp_M_ptr++);
		s2 = sqrt(*tmp_M_ptr++);
		s3 = sqrt(*tmp_M_ptr++);
		s4 = sqrt(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a5*a5;
		double t5 = a10*a10;
		double t7 = a15*a15;
		double t9 = a20*a20;
		double t12 = a0*s0;
		double t14 = a5*s1;
		double t16 = a10*s2;
		double t18 = a15*s3;
		double t20 = a20*s4;
		double t22 = t12*a1+t14*a6+t16*a11+t18*a16+t20*a21;
		double t28 = t12*a2+t14*a7+t16*a12+t18*a17+t20*a22;
		double t34 = t12*a3+t14*a8+t16*a13+t18*a18+t20*a23;
		double t40 = t12*a4+t14*a9+t16*a14+t18*a19+t20*a24;
		double t41 = a1*a1;
		double t43 = a6*a6;
		double t45 = a11*a11;
		double t47 = a16*a16;
		double t49 = a21*a21;
		double t52 = a1*s0;
		double t54 = a6*s1;
		double t56 = a11*s2;
		double t58 = a16*s3;
		double t60 = a21*s4;
		double t62 = t52*a2+t54*a7+t56*a12+t58*a17+t60*a22;
		double t68 = t52*a3+t54*a8+t56*a13+t58*a18+t60*a23;
		double t74 = t52*a4+t54*a9+t56*a14+t58*a19+t60*a24;
		double t75 = a2*a2;
		double t77 = a7*a7;
		double t79 = a12*a12;
		double t81 = a17*a17;
		double t83 = a22*a22;
		double t86 = a2*s0;
		double t88 = a7*s1;
		double t90 = a12*s2;
		double t92 = a17*s3;
		double t94 = a22*s4;
		double t96 = t86*a3+t88*a8+t90*a13+t92*a18+t94*a23;
		double t102 = t86*a4+t88*a9+t90*a14+t92*a19+t94*a24;
		double t103 = a3*a3;
		double t105 = a8*a8;
		double t107 = a13*a13;
		double t109 = a18*a18;
		double t111 = a23*a23;
		double t124 = a3*s0*a4+a8*s1*a9+a13*s2*a14+a18*s3*a19+a23*s4*a24;
		double t125 = a4*a4;
		double t127 = a9*a9;
		double t129 = a14*a14;
		double t131 = a19*a19;
		double t133 = a24*a24;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4;
		*tmp_M_ptr++  = t22;
		*tmp_M_ptr++  = t28;
		*tmp_M_ptr++  = t34;
		*tmp_M_ptr++  = t40;
		*tmp_M_ptr++  = t41*s0+t43*s1+t45*s2+t47*s3+t49*s4;
		*tmp_M_ptr++  = t62;
		*tmp_M_ptr++  = t68;
		*tmp_M_ptr++  = t74;
		*tmp_M_ptr++  = t75*s0+t77*s1+t79*s2+t81*s3+t83*s4;
		*tmp_M_ptr++  = t96;
		*tmp_M_ptr++  = t102;
		*tmp_M_ptr++  = t103*s0+t105*s1+t107*s2+t109*s3+t111*s4;
		*tmp_M_ptr++  = t124;
		*tmp_M_ptr++  = t125*s0+t127*s1+t129*s2+t131*s3+t133*s4;
	}
	else if ( this->M_size == 6 ) {
		double D[6];
		double U[36];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 6);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35;
		double s0,s1,s2,s3,s4,s5;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = sqrt(*tmp_M_ptr++);
		s1 = sqrt(*tmp_M_ptr++);
		s2 = sqrt(*tmp_M_ptr++);
		s3 = sqrt(*tmp_M_ptr++);
		s4 = sqrt(*tmp_M_ptr++);
		s5 = sqrt(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a6*a6;
		double t5 = a12*a12;
		double t7 = a18*a18;
		double t9 = a24*a24;
		double t11 = a30*a30;
		double t14 = a0*s0;
		double t16 = a6*s1;
		double t18 = a12*s2;
		double t20 = a18*s3;
		double t22 = a24*s4;
		double t24 = a30*s5;
		double t26 = t14*a1+t16*a7+t18*a13+t20*a19+t22*a25+t24*a31;
		double t33 = t14*a2+t16*a8+t18*a14+t20*a20+t22*a26+t24*a32;
		double t40 = t14*a3+t16*a9+t18*a15+t20*a21+t22*a27+t24*a33;
		double t47 = t14*a4+10.0*t16+t18*a16+t20*a22+t22*a28+t24*a34;
		double t54 = t14*a5+t16*a11+t18*a17+t20*a23+t22*a29+t24*a35;
		double t55 = a1*a1;
		double t57 = a7*a7;
		double t59 = a13*a13;
		double t61 = a19*a19;
		double t63 = a25*a25;
		double t65 = a31*a31;
		double t68 = a1*s0;
		double t70 = a7*s1;
		double t72 = a13*s2;
		double t74 = a19*s3;
		double t76 = a25*s4;
		double t78 = a31*s5;
		double t80 = t68*a2+t70*a8+t72*a14+t74*a20+t76*a26+t78*a32;
		double t87 = t68*a3+t70*a9+t72*a15+t74*a21+t76*a27+t78*a33;
		double t94 = t68*a4+10.0*t70+t72*a16+t74*a22+t76*a28+t78*a34;
		double t101 = t68*a5+t70*a11+t72*a17+t74*a23+t76*a29+t78*a35;
		double t102 = a2*a2;
		double t104 = a8*a8;
		double t106 = a14*a14;
		double t108 = a20*a20;
		double t110 = a26*a26;
		double t112 = a32*a32;
		double t115 = a2*s0;
		double t117 = a8*s1;
		double t119 = a14*s2;
		double t121 = a20*s3;
		double t123 = a26*s4;
		double t125 = a32*s5;
		double t127 = t115*a3+t117*a9+t119*a15+t121*a21+t123*a27+t125*a33;
		double t134 = t115*a4+10.0*t117+t119*a16+t121*a22+t123*a28+t125*a34;
		double t141 = t115*a5+t117*a11+t119*a17+t121*a23+t123*a29+t125*a35;
		double t142 = a3*a3;
		double t144 = a9*a9;
		double t146 = a15*a15;
		double t148 = a21*a21;
		double t150 = a27*a27;
		double t152 = a33*a33;
		double t155 = a3*s0;
		double t157 = a9*s1;
		double t159 = a15*s2;
		double t161 = a21*s3;
		double t163 = a27*s4;
		double t165 = a33*s5;
		double t167 = t155*a4+10.0*t157+t159*a16+t161*a22+t163*a28+t165*a34;
		double t174 = t155*a5+t157*a11+t159*a17+t161*a23+t163*a29+t165*a35;
		double t175 = a4*a4;
		double t178 = a16*a16;
		double t180 = a22*a22;
		double t182 = a28*a28;
		double t184 = a34*a34;
		double t199 = a4*s0*a5+10.0*a11*s1+a16*s2*a17+a22*s3*a23+a28*s4*a29+a34*s5*a35;
		double t200 = a5*a5;
		double t202 = a11*a11;
		double t204 = a17*a17;
		double t206 = a23*a23;
		double t208 = a29*a29;
		double t210 = a35*a35;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5;
		*tmp_M_ptr++ = t26;
		*tmp_M_ptr++ = t33;
		*tmp_M_ptr++ = t40;
		*tmp_M_ptr++ = t47;
		*tmp_M_ptr++ = t54;
		*tmp_M_ptr++ = t55*s0+t57*s1+t59*s2+t61*s3+t63*s4+t65*s5;
		*tmp_M_ptr++ = t80;
		*tmp_M_ptr++ = t87;
		*tmp_M_ptr++ = t94;
		*tmp_M_ptr++ = t101;
		*tmp_M_ptr++ = t102*s0+t104*s1+t106*s2+t108*s3+t110*s4+t112*s5;
		*tmp_M_ptr++ = t127;
		*tmp_M_ptr++ = t134;
		*tmp_M_ptr++ = t141;
		*tmp_M_ptr++ = t142*s0+t144*s1+t146*s2+t148*s3+t150*s4+t152*s5;
		*tmp_M_ptr++ = t167;
		*tmp_M_ptr++ = t174;
		*tmp_M_ptr++ = t175*s0+100.0*s1+t178*s2+t180*s3+t182*s4+t184*s5;
		*tmp_M_ptr++ = t199;
		*tmp_M_ptr++ = t200*s0+t202*s1+t204*s2+t206*s3+t208*s4+t210*s5;
	}
	else if ( this->M_size == 7 ) {
		double D[7];
		double U[49];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 7);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double s0,s1,s2,s3,s4,s5,s6;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = sqrt(*tmp_M_ptr++);
		s1 = sqrt(*tmp_M_ptr++);
		s2 = sqrt(*tmp_M_ptr++);
		s3 = sqrt(*tmp_M_ptr++);
		s4 = sqrt(*tmp_M_ptr++);
		s5 = sqrt(*tmp_M_ptr++);
		s6 = sqrt(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a7*a7;
		double t5 = a14*a14;
		double t7 = a21*a21;
		double t9 = a28*a28;
		double t11 = a35*a35;
		double t13 = a42*a42;
		double t16 = a0*s0;
		double t18 = a7*s1;
		double t20 = a14*s2;
		double t22 = a21*s3;
		double t24 = a28*s4;
		double t26 = a35*s5;
		double t28 = a42*s6;
		double t30 = t16*a1+t18*a8+t20*a15+t22*a22+t24*a29+t26*a36+t28*a43;
		double t38 = t16*a2+t18*a9+t20*a16+t22*a23+t24*a30+t26*a37+t28*a44;
		double t46 = t16*a3+t18*a10+t20*a17+t22*a24+t24*a31+t26*a38+t28*a45;
		double t54 = t16*a4+t18*a11+t20*a18+t22*a25+t24*a32+t26*a39+t28*a46;
		double t62 = t16*a5+t18*a12+t20*a19+t22*a26+t24*a33+t26*a40+t28*a47;
		double t70 = t16*a6+t18*a13+t20*a20+t22*a27+t24*a34+t26*a41+t28*a48;
		double t71 = a1*a1;
		double t73 = a8*a8;
		double t75 = a15*a15;
		double t77 = a22*a22;
		double t79 = a29*a29;
		double t81 = a36*a36;
		double t83 = a43*a43;
		double t86 = a1*s0;
		double t88 = a8*s1;
		double t90 = a15*s2;
		double t92 = a22*s3;
		double t94 = a29*s4;
		double t96 = a36*s5;
		double t98 = a43*s6;
		double t100 = t86*a2+t88*a9+t90*a16+t92*a23+t94*a30+t96*a37+t98*a44;
		double t108 = t86*a3+t88*a10+t90*a17+t92*a24+t94*a31+t96*a38+t98*a45;
		double t116 = t86*a4+t88*a11+t90*a18+t92*a25+t94*a32+t96*a39+t98*a46;
		double t124 = t86*a5+t88*a12+t90*a19+t92*a26+t94*a33+t96*a40+t98*a47;
		double t132 = t86*a6+t88*a13+t90*a20+t92*a27+t94*a34+t96*a41+t98*a48;
		double t133 = a2*a2;
		double t135 = a9*a9;
		double t137 = a16*a16;
		double t139 = a23*a23;
		double t141 = a30*a30;
		double t143 = a37*a37;
		double t145 = a44*a44;
		double t148 = a2*s0;
		double t150 = a9*s1;
		double t152 = a16*s2;
		double t154 = a23*s3;
		double t156 = a30*s4;
		double t158 = a37*s5;
		double t160 = a44*s6;
		double t162 = t148*a3+t150*a10+t152*a17+t154*a24+t156*a31+t158*a38+t160*a45;
		double t170 = t148*a4+t150*a11+t152*a18+t154*a25+t156*a32+t158*a39+t160*a46;
		double t178 = t148*a5+t150*a12+t152*a19+t154*a26+t156*a33+t158*a40+t160*a47;
		double t186 = t148*a6+t150*a13+t152*a20+t154*a27+t156*a34+t158*a41+t160*a48;
		double t187 = a3*a3;
		double t189 = a10*a10;
		double t191 = a17*a17;
		double t193 = a24*a24;
		double t195 = a31*a31;
		double t197 = a38*a38;
		double t199 = a45*a45;
		double t202 = a3*s0;
		double t204 = a10*s1;
		double t206 = a17*s2;
		double t208 = a24*s3;
		double t210 = a31*s4;
		double t212 = a38*s5;
		double t214 = a45*s6;
		double t216 = t202*a4+t204*a11+t206*a18+t208*a25+t210*a32+t212*a39+t214*a46;
		double t224 = t202*a5+t204*a12+t206*a19+t208*a26+t210*a33+t212*a40+t214*a47;
		double t232 = t202*a6+t204*a13+t206*a20+t208*a27+t210*a34+t212*a41+t214*a48;
		double t233 = a4*a4;
		double t235 = a11*a11;
		double t237 = a18*a18;
		double t239 = a25*a25;
		double t241 = a32*a32;
		double t243 = a39*a39;
		double t245 = a46*a46;
		double t248 = a4*s0;
		double t250 = a11*s1;
		double t252 = a18*s2;
		double t254 = a25*s3;
		double t256 = a32*s4;
		double t258 = a39*s5;
		double t260 = a46*s6;
		double t262 = t248*a5+t250*a12+t252*a19+t254*a26+t256*a33+t258*a40+t260*a47;
		double t270 = t248*a6+t250*a13+t252*a20+t254*a27+t256*a34+t258*a41+t260*a48;
		double t271 = a5*a5;
		double t273 = a12*a12;
		double t275 = a19*a19;
		double t277 = a26*a26;
		double t279 = a33*a33;
		double t281 = a40*a40;
		double t283 = a47*a47;
		double t300 = a5*s0*a6+a12*s1*a13+a19*s2*a20+a26*s3*a27+a33*s4*a34+a40*s5*a41+a47*s6*a48;
		double t301 = a6*a6;
		double t303 = a13*a13;
		double t305 = a20*a20;
		double t307 = a27*a27;
		double t309 = a34*a34;
		double t311 = a41*a41;
		double t313 = a48*a48;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6;
		*tmp_M_ptr++ = t30;
		*tmp_M_ptr++ = t38;
		*tmp_M_ptr++ = t46;
		*tmp_M_ptr++ = t54;
		*tmp_M_ptr++ = t62;
		*tmp_M_ptr++ = t70;
		*tmp_M_ptr++ = t71*s0+t73*s1+t75*s2+t77*s3+t79*s4+t81*s5+t83*s6;
		*tmp_M_ptr++ = t100;
		*tmp_M_ptr++ = t108;
		*tmp_M_ptr++ = t116;
		*tmp_M_ptr++ = t124;
		*tmp_M_ptr++ = t132;
		*tmp_M_ptr++ = t133*s0+t135*s1+t137*s2+t139*s3+t141*s4+t143*s5+t145*s6;
		*tmp_M_ptr++ = t162;
		*tmp_M_ptr++ = t170;
		*tmp_M_ptr++ = t178;
		*tmp_M_ptr++ = t186;
		*tmp_M_ptr++ = t187*s0+t189*s1+t191*s2+t193*s3+t195*s4+t197*s5+t199*s6;
		*tmp_M_ptr++ = t216;
		*tmp_M_ptr++ = t224;
		*tmp_M_ptr++ = t232;
		*tmp_M_ptr++ = t233*s0+t235*s1+t237*s2+t239*s3+t241*s4+t243*s5+t245*s6;
		*tmp_M_ptr++ = t262;
		*tmp_M_ptr++ = t270;
		*tmp_M_ptr++ = t271*s0+t273*s1+t275*s2+t277*s3+t279*s4+t281*s5+t283*s6;
		*tmp_M_ptr++ = t300;
		*tmp_M_ptr++ = t301*s0+t303*s1+t305*s2+t307*s3+t309*s4+t311*s5+t313*s6;

	}
	else if ( this->M_size == 8 ) {
		double D[8];
		double U[64];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 8);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63;
		double s0,s1,s2,s3,s4,s5,s6,s7;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;
		a49 = *tmp_M_ptr++;
		a50 = *tmp_M_ptr++;
		a51 = *tmp_M_ptr++;
		a52 = *tmp_M_ptr++;
		a53 = *tmp_M_ptr++;
		a54 = *tmp_M_ptr++;
		a55 = *tmp_M_ptr++;
		a56 = *tmp_M_ptr++;
		a57 = *tmp_M_ptr++;
		a58 = *tmp_M_ptr++;
		a59 = *tmp_M_ptr++;
		a60 = *tmp_M_ptr++;
		a61 = *tmp_M_ptr++;
		a62 = *tmp_M_ptr++;
		a63 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = sqrt(*tmp_M_ptr++);
		s1 = sqrt(*tmp_M_ptr++);
		s2 = sqrt(*tmp_M_ptr++);
		s3 = sqrt(*tmp_M_ptr++);
		s4 = sqrt(*tmp_M_ptr++);
		s5 = sqrt(*tmp_M_ptr++);
		s6 = sqrt(*tmp_M_ptr++);
		s7 = sqrt(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a8*a8;
		double t5 = a16*a16;
		double t7 = a24*a24;
		double t9 = a32*a32;
		double t11 = a40*a40;
		double t13 = a48*a48;
		double t15 = a56*a56;
		double t18 = a0*s0;
		double t20 = a8*s1;
		double t22 = a16*s2;
		double t24 = a24*s3;
		double t26 = a32*s4;
		double t28 = a40*s5;
		double t30 = a48*s6;
		double t32 = a56*s7;
		double t34 = t18*a1+t20*a9+t22*a17+t24*a25+t26*a33+t28*a41+t30*a49+t32*a57;
		double t43 = t18*a2+t20*a10+t22*a18+t24*a26+t26*a34+t28*a42+t30*a50+t32*a58;
		double t52 = t18*a3+t20*a11+t22*a19+t24*a27+t26*a35+t28*a43+t30*a51+t32*a59;
		double t61 = t18*a4+t20*a12+t22*a20+t24*a28+t26*a36+t28*a44+t30*a52+t32*a60;
		double t70 = t18*a5+t20*a13+t22*a21+t24*a29+t26*a37+t28*a45+t30*a53+t32*a61;
		double t79 = t18*a6+t20*a14+t22*a22+t24*a30+t26*a38+t28*a46+t30*a54+t32*a62;
		double t88 = t18*a7+t20*a15+t22*a23+t24*a31+t26*a39+t28*a47+t30*a55+t32*a63;
		double t89 = a1*a1;
		double t91 = a9*a9;
		double t93 = a17*a17;
		double t95 = a25*a25;
		double t97 = a33*a33;
		double t99 = a41*a41;
		double t101 = a49*a49;
		double t103 = a57*a57;
		double t106 = a1*s0;
		double t108 = a9*s1;
		double t110 = a17*s2;
		double t112 = a25*s3;
		double t114 = a33*s4;
		double t116 = a41*s5;
		double t118 = a49*s6;
		double t120 = a57*s7;
		double t122 = t106*a2+t108*a10+t110*a18+t112*a26+t114*a34+t116*a42+t118*a50+t120*a58;
		double t131 = t106*a3+t108*a11+t110*a19+t112*a27+t114*a35+t116*a43+t118*a51+t120*a59;
		double t140 = t106*a4+t108*a12+t110*a20+t112*a28+t114*a36+t116*a44+t118*a52+t120*a60;
		double t149 = t106*a5+t108*a13+t110*a21+t112*a29+t114*a37+t116*a45+t118*a53+t120*a61;
		double t158 = t106*a6+t108*a14+t110*a22+t112*a30+t114*a38+t116*a46+t118*a54+t120*a62;
		double t167 = t106*a7+t108*a15+t110*a23+t112*a31+t114*a39+t116*a47+t118*a55+t120*a63;
		double t168 = a2*a2;
		double t170 = a10*a10;
		double t172 = a18*a18;
		double t174 = a26*a26;
		double t176 = a34*a34;
		double t178 = a42*a42;
		double t180 = a50*a50;
		double t182 = a58*a58;
		double t185 = a2*s0;
		double t187 = a10*s1;
		double t189 = a18*s2;
		double t191 = a26*s3;
		double t193 = a34*s4;
		double t195 = a42*s5;
		double t197 = a50*s6;
		double t199 = a58*s7;
		double t201 = t185*a3+t187*a11+t189*a19+t191*a27+t193*a35+t195*a43+t197*a51+t199*a59;
		double t210 = t185*a4+t187*a12+t189*a20+t191*a28+t193*a36+t195*a44+t197*a52+t199*a60;
		double t219 = t185*a5+t187*a13+t189*a21+t191*a29+t193*a37+t195*a45+t197*a53+t199*a61;
		double t228 = t185*a6+t187*a14+t189*a22+t191*a30+t193*a38+t195*a46+t197*a54+t199*a62;
		double t237 = t185*a7+t187*a15+t189*a23+t191*a31+t193*a39+t195*a47+t197*a55+t199*a63;
		double t238 = a3*a3;
		double t240 = a11*a11;
		double t242 = a19*a19;
		double t244 = a27*a27;
		double t246 = a35*a35;
		double t248 = a43*a43;
		double t250 = a51*a51;
		double t252 = a59*a59;
		double t255 = a3*s0;
		double t257 = a11*s1;
		double t259 = a19*s2;
		double t261 = a27*s3;
		double t263 = a35*s4;
		double t265 = a43*s5;
		double t267 = a51*s6;
		double t269 = a59*s7;
		double t271 = t255*a4+t257*a12+t259*a20+t261*a28+t263*a36+t265*a44+t267*a52+t269*a60;
		double t280 = t255*a5+t257*a13+t259*a21+t261*a29+t263*a37+t265*a45+t267*a53+t269*a61;
		double t289 = t255*a6+t257*a14+t259*a22+t261*a30+t263*a38+t265*a46+t267*a54+t269*a62;
		double t298 = t255*a7+t257*a15+t259*a23+t261*a31+t263*a39+t265*a47+t267*a55+t269*a63;
		double t299 = a4*a4;
		double t301 = a12*a12;
		double t303 = a20*a20;
		double t305 = a28*a28;
		double t307 = a36*a36;
		double t309 = a44*a44;
		double t311 = a52*a52;
		double t313 = a60*a60;
		double t316 = a4*s0;
		double t318 = a12*s1;
		double t320 = a20*s2;
		double t322 = a28*s3;
		double t324 = a36*s4;
		double t326 = a44*s5;
		double t328 = a52*s6;
		double t330 = a60*s7;
		double t332 = t316*a5+t318*a13+t320*a21+t322*a29+t324*a37+t326*a45+t328*a53+t330*a61;
		double t341 = t316*a6+t318*a14+t320*a22+t322*a30+t324*a38+t326*a46+t328*a54+t330*a62;
		double t350 = t316*a7+t318*a15+t320*a23+t322*a31+t324*a39+t326*a47+t328*a55+t330*a63;
		double t351 = a5*a5;
		double t353 = a13*a13;
		double t355 = a21*a21;
		double t357 = a29*a29;
		double t359 = a37*a37;
		double t361 = a45*a45;
		double t363 = a53*a53;
		double t365 = a61*a61;
		double t368 = a5*s0;
		double t370 = a13*s1;
		double t372 = a21*s2;
		double t374 = a29*s3;
		double t376 = a37*s4;
		double t378 = a45*s5;
		double t380 = a53*s6;
		double t382 = a61*s7;
		double t384 = t368*a6+t370*a14+t372*a22+t374*a30+t376*a38+t378*a46+t380*a54+t382*a62;
		double t393 = t368*a7+t370*a15+t372*a23+t374*a31+t376*a39+t378*a47+t380*a55+t382*a63;
		double t394 = a6*a6;
		double t396 = a14*a14;
		double t398 = a22*a22;
		double t400 = a30*a30;
		double t402 = a38*a38;
		double t404 = a46*a46;
		double t406 = a54*a54;
		double t408 = a62*a62;
		double t427 = a6*s0*a7+a14*s1*a15+a22*s2*a23+a30*s3*a31+a38*s4*a39+a46*s5*a47+a54*s6*a55+a62*s7*a63;
		double t428 = a7*a7;
		double t430 = a15*a15;
		double t432 = a23*a23;
		double t434 = a31*a31;
		double t436 = a39*a39;
		double t438 = a47*a47;
		double t440 = a55*a55;
		double t442 = a63*a63;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6+t15*s7;
		*tmp_M_ptr++ = t34;
		*tmp_M_ptr++ = t43;
		*tmp_M_ptr++ = t52;
		*tmp_M_ptr++ = t61;
		*tmp_M_ptr++ = t70;
		*tmp_M_ptr++ = t79;
		*tmp_M_ptr++ = t88;
		*tmp_M_ptr++ = t89*s0+t91*s1+t93*s2+t95*s3+t97*s4+t99*s5+t101*s6+t103*s7;
		*tmp_M_ptr++ = t122;
		*tmp_M_ptr++ = t131;
		*tmp_M_ptr++ = t140;
		*tmp_M_ptr++ = t149;
		*tmp_M_ptr++ = t158;
		*tmp_M_ptr++ = t167;
		*tmp_M_ptr++ = t168*s0+t170*s1+t172*s2+t174*s3+t176*s4+t178*s5+t180*s6+t182*s7;
		*tmp_M_ptr++ = t201;
		*tmp_M_ptr++ = t210;
		*tmp_M_ptr++ = t219;
		*tmp_M_ptr++ = t228;
		*tmp_M_ptr++ = t237;
		*tmp_M_ptr++ = t238*s0+t240*s1+t242*s2+t244*s3+t246*s4+t248*s5+t250*s6+t252*s7;
		*tmp_M_ptr++ = t271;
		*tmp_M_ptr++ = t280;
		*tmp_M_ptr++ = t289;
		*tmp_M_ptr++ = t298;
		*tmp_M_ptr++ = t299*s0+t301*s1+t303*s2+t305*s3+t307*s4+t309*s5+t311*s6+t313*s7;
		*tmp_M_ptr++ = t332;
		*tmp_M_ptr++ = t341;
		*tmp_M_ptr++ = t350;
		*tmp_M_ptr++ = t351*s0+t353*s1+t355*s2+t357*s3+t359*s4+t361*s5+t363*s6+t365*s7;
		*tmp_M_ptr++ = t384;
		*tmp_M_ptr++ = t393;
		*tmp_M_ptr++ = t394*s0+t396*s1+t398*s2+t400*s3+t402*s4+t404*s5+t406*s6+t408*s7;
		*tmp_M_ptr++ = t427;
		*tmp_M_ptr++ = t428*s0+t430*s1+t432*s2+t434*s3+t436*s4+t438*s5+t440*s6+t442*s7;
	}
	else if ( this->M_size == 9 ) {
		double D[9];
		double U[81];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 9);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80;
		double s0,s1,s2,s3,s4,s5,s6,s7,s8;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;
		a49 = *tmp_M_ptr++;
		a50 = *tmp_M_ptr++;
		a51 = *tmp_M_ptr++;
		a52 = *tmp_M_ptr++;
		a53 = *tmp_M_ptr++;
		a54 = *tmp_M_ptr++;
		a55 = *tmp_M_ptr++;
		a56 = *tmp_M_ptr++;
		a57 = *tmp_M_ptr++;
		a58 = *tmp_M_ptr++;
		a59 = *tmp_M_ptr++;
		a60 = *tmp_M_ptr++;
		a61 = *tmp_M_ptr++;
		a62 = *tmp_M_ptr++;
		a63 = *tmp_M_ptr++;
		a64 = *tmp_M_ptr++;
		a65 = *tmp_M_ptr++;
		a66 = *tmp_M_ptr++;
		a67 = *tmp_M_ptr++;
		a68 = *tmp_M_ptr++;
		a69 = *tmp_M_ptr++;
		a70 = *tmp_M_ptr++;
		a71 = *tmp_M_ptr++;
		a72 = *tmp_M_ptr++;
		a73 = *tmp_M_ptr++;
		a74 = *tmp_M_ptr++;
		a75 = *tmp_M_ptr++;
		a76 = *tmp_M_ptr++;
		a77 = *tmp_M_ptr++;
		a78 = *tmp_M_ptr++;
		a79 = *tmp_M_ptr++;
		a80 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = sqrt(*tmp_M_ptr++);
		s1 = sqrt(*tmp_M_ptr++);
		s2 = sqrt(*tmp_M_ptr++);
		s3 = sqrt(*tmp_M_ptr++);
		s4 = sqrt(*tmp_M_ptr++);
		s5 = sqrt(*tmp_M_ptr++);
		s6 = sqrt(*tmp_M_ptr++);
		s7 = sqrt(*tmp_M_ptr++);
		s8 = sqrt(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a9*a9;
		double t5 = a18*a18;
		double t7 = a27*a27;
		double t9 = a36*a36;
		double t11 = a45*a45;
		double t13 = a54*a54;
		double t15 = a63*a63;
		double t17 = a72*a72;
		double t20 = a0*s0;
		double t22 = a9*s1;
		double t24 = a18*s2;
		double t26 = a27*s3;
		double t28 = a36*s4;
		double t30 = a45*s5;
		double t32 = a54*s6;
		double t34 = a63*s7;
		double t36 = a72*s8;
		double t38 = t20*a1+t22*a10+t24*a19+t26*a28+t28*a37+t30*a46+t32*a55+t34*a64+t36*a73;
		double t48 = t20*a2+t22*a11+t24*a20+t26*a29+t28*a38+t30*a47+t32*a56+t34*a65+t36*a74;
		double t58 = t20*a3+t22*a12+t24*a21+t26*a30+t28*a39+t30*a48+t32*a57+t34*a66+t36*a75;
		double t68 = t20*a4+t22*a13+t24*a22+t26*a31+t28*a40+t30*a49+t32*a58+t34*a67+t36*a76;
		double t78 = t20*a5+t22*a14+t24*a23+t26*a32+t28*a41+t30*a50+t32*a59+t34*a68+t36*a77;
		double t88 = t20*a6+t22*a15+t24*a24+t26*a33+t28*a42+t30*a51+t32*a60+t34*a69+t36*a78;
		double t98 = t20*a7+t22*a16+t24*a25+t26*a34+t28*a43+t30*a52+t32*a61+t34*a70+t36*a79;
		double t108 = t20*a8+t22*a17+t24*a26+t26*a35+t28*a44+t30*a53+t32*a62+t34*a71+t36*a80;
		double t109 = a1*a1;
		double t111 = a10*a10;
		double t113 = a19*a19;
		double t115 = a28*a28;
		double t117 = a37*a37;
		double t119 = a46*a46;
		double t121 = a55*a55;
		double t123 = a64*a64;
		double t125 = a73*a73;
		double t128 = a1*s0;
		double t130 = a10*s1;
		double t132 = a19*s2;
		double t134 = a28*s3;
		double t136 = a37*s4;
		double t138 = a46*s5;
		double t140 = a55*s6;
		double t142 = a64*s7;
		double t144 = a73*s8;
		double t146 = t128*a2+t130*a11+t132*a20+t134*a29+t136*a38+t138*a47+t140*a56+t142*a65+t144*a74;
		double t156 = t128*a3+t130*a12+t132*a21+t134*a30+t136*a39+t138*a48+t140*a57+t142*a66+t144*a75;
		double t166 = t128*a4+t130*a13+t132*a22+t134*a31+t136*a40+t138*a49+t140*a58+t142*a67+t144*a76;
		double t176 = t128*a5+t130*a14+t132*a23+t134*a32+t136*a41+t138*a50+t140*a59+t142*a68+t144*a77;
		double t186 = t128*a6+t130*a15+t132*a24+t134*a33+t136*a42+t138*a51+t140*a60+t142*a69+t144*a78;
		double t196 = t128*a7+t130*a16+t132*a25+t134*a34+t136*a43+t138*a52+t140*a61+t142*a70+t144*a79;
		double t206 = t128*a8+t130*a17+t132*a26+t134*a35+t136*a44+t138*a53+t140*a62+t142*a71+t144*a80;
		double t207 = a2*a2;
		double t209 = a11*a11;
		double t211 = a20*a20;
		double t213 = a29*a29;
		double t215 = a38*a38;
		double t217 = a47*a47;
		double t219 = a56*a56;
		double t221 = a65*a65;
		double t223 = a74*a74;
		double t226 = a2*s0;
		double t228 = a11*s1;
		double t230 = a20*s2;
		double t232 = a29*s3;
		double t234 = a38*s4;
		double t236 = a47*s5;
		double t238 = a56*s6;
		double t240 = a65*s7;
		double t242 = a74*s8;
		double t244 = t226*a3+t228*a12+t230*a21+t232*a30+t234*a39+t236*a48+t238*a57+t240*a66+t242*a75;
		double t254 = t226*a4+t228*a13+t230*a22+t232*a31+t234*a40+t236*a49+t238*a58+t240*a67+t242*a76;
		double t264 = t226*a5+t228*a14+t230*a23+t232*a32+t234*a41+t236*a50+t238*a59+t240*a68+t242*a77;
		double t274 = t226*a6+t228*a15+t230*a24+t232*a33+t234*a42+t236*a51+t238*a60+t240*a69+t242*a78;
		double t284 = t226*a7+t228*a16+t230*a25+t232*a34+t234*a43+t236*a52+t238*a61+t240*a70+t242*a79;
		double t294 = t226*a8+t228*a17+t230*a26+t232*a35+t234*a44+t236*a53+t238*a62+t240*a71+t242*a80;
		double t295 = a3*a3;
		double t297 = a12*a12;
		double t299 = a21*a21;
		double t301 = a30*a30;
		double t303 = a39*a39;
		double t305 = a48*a48;
		double t307 = a57*a57;
		double t309 = a66*a66;
		double t311 = a75*a75;
		double t314 = a3*s0;
		double t316 = a12*s1;
		double t318 = a21*s2;
		double t320 = a30*s3;
		double t322 = a39*s4;
		double t324 = a48*s5;
		double t326 = a57*s6;
		double t328 = a66*s7;
		double t330 = a75*s8;
		double t332 = t314*a4+t316*a13+t318*a22+t320*a31+t322*a40+t324*a49+t326*a58+t328*a67+t330*a76;
		double t342 = t314*a5+t316*a14+t318*a23+t320*a32+t322*a41+t324*a50+t326*a59+t328*a68+t330*a77;
		double t352 = t314*a6+t316*a15+t318*a24+t320*a33+t322*a42+t324*a51+t326*a60+t328*a69+t330*a78;
		double t362 = t314*a7+t316*a16+t318*a25+t320*a34+t322*a43+t324*a52+t326*a61+t328*a70+t330*a79;
		double t372 = t314*a8+t316*a17+t318*a26+t320*a35+t322*a44+t324*a53+t326*a62+t328*a71+t330*a80;
		double t373 = a4*a4;
		double t375 = a13*a13;
		double t377 = a22*a22;
		double t379 = a31*a31;
		double t381 = a40*a40;
		double t383 = a49*a49;
		double t385 = a58*a58;
		double t387 = a67*a67;
		double t389 = a76*a76;
		double t392 = a4*s0;
		double t394 = a13*s1;
		double t396 = a22*s2;
		double t398 = a31*s3;
		double t400 = a40*s4;
		double t402 = a49*s5;
		double t404 = a58*s6;
		double t406 = a67*s7;
		double t408 = a76*s8;
		double t410 = t392*a5+t394*a14+t396*a23+t398*a32+t400*a41+t402*a50+t404*a59+t406*a68+t408*a77;
		double t420 = t392*a6+t394*a15+t396*a24+t398*a33+t400*a42+t402*a51+t404*a60+t406*a69+t408*a78;
		double t430 = t392*a7+t394*a16+t396*a25+t398*a34+t400*a43+t402*a52+t404*a61+t406*a70+t408*a79;
		double t440 = t392*a8+t394*a17+t396*a26+t398*a35+t400*a44+t402*a53+t404*a62+t406*a71+t408*a80;
		double t441 = a5*a5;
		double t443 = a14*a14;
		double t445 = a23*a23;
		double t447 = a32*a32;
		double t449 = a41*a41;
		double t451 = a50*a50;
		double t453 = a59*a59;
		double t455 = a68*a68;
		double t457 = a77*a77;
		double t460 = a5*s0;
		double t462 = a14*s1;
		double t464 = a23*s2;
		double t466 = a32*s3;
		double t468 = a41*s4;
		double t470 = a50*s5;
		double t472 = a59*s6;
		double t474 = a68*s7;
		double t476 = a77*s8;
		double t478 = t460*a6+t462*a15+t464*a24+t466*a33+t468*a42+t470*a51+t472*a60+t474*a69+t476*a78;
		double t488 = t460*a7+t462*a16+t464*a25+t466*a34+t468*a43+t470*a52+t472*a61+t474*a70+t476*a79;
		double t498 = t460*a8+t462*a17+t464*a26+t466*a35+t468*a44+t470*a53+t472*a62+t474*a71+t476*a80;
		double t499 = a6*a6;
		double t501 = a15*a15;
		double t503 = a24*a24;
		double t505 = a33*a33;
		double t507 = a42*a42;
		double t509 = a51*a51;
		double t511 = a60*a60;
		double t513 = a69*a69;
		double t515 = a78*a78;
		double t518 = a6*s0;
		double t520 = a15*s1;
		double t522 = a24*s2;
		double t524 = a33*s3;
		double t526 = a42*s4;
		double t528 = a51*s5;
		double t530 = a60*s6;
		double t532 = a69*s7;
		double t534 = a78*s8;
		double t536 = t518*a7+t520*a16+t522*a25+t524*a34+t526*a43+t528*a52+t530*a61+t532*a70+t534*a79;
		double t546 = t518*a8+t520*a17+t522*a26+t524*a35+t526*a44+t528*a53+t530*a62+t532*a71+t534*a80;
		double t547 = a7*a7;
		double t549 = a16*a16;
		double t551 = a25*a25;
		double t553 = a34*a34;
		double t555 = a43*a43;
		double t557 = a52*a52;
		double t559 = a61*a61;
		double t561 = a70*a70;
		double t563 = a79*a79;
		double t584 = a7*s0*a8+a16*s1*a17+a25*s2*a26+a34*s3*a35+a43*s4*a44+a52*s5*a53+a61*s6*a62+a70*s7*a71+a79*s8*a80;
		double t585 = a8*a8;
		double t587 = a17*a17;
		double t589 = a26*a26;
		double t591 = a35*a35;
		double t593 = a44*a44;
		double t595 = a53*a53;
		double t597 = a62*a62;
		double t599 = a71*a71;
		double t601 = a80*a80;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6+t15*s7+t17*s8;
		*tmp_M_ptr++ = t38;
		*tmp_M_ptr++ = t48;
		*tmp_M_ptr++ = t58;
		*tmp_M_ptr++ = t68;
		*tmp_M_ptr++ = t78;
		*tmp_M_ptr++ = t88;
		*tmp_M_ptr++ = t98;
		*tmp_M_ptr++ = t108;
		*tmp_M_ptr++ = t109*s0+t111*s1+t113*s2+t115*s3+t117*s4+t119*s5+t121*s6+t123*s7+t125*s8;
		*tmp_M_ptr++ = t146;
		*tmp_M_ptr++ = t156;
		*tmp_M_ptr++ = t166;
		*tmp_M_ptr++ = t176;
		*tmp_M_ptr++ = t186;
		*tmp_M_ptr++ = t196;
		*tmp_M_ptr++ = t206;
		*tmp_M_ptr++ = t207*s0+t209*s1+t211*s2+t213*s3+t215*s4+t217*s5+t219*s6+t221*s7+t223*s8;
		*tmp_M_ptr++ = t244;
		*tmp_M_ptr++ = t254;
		*tmp_M_ptr++ = t264;
		*tmp_M_ptr++ = t274;
		*tmp_M_ptr++ = t284;
		*tmp_M_ptr++ = t294;
		*tmp_M_ptr++ = t295*s0+t297*s1+t299*s2+t301*s3+t303*s4+t305*s5+t307*s6+t309*s7+t311*s8;
		*tmp_M_ptr++ = t332;
		*tmp_M_ptr++ = t342;
		*tmp_M_ptr++ = t352;
		*tmp_M_ptr++ = t362;
		*tmp_M_ptr++ = t372;
		*tmp_M_ptr++ = t373*s0+t375*s1+t377*s2+t379*s3+t381*s4+t383*s5+t385*s6+t387*s7+t389*s8;
		*tmp_M_ptr++ = t410;
		*tmp_M_ptr++ = t420;
		*tmp_M_ptr++ = t430;
		*tmp_M_ptr++ = t440;
		*tmp_M_ptr++ = t441*s0+t443*s1+t445*s2+t447*s3+t449*s4+t451*s5+t453*s6+t455*s7+t457*s8;
		*tmp_M_ptr++ = t478;
		*tmp_M_ptr++ = t488;
		*tmp_M_ptr++ = t498;
		*tmp_M_ptr++ = t499*s0+t501*s1+t503*s2+t505*s3+t507*s4+t509*s5+t511*s6+t513*s7+t515*s8;
		*tmp_M_ptr++ = t536;
		*tmp_M_ptr++ = t546;
		*tmp_M_ptr++ = t547*s0+t549*s1+t551*s2+t553*s3+t555*s4+t557*s5+t559*s6+t561*s7+t563*s8;
		*tmp_M_ptr++ = t584;
		*tmp_M_ptr++ = t585*s0+t587*s1+t589*s2+t591*s3+t593*s4+t595*s5+t597*s6+t599*s7+t601*s8;
	}
	else if ( this->M_size == 10 ) {
		double D[10];
		double U[100];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 10);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80,a81,a82,a83,a84,a85,a86,a87,a88,a89,a90,a91,a92,a93,a94,a95,a96,a97,a98,a99;
		double s0,s1,s2,s3,s4,s5,s6,s7,s8,s9;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;
		a49 = *tmp_M_ptr++;
		a50 = *tmp_M_ptr++;
		a51 = *tmp_M_ptr++;
		a52 = *tmp_M_ptr++;
		a53 = *tmp_M_ptr++;
		a54 = *tmp_M_ptr++;
		a55 = *tmp_M_ptr++;
		a56 = *tmp_M_ptr++;
		a57 = *tmp_M_ptr++;
		a58 = *tmp_M_ptr++;
		a59 = *tmp_M_ptr++;
		a60 = *tmp_M_ptr++;
		a61 = *tmp_M_ptr++;
		a62 = *tmp_M_ptr++;
		a63 = *tmp_M_ptr++;
		a64 = *tmp_M_ptr++;
		a65 = *tmp_M_ptr++;
		a66 = *tmp_M_ptr++;
		a67 = *tmp_M_ptr++;
		a68 = *tmp_M_ptr++;
		a69 = *tmp_M_ptr++;
		a70 = *tmp_M_ptr++;
		a71 = *tmp_M_ptr++;
		a72 = *tmp_M_ptr++;
		a73 = *tmp_M_ptr++;
		a74 = *tmp_M_ptr++;
		a75 = *tmp_M_ptr++;
		a76 = *tmp_M_ptr++;
		a77 = *tmp_M_ptr++;
		a78 = *tmp_M_ptr++;
		a79 = *tmp_M_ptr++;
		a80 = *tmp_M_ptr++;
		a81 = *tmp_M_ptr++;
		a82 = *tmp_M_ptr++;
		a83 = *tmp_M_ptr++;
		a84 = *tmp_M_ptr++;
		a85 = *tmp_M_ptr++;
		a86 = *tmp_M_ptr++;
		a87 = *tmp_M_ptr++;
		a88 = *tmp_M_ptr++;
		a89 = *tmp_M_ptr++;
		a90 = *tmp_M_ptr++;
		a91 = *tmp_M_ptr++;
		a92 = *tmp_M_ptr++;
		a93 = *tmp_M_ptr++;
		a94 = *tmp_M_ptr++;
		a95 = *tmp_M_ptr++;
		a96 = *tmp_M_ptr++;
		a97 = *tmp_M_ptr++;
		a98 = *tmp_M_ptr++;
		a99 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = sqrt(*tmp_M_ptr++);
		s1 = sqrt(*tmp_M_ptr++);
		s2 = sqrt(*tmp_M_ptr++);
		s3 = sqrt(*tmp_M_ptr++);
		s4 = sqrt(*tmp_M_ptr++);
		s5 = sqrt(*tmp_M_ptr++);
		s6 = sqrt(*tmp_M_ptr++);
		s7 = sqrt(*tmp_M_ptr++);
		s8 = sqrt(*tmp_M_ptr++);
		s9 = sqrt(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a10*a10;
		double t5 = a20*a20;
		double t7 = a30*a30;
		double t9 = a40*a40;
		double t11 = a50*a50;
		double t13 = a60*a60;
		double t15 = a70*a70;
		double t17 = a80*a80;
		double t19 = a90*a90;
		double t22 = a0*s0;
		double t24 = a10*s1;
		double t26 = a20*s2;
		double t28 = a30*s3;
		double t30 = a40*s4;
		double t32 = a50*s5;
		double t34 = a60*s6;
		double t36 = a70*s7;
		double t38 = a80*s8;
		double t40 = a90*s9;
		double t42 = t22*a1+t24*a11+t26*a21+t28*a31+t30*a41+t32*a51+t34*a61+t36*a71+t38*a81+t40*a91;
		double t53 = t22*a2+t24*a12+t26*a22+t28*a32+t30*a42+t32*a52+t34*a62+t36*a72+t38*a82+t40*a92;
		double t64 = t22*a3+t24*a13+t26*a23+t28*a33+t30*a43+t32*a53+t34*a63+t36*a73+t38*a83+t40*a93;
		double t75 = t22*a4+t24*a14+t26*a24+t28*a34+t30*a44+t32*a54+t34*a64+t36*a74+t38*a84+t40*a94;
		double t86 = t22*a5+t24*a15+t26*a25+t28*a35+t30*a45+t32*a55+t34*a65+t36*a75+t38*a85+t40*a95;
		double t97 = t22*a6+t24*a16+t26*a26+t28*a36+t30*a46+t32*a56+t34*a66+t36*a76+t38*a86+t40*a96;
		double t108 = t22*a7+t24*a17+t26*a27+t28*a37+t30*a47+t32*a57+t34*a67+t36*a77+t38*a87+t40*a97;
		double t119 = t22*a8+t24*a18+t26*a28+t28*a38+t30*a48+t32*a58+t34*a68+t36*a78+t38*a88+t40*a98;
		double t130 = t22*a9+t24*a19+t26*a29+t28*a39+t30*a49+t32*a59+t34*a69+t36*a79+t38*a89+t40*a99;
		double t131 = a1*a1;
		double t133 = a11*a11;
		double t135 = a21*a21;
		double t137 = a31*a31;
		double t139 = a41*a41;
		double t141 = a51*a51;
		double t143 = a61*a61;
		double t145 = a71*a71;
		double t147 = a81*a81;
		double t149 = a91*a91;
		double t152 = a1*s0;
		double t154 = a11*s1;
		double t156 = a21*s2;
		double t158 = a31*s3;
		double t160 = a41*s4;
		double t162 = a51*s5;
		double t164 = a61*s6;
		double t166 = a71*s7;
		double t168 = a81*s8;
		double t170 = a91*s9;
		double t172 = t152*a2+t154*a12+t156*a22+t158*a32+t160*a42+t162*a52+t164*a62+t166*a72+t168*a82+t170*a92;
		double t183 = t152*a3+t154*a13+t156*a23+t158*a33+t160*a43+t162*a53+t164*a63+t166*a73+t168*a83+t170*a93;
		double t194 = t152*a4+t154*a14+t156*a24+t158*a34+t160*a44+t162*a54+t164*a64+t166*a74+t168*a84+t170*a94;
		double t205 = t152*a5+t154*a15+t156*a25+t158*a35+t160*a45+t162*a55+t164*a65+t166*a75+t168*a85+t170*a95;
		double t216 = t152*a6+t154*a16+t156*a26+t158*a36+t160*a46+t162*a56+t164*a66+t166*a76+t168*a86+t170*a96;
		double t227 = t152*a7+t154*a17+t156*a27+t158*a37+t160*a47+t162*a57+t164*a67+t166*a77+t168*a87+t170*a97;
		double t238 = t152*a8+t154*a18+t156*a28+t158*a38+t160*a48+t162*a58+t164*a68+t166*a78+t168*a88+t170*a98;
		double t249 = t152*a9+t154*a19+t156*a29+t158*a39+t160*a49+t162*a59+t164*a69+t166*a79+t168*a89+t170*a99;
		double t250 = a2*a2;
		double t252 = a12*a12;
		double t254 = a22*a22;
		double t256 = a32*a32;
		double t258 = a42*a42;
		double t260 = a52*a52;
		double t262 = a62*a62;
		double t264 = a72*a72;
		double t266 = a82*a82;
		double t268 = a92*a92;
		double t271 = a2*s0;
		double t273 = a12*s1;
		double t275 = a22*s2;
		double t277 = a32*s3;
		double t279 = a42*s4;
		double t281 = a52*s5;
		double t283 = a62*s6;
		double t285 = a72*s7;
		double t287 = a82*s8;
		double t289 = a92*s9;
		double t291 = t271*a3+t273*a13+t275*a23+t277*a33+t279*a43+t281*a53+t283*a63+t285*a73+t287*a83+t289*a93;
		double t302 = t271*a4+t273*a14+t275*a24+t277*a34+t279*a44+t281*a54+t283*a64+t285*a74+t287*a84+t289*a94;
		double t313 = t271*a5+t273*a15+t275*a25+t277*a35+t279*a45+t281*a55+t283*a65+t285*a75+t287*a85+t289*a95;
		double t324 = t271*a6+t273*a16+t275*a26+t277*a36+t279*a46+t281*a56+t283*a66+t285*a76+t287*a86+t289*a96;
		double t335 = t271*a7+t273*a17+t275*a27+t277*a37+t279*a47+t281*a57+t283*a67+t285*a77+t287*a87+t289*a97;
		double t346 = t271*a8+t273*a18+t275*a28+t277*a38+t279*a48+t281*a58+t283*a68+t285*a78+t287*a88+t289*a98;
		double t357 = t271*a9+t273*a19+t275*a29+t277*a39+t279*a49+t281*a59+t283*a69+t285*a79+t287*a89+t289*a99;
		double t358 = a3*a3;
		double t360 = a13*a13;
		double t362 = a23*a23;
		double t364 = a33*a33;
		double t366 = a43*a43;
		double t368 = a53*a53;
		double t370 = a63*a63;
		double t372 = a73*a73;
		double t374 = a83*a83;
		double t376 = a93*a93;
		double t379 = a3*s0;
		double t381 = a13*s1;
		double t383 = a23*s2;
		double t385 = a33*s3;
		double t387 = a43*s4;
		double t389 = a53*s5;
		double t391 = a63*s6;
		double t393 = a73*s7;
		double t395 = a83*s8;
		double t397 = a93*s9;
		double t399 = t379*a4+t381*a14+t383*a24+t385*a34+t387*a44+t389*a54+t391*a64+t393*a74+t395*a84+t397*a94;
		double t410 = t379*a5+t381*a15+t383*a25+t385*a35+t387*a45+t389*a55+t391*a65+t393*a75+t395*a85+t397*a95;
		double t421 = t379*a6+t381*a16+t383*a26+t385*a36+t387*a46+t389*a56+t391*a66+t393*a76+t395*a86+t397*a96;
		double t432 = t379*a7+t381*a17+t383*a27+t385*a37+t387*a47+t389*a57+t391*a67+t393*a77+t395*a87+t397*a97;
		double t443 = t379*a8+t381*a18+t383*a28+t385*a38+t387*a48+t389*a58+t391*a68+t393*a78+t395*a88+t397*a98;
		double t454 = t379*a9+t381*a19+t383*a29+t385*a39+t387*a49+t389*a59+t391*a69+t393*a79+t395*a89+t397*a99;
		double t455 = a4*a4;
		double t457 = a14*a14;
		double t459 = a24*a24;
		double t461 = a34*a34;
		double t463 = a44*a44;
		double t465 = a54*a54;
		double t467 = a64*a64;
		double t469 = a74*a74;
		double t471 = a84*a84;
		double t473 = a94*a94;
		double t476 = a4*s0;
		double t478 = a14*s1;
		double t480 = a24*s2;
		double t482 = a34*s3;
		double t484 = a44*s4;
		double t486 = a54*s5;
		double t488 = a64*s6;
		double t490 = a74*s7;
		double t492 = a84*s8;
		double t494 = a94*s9;
		double t496 = t476*a5+t478*a15+t480*a25+t482*a35+t484*a45+t486*a55+t488*a65+t490*a75+t492*a85+t494*a95;
		double t507 = t476*a6+t478*a16+t480*a26+t482*a36+t484*a46+t486*a56+t488*a66+t490*a76+t492*a86+t494*a96;
		double t518 = t476*a7+t478*a17+t480*a27+t482*a37+t484*a47+t486*a57+t488*a67+t490*a77+t492*a87+t494*a97;
		double t529 = t476*a8+t478*a18+t480*a28+t482*a38+t484*a48+t486*a58+t488*a68+t490*a78+t492*a88+t494*a98;
		double t540 = t476*a9+t478*a19+t480*a29+t482*a39+t484*a49+t486*a59+t488*a69+t490*a79+t492*a89+t494*a99;
		double t541 = a5*a5;
		double t543 = a15*a15;
		double t545 = a25*a25;
		double t547 = a35*a35;
		double t549 = a45*a45;
		double t551 = a55*a55;
		double t553 = a65*a65;
		double t555 = a75*a75;
		double t557 = a85*a85;
		double t559 = a95*a95;
		double t562 = a5*s0;
		double t564 = a15*s1;
		double t566 = a25*s2;
		double t568 = a35*s3;
		double t570 = a45*s4;
		double t572 = a55*s5;
		double t574 = a65*s6;
		double t576 = a75*s7;
		double t578 = a85*s8;
		double t580 = a95*s9;
		double t582 = t562*a6+t564*a16+t566*a26+t568*a36+t570*a46+t572*a56+t574*a66+t576*a76+t578*a86+t580*a96;
		double t593 = t562*a7+t564*a17+t566*a27+t568*a37+t570*a47+t572*a57+t574*a67+t576*a77+t578*a87+t580*a97;
		double t604 = t562*a8+t564*a18+t566*a28+t568*a38+t570*a48+t572*a58+t574*a68+t576*a78+t578*a88+t580*a98;
		double t615 = t562*a9+t564*a19+t566*a29+t568*a39+t570*a49+t572*a59+t574*a69+t576*a79+t578*a89+t580*a99;
		double t616 = a6*a6;
		double t618 = a16*a16;
		double t620 = a26*a26;
		double t622 = a36*a36;
		double t624 = a46*a46;
		double t626 = a56*a56;
		double t628 = a66*a66;
		double t630 = a76*a76;
		double t632 = a86*a86;
		double t634 = a96*a96;
		double t637 = a6*s0;
		double t639 = a16*s1;
		double t641 = a26*s2;
		double t643 = a36*s3;
		double t645 = a46*s4;
		double t647 = a56*s5;
		double t649 = a66*s6;
		double t651 = a76*s7;
		double t653 = a86*s8;
		double t655 = a96*s9;
		double t657 = t637*a7+t639*a17+t641*a27+t643*a37+t645*a47+t647*a57+t649*a67+t651*a77+t653*a87+t655*a97;
		double t668 = t637*a8+t639*a18+t641*a28+t643*a38+t645*a48+t647*a58+t649*a68+t651*a78+t653*a88+t655*a98;
		double t679 = t637*a9+t639*a19+t641*a29+t643*a39+t645*a49+t647*a59+t649*a69+t651*a79+t653*a89+t655*a99;
		double t680 = a7*a7;
		double t682 = a17*a17;
		double t684 = a27*a27;
		double t686 = a37*a37;
		double t688 = a47*a47;
		double t690 = a57*a57;
		double t692 = a67*a67;
		double t694 = a77*a77;
		double t696 = a87*a87;
		double t698 = a97*a97;
		double t701 = a7*s0;
		double t703 = a17*s1;
		double t705 = a27*s2;
		double t707 = a37*s3;
		double t709 = a47*s4;
		double t711 = a57*s5;
		double t713 = a67*s6;
		double t715 = a77*s7;
		double t717 = a87*s8;
		double t719 = a97*s9;
		double t721 = t701*a8+t703*a18+t705*a28+t707*a38+t709*a48+t711*a58+t713*a68+t715*a78+t717*a88+t719*a98;
		double t732 = t701*a9+t703*a19+t705*a29+t707*a39+t709*a49+t711*a59+t713*a69+t715*a79+t717*a89+t719*a99;
		double t733 = a8*a8;
		double t735 = a18*a18;
		double t737 = a28*a28;
		double t739 = a38*a38;
		double t741 = a48*a48;
		double t743 = a58*a58;
		double t745 = a68*a68;
		double t747 = a78*a78;
		double t749 = a88*a88;
		double t751 = a98*a98;
		double t774 = a8*s0*a9+a18*s1*a19+a28*s2*a29+a38*s3*a39+a48*s4*a49+a58*s5*a59+a68*s6*a69+a78*s7*a79+a88*s8*a89+a98*s9*a99;
		double t775 = a9*a9;
		double t777 = a19*a19;
		double t779 = a29*a29;
		double t781 = a39*a39;
		double t783 = a49*a49;
		double t785 = a59*a59;
		double t787 = a69*a69;
		double t789 = a79*a79;
		double t791 = a89*a89;
		double t793 = a99*a99;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6+t15*s7+t17*s8+t19*s9;
		*tmp_M_ptr++ = t42;
		*tmp_M_ptr++ = t53;
		*tmp_M_ptr++ = t64;
		*tmp_M_ptr++ = t75;
		*tmp_M_ptr++ = t86;
		*tmp_M_ptr++ = t97;
		*tmp_M_ptr++ = t108;
		*tmp_M_ptr++ = t119;
		*tmp_M_ptr++ = t130;
		*tmp_M_ptr++ = t131*s0+t133*s1+t135*s2+t137*s3+t139*s4+t141*s5+t143*s6+t145*s7+t147*s8+t149*s9;
		*tmp_M_ptr++ = t172;
		*tmp_M_ptr++ = t183;
		*tmp_M_ptr++ = t194;
		*tmp_M_ptr++ = t205;
		*tmp_M_ptr++ = t216;
		*tmp_M_ptr++ = t227;
		*tmp_M_ptr++ = t238;
		*tmp_M_ptr++ = t249;
		*tmp_M_ptr++ = t250*s0+t252*s1+t254*s2+t256*s3+t258*s4+t260*s5+t262*s6+t264*s7+t266*s8+t268*s9;
		*tmp_M_ptr++ = t291;
		*tmp_M_ptr++ = t302;
		*tmp_M_ptr++ = t313;
		*tmp_M_ptr++ = t324;
		*tmp_M_ptr++ = t335;
		*tmp_M_ptr++ = t346;
		*tmp_M_ptr++ = t357;
		*tmp_M_ptr++ = t358*s0+t360*s1+t362*s2+t364*s3+t366*s4+t368*s5+t370*s6+t372*s7+t374*s8+t376*s9;
		*tmp_M_ptr++ = t399;
		*tmp_M_ptr++ = t410;
		*tmp_M_ptr++ = t421;
		*tmp_M_ptr++ = t432;
		*tmp_M_ptr++ = t443;
		*tmp_M_ptr++ = t454;
		*tmp_M_ptr++ = t455*s0+t457*s1+t459*s2+t461*s3+t463*s4+t465*s5+t467*s6+t469*s7+t471*s8+t473*s9;
		*tmp_M_ptr++ = t496;
		*tmp_M_ptr++ = t507;
		*tmp_M_ptr++ = t518;
		*tmp_M_ptr++ = t529;
		*tmp_M_ptr++ = t540;
		*tmp_M_ptr++ = t541*s0+t543*s1+t545*s2+t547*s3+t549*s4+t551*s5+t553*s6+t555*s7+t557*s8+t559*s9;
		*tmp_M_ptr++ = t582;
		*tmp_M_ptr++ = t593;
		*tmp_M_ptr++ = t604;
		*tmp_M_ptr++ = t615;
		*tmp_M_ptr++ = t616*s0+t618*s1+t620*s2+t622*s3+t624*s4+t626*s5+t628*s6+t630*s7+t632*s8+t634*s9;
		*tmp_M_ptr++ = t657;
		*tmp_M_ptr++ = t668;
		*tmp_M_ptr++ = t679;
		*tmp_M_ptr++ = t680*s0+t682*s1+t684*s2+t686*s3+t688*s4+t690*s5+t692*s6+t694*s7+t696*s8+t698*s9;
		*tmp_M_ptr++ = t721;
		*tmp_M_ptr++ = t732;
		*tmp_M_ptr++ = t733*s0+t735*s1+t737*s2+t739*s3+t741*s4+t743*s5+t745*s6+t747*s7+t749*s8+t751*s9;
		*tmp_M_ptr++ = t774;
		*tmp_M_ptr++ = t775*s0+t777*s1+t779*s2+t781*s3+t783*s4+t785*s5+t787*s6+t789*s7+t791*s8+t793*s9;
	}
	else {
		CvMat* D = cvCreateMat(M_size, 1, CV_64F);
		CvMat* U = cvCreateMat(M_size, M_size, CV_64F);

		this->SVD(D, U);

		CvMat* U_t = cvCreateMat(M_size, M_size, CV_64F);
		CvMat* U_tU = cvCreateMat(M_size, M_size, CV_64F);

		cvTranspose(U, U_t);

		/* D = sqrt(D) */
		Sqrt_Diag(D);

		/* U_t = U_t * D */
		Mul_Diag(U_t, D);

		/* U_tU = U_t * U */
		cvMatMul(U_t, U, U_tU);

		dst_mat->SetData((double*)U_tU->data.ptr, M_size2);

		cvReleaseMat(&U_t);
		cvReleaseMat(&U_tU);

		cvReleaseMat(&D);
		cvReleaseMat(&U);
	}
}

void CCovarianceMatrix::Log(CCovarianceMatrix *log_mat)
{
	CCovarianceMatrix *dst_mat;

	if ( log_mat )
		dst_mat = log_mat;
	else
		dst_mat = this;

	if ( this->M_size == 1 ) {
		if ( SM_ptr[0] <= 0 ) {
			printf("SM_ptr[0] = %.20f - log error in Log(CCovarianceMatrix *log_mat)!\n", SM_ptr[0]);
			SM_ptr[0] = MIN_DOUBLE_NUM;
		}
		dst_mat->SM_ptr[0] = log(SM_ptr[0]);
	}
	else if ( this->M_size == 2 ) {
		double D[2];
		double U[4];

		eigen2_decomposition(this->SM_ptr, U, D);

		double a0,a1,a2,a3;
		double s0,s1;

		a0 = U[0];
		a1 = U[1];
		a2 = U[2];
		a3 = U[3];

		if ( D[0] <= 0 || D[1] <= 0 ) {
			printf("SM_ptr[0] = %.20f\tSM_ptr[1] = %.20f - log error in Log(CCovarianceMatrix *log_mat)!\n", D[0], D[1]);
			D[0] = MAX(D[0], MIN_DOUBLE_NUM);
			D[1] = MAX(D[1], MIN_DOUBLE_NUM);
		}

		s0 = log(D[0]);
		s1 = log(D[1]);

		double t1 = a0*a0;
		double t3 = a2*a2;
		double t10 = a0*s0*a1+a2*s1*a3;
		double t11 = a1*a1;
		double t13 = a3*a3;

		double* tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1;
		*tmp_M_ptr++ = t10;
		*tmp_M_ptr++ = t11*s0+t13*s1;
	}
	else if ( this->M_size == 3 ) {
		double D[3];
		double U[9];

		eigen3_decomposition(this->SM_ptr, U, D);

		double a0,a1,a2,a3,a4,a5,a6,a7,a8;
		double s0,s1,s2;

		a0 = U[0];
		a1 = U[1];
		a2 = U[2];
		a3 = U[3];
		a4 = U[4];
		a5 = U[5];
		a6 = U[6];
		a7 = U[7];
		a8 = U[8];

		if ( D[0] <= 0 || D[1] <= 0 || D[2] <= 0 ) {
			printf("SM_ptr[0] = %.20f\tSM_ptr[1] = %.20f\tSM_ptr[2] = %.20f - log error in Log(CCovarianceMatrix *log_mat)!\n", D[0], D[1], D[2]);
			D[0] = MAX(D[0], MIN_DOUBLE_NUM);
			D[1] = MAX(D[1], MIN_DOUBLE_NUM);
			D[2] = MAX(D[2], MIN_DOUBLE_NUM);
		}

		s0 = log(D[0]);
		s1 = log(D[1]);
		s2 = log(D[2]);

		double t1 = a0*a0;
		double t3 = a3*a3;
		double t5 = a6*a6;
		double t8 = a0*s0;
		double t10 = a3*s1;
		double t12 = a6*s2;
		double t14 = t8*a1+t10*a4+t12*a7;
		double t18 = t8*a2+t10*a5+t12*a8;
		double t19 = a1*a1;
		double t21 = a4*a4;
		double t23 = a7*a7;
		double t32 = a1*s0*a2+a4*s1*a5+a7*s2*a8;
		double t33 = a2*a2;
		double t35 = a5*a5;
		double t37 = a8*a8;

		double *tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2;
		*tmp_M_ptr++ = t14;
		*tmp_M_ptr++ = t18;
		*tmp_M_ptr++ = t19*s0+t21*s1+t23*s2;
		*tmp_M_ptr++ = t32;
		*tmp_M_ptr++ = t33*s0+t35*s1+t37*s2;
	}
	else if ( this->M_size == 4 ) {
		double D[4];
		double U[16];

		eigen4_decomposition(this->SM_ptr, U, D);

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
		double s0,s1,s2,s3;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = log(*tmp_M_ptr++);
		s1 = log(*tmp_M_ptr++);
		s2 = log(*tmp_M_ptr++);
		s3 = log(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a4*a4;
		double t5 = a8*a8;
		double t7 = a12*a12;
		double t10 = a0*s0;
		double t12 = a4*s1;
		double t14 = a8*s2;
		double t16 = a12*s3;
		double t18 = t10*a1+t12*a5+t14*a9+t16*a13;
		double t23 = t10*a2+t12*a6+t14*a10+t16*a14;
		double t28 = t10*a3+t12*a7+t14*a11+t16*a15;
		double t29 = a1*a1;
		double t31 = a5*a5;
		double t33 = a9*a9;
		double t35 = a13*a13;
		double t38 = a1*s0;
		double t40 = a5*s1;
		double t42 = a9*s2;
		double t44 = a13*s3;
		double t46 = t38*a2+t40*a6+t42*a10+t44*a14;
		double t51 = t38*a3+t40*a7+t42*a11+t44*a15;
		double t52 = a2*a2;
		double t54 = a6*a6;
		double t56 = a10*a10;
		double t58 = a14*a14;
		double t69 = a2*s0*a3+a6*s1*a7+a10*s2*a11+a14*s3*a15;
		double t70 = a3*a3;
		double t72 = a7*a7;
		double t74 = a11*a11;
		double t76 = a15*a15;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3;
		*tmp_M_ptr++ = t18;
		*tmp_M_ptr++ = t23;
		*tmp_M_ptr++ = t28;
		*tmp_M_ptr++ = t29*s0+t31*s1+t33*s2+t35*s3;
		*tmp_M_ptr++ = t46;
		*tmp_M_ptr++ = t51;
		*tmp_M_ptr++ = t52*s0+t54*s1+t56*s2+t58*s3;
		*tmp_M_ptr++ = t69;
		*tmp_M_ptr++ = t70*s0+t72*s1+t74*s2+t76*s3;
	}
	else if ( this->M_size == 5 ) {
		double D[5];
		double U[25];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 5);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24;
		double s0,s1,s2,s3,s4;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = log(*tmp_M_ptr++);
		s1 = log(*tmp_M_ptr++);
		s2 = log(*tmp_M_ptr++);
		s3 = log(*tmp_M_ptr++);
		s4 = log(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a5*a5;
		double t5 = a10*a10;
		double t7 = a15*a15;
		double t9 = a20*a20;
		double t12 = a0*s0;
		double t14 = a5*s1;
		double t16 = a10*s2;
		double t18 = a15*s3;
		double t20 = a20*s4;
		double t22 = t12*a1+t14*a6+t16*a11+t18*a16+t20*a21;
		double t28 = t12*a2+t14*a7+t16*a12+t18*a17+t20*a22;
		double t34 = t12*a3+t14*a8+t16*a13+t18*a18+t20*a23;
		double t40 = t12*a4+t14*a9+t16*a14+t18*a19+t20*a24;
		double t41 = a1*a1;
		double t43 = a6*a6;
		double t45 = a11*a11;
		double t47 = a16*a16;
		double t49 = a21*a21;
		double t52 = a1*s0;
		double t54 = a6*s1;
		double t56 = a11*s2;
		double t58 = a16*s3;
		double t60 = a21*s4;
		double t62 = t52*a2+t54*a7+t56*a12+t58*a17+t60*a22;
		double t68 = t52*a3+t54*a8+t56*a13+t58*a18+t60*a23;
		double t74 = t52*a4+t54*a9+t56*a14+t58*a19+t60*a24;
		double t75 = a2*a2;
		double t77 = a7*a7;
		double t79 = a12*a12;
		double t81 = a17*a17;
		double t83 = a22*a22;
		double t86 = a2*s0;
		double t88 = a7*s1;
		double t90 = a12*s2;
		double t92 = a17*s3;
		double t94 = a22*s4;
		double t96 = t86*a3+t88*a8+t90*a13+t92*a18+t94*a23;
		double t102 = t86*a4+t88*a9+t90*a14+t92*a19+t94*a24;
		double t103 = a3*a3;
		double t105 = a8*a8;
		double t107 = a13*a13;
		double t109 = a18*a18;
		double t111 = a23*a23;
		double t124 = a3*s0*a4+a8*s1*a9+a13*s2*a14+a18*s3*a19+a23*s4*a24;
		double t125 = a4*a4;
		double t127 = a9*a9;
		double t129 = a14*a14;
		double t131 = a19*a19;
		double t133 = a24*a24;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4;
		*tmp_M_ptr++  = t22;
		*tmp_M_ptr++  = t28;
		*tmp_M_ptr++  = t34;
		*tmp_M_ptr++  = t40;
		*tmp_M_ptr++  = t41*s0+t43*s1+t45*s2+t47*s3+t49*s4;
		*tmp_M_ptr++  = t62;
		*tmp_M_ptr++  = t68;
		*tmp_M_ptr++  = t74;
		*tmp_M_ptr++  = t75*s0+t77*s1+t79*s2+t81*s3+t83*s4;
		*tmp_M_ptr++  = t96;
		*tmp_M_ptr++  = t102;
		*tmp_M_ptr++  = t103*s0+t105*s1+t107*s2+t109*s3+t111*s4;
		*tmp_M_ptr++  = t124;
		*tmp_M_ptr++  = t125*s0+t127*s1+t129*s2+t131*s3+t133*s4;
	}
	else if ( this->M_size == 6 ) {
		double D[6];
		double U[36];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 6);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35;
		double s0,s1,s2,s3,s4,s5;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = log(*tmp_M_ptr++);
		s1 = log(*tmp_M_ptr++);
		s2 = log(*tmp_M_ptr++);
		s3 = log(*tmp_M_ptr++);
		s4 = log(*tmp_M_ptr++);
		s5 = log(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a6*a6;
		double t5 = a12*a12;
		double t7 = a18*a18;
		double t9 = a24*a24;
		double t11 = a30*a30;
		double t14 = a0*s0;
		double t16 = a6*s1;
		double t18 = a12*s2;
		double t20 = a18*s3;
		double t22 = a24*s4;
		double t24 = a30*s5;
		double t26 = t14*a1+t16*a7+t18*a13+t20*a19+t22*a25+t24*a31;
		double t33 = t14*a2+t16*a8+t18*a14+t20*a20+t22*a26+t24*a32;
		double t40 = t14*a3+t16*a9+t18*a15+t20*a21+t22*a27+t24*a33;
		double t47 = t14*a4+10.0*t16+t18*a16+t20*a22+t22*a28+t24*a34;
		double t54 = t14*a5+t16*a11+t18*a17+t20*a23+t22*a29+t24*a35;
		double t55 = a1*a1;
		double t57 = a7*a7;
		double t59 = a13*a13;
		double t61 = a19*a19;
		double t63 = a25*a25;
		double t65 = a31*a31;
		double t68 = a1*s0;
		double t70 = a7*s1;
		double t72 = a13*s2;
		double t74 = a19*s3;
		double t76 = a25*s4;
		double t78 = a31*s5;
		double t80 = t68*a2+t70*a8+t72*a14+t74*a20+t76*a26+t78*a32;
		double t87 = t68*a3+t70*a9+t72*a15+t74*a21+t76*a27+t78*a33;
		double t94 = t68*a4+10.0*t70+t72*a16+t74*a22+t76*a28+t78*a34;
		double t101 = t68*a5+t70*a11+t72*a17+t74*a23+t76*a29+t78*a35;
		double t102 = a2*a2;
		double t104 = a8*a8;
		double t106 = a14*a14;
		double t108 = a20*a20;
		double t110 = a26*a26;
		double t112 = a32*a32;
		double t115 = a2*s0;
		double t117 = a8*s1;
		double t119 = a14*s2;
		double t121 = a20*s3;
		double t123 = a26*s4;
		double t125 = a32*s5;
		double t127 = t115*a3+t117*a9+t119*a15+t121*a21+t123*a27+t125*a33;
		double t134 = t115*a4+10.0*t117+t119*a16+t121*a22+t123*a28+t125*a34;
		double t141 = t115*a5+t117*a11+t119*a17+t121*a23+t123*a29+t125*a35;
		double t142 = a3*a3;
		double t144 = a9*a9;
		double t146 = a15*a15;
		double t148 = a21*a21;
		double t150 = a27*a27;
		double t152 = a33*a33;
		double t155 = a3*s0;
		double t157 = a9*s1;
		double t159 = a15*s2;
		double t161 = a21*s3;
		double t163 = a27*s4;
		double t165 = a33*s5;
		double t167 = t155*a4+10.0*t157+t159*a16+t161*a22+t163*a28+t165*a34;
		double t174 = t155*a5+t157*a11+t159*a17+t161*a23+t163*a29+t165*a35;
		double t175 = a4*a4;
		double t178 = a16*a16;
		double t180 = a22*a22;
		double t182 = a28*a28;
		double t184 = a34*a34;
		double t199 = a4*s0*a5+10.0*a11*s1+a16*s2*a17+a22*s3*a23+a28*s4*a29+a34*s5*a35;
		double t200 = a5*a5;
		double t202 = a11*a11;
		double t204 = a17*a17;
		double t206 = a23*a23;
		double t208 = a29*a29;
		double t210 = a35*a35;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5;
		*tmp_M_ptr++ = t26;
		*tmp_M_ptr++ = t33;
		*tmp_M_ptr++ = t40;
		*tmp_M_ptr++ = t47;
		*tmp_M_ptr++ = t54;
		*tmp_M_ptr++ = t55*s0+t57*s1+t59*s2+t61*s3+t63*s4+t65*s5;
		*tmp_M_ptr++ = t80;
		*tmp_M_ptr++ = t87;
		*tmp_M_ptr++ = t94;
		*tmp_M_ptr++ = t101;
		*tmp_M_ptr++ = t102*s0+t104*s1+t106*s2+t108*s3+t110*s4+t112*s5;
		*tmp_M_ptr++ = t127;
		*tmp_M_ptr++ = t134;
		*tmp_M_ptr++ = t141;
		*tmp_M_ptr++ = t142*s0+t144*s1+t146*s2+t148*s3+t150*s4+t152*s5;
		*tmp_M_ptr++ = t167;
		*tmp_M_ptr++ = t174;
		*tmp_M_ptr++ = t175*s0+100.0*s1+t178*s2+t180*s3+t182*s4+t184*s5;
		*tmp_M_ptr++ = t199;
		*tmp_M_ptr++ = t200*s0+t202*s1+t204*s2+t206*s3+t208*s4+t210*s5;
	}
	else if ( this->M_size == 7 ) {
		double D[7];
		double U[49];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 7);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double s0,s1,s2,s3,s4,s5,s6;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = log(*tmp_M_ptr++);
		s1 = log(*tmp_M_ptr++);
		s2 = log(*tmp_M_ptr++);
		s3 = log(*tmp_M_ptr++);
		s4 = log(*tmp_M_ptr++);
		s5 = log(*tmp_M_ptr++);
		s6 = log(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a7*a7;
		double t5 = a14*a14;
		double t7 = a21*a21;
		double t9 = a28*a28;
		double t11 = a35*a35;
		double t13 = a42*a42;
		double t16 = a0*s0;
		double t18 = a7*s1;
		double t20 = a14*s2;
		double t22 = a21*s3;
		double t24 = a28*s4;
		double t26 = a35*s5;
		double t28 = a42*s6;
		double t30 = t16*a1+t18*a8+t20*a15+t22*a22+t24*a29+t26*a36+t28*a43;
		double t38 = t16*a2+t18*a9+t20*a16+t22*a23+t24*a30+t26*a37+t28*a44;
		double t46 = t16*a3+t18*a10+t20*a17+t22*a24+t24*a31+t26*a38+t28*a45;
		double t54 = t16*a4+t18*a11+t20*a18+t22*a25+t24*a32+t26*a39+t28*a46;
		double t62 = t16*a5+t18*a12+t20*a19+t22*a26+t24*a33+t26*a40+t28*a47;
		double t70 = t16*a6+t18*a13+t20*a20+t22*a27+t24*a34+t26*a41+t28*a48;
		double t71 = a1*a1;
		double t73 = a8*a8;
		double t75 = a15*a15;
		double t77 = a22*a22;
		double t79 = a29*a29;
		double t81 = a36*a36;
		double t83 = a43*a43;
		double t86 = a1*s0;
		double t88 = a8*s1;
		double t90 = a15*s2;
		double t92 = a22*s3;
		double t94 = a29*s4;
		double t96 = a36*s5;
		double t98 = a43*s6;
		double t100 = t86*a2+t88*a9+t90*a16+t92*a23+t94*a30+t96*a37+t98*a44;
		double t108 = t86*a3+t88*a10+t90*a17+t92*a24+t94*a31+t96*a38+t98*a45;
		double t116 = t86*a4+t88*a11+t90*a18+t92*a25+t94*a32+t96*a39+t98*a46;
		double t124 = t86*a5+t88*a12+t90*a19+t92*a26+t94*a33+t96*a40+t98*a47;
		double t132 = t86*a6+t88*a13+t90*a20+t92*a27+t94*a34+t96*a41+t98*a48;
		double t133 = a2*a2;
		double t135 = a9*a9;
		double t137 = a16*a16;
		double t139 = a23*a23;
		double t141 = a30*a30;
		double t143 = a37*a37;
		double t145 = a44*a44;
		double t148 = a2*s0;
		double t150 = a9*s1;
		double t152 = a16*s2;
		double t154 = a23*s3;
		double t156 = a30*s4;
		double t158 = a37*s5;
		double t160 = a44*s6;
		double t162 = t148*a3+t150*a10+t152*a17+t154*a24+t156*a31+t158*a38+t160*a45;
		double t170 = t148*a4+t150*a11+t152*a18+t154*a25+t156*a32+t158*a39+t160*a46;
		double t178 = t148*a5+t150*a12+t152*a19+t154*a26+t156*a33+t158*a40+t160*a47;
		double t186 = t148*a6+t150*a13+t152*a20+t154*a27+t156*a34+t158*a41+t160*a48;
		double t187 = a3*a3;
		double t189 = a10*a10;
		double t191 = a17*a17;
		double t193 = a24*a24;
		double t195 = a31*a31;
		double t197 = a38*a38;
		double t199 = a45*a45;
		double t202 = a3*s0;
		double t204 = a10*s1;
		double t206 = a17*s2;
		double t208 = a24*s3;
		double t210 = a31*s4;
		double t212 = a38*s5;
		double t214 = a45*s6;
		double t216 = t202*a4+t204*a11+t206*a18+t208*a25+t210*a32+t212*a39+t214*a46;
		double t224 = t202*a5+t204*a12+t206*a19+t208*a26+t210*a33+t212*a40+t214*a47;
		double t232 = t202*a6+t204*a13+t206*a20+t208*a27+t210*a34+t212*a41+t214*a48;
		double t233 = a4*a4;
		double t235 = a11*a11;
		double t237 = a18*a18;
		double t239 = a25*a25;
		double t241 = a32*a32;
		double t243 = a39*a39;
		double t245 = a46*a46;
		double t248 = a4*s0;
		double t250 = a11*s1;
		double t252 = a18*s2;
		double t254 = a25*s3;
		double t256 = a32*s4;
		double t258 = a39*s5;
		double t260 = a46*s6;
		double t262 = t248*a5+t250*a12+t252*a19+t254*a26+t256*a33+t258*a40+t260*a47;
		double t270 = t248*a6+t250*a13+t252*a20+t254*a27+t256*a34+t258*a41+t260*a48;
		double t271 = a5*a5;
		double t273 = a12*a12;
		double t275 = a19*a19;
		double t277 = a26*a26;
		double t279 = a33*a33;
		double t281 = a40*a40;
		double t283 = a47*a47;
		double t300 = a5*s0*a6+a12*s1*a13+a19*s2*a20+a26*s3*a27+a33*s4*a34+a40*s5*a41+a47*s6*a48;
		double t301 = a6*a6;
		double t303 = a13*a13;
		double t305 = a20*a20;
		double t307 = a27*a27;
		double t309 = a34*a34;
		double t311 = a41*a41;
		double t313 = a48*a48;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6;
		*tmp_M_ptr++ = t30;
		*tmp_M_ptr++ = t38;
		*tmp_M_ptr++ = t46;
		*tmp_M_ptr++ = t54;
		*tmp_M_ptr++ = t62;
		*tmp_M_ptr++ = t70;
		*tmp_M_ptr++ = t71*s0+t73*s1+t75*s2+t77*s3+t79*s4+t81*s5+t83*s6;
		*tmp_M_ptr++ = t100;
		*tmp_M_ptr++ = t108;
		*tmp_M_ptr++ = t116;
		*tmp_M_ptr++ = t124;
		*tmp_M_ptr++ = t132;
		*tmp_M_ptr++ = t133*s0+t135*s1+t137*s2+t139*s3+t141*s4+t143*s5+t145*s6;
		*tmp_M_ptr++ = t162;
		*tmp_M_ptr++ = t170;
		*tmp_M_ptr++ = t178;
		*tmp_M_ptr++ = t186;
		*tmp_M_ptr++ = t187*s0+t189*s1+t191*s2+t193*s3+t195*s4+t197*s5+t199*s6;
		*tmp_M_ptr++ = t216;
		*tmp_M_ptr++ = t224;
		*tmp_M_ptr++ = t232;
		*tmp_M_ptr++ = t233*s0+t235*s1+t237*s2+t239*s3+t241*s4+t243*s5+t245*s6;
		*tmp_M_ptr++ = t262;
		*tmp_M_ptr++ = t270;
		*tmp_M_ptr++ = t271*s0+t273*s1+t275*s2+t277*s3+t279*s4+t281*s5+t283*s6;
		*tmp_M_ptr++ = t300;
		*tmp_M_ptr++ = t301*s0+t303*s1+t305*s2+t307*s3+t309*s4+t311*s5+t313*s6;

	}
	else if ( this->M_size == 8 ) {
		double D[8];
		double U[64];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 8);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63;
		double s0,s1,s2,s3,s4,s5,s6,s7;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;
		a49 = *tmp_M_ptr++;
		a50 = *tmp_M_ptr++;
		a51 = *tmp_M_ptr++;
		a52 = *tmp_M_ptr++;
		a53 = *tmp_M_ptr++;
		a54 = *tmp_M_ptr++;
		a55 = *tmp_M_ptr++;
		a56 = *tmp_M_ptr++;
		a57 = *tmp_M_ptr++;
		a58 = *tmp_M_ptr++;
		a59 = *tmp_M_ptr++;
		a60 = *tmp_M_ptr++;
		a61 = *tmp_M_ptr++;
		a62 = *tmp_M_ptr++;
		a63 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = log(*tmp_M_ptr++);
		s1 = log(*tmp_M_ptr++);
		s2 = log(*tmp_M_ptr++);
		s3 = log(*tmp_M_ptr++);
		s4 = log(*tmp_M_ptr++);
		s5 = log(*tmp_M_ptr++);
		s6 = log(*tmp_M_ptr++);
		s7 = log(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a8*a8;
		double t5 = a16*a16;
		double t7 = a24*a24;
		double t9 = a32*a32;
		double t11 = a40*a40;
		double t13 = a48*a48;
		double t15 = a56*a56;
		double t18 = a0*s0;
		double t20 = a8*s1;
		double t22 = a16*s2;
		double t24 = a24*s3;
		double t26 = a32*s4;
		double t28 = a40*s5;
		double t30 = a48*s6;
		double t32 = a56*s7;
		double t34 = t18*a1+t20*a9+t22*a17+t24*a25+t26*a33+t28*a41+t30*a49+t32*a57;
		double t43 = t18*a2+t20*a10+t22*a18+t24*a26+t26*a34+t28*a42+t30*a50+t32*a58;
		double t52 = t18*a3+t20*a11+t22*a19+t24*a27+t26*a35+t28*a43+t30*a51+t32*a59;
		double t61 = t18*a4+t20*a12+t22*a20+t24*a28+t26*a36+t28*a44+t30*a52+t32*a60;
		double t70 = t18*a5+t20*a13+t22*a21+t24*a29+t26*a37+t28*a45+t30*a53+t32*a61;
		double t79 = t18*a6+t20*a14+t22*a22+t24*a30+t26*a38+t28*a46+t30*a54+t32*a62;
		double t88 = t18*a7+t20*a15+t22*a23+t24*a31+t26*a39+t28*a47+t30*a55+t32*a63;
		double t89 = a1*a1;
		double t91 = a9*a9;
		double t93 = a17*a17;
		double t95 = a25*a25;
		double t97 = a33*a33;
		double t99 = a41*a41;
		double t101 = a49*a49;
		double t103 = a57*a57;
		double t106 = a1*s0;
		double t108 = a9*s1;
		double t110 = a17*s2;
		double t112 = a25*s3;
		double t114 = a33*s4;
		double t116 = a41*s5;
		double t118 = a49*s6;
		double t120 = a57*s7;
		double t122 = t106*a2+t108*a10+t110*a18+t112*a26+t114*a34+t116*a42+t118*a50+t120*a58;
		double t131 = t106*a3+t108*a11+t110*a19+t112*a27+t114*a35+t116*a43+t118*a51+t120*a59;
		double t140 = t106*a4+t108*a12+t110*a20+t112*a28+t114*a36+t116*a44+t118*a52+t120*a60;
		double t149 = t106*a5+t108*a13+t110*a21+t112*a29+t114*a37+t116*a45+t118*a53+t120*a61;
		double t158 = t106*a6+t108*a14+t110*a22+t112*a30+t114*a38+t116*a46+t118*a54+t120*a62;
		double t167 = t106*a7+t108*a15+t110*a23+t112*a31+t114*a39+t116*a47+t118*a55+t120*a63;
		double t168 = a2*a2;
		double t170 = a10*a10;
		double t172 = a18*a18;
		double t174 = a26*a26;
		double t176 = a34*a34;
		double t178 = a42*a42;
		double t180 = a50*a50;
		double t182 = a58*a58;
		double t185 = a2*s0;
		double t187 = a10*s1;
		double t189 = a18*s2;
		double t191 = a26*s3;
		double t193 = a34*s4;
		double t195 = a42*s5;
		double t197 = a50*s6;
		double t199 = a58*s7;
		double t201 = t185*a3+t187*a11+t189*a19+t191*a27+t193*a35+t195*a43+t197*a51+t199*a59;
		double t210 = t185*a4+t187*a12+t189*a20+t191*a28+t193*a36+t195*a44+t197*a52+t199*a60;
		double t219 = t185*a5+t187*a13+t189*a21+t191*a29+t193*a37+t195*a45+t197*a53+t199*a61;
		double t228 = t185*a6+t187*a14+t189*a22+t191*a30+t193*a38+t195*a46+t197*a54+t199*a62;
		double t237 = t185*a7+t187*a15+t189*a23+t191*a31+t193*a39+t195*a47+t197*a55+t199*a63;
		double t238 = a3*a3;
		double t240 = a11*a11;
		double t242 = a19*a19;
		double t244 = a27*a27;
		double t246 = a35*a35;
		double t248 = a43*a43;
		double t250 = a51*a51;
		double t252 = a59*a59;
		double t255 = a3*s0;
		double t257 = a11*s1;
		double t259 = a19*s2;
		double t261 = a27*s3;
		double t263 = a35*s4;
		double t265 = a43*s5;
		double t267 = a51*s6;
		double t269 = a59*s7;
		double t271 = t255*a4+t257*a12+t259*a20+t261*a28+t263*a36+t265*a44+t267*a52+t269*a60;
		double t280 = t255*a5+t257*a13+t259*a21+t261*a29+t263*a37+t265*a45+t267*a53+t269*a61;
		double t289 = t255*a6+t257*a14+t259*a22+t261*a30+t263*a38+t265*a46+t267*a54+t269*a62;
		double t298 = t255*a7+t257*a15+t259*a23+t261*a31+t263*a39+t265*a47+t267*a55+t269*a63;
		double t299 = a4*a4;
		double t301 = a12*a12;
		double t303 = a20*a20;
		double t305 = a28*a28;
		double t307 = a36*a36;
		double t309 = a44*a44;
		double t311 = a52*a52;
		double t313 = a60*a60;
		double t316 = a4*s0;
		double t318 = a12*s1;
		double t320 = a20*s2;
		double t322 = a28*s3;
		double t324 = a36*s4;
		double t326 = a44*s5;
		double t328 = a52*s6;
		double t330 = a60*s7;
		double t332 = t316*a5+t318*a13+t320*a21+t322*a29+t324*a37+t326*a45+t328*a53+t330*a61;
		double t341 = t316*a6+t318*a14+t320*a22+t322*a30+t324*a38+t326*a46+t328*a54+t330*a62;
		double t350 = t316*a7+t318*a15+t320*a23+t322*a31+t324*a39+t326*a47+t328*a55+t330*a63;
		double t351 = a5*a5;
		double t353 = a13*a13;
		double t355 = a21*a21;
		double t357 = a29*a29;
		double t359 = a37*a37;
		double t361 = a45*a45;
		double t363 = a53*a53;
		double t365 = a61*a61;
		double t368 = a5*s0;
		double t370 = a13*s1;
		double t372 = a21*s2;
		double t374 = a29*s3;
		double t376 = a37*s4;
		double t378 = a45*s5;
		double t380 = a53*s6;
		double t382 = a61*s7;
		double t384 = t368*a6+t370*a14+t372*a22+t374*a30+t376*a38+t378*a46+t380*a54+t382*a62;
		double t393 = t368*a7+t370*a15+t372*a23+t374*a31+t376*a39+t378*a47+t380*a55+t382*a63;
		double t394 = a6*a6;
		double t396 = a14*a14;
		double t398 = a22*a22;
		double t400 = a30*a30;
		double t402 = a38*a38;
		double t404 = a46*a46;
		double t406 = a54*a54;
		double t408 = a62*a62;
		double t427 = a6*s0*a7+a14*s1*a15+a22*s2*a23+a30*s3*a31+a38*s4*a39+a46*s5*a47+a54*s6*a55+a62*s7*a63;
		double t428 = a7*a7;
		double t430 = a15*a15;
		double t432 = a23*a23;
		double t434 = a31*a31;
		double t436 = a39*a39;
		double t438 = a47*a47;
		double t440 = a55*a55;
		double t442 = a63*a63;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6+t15*s7;
		*tmp_M_ptr++ = t34;
		*tmp_M_ptr++ = t43;
		*tmp_M_ptr++ = t52;
		*tmp_M_ptr++ = t61;
		*tmp_M_ptr++ = t70;
		*tmp_M_ptr++ = t79;
		*tmp_M_ptr++ = t88;
		*tmp_M_ptr++ = t89*s0+t91*s1+t93*s2+t95*s3+t97*s4+t99*s5+t101*s6+t103*s7;
		*tmp_M_ptr++ = t122;
		*tmp_M_ptr++ = t131;
		*tmp_M_ptr++ = t140;
		*tmp_M_ptr++ = t149;
		*tmp_M_ptr++ = t158;
		*tmp_M_ptr++ = t167;
		*tmp_M_ptr++ = t168*s0+t170*s1+t172*s2+t174*s3+t176*s4+t178*s5+t180*s6+t182*s7;
		*tmp_M_ptr++ = t201;
		*tmp_M_ptr++ = t210;
		*tmp_M_ptr++ = t219;
		*tmp_M_ptr++ = t228;
		*tmp_M_ptr++ = t237;
		*tmp_M_ptr++ = t238*s0+t240*s1+t242*s2+t244*s3+t246*s4+t248*s5+t250*s6+t252*s7;
		*tmp_M_ptr++ = t271;
		*tmp_M_ptr++ = t280;
		*tmp_M_ptr++ = t289;
		*tmp_M_ptr++ = t298;
		*tmp_M_ptr++ = t299*s0+t301*s1+t303*s2+t305*s3+t307*s4+t309*s5+t311*s6+t313*s7;
		*tmp_M_ptr++ = t332;
		*tmp_M_ptr++ = t341;
		*tmp_M_ptr++ = t350;
		*tmp_M_ptr++ = t351*s0+t353*s1+t355*s2+t357*s3+t359*s4+t361*s5+t363*s6+t365*s7;
		*tmp_M_ptr++ = t384;
		*tmp_M_ptr++ = t393;
		*tmp_M_ptr++ = t394*s0+t396*s1+t398*s2+t400*s3+t402*s4+t404*s5+t406*s6+t408*s7;
		*tmp_M_ptr++ = t427;
		*tmp_M_ptr++ = t428*s0+t430*s1+t432*s2+t434*s3+t436*s4+t438*s5+t440*s6+t442*s7;
	}
	else if ( this->M_size == 9 ) {
		double D[9];
		double U[81];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 9);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80;
		double s0,s1,s2,s3,s4,s5,s6,s7,s8;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;
		a49 = *tmp_M_ptr++;
		a50 = *tmp_M_ptr++;
		a51 = *tmp_M_ptr++;
		a52 = *tmp_M_ptr++;
		a53 = *tmp_M_ptr++;
		a54 = *tmp_M_ptr++;
		a55 = *tmp_M_ptr++;
		a56 = *tmp_M_ptr++;
		a57 = *tmp_M_ptr++;
		a58 = *tmp_M_ptr++;
		a59 = *tmp_M_ptr++;
		a60 = *tmp_M_ptr++;
		a61 = *tmp_M_ptr++;
		a62 = *tmp_M_ptr++;
		a63 = *tmp_M_ptr++;
		a64 = *tmp_M_ptr++;
		a65 = *tmp_M_ptr++;
		a66 = *tmp_M_ptr++;
		a67 = *tmp_M_ptr++;
		a68 = *tmp_M_ptr++;
		a69 = *tmp_M_ptr++;
		a70 = *tmp_M_ptr++;
		a71 = *tmp_M_ptr++;
		a72 = *tmp_M_ptr++;
		a73 = *tmp_M_ptr++;
		a74 = *tmp_M_ptr++;
		a75 = *tmp_M_ptr++;
		a76 = *tmp_M_ptr++;
		a77 = *tmp_M_ptr++;
		a78 = *tmp_M_ptr++;
		a79 = *tmp_M_ptr++;
		a80 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = log(*tmp_M_ptr++);
		s1 = log(*tmp_M_ptr++);
		s2 = log(*tmp_M_ptr++);
		s3 = log(*tmp_M_ptr++);
		s4 = log(*tmp_M_ptr++);
		s5 = log(*tmp_M_ptr++);
		s6 = log(*tmp_M_ptr++);
		s7 = log(*tmp_M_ptr++);
		s8 = log(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a9*a9;
		double t5 = a18*a18;
		double t7 = a27*a27;
		double t9 = a36*a36;
		double t11 = a45*a45;
		double t13 = a54*a54;
		double t15 = a63*a63;
		double t17 = a72*a72;
		double t20 = a0*s0;
		double t22 = a9*s1;
		double t24 = a18*s2;
		double t26 = a27*s3;
		double t28 = a36*s4;
		double t30 = a45*s5;
		double t32 = a54*s6;
		double t34 = a63*s7;
		double t36 = a72*s8;
		double t38 = t20*a1+t22*a10+t24*a19+t26*a28+t28*a37+t30*a46+t32*a55+t34*a64+t36*a73;
		double t48 = t20*a2+t22*a11+t24*a20+t26*a29+t28*a38+t30*a47+t32*a56+t34*a65+t36*a74;
		double t58 = t20*a3+t22*a12+t24*a21+t26*a30+t28*a39+t30*a48+t32*a57+t34*a66+t36*a75;
		double t68 = t20*a4+t22*a13+t24*a22+t26*a31+t28*a40+t30*a49+t32*a58+t34*a67+t36*a76;
		double t78 = t20*a5+t22*a14+t24*a23+t26*a32+t28*a41+t30*a50+t32*a59+t34*a68+t36*a77;
		double t88 = t20*a6+t22*a15+t24*a24+t26*a33+t28*a42+t30*a51+t32*a60+t34*a69+t36*a78;
		double t98 = t20*a7+t22*a16+t24*a25+t26*a34+t28*a43+t30*a52+t32*a61+t34*a70+t36*a79;
		double t108 = t20*a8+t22*a17+t24*a26+t26*a35+t28*a44+t30*a53+t32*a62+t34*a71+t36*a80;
		double t109 = a1*a1;
		double t111 = a10*a10;
		double t113 = a19*a19;
		double t115 = a28*a28;
		double t117 = a37*a37;
		double t119 = a46*a46;
		double t121 = a55*a55;
		double t123 = a64*a64;
		double t125 = a73*a73;
		double t128 = a1*s0;
		double t130 = a10*s1;
		double t132 = a19*s2;
		double t134 = a28*s3;
		double t136 = a37*s4;
		double t138 = a46*s5;
		double t140 = a55*s6;
		double t142 = a64*s7;
		double t144 = a73*s8;
		double t146 = t128*a2+t130*a11+t132*a20+t134*a29+t136*a38+t138*a47+t140*a56+t142*a65+t144*a74;
		double t156 = t128*a3+t130*a12+t132*a21+t134*a30+t136*a39+t138*a48+t140*a57+t142*a66+t144*a75;
		double t166 = t128*a4+t130*a13+t132*a22+t134*a31+t136*a40+t138*a49+t140*a58+t142*a67+t144*a76;
		double t176 = t128*a5+t130*a14+t132*a23+t134*a32+t136*a41+t138*a50+t140*a59+t142*a68+t144*a77;
		double t186 = t128*a6+t130*a15+t132*a24+t134*a33+t136*a42+t138*a51+t140*a60+t142*a69+t144*a78;
		double t196 = t128*a7+t130*a16+t132*a25+t134*a34+t136*a43+t138*a52+t140*a61+t142*a70+t144*a79;
		double t206 = t128*a8+t130*a17+t132*a26+t134*a35+t136*a44+t138*a53+t140*a62+t142*a71+t144*a80;
		double t207 = a2*a2;
		double t209 = a11*a11;
		double t211 = a20*a20;
		double t213 = a29*a29;
		double t215 = a38*a38;
		double t217 = a47*a47;
		double t219 = a56*a56;
		double t221 = a65*a65;
		double t223 = a74*a74;
		double t226 = a2*s0;
		double t228 = a11*s1;
		double t230 = a20*s2;
		double t232 = a29*s3;
		double t234 = a38*s4;
		double t236 = a47*s5;
		double t238 = a56*s6;
		double t240 = a65*s7;
		double t242 = a74*s8;
		double t244 = t226*a3+t228*a12+t230*a21+t232*a30+t234*a39+t236*a48+t238*a57+t240*a66+t242*a75;
		double t254 = t226*a4+t228*a13+t230*a22+t232*a31+t234*a40+t236*a49+t238*a58+t240*a67+t242*a76;
		double t264 = t226*a5+t228*a14+t230*a23+t232*a32+t234*a41+t236*a50+t238*a59+t240*a68+t242*a77;
		double t274 = t226*a6+t228*a15+t230*a24+t232*a33+t234*a42+t236*a51+t238*a60+t240*a69+t242*a78;
		double t284 = t226*a7+t228*a16+t230*a25+t232*a34+t234*a43+t236*a52+t238*a61+t240*a70+t242*a79;
		double t294 = t226*a8+t228*a17+t230*a26+t232*a35+t234*a44+t236*a53+t238*a62+t240*a71+t242*a80;
		double t295 = a3*a3;
		double t297 = a12*a12;
		double t299 = a21*a21;
		double t301 = a30*a30;
		double t303 = a39*a39;
		double t305 = a48*a48;
		double t307 = a57*a57;
		double t309 = a66*a66;
		double t311 = a75*a75;
		double t314 = a3*s0;
		double t316 = a12*s1;
		double t318 = a21*s2;
		double t320 = a30*s3;
		double t322 = a39*s4;
		double t324 = a48*s5;
		double t326 = a57*s6;
		double t328 = a66*s7;
		double t330 = a75*s8;
		double t332 = t314*a4+t316*a13+t318*a22+t320*a31+t322*a40+t324*a49+t326*a58+t328*a67+t330*a76;
		double t342 = t314*a5+t316*a14+t318*a23+t320*a32+t322*a41+t324*a50+t326*a59+t328*a68+t330*a77;
		double t352 = t314*a6+t316*a15+t318*a24+t320*a33+t322*a42+t324*a51+t326*a60+t328*a69+t330*a78;
		double t362 = t314*a7+t316*a16+t318*a25+t320*a34+t322*a43+t324*a52+t326*a61+t328*a70+t330*a79;
		double t372 = t314*a8+t316*a17+t318*a26+t320*a35+t322*a44+t324*a53+t326*a62+t328*a71+t330*a80;
		double t373 = a4*a4;
		double t375 = a13*a13;
		double t377 = a22*a22;
		double t379 = a31*a31;
		double t381 = a40*a40;
		double t383 = a49*a49;
		double t385 = a58*a58;
		double t387 = a67*a67;
		double t389 = a76*a76;
		double t392 = a4*s0;
		double t394 = a13*s1;
		double t396 = a22*s2;
		double t398 = a31*s3;
		double t400 = a40*s4;
		double t402 = a49*s5;
		double t404 = a58*s6;
		double t406 = a67*s7;
		double t408 = a76*s8;
		double t410 = t392*a5+t394*a14+t396*a23+t398*a32+t400*a41+t402*a50+t404*a59+t406*a68+t408*a77;
		double t420 = t392*a6+t394*a15+t396*a24+t398*a33+t400*a42+t402*a51+t404*a60+t406*a69+t408*a78;
		double t430 = t392*a7+t394*a16+t396*a25+t398*a34+t400*a43+t402*a52+t404*a61+t406*a70+t408*a79;
		double t440 = t392*a8+t394*a17+t396*a26+t398*a35+t400*a44+t402*a53+t404*a62+t406*a71+t408*a80;
		double t441 = a5*a5;
		double t443 = a14*a14;
		double t445 = a23*a23;
		double t447 = a32*a32;
		double t449 = a41*a41;
		double t451 = a50*a50;
		double t453 = a59*a59;
		double t455 = a68*a68;
		double t457 = a77*a77;
		double t460 = a5*s0;
		double t462 = a14*s1;
		double t464 = a23*s2;
		double t466 = a32*s3;
		double t468 = a41*s4;
		double t470 = a50*s5;
		double t472 = a59*s6;
		double t474 = a68*s7;
		double t476 = a77*s8;
		double t478 = t460*a6+t462*a15+t464*a24+t466*a33+t468*a42+t470*a51+t472*a60+t474*a69+t476*a78;
		double t488 = t460*a7+t462*a16+t464*a25+t466*a34+t468*a43+t470*a52+t472*a61+t474*a70+t476*a79;
		double t498 = t460*a8+t462*a17+t464*a26+t466*a35+t468*a44+t470*a53+t472*a62+t474*a71+t476*a80;
		double t499 = a6*a6;
		double t501 = a15*a15;
		double t503 = a24*a24;
		double t505 = a33*a33;
		double t507 = a42*a42;
		double t509 = a51*a51;
		double t511 = a60*a60;
		double t513 = a69*a69;
		double t515 = a78*a78;
		double t518 = a6*s0;
		double t520 = a15*s1;
		double t522 = a24*s2;
		double t524 = a33*s3;
		double t526 = a42*s4;
		double t528 = a51*s5;
		double t530 = a60*s6;
		double t532 = a69*s7;
		double t534 = a78*s8;
		double t536 = t518*a7+t520*a16+t522*a25+t524*a34+t526*a43+t528*a52+t530*a61+t532*a70+t534*a79;
		double t546 = t518*a8+t520*a17+t522*a26+t524*a35+t526*a44+t528*a53+t530*a62+t532*a71+t534*a80;
		double t547 = a7*a7;
		double t549 = a16*a16;
		double t551 = a25*a25;
		double t553 = a34*a34;
		double t555 = a43*a43;
		double t557 = a52*a52;
		double t559 = a61*a61;
		double t561 = a70*a70;
		double t563 = a79*a79;
		double t584 = a7*s0*a8+a16*s1*a17+a25*s2*a26+a34*s3*a35+a43*s4*a44+a52*s5*a53+a61*s6*a62+a70*s7*a71+a79*s8*a80;
		double t585 = a8*a8;
		double t587 = a17*a17;
		double t589 = a26*a26;
		double t591 = a35*a35;
		double t593 = a44*a44;
		double t595 = a53*a53;
		double t597 = a62*a62;
		double t599 = a71*a71;
		double t601 = a80*a80;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6+t15*s7+t17*s8;
		*tmp_M_ptr++ = t38;
		*tmp_M_ptr++ = t48;
		*tmp_M_ptr++ = t58;
		*tmp_M_ptr++ = t68;
		*tmp_M_ptr++ = t78;
		*tmp_M_ptr++ = t88;
		*tmp_M_ptr++ = t98;
		*tmp_M_ptr++ = t108;
		*tmp_M_ptr++ = t109*s0+t111*s1+t113*s2+t115*s3+t117*s4+t119*s5+t121*s6+t123*s7+t125*s8;
		*tmp_M_ptr++ = t146;
		*tmp_M_ptr++ = t156;
		*tmp_M_ptr++ = t166;
		*tmp_M_ptr++ = t176;
		*tmp_M_ptr++ = t186;
		*tmp_M_ptr++ = t196;
		*tmp_M_ptr++ = t206;
		*tmp_M_ptr++ = t207*s0+t209*s1+t211*s2+t213*s3+t215*s4+t217*s5+t219*s6+t221*s7+t223*s8;
		*tmp_M_ptr++ = t244;
		*tmp_M_ptr++ = t254;
		*tmp_M_ptr++ = t264;
		*tmp_M_ptr++ = t274;
		*tmp_M_ptr++ = t284;
		*tmp_M_ptr++ = t294;
		*tmp_M_ptr++ = t295*s0+t297*s1+t299*s2+t301*s3+t303*s4+t305*s5+t307*s6+t309*s7+t311*s8;
		*tmp_M_ptr++ = t332;
		*tmp_M_ptr++ = t342;
		*tmp_M_ptr++ = t352;
		*tmp_M_ptr++ = t362;
		*tmp_M_ptr++ = t372;
		*tmp_M_ptr++ = t373*s0+t375*s1+t377*s2+t379*s3+t381*s4+t383*s5+t385*s6+t387*s7+t389*s8;
		*tmp_M_ptr++ = t410;
		*tmp_M_ptr++ = t420;
		*tmp_M_ptr++ = t430;
		*tmp_M_ptr++ = t440;
		*tmp_M_ptr++ = t441*s0+t443*s1+t445*s2+t447*s3+t449*s4+t451*s5+t453*s6+t455*s7+t457*s8;
		*tmp_M_ptr++ = t478;
		*tmp_M_ptr++ = t488;
		*tmp_M_ptr++ = t498;
		*tmp_M_ptr++ = t499*s0+t501*s1+t503*s2+t505*s3+t507*s4+t509*s5+t511*s6+t513*s7+t515*s8;
		*tmp_M_ptr++ = t536;
		*tmp_M_ptr++ = t546;
		*tmp_M_ptr++ = t547*s0+t549*s1+t551*s2+t553*s3+t555*s4+t557*s5+t559*s6+t561*s7+t563*s8;
		*tmp_M_ptr++ = t584;
		*tmp_M_ptr++ = t585*s0+t587*s1+t589*s2+t591*s3+t593*s4+t595*s5+t597*s6+t599*s7+t601*s8;
	}
	else if ( this->M_size == 10 ) {
		double D[10];
		double U[100];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 10);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80,a81,a82,a83,a84,a85,a86,a87,a88,a89,a90,a91,a92,a93,a94,a95,a96,a97,a98,a99;
		double s0,s1,s2,s3,s4,s5,s6,s7,s8,s9;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;
		a49 = *tmp_M_ptr++;
		a50 = *tmp_M_ptr++;
		a51 = *tmp_M_ptr++;
		a52 = *tmp_M_ptr++;
		a53 = *tmp_M_ptr++;
		a54 = *tmp_M_ptr++;
		a55 = *tmp_M_ptr++;
		a56 = *tmp_M_ptr++;
		a57 = *tmp_M_ptr++;
		a58 = *tmp_M_ptr++;
		a59 = *tmp_M_ptr++;
		a60 = *tmp_M_ptr++;
		a61 = *tmp_M_ptr++;
		a62 = *tmp_M_ptr++;
		a63 = *tmp_M_ptr++;
		a64 = *tmp_M_ptr++;
		a65 = *tmp_M_ptr++;
		a66 = *tmp_M_ptr++;
		a67 = *tmp_M_ptr++;
		a68 = *tmp_M_ptr++;
		a69 = *tmp_M_ptr++;
		a70 = *tmp_M_ptr++;
		a71 = *tmp_M_ptr++;
		a72 = *tmp_M_ptr++;
		a73 = *tmp_M_ptr++;
		a74 = *tmp_M_ptr++;
		a75 = *tmp_M_ptr++;
		a76 = *tmp_M_ptr++;
		a77 = *tmp_M_ptr++;
		a78 = *tmp_M_ptr++;
		a79 = *tmp_M_ptr++;
		a80 = *tmp_M_ptr++;
		a81 = *tmp_M_ptr++;
		a82 = *tmp_M_ptr++;
		a83 = *tmp_M_ptr++;
		a84 = *tmp_M_ptr++;
		a85 = *tmp_M_ptr++;
		a86 = *tmp_M_ptr++;
		a87 = *tmp_M_ptr++;
		a88 = *tmp_M_ptr++;
		a89 = *tmp_M_ptr++;
		a90 = *tmp_M_ptr++;
		a91 = *tmp_M_ptr++;
		a92 = *tmp_M_ptr++;
		a93 = *tmp_M_ptr++;
		a94 = *tmp_M_ptr++;
		a95 = *tmp_M_ptr++;
		a96 = *tmp_M_ptr++;
		a97 = *tmp_M_ptr++;
		a98 = *tmp_M_ptr++;
		a99 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = log(*tmp_M_ptr++);
		s1 = log(*tmp_M_ptr++);
		s2 = log(*tmp_M_ptr++);
		s3 = log(*tmp_M_ptr++);
		s4 = log(*tmp_M_ptr++);
		s5 = log(*tmp_M_ptr++);
		s6 = log(*tmp_M_ptr++);
		s7 = log(*tmp_M_ptr++);
		s8 = log(*tmp_M_ptr++);
		s9 = log(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a10*a10;
		double t5 = a20*a20;
		double t7 = a30*a30;
		double t9 = a40*a40;
		double t11 = a50*a50;
		double t13 = a60*a60;
		double t15 = a70*a70;
		double t17 = a80*a80;
		double t19 = a90*a90;
		double t22 = a0*s0;
		double t24 = a10*s1;
		double t26 = a20*s2;
		double t28 = a30*s3;
		double t30 = a40*s4;
		double t32 = a50*s5;
		double t34 = a60*s6;
		double t36 = a70*s7;
		double t38 = a80*s8;
		double t40 = a90*s9;
		double t42 = t22*a1+t24*a11+t26*a21+t28*a31+t30*a41+t32*a51+t34*a61+t36*a71+t38*a81+t40*a91;
		double t53 = t22*a2+t24*a12+t26*a22+t28*a32+t30*a42+t32*a52+t34*a62+t36*a72+t38*a82+t40*a92;
		double t64 = t22*a3+t24*a13+t26*a23+t28*a33+t30*a43+t32*a53+t34*a63+t36*a73+t38*a83+t40*a93;
		double t75 = t22*a4+t24*a14+t26*a24+t28*a34+t30*a44+t32*a54+t34*a64+t36*a74+t38*a84+t40*a94;
		double t86 = t22*a5+t24*a15+t26*a25+t28*a35+t30*a45+t32*a55+t34*a65+t36*a75+t38*a85+t40*a95;
		double t97 = t22*a6+t24*a16+t26*a26+t28*a36+t30*a46+t32*a56+t34*a66+t36*a76+t38*a86+t40*a96;
		double t108 = t22*a7+t24*a17+t26*a27+t28*a37+t30*a47+t32*a57+t34*a67+t36*a77+t38*a87+t40*a97;
		double t119 = t22*a8+t24*a18+t26*a28+t28*a38+t30*a48+t32*a58+t34*a68+t36*a78+t38*a88+t40*a98;
		double t130 = t22*a9+t24*a19+t26*a29+t28*a39+t30*a49+t32*a59+t34*a69+t36*a79+t38*a89+t40*a99;
		double t131 = a1*a1;
		double t133 = a11*a11;
		double t135 = a21*a21;
		double t137 = a31*a31;
		double t139 = a41*a41;
		double t141 = a51*a51;
		double t143 = a61*a61;
		double t145 = a71*a71;
		double t147 = a81*a81;
		double t149 = a91*a91;
		double t152 = a1*s0;
		double t154 = a11*s1;
		double t156 = a21*s2;
		double t158 = a31*s3;
		double t160 = a41*s4;
		double t162 = a51*s5;
		double t164 = a61*s6;
		double t166 = a71*s7;
		double t168 = a81*s8;
		double t170 = a91*s9;
		double t172 = t152*a2+t154*a12+t156*a22+t158*a32+t160*a42+t162*a52+t164*a62+t166*a72+t168*a82+t170*a92;
		double t183 = t152*a3+t154*a13+t156*a23+t158*a33+t160*a43+t162*a53+t164*a63+t166*a73+t168*a83+t170*a93;
		double t194 = t152*a4+t154*a14+t156*a24+t158*a34+t160*a44+t162*a54+t164*a64+t166*a74+t168*a84+t170*a94;
		double t205 = t152*a5+t154*a15+t156*a25+t158*a35+t160*a45+t162*a55+t164*a65+t166*a75+t168*a85+t170*a95;
		double t216 = t152*a6+t154*a16+t156*a26+t158*a36+t160*a46+t162*a56+t164*a66+t166*a76+t168*a86+t170*a96;
		double t227 = t152*a7+t154*a17+t156*a27+t158*a37+t160*a47+t162*a57+t164*a67+t166*a77+t168*a87+t170*a97;
		double t238 = t152*a8+t154*a18+t156*a28+t158*a38+t160*a48+t162*a58+t164*a68+t166*a78+t168*a88+t170*a98;
		double t249 = t152*a9+t154*a19+t156*a29+t158*a39+t160*a49+t162*a59+t164*a69+t166*a79+t168*a89+t170*a99;
		double t250 = a2*a2;
		double t252 = a12*a12;
		double t254 = a22*a22;
		double t256 = a32*a32;
		double t258 = a42*a42;
		double t260 = a52*a52;
		double t262 = a62*a62;
		double t264 = a72*a72;
		double t266 = a82*a82;
		double t268 = a92*a92;
		double t271 = a2*s0;
		double t273 = a12*s1;
		double t275 = a22*s2;
		double t277 = a32*s3;
		double t279 = a42*s4;
		double t281 = a52*s5;
		double t283 = a62*s6;
		double t285 = a72*s7;
		double t287 = a82*s8;
		double t289 = a92*s9;
		double t291 = t271*a3+t273*a13+t275*a23+t277*a33+t279*a43+t281*a53+t283*a63+t285*a73+t287*a83+t289*a93;
		double t302 = t271*a4+t273*a14+t275*a24+t277*a34+t279*a44+t281*a54+t283*a64+t285*a74+t287*a84+t289*a94;
		double t313 = t271*a5+t273*a15+t275*a25+t277*a35+t279*a45+t281*a55+t283*a65+t285*a75+t287*a85+t289*a95;
		double t324 = t271*a6+t273*a16+t275*a26+t277*a36+t279*a46+t281*a56+t283*a66+t285*a76+t287*a86+t289*a96;
		double t335 = t271*a7+t273*a17+t275*a27+t277*a37+t279*a47+t281*a57+t283*a67+t285*a77+t287*a87+t289*a97;
		double t346 = t271*a8+t273*a18+t275*a28+t277*a38+t279*a48+t281*a58+t283*a68+t285*a78+t287*a88+t289*a98;
		double t357 = t271*a9+t273*a19+t275*a29+t277*a39+t279*a49+t281*a59+t283*a69+t285*a79+t287*a89+t289*a99;
		double t358 = a3*a3;
		double t360 = a13*a13;
		double t362 = a23*a23;
		double t364 = a33*a33;
		double t366 = a43*a43;
		double t368 = a53*a53;
		double t370 = a63*a63;
		double t372 = a73*a73;
		double t374 = a83*a83;
		double t376 = a93*a93;
		double t379 = a3*s0;
		double t381 = a13*s1;
		double t383 = a23*s2;
		double t385 = a33*s3;
		double t387 = a43*s4;
		double t389 = a53*s5;
		double t391 = a63*s6;
		double t393 = a73*s7;
		double t395 = a83*s8;
		double t397 = a93*s9;
		double t399 = t379*a4+t381*a14+t383*a24+t385*a34+t387*a44+t389*a54+t391*a64+t393*a74+t395*a84+t397*a94;
		double t410 = t379*a5+t381*a15+t383*a25+t385*a35+t387*a45+t389*a55+t391*a65+t393*a75+t395*a85+t397*a95;
		double t421 = t379*a6+t381*a16+t383*a26+t385*a36+t387*a46+t389*a56+t391*a66+t393*a76+t395*a86+t397*a96;
		double t432 = t379*a7+t381*a17+t383*a27+t385*a37+t387*a47+t389*a57+t391*a67+t393*a77+t395*a87+t397*a97;
		double t443 = t379*a8+t381*a18+t383*a28+t385*a38+t387*a48+t389*a58+t391*a68+t393*a78+t395*a88+t397*a98;
		double t454 = t379*a9+t381*a19+t383*a29+t385*a39+t387*a49+t389*a59+t391*a69+t393*a79+t395*a89+t397*a99;
		double t455 = a4*a4;
		double t457 = a14*a14;
		double t459 = a24*a24;
		double t461 = a34*a34;
		double t463 = a44*a44;
		double t465 = a54*a54;
		double t467 = a64*a64;
		double t469 = a74*a74;
		double t471 = a84*a84;
		double t473 = a94*a94;
		double t476 = a4*s0;
		double t478 = a14*s1;
		double t480 = a24*s2;
		double t482 = a34*s3;
		double t484 = a44*s4;
		double t486 = a54*s5;
		double t488 = a64*s6;
		double t490 = a74*s7;
		double t492 = a84*s8;
		double t494 = a94*s9;
		double t496 = t476*a5+t478*a15+t480*a25+t482*a35+t484*a45+t486*a55+t488*a65+t490*a75+t492*a85+t494*a95;
		double t507 = t476*a6+t478*a16+t480*a26+t482*a36+t484*a46+t486*a56+t488*a66+t490*a76+t492*a86+t494*a96;
		double t518 = t476*a7+t478*a17+t480*a27+t482*a37+t484*a47+t486*a57+t488*a67+t490*a77+t492*a87+t494*a97;
		double t529 = t476*a8+t478*a18+t480*a28+t482*a38+t484*a48+t486*a58+t488*a68+t490*a78+t492*a88+t494*a98;
		double t540 = t476*a9+t478*a19+t480*a29+t482*a39+t484*a49+t486*a59+t488*a69+t490*a79+t492*a89+t494*a99;
		double t541 = a5*a5;
		double t543 = a15*a15;
		double t545 = a25*a25;
		double t547 = a35*a35;
		double t549 = a45*a45;
		double t551 = a55*a55;
		double t553 = a65*a65;
		double t555 = a75*a75;
		double t557 = a85*a85;
		double t559 = a95*a95;
		double t562 = a5*s0;
		double t564 = a15*s1;
		double t566 = a25*s2;
		double t568 = a35*s3;
		double t570 = a45*s4;
		double t572 = a55*s5;
		double t574 = a65*s6;
		double t576 = a75*s7;
		double t578 = a85*s8;
		double t580 = a95*s9;
		double t582 = t562*a6+t564*a16+t566*a26+t568*a36+t570*a46+t572*a56+t574*a66+t576*a76+t578*a86+t580*a96;
		double t593 = t562*a7+t564*a17+t566*a27+t568*a37+t570*a47+t572*a57+t574*a67+t576*a77+t578*a87+t580*a97;
		double t604 = t562*a8+t564*a18+t566*a28+t568*a38+t570*a48+t572*a58+t574*a68+t576*a78+t578*a88+t580*a98;
		double t615 = t562*a9+t564*a19+t566*a29+t568*a39+t570*a49+t572*a59+t574*a69+t576*a79+t578*a89+t580*a99;
		double t616 = a6*a6;
		double t618 = a16*a16;
		double t620 = a26*a26;
		double t622 = a36*a36;
		double t624 = a46*a46;
		double t626 = a56*a56;
		double t628 = a66*a66;
		double t630 = a76*a76;
		double t632 = a86*a86;
		double t634 = a96*a96;
		double t637 = a6*s0;
		double t639 = a16*s1;
		double t641 = a26*s2;
		double t643 = a36*s3;
		double t645 = a46*s4;
		double t647 = a56*s5;
		double t649 = a66*s6;
		double t651 = a76*s7;
		double t653 = a86*s8;
		double t655 = a96*s9;
		double t657 = t637*a7+t639*a17+t641*a27+t643*a37+t645*a47+t647*a57+t649*a67+t651*a77+t653*a87+t655*a97;
		double t668 = t637*a8+t639*a18+t641*a28+t643*a38+t645*a48+t647*a58+t649*a68+t651*a78+t653*a88+t655*a98;
		double t679 = t637*a9+t639*a19+t641*a29+t643*a39+t645*a49+t647*a59+t649*a69+t651*a79+t653*a89+t655*a99;
		double t680 = a7*a7;
		double t682 = a17*a17;
		double t684 = a27*a27;
		double t686 = a37*a37;
		double t688 = a47*a47;
		double t690 = a57*a57;
		double t692 = a67*a67;
		double t694 = a77*a77;
		double t696 = a87*a87;
		double t698 = a97*a97;
		double t701 = a7*s0;
		double t703 = a17*s1;
		double t705 = a27*s2;
		double t707 = a37*s3;
		double t709 = a47*s4;
		double t711 = a57*s5;
		double t713 = a67*s6;
		double t715 = a77*s7;
		double t717 = a87*s8;
		double t719 = a97*s9;
		double t721 = t701*a8+t703*a18+t705*a28+t707*a38+t709*a48+t711*a58+t713*a68+t715*a78+t717*a88+t719*a98;
		double t732 = t701*a9+t703*a19+t705*a29+t707*a39+t709*a49+t711*a59+t713*a69+t715*a79+t717*a89+t719*a99;
		double t733 = a8*a8;
		double t735 = a18*a18;
		double t737 = a28*a28;
		double t739 = a38*a38;
		double t741 = a48*a48;
		double t743 = a58*a58;
		double t745 = a68*a68;
		double t747 = a78*a78;
		double t749 = a88*a88;
		double t751 = a98*a98;
		double t774 = a8*s0*a9+a18*s1*a19+a28*s2*a29+a38*s3*a39+a48*s4*a49+a58*s5*a59+a68*s6*a69+a78*s7*a79+a88*s8*a89+a98*s9*a99;
		double t775 = a9*a9;
		double t777 = a19*a19;
		double t779 = a29*a29;
		double t781 = a39*a39;
		double t783 = a49*a49;
		double t785 = a59*a59;
		double t787 = a69*a69;
		double t789 = a79*a79;
		double t791 = a89*a89;
		double t793 = a99*a99;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6+t15*s7+t17*s8+t19*s9;
		*tmp_M_ptr++ = t42;
		*tmp_M_ptr++ = t53;
		*tmp_M_ptr++ = t64;
		*tmp_M_ptr++ = t75;
		*tmp_M_ptr++ = t86;
		*tmp_M_ptr++ = t97;
		*tmp_M_ptr++ = t108;
		*tmp_M_ptr++ = t119;
		*tmp_M_ptr++ = t130;
		*tmp_M_ptr++ = t131*s0+t133*s1+t135*s2+t137*s3+t139*s4+t141*s5+t143*s6+t145*s7+t147*s8+t149*s9;
		*tmp_M_ptr++ = t172;
		*tmp_M_ptr++ = t183;
		*tmp_M_ptr++ = t194;
		*tmp_M_ptr++ = t205;
		*tmp_M_ptr++ = t216;
		*tmp_M_ptr++ = t227;
		*tmp_M_ptr++ = t238;
		*tmp_M_ptr++ = t249;
		*tmp_M_ptr++ = t250*s0+t252*s1+t254*s2+t256*s3+t258*s4+t260*s5+t262*s6+t264*s7+t266*s8+t268*s9;
		*tmp_M_ptr++ = t291;
		*tmp_M_ptr++ = t302;
		*tmp_M_ptr++ = t313;
		*tmp_M_ptr++ = t324;
		*tmp_M_ptr++ = t335;
		*tmp_M_ptr++ = t346;
		*tmp_M_ptr++ = t357;
		*tmp_M_ptr++ = t358*s0+t360*s1+t362*s2+t364*s3+t366*s4+t368*s5+t370*s6+t372*s7+t374*s8+t376*s9;
		*tmp_M_ptr++ = t399;
		*tmp_M_ptr++ = t410;
		*tmp_M_ptr++ = t421;
		*tmp_M_ptr++ = t432;
		*tmp_M_ptr++ = t443;
		*tmp_M_ptr++ = t454;
		*tmp_M_ptr++ = t455*s0+t457*s1+t459*s2+t461*s3+t463*s4+t465*s5+t467*s6+t469*s7+t471*s8+t473*s9;
		*tmp_M_ptr++ = t496;
		*tmp_M_ptr++ = t507;
		*tmp_M_ptr++ = t518;
		*tmp_M_ptr++ = t529;
		*tmp_M_ptr++ = t540;
		*tmp_M_ptr++ = t541*s0+t543*s1+t545*s2+t547*s3+t549*s4+t551*s5+t553*s6+t555*s7+t557*s8+t559*s9;
		*tmp_M_ptr++ = t582;
		*tmp_M_ptr++ = t593;
		*tmp_M_ptr++ = t604;
		*tmp_M_ptr++ = t615;
		*tmp_M_ptr++ = t616*s0+t618*s1+t620*s2+t622*s3+t624*s4+t626*s5+t628*s6+t630*s7+t632*s8+t634*s9;
		*tmp_M_ptr++ = t657;
		*tmp_M_ptr++ = t668;
		*tmp_M_ptr++ = t679;
		*tmp_M_ptr++ = t680*s0+t682*s1+t684*s2+t686*s3+t688*s4+t690*s5+t692*s6+t694*s7+t696*s8+t698*s9;
		*tmp_M_ptr++ = t721;
		*tmp_M_ptr++ = t732;
		*tmp_M_ptr++ = t733*s0+t735*s1+t737*s2+t739*s3+t741*s4+t743*s5+t745*s6+t747*s7+t749*s8+t751*s9;
		*tmp_M_ptr++ = t774;
		*tmp_M_ptr++ = t775*s0+t777*s1+t779*s2+t781*s3+t783*s4+t785*s5+t787*s6+t789*s7+t791*s8+t793*s9;
	}
	else {
		CvMat* D = cvCreateMat(M_size, 1, CV_64F);
		CvMat* U = cvCreateMat(M_size, M_size, CV_64F);

		this->SVD(D, U);

		CvMat* U_t = cvCreateMat(M_size, M_size, CV_64F);
		CvMat* U_tU = cvCreateMat(M_size, M_size, CV_64F);

		cvTranspose(U, U_t);

		/* D = log(D) */
		Log_Diag(D);

		/* U_t = U_t * D */
		Mul_Diag(U_t, D);

		/* U_tU = U_t * U */
		cvMatMul(U_t, U, U_tU);

		dst_mat->SetData((double*)U_tU->data.ptr, M_size2);

		cvReleaseMat(&U_t);
		cvReleaseMat(&U_tU);

		cvReleaseMat(&D);
		cvReleaseMat(&U);
	}
}

void CCovarianceMatrix::Exp(CCovarianceMatrix *exp_mat)
{
	CCovarianceMatrix *dst_mat;

	if ( exp_mat )
		dst_mat = exp_mat;
	else
		dst_mat = this;

	if ( this->M_size == 1 ) {
		dst_mat->SM_ptr[0] = exp(SM_ptr[0]);
	}
	else if ( this->M_size == 2 ) {
		double D[2];
		double U[4];

		eigen2_decomposition(this->SM_ptr, U, D);

		double a0,a1,a2,a3;
		double s0,s1;

		a0 = U[0];
		a1 = U[1];
		a2 = U[2];
		a3 = U[3];

		s0 = exp(D[0]);
		s1 = exp(D[1]);

		double t1 = a0*a0;
		double t3 = a2*a2;
		double t10 = a0*s0*a1+a2*s1*a3;
		double t11 = a1*a1;
		double t13 = a3*a3;

		double* tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1;
		*tmp_M_ptr++ = t10;
		*tmp_M_ptr++ = t11*s0+t13*s1;
	}
	else if ( this->M_size == 3 ) {
		double D[3];
		double U[9];

		eigen3_decomposition(this->SM_ptr, U, D);

		double a0,a1,a2,a3,a4,a5,a6,a7,a8;
		double s0,s1,s2;

		a0 = U[0];
		a1 = U[1];
		a2 = U[2];
		a3 = U[3];
		a4 = U[4];
		a5 = U[5];
		a6 = U[6];
		a7 = U[7];
		a8 = U[8];

		s0 = exp(D[0]);
		s1 = exp(D[1]);
		s2 = exp(D[2]);

		double t1 = a0*a0;
		double t3 = a3*a3;
		double t5 = a6*a6;
		double t8 = a0*s0;
		double t10 = a3*s1;
		double t12 = a6*s2;
		double t14 = t8*a1+t10*a4+t12*a7;
		double t18 = t8*a2+t10*a5+t12*a8;
		double t19 = a1*a1;
		double t21 = a4*a4;
		double t23 = a7*a7;
		double t32 = a1*s0*a2+a4*s1*a5+a7*s2*a8;
		double t33 = a2*a2;
		double t35 = a5*a5;
		double t37 = a8*a8;

		double *tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2;
		*tmp_M_ptr++ = t14;
		*tmp_M_ptr++ = t18;
		*tmp_M_ptr++ = t19*s0+t21*s1+t23*s2;
		*tmp_M_ptr++ = t32;
		*tmp_M_ptr++ = t33*s0+t35*s1+t37*s2;
	}
	else if ( this->M_size == 4 ) {
		double D[4];
		double U[16];

		eigen4_decomposition(this->SM_ptr, U, D);

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
		double s0,s1,s2,s3;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = exp(*tmp_M_ptr++);
		s1 = exp(*tmp_M_ptr++);
		s2 = exp(*tmp_M_ptr++);
		s3 = exp(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a4*a4;
		double t5 = a8*a8;
		double t7 = a12*a12;
		double t10 = a0*s0;
		double t12 = a4*s1;
		double t14 = a8*s2;
		double t16 = a12*s3;
		double t18 = t10*a1+t12*a5+t14*a9+t16*a13;
		double t23 = t10*a2+t12*a6+t14*a10+t16*a14;
		double t28 = t10*a3+t12*a7+t14*a11+t16*a15;
		double t29 = a1*a1;
		double t31 = a5*a5;
		double t33 = a9*a9;
		double t35 = a13*a13;
		double t38 = a1*s0;
		double t40 = a5*s1;
		double t42 = a9*s2;
		double t44 = a13*s3;
		double t46 = t38*a2+t40*a6+t42*a10+t44*a14;
		double t51 = t38*a3+t40*a7+t42*a11+t44*a15;
		double t52 = a2*a2;
		double t54 = a6*a6;
		double t56 = a10*a10;
		double t58 = a14*a14;
		double t69 = a2*s0*a3+a6*s1*a7+a10*s2*a11+a14*s3*a15;
		double t70 = a3*a3;
		double t72 = a7*a7;
		double t74 = a11*a11;
		double t76 = a15*a15;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3;
		*tmp_M_ptr++ = t18;
		*tmp_M_ptr++ = t23;
		*tmp_M_ptr++ = t28;
		*tmp_M_ptr++ = t29*s0+t31*s1+t33*s2+t35*s3;
		*tmp_M_ptr++ = t46;
		*tmp_M_ptr++ = t51;
		*tmp_M_ptr++ = t52*s0+t54*s1+t56*s2+t58*s3;
		*tmp_M_ptr++ = t69;
		*tmp_M_ptr++ = t70*s0+t72*s1+t74*s2+t76*s3;
	}
	else if ( this->M_size == 5 ) {
		double D[5];
		double U[25];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 5);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24;
		double s0,s1,s2,s3,s4;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = exp(*tmp_M_ptr++);
		s1 = exp(*tmp_M_ptr++);
		s2 = exp(*tmp_M_ptr++);
		s3 = exp(*tmp_M_ptr++);
		s4 = exp(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a5*a5;
		double t5 = a10*a10;
		double t7 = a15*a15;
		double t9 = a20*a20;
		double t12 = a0*s0;
		double t14 = a5*s1;
		double t16 = a10*s2;
		double t18 = a15*s3;
		double t20 = a20*s4;
		double t22 = t12*a1+t14*a6+t16*a11+t18*a16+t20*a21;
		double t28 = t12*a2+t14*a7+t16*a12+t18*a17+t20*a22;
		double t34 = t12*a3+t14*a8+t16*a13+t18*a18+t20*a23;
		double t40 = t12*a4+t14*a9+t16*a14+t18*a19+t20*a24;
		double t41 = a1*a1;
		double t43 = a6*a6;
		double t45 = a11*a11;
		double t47 = a16*a16;
		double t49 = a21*a21;
		double t52 = a1*s0;
		double t54 = a6*s1;
		double t56 = a11*s2;
		double t58 = a16*s3;
		double t60 = a21*s4;
		double t62 = t52*a2+t54*a7+t56*a12+t58*a17+t60*a22;
		double t68 = t52*a3+t54*a8+t56*a13+t58*a18+t60*a23;
		double t74 = t52*a4+t54*a9+t56*a14+t58*a19+t60*a24;
		double t75 = a2*a2;
		double t77 = a7*a7;
		double t79 = a12*a12;
		double t81 = a17*a17;
		double t83 = a22*a22;
		double t86 = a2*s0;
		double t88 = a7*s1;
		double t90 = a12*s2;
		double t92 = a17*s3;
		double t94 = a22*s4;
		double t96 = t86*a3+t88*a8+t90*a13+t92*a18+t94*a23;
		double t102 = t86*a4+t88*a9+t90*a14+t92*a19+t94*a24;
		double t103 = a3*a3;
		double t105 = a8*a8;
		double t107 = a13*a13;
		double t109 = a18*a18;
		double t111 = a23*a23;
		double t124 = a3*s0*a4+a8*s1*a9+a13*s2*a14+a18*s3*a19+a23*s4*a24;
		double t125 = a4*a4;
		double t127 = a9*a9;
		double t129 = a14*a14;
		double t131 = a19*a19;
		double t133 = a24*a24;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4;
		*tmp_M_ptr++  = t22;
		*tmp_M_ptr++  = t28;
		*tmp_M_ptr++  = t34;
		*tmp_M_ptr++  = t40;
		*tmp_M_ptr++  = t41*s0+t43*s1+t45*s2+t47*s3+t49*s4;
		*tmp_M_ptr++  = t62;
		*tmp_M_ptr++  = t68;
		*tmp_M_ptr++  = t74;
		*tmp_M_ptr++  = t75*s0+t77*s1+t79*s2+t81*s3+t83*s4;
		*tmp_M_ptr++  = t96;
		*tmp_M_ptr++  = t102;
		*tmp_M_ptr++  = t103*s0+t105*s1+t107*s2+t109*s3+t111*s4;
		*tmp_M_ptr++  = t124;
		*tmp_M_ptr++  = t125*s0+t127*s1+t129*s2+t131*s3+t133*s4;
	}
	else if ( this->M_size == 6 ) {
		double D[6];
		double U[36];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 6);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35;
		double s0,s1,s2,s3,s4,s5;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = exp(*tmp_M_ptr++);
		s1 = exp(*tmp_M_ptr++);
		s2 = exp(*tmp_M_ptr++);
		s3 = exp(*tmp_M_ptr++);
		s4 = exp(*tmp_M_ptr++);
		s5 = exp(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a6*a6;
		double t5 = a12*a12;
		double t7 = a18*a18;
		double t9 = a24*a24;
		double t11 = a30*a30;
		double t14 = a0*s0;
		double t16 = a6*s1;
		double t18 = a12*s2;
		double t20 = a18*s3;
		double t22 = a24*s4;
		double t24 = a30*s5;
		double t26 = t14*a1+t16*a7+t18*a13+t20*a19+t22*a25+t24*a31;
		double t33 = t14*a2+t16*a8+t18*a14+t20*a20+t22*a26+t24*a32;
		double t40 = t14*a3+t16*a9+t18*a15+t20*a21+t22*a27+t24*a33;
		double t47 = t14*a4+10.0*t16+t18*a16+t20*a22+t22*a28+t24*a34;
		double t54 = t14*a5+t16*a11+t18*a17+t20*a23+t22*a29+t24*a35;
		double t55 = a1*a1;
		double t57 = a7*a7;
		double t59 = a13*a13;
		double t61 = a19*a19;
		double t63 = a25*a25;
		double t65 = a31*a31;
		double t68 = a1*s0;
		double t70 = a7*s1;
		double t72 = a13*s2;
		double t74 = a19*s3;
		double t76 = a25*s4;
		double t78 = a31*s5;
		double t80 = t68*a2+t70*a8+t72*a14+t74*a20+t76*a26+t78*a32;
		double t87 = t68*a3+t70*a9+t72*a15+t74*a21+t76*a27+t78*a33;
		double t94 = t68*a4+10.0*t70+t72*a16+t74*a22+t76*a28+t78*a34;
		double t101 = t68*a5+t70*a11+t72*a17+t74*a23+t76*a29+t78*a35;
		double t102 = a2*a2;
		double t104 = a8*a8;
		double t106 = a14*a14;
		double t108 = a20*a20;
		double t110 = a26*a26;
		double t112 = a32*a32;
		double t115 = a2*s0;
		double t117 = a8*s1;
		double t119 = a14*s2;
		double t121 = a20*s3;
		double t123 = a26*s4;
		double t125 = a32*s5;
		double t127 = t115*a3+t117*a9+t119*a15+t121*a21+t123*a27+t125*a33;
		double t134 = t115*a4+10.0*t117+t119*a16+t121*a22+t123*a28+t125*a34;
		double t141 = t115*a5+t117*a11+t119*a17+t121*a23+t123*a29+t125*a35;
		double t142 = a3*a3;
		double t144 = a9*a9;
		double t146 = a15*a15;
		double t148 = a21*a21;
		double t150 = a27*a27;
		double t152 = a33*a33;
		double t155 = a3*s0;
		double t157 = a9*s1;
		double t159 = a15*s2;
		double t161 = a21*s3;
		double t163 = a27*s4;
		double t165 = a33*s5;
		double t167 = t155*a4+10.0*t157+t159*a16+t161*a22+t163*a28+t165*a34;
		double t174 = t155*a5+t157*a11+t159*a17+t161*a23+t163*a29+t165*a35;
		double t175 = a4*a4;
		double t178 = a16*a16;
		double t180 = a22*a22;
		double t182 = a28*a28;
		double t184 = a34*a34;
		double t199 = a4*s0*a5+10.0*a11*s1+a16*s2*a17+a22*s3*a23+a28*s4*a29+a34*s5*a35;
		double t200 = a5*a5;
		double t202 = a11*a11;
		double t204 = a17*a17;
		double t206 = a23*a23;
		double t208 = a29*a29;
		double t210 = a35*a35;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5;
		*tmp_M_ptr++ = t26;
		*tmp_M_ptr++ = t33;
		*tmp_M_ptr++ = t40;
		*tmp_M_ptr++ = t47;
		*tmp_M_ptr++ = t54;
		*tmp_M_ptr++ = t55*s0+t57*s1+t59*s2+t61*s3+t63*s4+t65*s5;
		*tmp_M_ptr++ = t80;
		*tmp_M_ptr++ = t87;
		*tmp_M_ptr++ = t94;
		*tmp_M_ptr++ = t101;
		*tmp_M_ptr++ = t102*s0+t104*s1+t106*s2+t108*s3+t110*s4+t112*s5;
		*tmp_M_ptr++ = t127;
		*tmp_M_ptr++ = t134;
		*tmp_M_ptr++ = t141;
		*tmp_M_ptr++ = t142*s0+t144*s1+t146*s2+t148*s3+t150*s4+t152*s5;
		*tmp_M_ptr++ = t167;
		*tmp_M_ptr++ = t174;
		*tmp_M_ptr++ = t175*s0+100.0*s1+t178*s2+t180*s3+t182*s4+t184*s5;
		*tmp_M_ptr++ = t199;
		*tmp_M_ptr++ = t200*s0+t202*s1+t204*s2+t206*s3+t208*s4+t210*s5;
	}
	else if ( this->M_size == 7 ) {
		double D[7];
		double U[49];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 7);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double s0,s1,s2,s3,s4,s5,s6;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = exp(*tmp_M_ptr++);
		s1 = exp(*tmp_M_ptr++);
		s2 = exp(*tmp_M_ptr++);
		s3 = exp(*tmp_M_ptr++);
		s4 = exp(*tmp_M_ptr++);
		s5 = exp(*tmp_M_ptr++);
		s6 = exp(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a7*a7;
		double t5 = a14*a14;
		double t7 = a21*a21;
		double t9 = a28*a28;
		double t11 = a35*a35;
		double t13 = a42*a42;
		double t16 = a0*s0;
		double t18 = a7*s1;
		double t20 = a14*s2;
		double t22 = a21*s3;
		double t24 = a28*s4;
		double t26 = a35*s5;
		double t28 = a42*s6;
		double t30 = t16*a1+t18*a8+t20*a15+t22*a22+t24*a29+t26*a36+t28*a43;
		double t38 = t16*a2+t18*a9+t20*a16+t22*a23+t24*a30+t26*a37+t28*a44;
		double t46 = t16*a3+t18*a10+t20*a17+t22*a24+t24*a31+t26*a38+t28*a45;
		double t54 = t16*a4+t18*a11+t20*a18+t22*a25+t24*a32+t26*a39+t28*a46;
		double t62 = t16*a5+t18*a12+t20*a19+t22*a26+t24*a33+t26*a40+t28*a47;
		double t70 = t16*a6+t18*a13+t20*a20+t22*a27+t24*a34+t26*a41+t28*a48;
		double t71 = a1*a1;
		double t73 = a8*a8;
		double t75 = a15*a15;
		double t77 = a22*a22;
		double t79 = a29*a29;
		double t81 = a36*a36;
		double t83 = a43*a43;
		double t86 = a1*s0;
		double t88 = a8*s1;
		double t90 = a15*s2;
		double t92 = a22*s3;
		double t94 = a29*s4;
		double t96 = a36*s5;
		double t98 = a43*s6;
		double t100 = t86*a2+t88*a9+t90*a16+t92*a23+t94*a30+t96*a37+t98*a44;
		double t108 = t86*a3+t88*a10+t90*a17+t92*a24+t94*a31+t96*a38+t98*a45;
		double t116 = t86*a4+t88*a11+t90*a18+t92*a25+t94*a32+t96*a39+t98*a46;
		double t124 = t86*a5+t88*a12+t90*a19+t92*a26+t94*a33+t96*a40+t98*a47;
		double t132 = t86*a6+t88*a13+t90*a20+t92*a27+t94*a34+t96*a41+t98*a48;
		double t133 = a2*a2;
		double t135 = a9*a9;
		double t137 = a16*a16;
		double t139 = a23*a23;
		double t141 = a30*a30;
		double t143 = a37*a37;
		double t145 = a44*a44;
		double t148 = a2*s0;
		double t150 = a9*s1;
		double t152 = a16*s2;
		double t154 = a23*s3;
		double t156 = a30*s4;
		double t158 = a37*s5;
		double t160 = a44*s6;
		double t162 = t148*a3+t150*a10+t152*a17+t154*a24+t156*a31+t158*a38+t160*a45;
		double t170 = t148*a4+t150*a11+t152*a18+t154*a25+t156*a32+t158*a39+t160*a46;
		double t178 = t148*a5+t150*a12+t152*a19+t154*a26+t156*a33+t158*a40+t160*a47;
		double t186 = t148*a6+t150*a13+t152*a20+t154*a27+t156*a34+t158*a41+t160*a48;
		double t187 = a3*a3;
		double t189 = a10*a10;
		double t191 = a17*a17;
		double t193 = a24*a24;
		double t195 = a31*a31;
		double t197 = a38*a38;
		double t199 = a45*a45;
		double t202 = a3*s0;
		double t204 = a10*s1;
		double t206 = a17*s2;
		double t208 = a24*s3;
		double t210 = a31*s4;
		double t212 = a38*s5;
		double t214 = a45*s6;
		double t216 = t202*a4+t204*a11+t206*a18+t208*a25+t210*a32+t212*a39+t214*a46;
		double t224 = t202*a5+t204*a12+t206*a19+t208*a26+t210*a33+t212*a40+t214*a47;
		double t232 = t202*a6+t204*a13+t206*a20+t208*a27+t210*a34+t212*a41+t214*a48;
		double t233 = a4*a4;
		double t235 = a11*a11;
		double t237 = a18*a18;
		double t239 = a25*a25;
		double t241 = a32*a32;
		double t243 = a39*a39;
		double t245 = a46*a46;
		double t248 = a4*s0;
		double t250 = a11*s1;
		double t252 = a18*s2;
		double t254 = a25*s3;
		double t256 = a32*s4;
		double t258 = a39*s5;
		double t260 = a46*s6;
		double t262 = t248*a5+t250*a12+t252*a19+t254*a26+t256*a33+t258*a40+t260*a47;
		double t270 = t248*a6+t250*a13+t252*a20+t254*a27+t256*a34+t258*a41+t260*a48;
		double t271 = a5*a5;
		double t273 = a12*a12;
		double t275 = a19*a19;
		double t277 = a26*a26;
		double t279 = a33*a33;
		double t281 = a40*a40;
		double t283 = a47*a47;
		double t300 = a5*s0*a6+a12*s1*a13+a19*s2*a20+a26*s3*a27+a33*s4*a34+a40*s5*a41+a47*s6*a48;
		double t301 = a6*a6;
		double t303 = a13*a13;
		double t305 = a20*a20;
		double t307 = a27*a27;
		double t309 = a34*a34;
		double t311 = a41*a41;
		double t313 = a48*a48;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6;
		*tmp_M_ptr++ = t30;
		*tmp_M_ptr++ = t38;
		*tmp_M_ptr++ = t46;
		*tmp_M_ptr++ = t54;
		*tmp_M_ptr++ = t62;
		*tmp_M_ptr++ = t70;
		*tmp_M_ptr++ = t71*s0+t73*s1+t75*s2+t77*s3+t79*s4+t81*s5+t83*s6;
		*tmp_M_ptr++ = t100;
		*tmp_M_ptr++ = t108;
		*tmp_M_ptr++ = t116;
		*tmp_M_ptr++ = t124;
		*tmp_M_ptr++ = t132;
		*tmp_M_ptr++ = t133*s0+t135*s1+t137*s2+t139*s3+t141*s4+t143*s5+t145*s6;
		*tmp_M_ptr++ = t162;
		*tmp_M_ptr++ = t170;
		*tmp_M_ptr++ = t178;
		*tmp_M_ptr++ = t186;
		*tmp_M_ptr++ = t187*s0+t189*s1+t191*s2+t193*s3+t195*s4+t197*s5+t199*s6;
		*tmp_M_ptr++ = t216;
		*tmp_M_ptr++ = t224;
		*tmp_M_ptr++ = t232;
		*tmp_M_ptr++ = t233*s0+t235*s1+t237*s2+t239*s3+t241*s4+t243*s5+t245*s6;
		*tmp_M_ptr++ = t262;
		*tmp_M_ptr++ = t270;
		*tmp_M_ptr++ = t271*s0+t273*s1+t275*s2+t277*s3+t279*s4+t281*s5+t283*s6;
		*tmp_M_ptr++ = t300;
		*tmp_M_ptr++ = t301*s0+t303*s1+t305*s2+t307*s3+t309*s4+t311*s5+t313*s6;

	}
	else if ( this->M_size == 8 ) {
		double D[8];
		double U[64];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 8);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63;
		double s0,s1,s2,s3,s4,s5,s6,s7;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;
		a49 = *tmp_M_ptr++;
		a50 = *tmp_M_ptr++;
		a51 = *tmp_M_ptr++;
		a52 = *tmp_M_ptr++;
		a53 = *tmp_M_ptr++;
		a54 = *tmp_M_ptr++;
		a55 = *tmp_M_ptr++;
		a56 = *tmp_M_ptr++;
		a57 = *tmp_M_ptr++;
		a58 = *tmp_M_ptr++;
		a59 = *tmp_M_ptr++;
		a60 = *tmp_M_ptr++;
		a61 = *tmp_M_ptr++;
		a62 = *tmp_M_ptr++;
		a63 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = exp(*tmp_M_ptr++);
		s1 = exp(*tmp_M_ptr++);
		s2 = exp(*tmp_M_ptr++);
		s3 = exp(*tmp_M_ptr++);
		s4 = exp(*tmp_M_ptr++);
		s5 = exp(*tmp_M_ptr++);
		s6 = exp(*tmp_M_ptr++);
		s7 = exp(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a8*a8;
		double t5 = a16*a16;
		double t7 = a24*a24;
		double t9 = a32*a32;
		double t11 = a40*a40;
		double t13 = a48*a48;
		double t15 = a56*a56;
		double t18 = a0*s0;
		double t20 = a8*s1;
		double t22 = a16*s2;
		double t24 = a24*s3;
		double t26 = a32*s4;
		double t28 = a40*s5;
		double t30 = a48*s6;
		double t32 = a56*s7;
		double t34 = t18*a1+t20*a9+t22*a17+t24*a25+t26*a33+t28*a41+t30*a49+t32*a57;
		double t43 = t18*a2+t20*a10+t22*a18+t24*a26+t26*a34+t28*a42+t30*a50+t32*a58;
		double t52 = t18*a3+t20*a11+t22*a19+t24*a27+t26*a35+t28*a43+t30*a51+t32*a59;
		double t61 = t18*a4+t20*a12+t22*a20+t24*a28+t26*a36+t28*a44+t30*a52+t32*a60;
		double t70 = t18*a5+t20*a13+t22*a21+t24*a29+t26*a37+t28*a45+t30*a53+t32*a61;
		double t79 = t18*a6+t20*a14+t22*a22+t24*a30+t26*a38+t28*a46+t30*a54+t32*a62;
		double t88 = t18*a7+t20*a15+t22*a23+t24*a31+t26*a39+t28*a47+t30*a55+t32*a63;
		double t89 = a1*a1;
		double t91 = a9*a9;
		double t93 = a17*a17;
		double t95 = a25*a25;
		double t97 = a33*a33;
		double t99 = a41*a41;
		double t101 = a49*a49;
		double t103 = a57*a57;
		double t106 = a1*s0;
		double t108 = a9*s1;
		double t110 = a17*s2;
		double t112 = a25*s3;
		double t114 = a33*s4;
		double t116 = a41*s5;
		double t118 = a49*s6;
		double t120 = a57*s7;
		double t122 = t106*a2+t108*a10+t110*a18+t112*a26+t114*a34+t116*a42+t118*a50+t120*a58;
		double t131 = t106*a3+t108*a11+t110*a19+t112*a27+t114*a35+t116*a43+t118*a51+t120*a59;
		double t140 = t106*a4+t108*a12+t110*a20+t112*a28+t114*a36+t116*a44+t118*a52+t120*a60;
		double t149 = t106*a5+t108*a13+t110*a21+t112*a29+t114*a37+t116*a45+t118*a53+t120*a61;
		double t158 = t106*a6+t108*a14+t110*a22+t112*a30+t114*a38+t116*a46+t118*a54+t120*a62;
		double t167 = t106*a7+t108*a15+t110*a23+t112*a31+t114*a39+t116*a47+t118*a55+t120*a63;
		double t168 = a2*a2;
		double t170 = a10*a10;
		double t172 = a18*a18;
		double t174 = a26*a26;
		double t176 = a34*a34;
		double t178 = a42*a42;
		double t180 = a50*a50;
		double t182 = a58*a58;
		double t185 = a2*s0;
		double t187 = a10*s1;
		double t189 = a18*s2;
		double t191 = a26*s3;
		double t193 = a34*s4;
		double t195 = a42*s5;
		double t197 = a50*s6;
		double t199 = a58*s7;
		double t201 = t185*a3+t187*a11+t189*a19+t191*a27+t193*a35+t195*a43+t197*a51+t199*a59;
		double t210 = t185*a4+t187*a12+t189*a20+t191*a28+t193*a36+t195*a44+t197*a52+t199*a60;
		double t219 = t185*a5+t187*a13+t189*a21+t191*a29+t193*a37+t195*a45+t197*a53+t199*a61;
		double t228 = t185*a6+t187*a14+t189*a22+t191*a30+t193*a38+t195*a46+t197*a54+t199*a62;
		double t237 = t185*a7+t187*a15+t189*a23+t191*a31+t193*a39+t195*a47+t197*a55+t199*a63;
		double t238 = a3*a3;
		double t240 = a11*a11;
		double t242 = a19*a19;
		double t244 = a27*a27;
		double t246 = a35*a35;
		double t248 = a43*a43;
		double t250 = a51*a51;
		double t252 = a59*a59;
		double t255 = a3*s0;
		double t257 = a11*s1;
		double t259 = a19*s2;
		double t261 = a27*s3;
		double t263 = a35*s4;
		double t265 = a43*s5;
		double t267 = a51*s6;
		double t269 = a59*s7;
		double t271 = t255*a4+t257*a12+t259*a20+t261*a28+t263*a36+t265*a44+t267*a52+t269*a60;
		double t280 = t255*a5+t257*a13+t259*a21+t261*a29+t263*a37+t265*a45+t267*a53+t269*a61;
		double t289 = t255*a6+t257*a14+t259*a22+t261*a30+t263*a38+t265*a46+t267*a54+t269*a62;
		double t298 = t255*a7+t257*a15+t259*a23+t261*a31+t263*a39+t265*a47+t267*a55+t269*a63;
		double t299 = a4*a4;
		double t301 = a12*a12;
		double t303 = a20*a20;
		double t305 = a28*a28;
		double t307 = a36*a36;
		double t309 = a44*a44;
		double t311 = a52*a52;
		double t313 = a60*a60;
		double t316 = a4*s0;
		double t318 = a12*s1;
		double t320 = a20*s2;
		double t322 = a28*s3;
		double t324 = a36*s4;
		double t326 = a44*s5;
		double t328 = a52*s6;
		double t330 = a60*s7;
		double t332 = t316*a5+t318*a13+t320*a21+t322*a29+t324*a37+t326*a45+t328*a53+t330*a61;
		double t341 = t316*a6+t318*a14+t320*a22+t322*a30+t324*a38+t326*a46+t328*a54+t330*a62;
		double t350 = t316*a7+t318*a15+t320*a23+t322*a31+t324*a39+t326*a47+t328*a55+t330*a63;
		double t351 = a5*a5;
		double t353 = a13*a13;
		double t355 = a21*a21;
		double t357 = a29*a29;
		double t359 = a37*a37;
		double t361 = a45*a45;
		double t363 = a53*a53;
		double t365 = a61*a61;
		double t368 = a5*s0;
		double t370 = a13*s1;
		double t372 = a21*s2;
		double t374 = a29*s3;
		double t376 = a37*s4;
		double t378 = a45*s5;
		double t380 = a53*s6;
		double t382 = a61*s7;
		double t384 = t368*a6+t370*a14+t372*a22+t374*a30+t376*a38+t378*a46+t380*a54+t382*a62;
		double t393 = t368*a7+t370*a15+t372*a23+t374*a31+t376*a39+t378*a47+t380*a55+t382*a63;
		double t394 = a6*a6;
		double t396 = a14*a14;
		double t398 = a22*a22;
		double t400 = a30*a30;
		double t402 = a38*a38;
		double t404 = a46*a46;
		double t406 = a54*a54;
		double t408 = a62*a62;
		double t427 = a6*s0*a7+a14*s1*a15+a22*s2*a23+a30*s3*a31+a38*s4*a39+a46*s5*a47+a54*s6*a55+a62*s7*a63;
		double t428 = a7*a7;
		double t430 = a15*a15;
		double t432 = a23*a23;
		double t434 = a31*a31;
		double t436 = a39*a39;
		double t438 = a47*a47;
		double t440 = a55*a55;
		double t442 = a63*a63;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6+t15*s7;
		*tmp_M_ptr++ = t34;
		*tmp_M_ptr++ = t43;
		*tmp_M_ptr++ = t52;
		*tmp_M_ptr++ = t61;
		*tmp_M_ptr++ = t70;
		*tmp_M_ptr++ = t79;
		*tmp_M_ptr++ = t88;
		*tmp_M_ptr++ = t89*s0+t91*s1+t93*s2+t95*s3+t97*s4+t99*s5+t101*s6+t103*s7;
		*tmp_M_ptr++ = t122;
		*tmp_M_ptr++ = t131;
		*tmp_M_ptr++ = t140;
		*tmp_M_ptr++ = t149;
		*tmp_M_ptr++ = t158;
		*tmp_M_ptr++ = t167;
		*tmp_M_ptr++ = t168*s0+t170*s1+t172*s2+t174*s3+t176*s4+t178*s5+t180*s6+t182*s7;
		*tmp_M_ptr++ = t201;
		*tmp_M_ptr++ = t210;
		*tmp_M_ptr++ = t219;
		*tmp_M_ptr++ = t228;
		*tmp_M_ptr++ = t237;
		*tmp_M_ptr++ = t238*s0+t240*s1+t242*s2+t244*s3+t246*s4+t248*s5+t250*s6+t252*s7;
		*tmp_M_ptr++ = t271;
		*tmp_M_ptr++ = t280;
		*tmp_M_ptr++ = t289;
		*tmp_M_ptr++ = t298;
		*tmp_M_ptr++ = t299*s0+t301*s1+t303*s2+t305*s3+t307*s4+t309*s5+t311*s6+t313*s7;
		*tmp_M_ptr++ = t332;
		*tmp_M_ptr++ = t341;
		*tmp_M_ptr++ = t350;
		*tmp_M_ptr++ = t351*s0+t353*s1+t355*s2+t357*s3+t359*s4+t361*s5+t363*s6+t365*s7;
		*tmp_M_ptr++ = t384;
		*tmp_M_ptr++ = t393;
		*tmp_M_ptr++ = t394*s0+t396*s1+t398*s2+t400*s3+t402*s4+t404*s5+t406*s6+t408*s7;
		*tmp_M_ptr++ = t427;
		*tmp_M_ptr++ = t428*s0+t430*s1+t432*s2+t434*s3+t436*s4+t438*s5+t440*s6+t442*s7;
	}
	else if ( this->M_size == 9 ) {
		double D[9];
		double U[81];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 9);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80;
		double s0,s1,s2,s3,s4,s5,s6,s7,s8;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;
		a49 = *tmp_M_ptr++;
		a50 = *tmp_M_ptr++;
		a51 = *tmp_M_ptr++;
		a52 = *tmp_M_ptr++;
		a53 = *tmp_M_ptr++;
		a54 = *tmp_M_ptr++;
		a55 = *tmp_M_ptr++;
		a56 = *tmp_M_ptr++;
		a57 = *tmp_M_ptr++;
		a58 = *tmp_M_ptr++;
		a59 = *tmp_M_ptr++;
		a60 = *tmp_M_ptr++;
		a61 = *tmp_M_ptr++;
		a62 = *tmp_M_ptr++;
		a63 = *tmp_M_ptr++;
		a64 = *tmp_M_ptr++;
		a65 = *tmp_M_ptr++;
		a66 = *tmp_M_ptr++;
		a67 = *tmp_M_ptr++;
		a68 = *tmp_M_ptr++;
		a69 = *tmp_M_ptr++;
		a70 = *tmp_M_ptr++;
		a71 = *tmp_M_ptr++;
		a72 = *tmp_M_ptr++;
		a73 = *tmp_M_ptr++;
		a74 = *tmp_M_ptr++;
		a75 = *tmp_M_ptr++;
		a76 = *tmp_M_ptr++;
		a77 = *tmp_M_ptr++;
		a78 = *tmp_M_ptr++;
		a79 = *tmp_M_ptr++;
		a80 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = exp(*tmp_M_ptr++);
		s1 = exp(*tmp_M_ptr++);
		s2 = exp(*tmp_M_ptr++);
		s3 = exp(*tmp_M_ptr++);
		s4 = exp(*tmp_M_ptr++);
		s5 = exp(*tmp_M_ptr++);
		s6 = exp(*tmp_M_ptr++);
		s7 = exp(*tmp_M_ptr++);
		s8 = exp(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a9*a9;
		double t5 = a18*a18;
		double t7 = a27*a27;
		double t9 = a36*a36;
		double t11 = a45*a45;
		double t13 = a54*a54;
		double t15 = a63*a63;
		double t17 = a72*a72;
		double t20 = a0*s0;
		double t22 = a9*s1;
		double t24 = a18*s2;
		double t26 = a27*s3;
		double t28 = a36*s4;
		double t30 = a45*s5;
		double t32 = a54*s6;
		double t34 = a63*s7;
		double t36 = a72*s8;
		double t38 = t20*a1+t22*a10+t24*a19+t26*a28+t28*a37+t30*a46+t32*a55+t34*a64+t36*a73;
		double t48 = t20*a2+t22*a11+t24*a20+t26*a29+t28*a38+t30*a47+t32*a56+t34*a65+t36*a74;
		double t58 = t20*a3+t22*a12+t24*a21+t26*a30+t28*a39+t30*a48+t32*a57+t34*a66+t36*a75;
		double t68 = t20*a4+t22*a13+t24*a22+t26*a31+t28*a40+t30*a49+t32*a58+t34*a67+t36*a76;
		double t78 = t20*a5+t22*a14+t24*a23+t26*a32+t28*a41+t30*a50+t32*a59+t34*a68+t36*a77;
		double t88 = t20*a6+t22*a15+t24*a24+t26*a33+t28*a42+t30*a51+t32*a60+t34*a69+t36*a78;
		double t98 = t20*a7+t22*a16+t24*a25+t26*a34+t28*a43+t30*a52+t32*a61+t34*a70+t36*a79;
		double t108 = t20*a8+t22*a17+t24*a26+t26*a35+t28*a44+t30*a53+t32*a62+t34*a71+t36*a80;
		double t109 = a1*a1;
		double t111 = a10*a10;
		double t113 = a19*a19;
		double t115 = a28*a28;
		double t117 = a37*a37;
		double t119 = a46*a46;
		double t121 = a55*a55;
		double t123 = a64*a64;
		double t125 = a73*a73;
		double t128 = a1*s0;
		double t130 = a10*s1;
		double t132 = a19*s2;
		double t134 = a28*s3;
		double t136 = a37*s4;
		double t138 = a46*s5;
		double t140 = a55*s6;
		double t142 = a64*s7;
		double t144 = a73*s8;
		double t146 = t128*a2+t130*a11+t132*a20+t134*a29+t136*a38+t138*a47+t140*a56+t142*a65+t144*a74;
		double t156 = t128*a3+t130*a12+t132*a21+t134*a30+t136*a39+t138*a48+t140*a57+t142*a66+t144*a75;
		double t166 = t128*a4+t130*a13+t132*a22+t134*a31+t136*a40+t138*a49+t140*a58+t142*a67+t144*a76;
		double t176 = t128*a5+t130*a14+t132*a23+t134*a32+t136*a41+t138*a50+t140*a59+t142*a68+t144*a77;
		double t186 = t128*a6+t130*a15+t132*a24+t134*a33+t136*a42+t138*a51+t140*a60+t142*a69+t144*a78;
		double t196 = t128*a7+t130*a16+t132*a25+t134*a34+t136*a43+t138*a52+t140*a61+t142*a70+t144*a79;
		double t206 = t128*a8+t130*a17+t132*a26+t134*a35+t136*a44+t138*a53+t140*a62+t142*a71+t144*a80;
		double t207 = a2*a2;
		double t209 = a11*a11;
		double t211 = a20*a20;
		double t213 = a29*a29;
		double t215 = a38*a38;
		double t217 = a47*a47;
		double t219 = a56*a56;
		double t221 = a65*a65;
		double t223 = a74*a74;
		double t226 = a2*s0;
		double t228 = a11*s1;
		double t230 = a20*s2;
		double t232 = a29*s3;
		double t234 = a38*s4;
		double t236 = a47*s5;
		double t238 = a56*s6;
		double t240 = a65*s7;
		double t242 = a74*s8;
		double t244 = t226*a3+t228*a12+t230*a21+t232*a30+t234*a39+t236*a48+t238*a57+t240*a66+t242*a75;
		double t254 = t226*a4+t228*a13+t230*a22+t232*a31+t234*a40+t236*a49+t238*a58+t240*a67+t242*a76;
		double t264 = t226*a5+t228*a14+t230*a23+t232*a32+t234*a41+t236*a50+t238*a59+t240*a68+t242*a77;
		double t274 = t226*a6+t228*a15+t230*a24+t232*a33+t234*a42+t236*a51+t238*a60+t240*a69+t242*a78;
		double t284 = t226*a7+t228*a16+t230*a25+t232*a34+t234*a43+t236*a52+t238*a61+t240*a70+t242*a79;
		double t294 = t226*a8+t228*a17+t230*a26+t232*a35+t234*a44+t236*a53+t238*a62+t240*a71+t242*a80;
		double t295 = a3*a3;
		double t297 = a12*a12;
		double t299 = a21*a21;
		double t301 = a30*a30;
		double t303 = a39*a39;
		double t305 = a48*a48;
		double t307 = a57*a57;
		double t309 = a66*a66;
		double t311 = a75*a75;
		double t314 = a3*s0;
		double t316 = a12*s1;
		double t318 = a21*s2;
		double t320 = a30*s3;
		double t322 = a39*s4;
		double t324 = a48*s5;
		double t326 = a57*s6;
		double t328 = a66*s7;
		double t330 = a75*s8;
		double t332 = t314*a4+t316*a13+t318*a22+t320*a31+t322*a40+t324*a49+t326*a58+t328*a67+t330*a76;
		double t342 = t314*a5+t316*a14+t318*a23+t320*a32+t322*a41+t324*a50+t326*a59+t328*a68+t330*a77;
		double t352 = t314*a6+t316*a15+t318*a24+t320*a33+t322*a42+t324*a51+t326*a60+t328*a69+t330*a78;
		double t362 = t314*a7+t316*a16+t318*a25+t320*a34+t322*a43+t324*a52+t326*a61+t328*a70+t330*a79;
		double t372 = t314*a8+t316*a17+t318*a26+t320*a35+t322*a44+t324*a53+t326*a62+t328*a71+t330*a80;
		double t373 = a4*a4;
		double t375 = a13*a13;
		double t377 = a22*a22;
		double t379 = a31*a31;
		double t381 = a40*a40;
		double t383 = a49*a49;
		double t385 = a58*a58;
		double t387 = a67*a67;
		double t389 = a76*a76;
		double t392 = a4*s0;
		double t394 = a13*s1;
		double t396 = a22*s2;
		double t398 = a31*s3;
		double t400 = a40*s4;
		double t402 = a49*s5;
		double t404 = a58*s6;
		double t406 = a67*s7;
		double t408 = a76*s8;
		double t410 = t392*a5+t394*a14+t396*a23+t398*a32+t400*a41+t402*a50+t404*a59+t406*a68+t408*a77;
		double t420 = t392*a6+t394*a15+t396*a24+t398*a33+t400*a42+t402*a51+t404*a60+t406*a69+t408*a78;
		double t430 = t392*a7+t394*a16+t396*a25+t398*a34+t400*a43+t402*a52+t404*a61+t406*a70+t408*a79;
		double t440 = t392*a8+t394*a17+t396*a26+t398*a35+t400*a44+t402*a53+t404*a62+t406*a71+t408*a80;
		double t441 = a5*a5;
		double t443 = a14*a14;
		double t445 = a23*a23;
		double t447 = a32*a32;
		double t449 = a41*a41;
		double t451 = a50*a50;
		double t453 = a59*a59;
		double t455 = a68*a68;
		double t457 = a77*a77;
		double t460 = a5*s0;
		double t462 = a14*s1;
		double t464 = a23*s2;
		double t466 = a32*s3;
		double t468 = a41*s4;
		double t470 = a50*s5;
		double t472 = a59*s6;
		double t474 = a68*s7;
		double t476 = a77*s8;
		double t478 = t460*a6+t462*a15+t464*a24+t466*a33+t468*a42+t470*a51+t472*a60+t474*a69+t476*a78;
		double t488 = t460*a7+t462*a16+t464*a25+t466*a34+t468*a43+t470*a52+t472*a61+t474*a70+t476*a79;
		double t498 = t460*a8+t462*a17+t464*a26+t466*a35+t468*a44+t470*a53+t472*a62+t474*a71+t476*a80;
		double t499 = a6*a6;
		double t501 = a15*a15;
		double t503 = a24*a24;
		double t505 = a33*a33;
		double t507 = a42*a42;
		double t509 = a51*a51;
		double t511 = a60*a60;
		double t513 = a69*a69;
		double t515 = a78*a78;
		double t518 = a6*s0;
		double t520 = a15*s1;
		double t522 = a24*s2;
		double t524 = a33*s3;
		double t526 = a42*s4;
		double t528 = a51*s5;
		double t530 = a60*s6;
		double t532 = a69*s7;
		double t534 = a78*s8;
		double t536 = t518*a7+t520*a16+t522*a25+t524*a34+t526*a43+t528*a52+t530*a61+t532*a70+t534*a79;
		double t546 = t518*a8+t520*a17+t522*a26+t524*a35+t526*a44+t528*a53+t530*a62+t532*a71+t534*a80;
		double t547 = a7*a7;
		double t549 = a16*a16;
		double t551 = a25*a25;
		double t553 = a34*a34;
		double t555 = a43*a43;
		double t557 = a52*a52;
		double t559 = a61*a61;
		double t561 = a70*a70;
		double t563 = a79*a79;
		double t584 = a7*s0*a8+a16*s1*a17+a25*s2*a26+a34*s3*a35+a43*s4*a44+a52*s5*a53+a61*s6*a62+a70*s7*a71+a79*s8*a80;
		double t585 = a8*a8;
		double t587 = a17*a17;
		double t589 = a26*a26;
		double t591 = a35*a35;
		double t593 = a44*a44;
		double t595 = a53*a53;
		double t597 = a62*a62;
		double t599 = a71*a71;
		double t601 = a80*a80;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6+t15*s7+t17*s8;
		*tmp_M_ptr++ = t38;
		*tmp_M_ptr++ = t48;
		*tmp_M_ptr++ = t58;
		*tmp_M_ptr++ = t68;
		*tmp_M_ptr++ = t78;
		*tmp_M_ptr++ = t88;
		*tmp_M_ptr++ = t98;
		*tmp_M_ptr++ = t108;
		*tmp_M_ptr++ = t109*s0+t111*s1+t113*s2+t115*s3+t117*s4+t119*s5+t121*s6+t123*s7+t125*s8;
		*tmp_M_ptr++ = t146;
		*tmp_M_ptr++ = t156;
		*tmp_M_ptr++ = t166;
		*tmp_M_ptr++ = t176;
		*tmp_M_ptr++ = t186;
		*tmp_M_ptr++ = t196;
		*tmp_M_ptr++ = t206;
		*tmp_M_ptr++ = t207*s0+t209*s1+t211*s2+t213*s3+t215*s4+t217*s5+t219*s6+t221*s7+t223*s8;
		*tmp_M_ptr++ = t244;
		*tmp_M_ptr++ = t254;
		*tmp_M_ptr++ = t264;
		*tmp_M_ptr++ = t274;
		*tmp_M_ptr++ = t284;
		*tmp_M_ptr++ = t294;
		*tmp_M_ptr++ = t295*s0+t297*s1+t299*s2+t301*s3+t303*s4+t305*s5+t307*s6+t309*s7+t311*s8;
		*tmp_M_ptr++ = t332;
		*tmp_M_ptr++ = t342;
		*tmp_M_ptr++ = t352;
		*tmp_M_ptr++ = t362;
		*tmp_M_ptr++ = t372;
		*tmp_M_ptr++ = t373*s0+t375*s1+t377*s2+t379*s3+t381*s4+t383*s5+t385*s6+t387*s7+t389*s8;
		*tmp_M_ptr++ = t410;
		*tmp_M_ptr++ = t420;
		*tmp_M_ptr++ = t430;
		*tmp_M_ptr++ = t440;
		*tmp_M_ptr++ = t441*s0+t443*s1+t445*s2+t447*s3+t449*s4+t451*s5+t453*s6+t455*s7+t457*s8;
		*tmp_M_ptr++ = t478;
		*tmp_M_ptr++ = t488;
		*tmp_M_ptr++ = t498;
		*tmp_M_ptr++ = t499*s0+t501*s1+t503*s2+t505*s3+t507*s4+t509*s5+t511*s6+t513*s7+t515*s8;
		*tmp_M_ptr++ = t536;
		*tmp_M_ptr++ = t546;
		*tmp_M_ptr++ = t547*s0+t549*s1+t551*s2+t553*s3+t555*s4+t557*s5+t559*s6+t561*s7+t563*s8;
		*tmp_M_ptr++ = t584;
		*tmp_M_ptr++ = t585*s0+t587*s1+t589*s2+t591*s3+t593*s4+t595*s5+t597*s6+t599*s7+t601*s8;
	}
	else if ( this->M_size == 10 ) {
		double D[10];
		double U[100];

		double *M_ptr = this->GetData(0);

		eigen_n_decomposition(M_ptr, U, D, 10);

		delete [] M_ptr;

		double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48;
		double a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80,a81,a82,a83,a84,a85,a86,a87,a88,a89,a90,a91,a92,a93,a94,a95,a96,a97,a98,a99;
		double s0,s1,s2,s3,s4,s5,s6,s7,s8,s9;

		double *tmp_M_ptr = U;

		a0 = *tmp_M_ptr++;
		a1 = *tmp_M_ptr++;
		a2 = *tmp_M_ptr++;
		a3 = *tmp_M_ptr++;
		a4 = *tmp_M_ptr++;
		a5 = *tmp_M_ptr++;
		a6 = *tmp_M_ptr++;
		a7 = *tmp_M_ptr++;
		a8 = *tmp_M_ptr++;
		a9 = *tmp_M_ptr++;
		a10 = *tmp_M_ptr++;
		a11 = *tmp_M_ptr++;
		a12 = *tmp_M_ptr++;
		a13 = *tmp_M_ptr++;
		a14 = *tmp_M_ptr++;
		a15 = *tmp_M_ptr++;
		a16 = *tmp_M_ptr++;
		a17 = *tmp_M_ptr++;
		a18 = *tmp_M_ptr++;
		a19 = *tmp_M_ptr++;
		a20 = *tmp_M_ptr++;
		a21 = *tmp_M_ptr++;
		a22 = *tmp_M_ptr++;
		a23 = *tmp_M_ptr++;
		a24 = *tmp_M_ptr++;
		a25 = *tmp_M_ptr++;
		a26 = *tmp_M_ptr++;
		a27 = *tmp_M_ptr++;
		a28 = *tmp_M_ptr++;
		a29 = *tmp_M_ptr++;
		a30 = *tmp_M_ptr++;
		a31 = *tmp_M_ptr++;
		a32 = *tmp_M_ptr++;
		a33 = *tmp_M_ptr++;
		a34 = *tmp_M_ptr++;
		a35 = *tmp_M_ptr++;
		a36 = *tmp_M_ptr++;
		a37 = *tmp_M_ptr++;
		a38 = *tmp_M_ptr++;
		a39 = *tmp_M_ptr++;
		a40 = *tmp_M_ptr++;
		a41 = *tmp_M_ptr++;
		a42 = *tmp_M_ptr++;
		a43 = *tmp_M_ptr++;
		a44 = *tmp_M_ptr++;
		a45 = *tmp_M_ptr++;
		a46 = *tmp_M_ptr++;
		a47 = *tmp_M_ptr++;
		a48 = *tmp_M_ptr++;
		a49 = *tmp_M_ptr++;
		a50 = *tmp_M_ptr++;
		a51 = *tmp_M_ptr++;
		a52 = *tmp_M_ptr++;
		a53 = *tmp_M_ptr++;
		a54 = *tmp_M_ptr++;
		a55 = *tmp_M_ptr++;
		a56 = *tmp_M_ptr++;
		a57 = *tmp_M_ptr++;
		a58 = *tmp_M_ptr++;
		a59 = *tmp_M_ptr++;
		a60 = *tmp_M_ptr++;
		a61 = *tmp_M_ptr++;
		a62 = *tmp_M_ptr++;
		a63 = *tmp_M_ptr++;
		a64 = *tmp_M_ptr++;
		a65 = *tmp_M_ptr++;
		a66 = *tmp_M_ptr++;
		a67 = *tmp_M_ptr++;
		a68 = *tmp_M_ptr++;
		a69 = *tmp_M_ptr++;
		a70 = *tmp_M_ptr++;
		a71 = *tmp_M_ptr++;
		a72 = *tmp_M_ptr++;
		a73 = *tmp_M_ptr++;
		a74 = *tmp_M_ptr++;
		a75 = *tmp_M_ptr++;
		a76 = *tmp_M_ptr++;
		a77 = *tmp_M_ptr++;
		a78 = *tmp_M_ptr++;
		a79 = *tmp_M_ptr++;
		a80 = *tmp_M_ptr++;
		a81 = *tmp_M_ptr++;
		a82 = *tmp_M_ptr++;
		a83 = *tmp_M_ptr++;
		a84 = *tmp_M_ptr++;
		a85 = *tmp_M_ptr++;
		a86 = *tmp_M_ptr++;
		a87 = *tmp_M_ptr++;
		a88 = *tmp_M_ptr++;
		a89 = *tmp_M_ptr++;
		a90 = *tmp_M_ptr++;
		a91 = *tmp_M_ptr++;
		a92 = *tmp_M_ptr++;
		a93 = *tmp_M_ptr++;
		a94 = *tmp_M_ptr++;
		a95 = *tmp_M_ptr++;
		a96 = *tmp_M_ptr++;
		a97 = *tmp_M_ptr++;
		a98 = *tmp_M_ptr++;
		a99 = *tmp_M_ptr++;

		tmp_M_ptr = D;

		s0 = exp(*tmp_M_ptr++);
		s1 = exp(*tmp_M_ptr++);
		s2 = exp(*tmp_M_ptr++);
		s3 = exp(*tmp_M_ptr++);
		s4 = exp(*tmp_M_ptr++);
		s5 = exp(*tmp_M_ptr++);
		s6 = exp(*tmp_M_ptr++);
		s7 = exp(*tmp_M_ptr++);
		s8 = exp(*tmp_M_ptr++);
		s9 = exp(*tmp_M_ptr++);

		double t1 = a0*a0;
		double t3 = a10*a10;
		double t5 = a20*a20;
		double t7 = a30*a30;
		double t9 = a40*a40;
		double t11 = a50*a50;
		double t13 = a60*a60;
		double t15 = a70*a70;
		double t17 = a80*a80;
		double t19 = a90*a90;
		double t22 = a0*s0;
		double t24 = a10*s1;
		double t26 = a20*s2;
		double t28 = a30*s3;
		double t30 = a40*s4;
		double t32 = a50*s5;
		double t34 = a60*s6;
		double t36 = a70*s7;
		double t38 = a80*s8;
		double t40 = a90*s9;
		double t42 = t22*a1+t24*a11+t26*a21+t28*a31+t30*a41+t32*a51+t34*a61+t36*a71+t38*a81+t40*a91;
		double t53 = t22*a2+t24*a12+t26*a22+t28*a32+t30*a42+t32*a52+t34*a62+t36*a72+t38*a82+t40*a92;
		double t64 = t22*a3+t24*a13+t26*a23+t28*a33+t30*a43+t32*a53+t34*a63+t36*a73+t38*a83+t40*a93;
		double t75 = t22*a4+t24*a14+t26*a24+t28*a34+t30*a44+t32*a54+t34*a64+t36*a74+t38*a84+t40*a94;
		double t86 = t22*a5+t24*a15+t26*a25+t28*a35+t30*a45+t32*a55+t34*a65+t36*a75+t38*a85+t40*a95;
		double t97 = t22*a6+t24*a16+t26*a26+t28*a36+t30*a46+t32*a56+t34*a66+t36*a76+t38*a86+t40*a96;
		double t108 = t22*a7+t24*a17+t26*a27+t28*a37+t30*a47+t32*a57+t34*a67+t36*a77+t38*a87+t40*a97;
		double t119 = t22*a8+t24*a18+t26*a28+t28*a38+t30*a48+t32*a58+t34*a68+t36*a78+t38*a88+t40*a98;
		double t130 = t22*a9+t24*a19+t26*a29+t28*a39+t30*a49+t32*a59+t34*a69+t36*a79+t38*a89+t40*a99;
		double t131 = a1*a1;
		double t133 = a11*a11;
		double t135 = a21*a21;
		double t137 = a31*a31;
		double t139 = a41*a41;
		double t141 = a51*a51;
		double t143 = a61*a61;
		double t145 = a71*a71;
		double t147 = a81*a81;
		double t149 = a91*a91;
		double t152 = a1*s0;
		double t154 = a11*s1;
		double t156 = a21*s2;
		double t158 = a31*s3;
		double t160 = a41*s4;
		double t162 = a51*s5;
		double t164 = a61*s6;
		double t166 = a71*s7;
		double t168 = a81*s8;
		double t170 = a91*s9;
		double t172 = t152*a2+t154*a12+t156*a22+t158*a32+t160*a42+t162*a52+t164*a62+t166*a72+t168*a82+t170*a92;
		double t183 = t152*a3+t154*a13+t156*a23+t158*a33+t160*a43+t162*a53+t164*a63+t166*a73+t168*a83+t170*a93;
		double t194 = t152*a4+t154*a14+t156*a24+t158*a34+t160*a44+t162*a54+t164*a64+t166*a74+t168*a84+t170*a94;
		double t205 = t152*a5+t154*a15+t156*a25+t158*a35+t160*a45+t162*a55+t164*a65+t166*a75+t168*a85+t170*a95;
		double t216 = t152*a6+t154*a16+t156*a26+t158*a36+t160*a46+t162*a56+t164*a66+t166*a76+t168*a86+t170*a96;
		double t227 = t152*a7+t154*a17+t156*a27+t158*a37+t160*a47+t162*a57+t164*a67+t166*a77+t168*a87+t170*a97;
		double t238 = t152*a8+t154*a18+t156*a28+t158*a38+t160*a48+t162*a58+t164*a68+t166*a78+t168*a88+t170*a98;
		double t249 = t152*a9+t154*a19+t156*a29+t158*a39+t160*a49+t162*a59+t164*a69+t166*a79+t168*a89+t170*a99;
		double t250 = a2*a2;
		double t252 = a12*a12;
		double t254 = a22*a22;
		double t256 = a32*a32;
		double t258 = a42*a42;
		double t260 = a52*a52;
		double t262 = a62*a62;
		double t264 = a72*a72;
		double t266 = a82*a82;
		double t268 = a92*a92;
		double t271 = a2*s0;
		double t273 = a12*s1;
		double t275 = a22*s2;
		double t277 = a32*s3;
		double t279 = a42*s4;
		double t281 = a52*s5;
		double t283 = a62*s6;
		double t285 = a72*s7;
		double t287 = a82*s8;
		double t289 = a92*s9;
		double t291 = t271*a3+t273*a13+t275*a23+t277*a33+t279*a43+t281*a53+t283*a63+t285*a73+t287*a83+t289*a93;
		double t302 = t271*a4+t273*a14+t275*a24+t277*a34+t279*a44+t281*a54+t283*a64+t285*a74+t287*a84+t289*a94;
		double t313 = t271*a5+t273*a15+t275*a25+t277*a35+t279*a45+t281*a55+t283*a65+t285*a75+t287*a85+t289*a95;
		double t324 = t271*a6+t273*a16+t275*a26+t277*a36+t279*a46+t281*a56+t283*a66+t285*a76+t287*a86+t289*a96;
		double t335 = t271*a7+t273*a17+t275*a27+t277*a37+t279*a47+t281*a57+t283*a67+t285*a77+t287*a87+t289*a97;
		double t346 = t271*a8+t273*a18+t275*a28+t277*a38+t279*a48+t281*a58+t283*a68+t285*a78+t287*a88+t289*a98;
		double t357 = t271*a9+t273*a19+t275*a29+t277*a39+t279*a49+t281*a59+t283*a69+t285*a79+t287*a89+t289*a99;
		double t358 = a3*a3;
		double t360 = a13*a13;
		double t362 = a23*a23;
		double t364 = a33*a33;
		double t366 = a43*a43;
		double t368 = a53*a53;
		double t370 = a63*a63;
		double t372 = a73*a73;
		double t374 = a83*a83;
		double t376 = a93*a93;
		double t379 = a3*s0;
		double t381 = a13*s1;
		double t383 = a23*s2;
		double t385 = a33*s3;
		double t387 = a43*s4;
		double t389 = a53*s5;
		double t391 = a63*s6;
		double t393 = a73*s7;
		double t395 = a83*s8;
		double t397 = a93*s9;
		double t399 = t379*a4+t381*a14+t383*a24+t385*a34+t387*a44+t389*a54+t391*a64+t393*a74+t395*a84+t397*a94;
		double t410 = t379*a5+t381*a15+t383*a25+t385*a35+t387*a45+t389*a55+t391*a65+t393*a75+t395*a85+t397*a95;
		double t421 = t379*a6+t381*a16+t383*a26+t385*a36+t387*a46+t389*a56+t391*a66+t393*a76+t395*a86+t397*a96;
		double t432 = t379*a7+t381*a17+t383*a27+t385*a37+t387*a47+t389*a57+t391*a67+t393*a77+t395*a87+t397*a97;
		double t443 = t379*a8+t381*a18+t383*a28+t385*a38+t387*a48+t389*a58+t391*a68+t393*a78+t395*a88+t397*a98;
		double t454 = t379*a9+t381*a19+t383*a29+t385*a39+t387*a49+t389*a59+t391*a69+t393*a79+t395*a89+t397*a99;
		double t455 = a4*a4;
		double t457 = a14*a14;
		double t459 = a24*a24;
		double t461 = a34*a34;
		double t463 = a44*a44;
		double t465 = a54*a54;
		double t467 = a64*a64;
		double t469 = a74*a74;
		double t471 = a84*a84;
		double t473 = a94*a94;
		double t476 = a4*s0;
		double t478 = a14*s1;
		double t480 = a24*s2;
		double t482 = a34*s3;
		double t484 = a44*s4;
		double t486 = a54*s5;
		double t488 = a64*s6;
		double t490 = a74*s7;
		double t492 = a84*s8;
		double t494 = a94*s9;
		double t496 = t476*a5+t478*a15+t480*a25+t482*a35+t484*a45+t486*a55+t488*a65+t490*a75+t492*a85+t494*a95;
		double t507 = t476*a6+t478*a16+t480*a26+t482*a36+t484*a46+t486*a56+t488*a66+t490*a76+t492*a86+t494*a96;
		double t518 = t476*a7+t478*a17+t480*a27+t482*a37+t484*a47+t486*a57+t488*a67+t490*a77+t492*a87+t494*a97;
		double t529 = t476*a8+t478*a18+t480*a28+t482*a38+t484*a48+t486*a58+t488*a68+t490*a78+t492*a88+t494*a98;
		double t540 = t476*a9+t478*a19+t480*a29+t482*a39+t484*a49+t486*a59+t488*a69+t490*a79+t492*a89+t494*a99;
		double t541 = a5*a5;
		double t543 = a15*a15;
		double t545 = a25*a25;
		double t547 = a35*a35;
		double t549 = a45*a45;
		double t551 = a55*a55;
		double t553 = a65*a65;
		double t555 = a75*a75;
		double t557 = a85*a85;
		double t559 = a95*a95;
		double t562 = a5*s0;
		double t564 = a15*s1;
		double t566 = a25*s2;
		double t568 = a35*s3;
		double t570 = a45*s4;
		double t572 = a55*s5;
		double t574 = a65*s6;
		double t576 = a75*s7;
		double t578 = a85*s8;
		double t580 = a95*s9;
		double t582 = t562*a6+t564*a16+t566*a26+t568*a36+t570*a46+t572*a56+t574*a66+t576*a76+t578*a86+t580*a96;
		double t593 = t562*a7+t564*a17+t566*a27+t568*a37+t570*a47+t572*a57+t574*a67+t576*a77+t578*a87+t580*a97;
		double t604 = t562*a8+t564*a18+t566*a28+t568*a38+t570*a48+t572*a58+t574*a68+t576*a78+t578*a88+t580*a98;
		double t615 = t562*a9+t564*a19+t566*a29+t568*a39+t570*a49+t572*a59+t574*a69+t576*a79+t578*a89+t580*a99;
		double t616 = a6*a6;
		double t618 = a16*a16;
		double t620 = a26*a26;
		double t622 = a36*a36;
		double t624 = a46*a46;
		double t626 = a56*a56;
		double t628 = a66*a66;
		double t630 = a76*a76;
		double t632 = a86*a86;
		double t634 = a96*a96;
		double t637 = a6*s0;
		double t639 = a16*s1;
		double t641 = a26*s2;
		double t643 = a36*s3;
		double t645 = a46*s4;
		double t647 = a56*s5;
		double t649 = a66*s6;
		double t651 = a76*s7;
		double t653 = a86*s8;
		double t655 = a96*s9;
		double t657 = t637*a7+t639*a17+t641*a27+t643*a37+t645*a47+t647*a57+t649*a67+t651*a77+t653*a87+t655*a97;
		double t668 = t637*a8+t639*a18+t641*a28+t643*a38+t645*a48+t647*a58+t649*a68+t651*a78+t653*a88+t655*a98;
		double t679 = t637*a9+t639*a19+t641*a29+t643*a39+t645*a49+t647*a59+t649*a69+t651*a79+t653*a89+t655*a99;
		double t680 = a7*a7;
		double t682 = a17*a17;
		double t684 = a27*a27;
		double t686 = a37*a37;
		double t688 = a47*a47;
		double t690 = a57*a57;
		double t692 = a67*a67;
		double t694 = a77*a77;
		double t696 = a87*a87;
		double t698 = a97*a97;
		double t701 = a7*s0;
		double t703 = a17*s1;
		double t705 = a27*s2;
		double t707 = a37*s3;
		double t709 = a47*s4;
		double t711 = a57*s5;
		double t713 = a67*s6;
		double t715 = a77*s7;
		double t717 = a87*s8;
		double t719 = a97*s9;
		double t721 = t701*a8+t703*a18+t705*a28+t707*a38+t709*a48+t711*a58+t713*a68+t715*a78+t717*a88+t719*a98;
		double t732 = t701*a9+t703*a19+t705*a29+t707*a39+t709*a49+t711*a59+t713*a69+t715*a79+t717*a89+t719*a99;
		double t733 = a8*a8;
		double t735 = a18*a18;
		double t737 = a28*a28;
		double t739 = a38*a38;
		double t741 = a48*a48;
		double t743 = a58*a58;
		double t745 = a68*a68;
		double t747 = a78*a78;
		double t749 = a88*a88;
		double t751 = a98*a98;
		double t774 = a8*s0*a9+a18*s1*a19+a28*s2*a29+a38*s3*a39+a48*s4*a49+a58*s5*a59+a68*s6*a69+a78*s7*a79+a88*s8*a89+a98*s9*a99;
		double t775 = a9*a9;
		double t777 = a19*a19;
		double t779 = a29*a29;
		double t781 = a39*a39;
		double t783 = a49*a49;
		double t785 = a59*a59;
		double t787 = a69*a69;
		double t789 = a79*a79;
		double t791 = a89*a89;
		double t793 = a99*a99;

		tmp_M_ptr = dst_mat->SM_ptr;

		*tmp_M_ptr++ = t1*s0+t3*s1+t5*s2+t7*s3+t9*s4+t11*s5+t13*s6+t15*s7+t17*s8+t19*s9;
		*tmp_M_ptr++ = t42;
		*tmp_M_ptr++ = t53;
		*tmp_M_ptr++ = t64;
		*tmp_M_ptr++ = t75;
		*tmp_M_ptr++ = t86;
		*tmp_M_ptr++ = t97;
		*tmp_M_ptr++ = t108;
		*tmp_M_ptr++ = t119;
		*tmp_M_ptr++ = t130;
		*tmp_M_ptr++ = t131*s0+t133*s1+t135*s2+t137*s3+t139*s4+t141*s5+t143*s6+t145*s7+t147*s8+t149*s9;
		*tmp_M_ptr++ = t172;
		*tmp_M_ptr++ = t183;
		*tmp_M_ptr++ = t194;
		*tmp_M_ptr++ = t205;
		*tmp_M_ptr++ = t216;
		*tmp_M_ptr++ = t227;
		*tmp_M_ptr++ = t238;
		*tmp_M_ptr++ = t249;
		*tmp_M_ptr++ = t250*s0+t252*s1+t254*s2+t256*s3+t258*s4+t260*s5+t262*s6+t264*s7+t266*s8+t268*s9;
		*tmp_M_ptr++ = t291;
		*tmp_M_ptr++ = t302;
		*tmp_M_ptr++ = t313;
		*tmp_M_ptr++ = t324;
		*tmp_M_ptr++ = t335;
		*tmp_M_ptr++ = t346;
		*tmp_M_ptr++ = t357;
		*tmp_M_ptr++ = t358*s0+t360*s1+t362*s2+t364*s3+t366*s4+t368*s5+t370*s6+t372*s7+t374*s8+t376*s9;
		*tmp_M_ptr++ = t399;
		*tmp_M_ptr++ = t410;
		*tmp_M_ptr++ = t421;
		*tmp_M_ptr++ = t432;
		*tmp_M_ptr++ = t443;
		*tmp_M_ptr++ = t454;
		*tmp_M_ptr++ = t455*s0+t457*s1+t459*s2+t461*s3+t463*s4+t465*s5+t467*s6+t469*s7+t471*s8+t473*s9;
		*tmp_M_ptr++ = t496;
		*tmp_M_ptr++ = t507;
		*tmp_M_ptr++ = t518;
		*tmp_M_ptr++ = t529;
		*tmp_M_ptr++ = t540;
		*tmp_M_ptr++ = t541*s0+t543*s1+t545*s2+t547*s3+t549*s4+t551*s5+t553*s6+t555*s7+t557*s8+t559*s9;
		*tmp_M_ptr++ = t582;
		*tmp_M_ptr++ = t593;
		*tmp_M_ptr++ = t604;
		*tmp_M_ptr++ = t615;
		*tmp_M_ptr++ = t616*s0+t618*s1+t620*s2+t622*s3+t624*s4+t626*s5+t628*s6+t630*s7+t632*s8+t634*s9;
		*tmp_M_ptr++ = t657;
		*tmp_M_ptr++ = t668;
		*tmp_M_ptr++ = t679;
		*tmp_M_ptr++ = t680*s0+t682*s1+t684*s2+t686*s3+t688*s4+t690*s5+t692*s6+t694*s7+t696*s8+t698*s9;
		*tmp_M_ptr++ = t721;
		*tmp_M_ptr++ = t732;
		*tmp_M_ptr++ = t733*s0+t735*s1+t737*s2+t739*s3+t741*s4+t743*s5+t745*s6+t747*s7+t749*s8+t751*s9;
		*tmp_M_ptr++ = t774;
		*tmp_M_ptr++ = t775*s0+t777*s1+t779*s2+t781*s3+t783*s4+t785*s5+t787*s6+t789*s7+t791*s8+t793*s9;
	}
	else {
		CvMat* D = cvCreateMat(M_size, 1, CV_64F);
		CvMat* U = cvCreateMat(M_size, M_size, CV_64F);

		this->SVD(D, U);

		CvMat* U_t = cvCreateMat(M_size, M_size, CV_64F);
		CvMat* U_tU = cvCreateMat(M_size, M_size, CV_64F);

		cvTranspose(U, U_t);

		/* D = exp(D) */
		Exp_Diag(D);

		/* U_t = U_t * D */
		Mul_Diag(U_t, D);

		/* U_tU = U_t * U */
		cvMatMul(U_t, U, U_tU);

		dst_mat->SetData((double*)U_tU->data.ptr, M_size2);

		cvReleaseMat(&U_t);
		cvReleaseMat(&U_tU);

		cvReleaseMat(&D);
		cvReleaseMat(&U);
	}
}

double CCovarianceMatrix::GeodesicDistance(CCovarianceMatrix *Y_mat, CCovarianceMatrix *this_mat_inv_sqrt)
{
	CCovarianceMatrix *log_Y = new CCovarianceMatrix(Y_mat->M_size);
	CCovarianceMatrix *log_2_mul_XYX = new CCovarianceMatrix(Y_mat->M_size);

	CCovarianceMatrix *inv_sqrt_X = this_mat_inv_sqrt;
	if ( !inv_sqrt_X ) {
		inv_sqrt_X = new CCovarianceMatrix(Y_mat->M_size);
		this->Inverse(inv_sqrt_X);
		inv_sqrt_X->Sqrt();
	}

	inv_sqrt_X->Mul2(Y_mat, log_2_mul_XYX);
	log_2_mul_XYX->Log();
	log_2_mul_XYX->Square();

	double trace = log_2_mul_XYX->Trace();

	delete log_Y;
	delete log_2_mul_XYX;

	if ( !this_mat_inv_sqrt )
		delete inv_sqrt_X;
	return trace;
}

void CCovarianceMatrix::SetData(double *data, int length, double *avg_data)
{
	if ( length == SM_length ) {
		memcpy(SM_ptr, data, sizeof(double)*SM_length);
	}
	else if ( length == M_size2 ) {
		if ( this->M_size == 2 ) {
			SM_ptr[0] = data[0];
			SM_ptr[1] = data[1];
			SM_ptr[2] = data[3];
		}
		else if ( this->M_size == 3 ) {
			SM_ptr[0] = data[0];
			SM_ptr[1] = data[1];
			SM_ptr[2] = data[2];
			SM_ptr[3] = data[4];
			SM_ptr[4] = data[5];
			SM_ptr[5] = data[8];
		}
		else {
			double *tmp_M_ptr = data;
			double *tmp_SM_ptr = SM_ptr;

			for ( int dy = 0 ; dy < M_size ; dy++ ) {
				tmp_M_ptr += dy;
				for ( int dx = dy ; dx < M_size ; dx++ )
					*tmp_SM_ptr++ = *tmp_M_ptr++;
			}
		}
	}
	else {
		printf("Set symmetric matrix data error!\n");
		exit(1);
	}

	if ( avg_data && AVG_ptr ) {
		memcpy(AVG_ptr, avg_data, sizeof(double)*M_size);
	}
}

double* CCovarianceMatrix::GetData(int one_upper_zero_full)
{
	double *data = NULL;
	if ( one_upper_zero_full ) {
		data = new double[SM_length];
		memcpy(data, SM_ptr, sizeof(double)*SM_length);
	}
	else {
		data = new double[M_size2];

		if ( this->M_size == 1 ) {
			data[0] = SM_ptr[0];
		}
		else if ( this->M_size == 2 ) {
			data[0] = SM_ptr[0];
			data[1] = data[2] = SM_ptr[1];
			data[3] = SM_ptr[2];
		}
		else if ( this->M_size == 3 ) {
			data[0] = SM_ptr[0];
			data[1] = data[3] = SM_ptr[1];
			data[2] = data[6] = SM_ptr[2];
			data[4] = SM_ptr[3];
			data[5] = data[7] = SM_ptr[4];
			data[8] = SM_ptr[5];
		}
		else {
			double *tmp_SM_ptr = SM_ptr;
			double *tmp_M_ptr1 = data;
			double *tmp_M_ptr2 = data;

			for ( int dy = 0 ; dy < M_size ; dy++ ) {
				tmp_M_ptr1 += dy;
				tmp_M_ptr2 = data + dy*M_size + dy;
				for ( int dx = dy ; dx < M_size ; dx++ ) {
					*tmp_M_ptr1++ = *tmp_SM_ptr;
					*tmp_M_ptr2 = *tmp_SM_ptr++;
					tmp_M_ptr2 += M_size;
				}
			}
		}
	}
	return data;
}

void CCovarianceMatrix::Log_Diag()
{
	double *tmp_SM_ptr = SM_ptr;
	for ( int dy = M_size ; dy > 0 ; dy-- ) {
		if ( *tmp_SM_ptr <= 0 ) {
			*tmp_SM_ptr = MIN_DOUBLE_NUM;
			printf("Log function error in Log_Diag()!\n");
			//exit(1);
		}
		*tmp_SM_ptr = log(*tmp_SM_ptr);
		tmp_SM_ptr += dy;
	}
}

void CCovarianceMatrix::Exp_Diag()
{
	double *tmp_SM_ptr = SM_ptr;
	for ( int dy = M_size ; dy > 0 ; dy-- ) {
		*tmp_SM_ptr = exp(*tmp_SM_ptr);
		tmp_SM_ptr += dy;
	}
}

void CCovarianceMatrix::Sqrt_Diag()
{
	double *tmp_SM_ptr = SM_ptr;
	for ( int dy = M_size ; dy > 0 ; dy-- ) {
		if ( *tmp_SM_ptr <= 0 ) {
			*tmp_SM_ptr = MIN_DOUBLE_NUM;
			printf("sqrt error in Sqrt_Diag!\n");
			//exit(1);
		}
		*tmp_SM_ptr = sqrt(*tmp_SM_ptr);
		tmp_SM_ptr += dy;
	}
}

void CCovarianceMatrix::Sqrt_Diag(CvMat *D)
{
	double *tmp_M_ptr = (double*)(D->data.ptr);
	int dy;
	for ( dy = 0 ; dy < M_size ; dy++ ) {
		if ( *tmp_M_ptr <= 0 ) {
			*tmp_M_ptr = MIN_DOUBLE_NUM;
			printf("sqrt error in Sqrt_Diag(CvMat *D)!\n");
			//exit(1);
		}
		*tmp_M_ptr = sqrt(*tmp_M_ptr);
		tmp_M_ptr++;
	}
}

void CCovarianceMatrix::Exp_Diag(CvMat *D)
{
	double *tmp_M_ptr = (double*)(D->data.ptr);
	int dy;
	for ( dy = 0 ; dy < M_size ; dy++ ) {
		*tmp_M_ptr = exp(*tmp_M_ptr);
		tmp_M_ptr++;
	}
}

void CCovarianceMatrix::Log_Diag(CvMat *D)
{
	double *tmp_M_ptr = (double*)(D->data.ptr);
	int dy;
	for ( dy = 0 ; dy < M_size ; dy++ ) {
		if ( *tmp_M_ptr <= 0 ) {
			*tmp_M_ptr = MIN_DOUBLE_NUM;
			printf("Log function error in Log_Diag(CvMat *D)!\n");
			//exit(1);
		}
		*tmp_M_ptr = log(*tmp_M_ptr);
		tmp_M_ptr++;
	}
}

void CCovarianceMatrix::Mul_Diag(CvMat *A, CvMat *D)
{
	int dy, dx;

	double *tmp_diag_vec;

	double *tmp_M_ptr = (double*)A->data.ptr;
	for ( dy = 0 ; dy < M_size ; dy++ ) {
		tmp_diag_vec = (double*)D->data.ptr;
		for ( dx = 0 ; dx < M_size ; dx++ ) {
			*tmp_M_ptr++ *= *tmp_diag_vec++;
		}
	}
}

void CCovarianceMatrix::Update(CCovarianceMatrix *Y_mat, double update_rate, CCovarianceMatrix *updated_mat)
{
	CCovarianceMatrix *log_X_X = new CCovarianceMatrix(M_size);
	CCovarianceMatrix *log_X_Y = new CCovarianceMatrix(M_size);

	this->Log_Y(this, log_X_X);
	this->Log_Y(Y_mat, log_X_Y);

	double one_update_rate = 1.0 - update_rate;

	int a;
	double *tmp_mat1 = log_X_X->SM_ptr;
	double *tmp_mat2 = log_X_Y->SM_ptr;
	for ( a = 0 ; a < SM_length ; a++ ) {
		*tmp_mat1 = (*tmp_mat1) * one_update_rate + (*tmp_mat2++) * update_rate;
		tmp_mat1++;
	}

	this->Exp_Y(log_X_X, updated_mat);

	delete log_X_X;
	delete log_X_Y;
}

double CCovarianceMatrix::CovDistance(CCovarianceMatrix *Y_mat)
{
	/* compute the generalized eigen-values */
	double *eig_v = this->GeneralizedEigenvalues(Y_mat);

	/* compute the covariance distance */
	double *_eig_v = eig_v;
	double cov_dist = 0;
	int d;
	double ln_eig_v;
	for ( d = 0 ; d < M_size; d++ ) {
		if ( *_eig_v > 0 ) {
			ln_eig_v = log(*_eig_v++);
			cov_dist += ln_eig_v*ln_eig_v;
		}
		else
			cov_dist += 99999.9;
	}
	cov_dist = sqrt(cov_dist);

	delete [] eig_v;

	return cov_dist;

}

double* CCovarianceMatrix::GeneralizedEigenvalues(CCovarianceMatrix *Y_mat)
{
	CvMat *F = cvCreateMat(M_size, M_size, CV_64F);

	/* do Cholesky factorization for Y, i.e. Y = F * F^t */
	Y_mat->CholeskyFactorization(F);

	CvMat *X = cvCreateMat(M_size, M_size, CV_64F);
	CvMat *Z = cvCreateMat(M_size, M_size, CV_64F);
	CvMat *invF = cvCreateMat(M_size, M_size, CV_64F);
	CvMat *invF_tran = cvCreateMat(M_size, M_size, CV_64F);
	CvMat *invFX = cvCreateMat(M_size, M_size, CV_64F);

	double *M_ptr = this->GetData(0);
	memcpy(X->data.ptr, M_ptr, sizeof(double)*M_size2);
	delete [] M_ptr;

	/* compute the inverse of F and the inverse's transpose of F */
	//cvInv(F, invF, CV_SVD_SYM);
	cvInv(F, invF);
	cvTranspose(invF, invF_tran);

	/* compute Z = F^{-1} * X * F^{-t} */
	cvMatMul(invF, X, invFX);
	cvMatMul(invFX, invF_tran, Z);

	/* compute the eigenvalues for the matrix Z */
	CvMat *diag = cvCreateMat(M_size, 1, CV_64F);
	//cvSVD(Z, diag, NULL, NULL, CV_SVD_MODIFY_A);
	//cvSVD(Z, diag);
	cvEigenVV(Z, NULL, diag);

	double *eig_v = new double[M_size];
	memcpy(eig_v, diag->data.ptr, sizeof(double)*M_size);

	/* release memories */
	cvReleaseMat(&F);
	cvReleaseMat(&X);
	cvReleaseMat(&Z);
	cvReleaseMat(&invF);
	cvReleaseMat(&invF_tran);
	cvReleaseMat(&invFX);
	cvReleaseMat(&diag);

	return eig_v;
}

void CCovarianceMatrix::CholeskyFactorization_2(CvMat *F)
{
	double *A_data = this->GetData(0);
	double *L_data = new double[M_size2];

	memset( L_data, 0, M_size2 * sizeof( double ));

	int i,j,k;
	double S;

	for( j=0 ; j<M_size ; j++ )
	{
		S=0;
		for( k=0 ; k<j ; k++ )
			S += powf(L_data[j*M_size+k],2);
		L_data[j*M_size+j] = sqrt(A_data[j*M_size+j] - S);
		for( i=j ; i<M_size ; i++ )
		{
			S = 0;
			for( k=0 ; k<j ; k++ )
				S += L_data[i*M_size+k] * L_data[j*M_size+k];

			if (L_data[j*M_size+j]==0.0)  {
			      printf("Cholesky Factorization failure!\n");
			      exit(1);
			}
			L_data[i*M_size+j] = (A_data[i*M_size+j] - S) / L_data[j*M_size+j];
		}
	}

	memcpy(F->data.ptr, L_data, sizeof(double)*M_size2);

	delete [] A_data;
	delete [] L_data;
}

void CCovarianceMatrix::CholeskyFactorization(CvMat *F)
{
	double *s = new double[SM_length];
	double *t = new double[SM_length];
	double *ti, *tj, *tk, d, sum;
	int i, j, k;

	double *M_ptr = this->GetData(0);
	double *s1 = s;
	double *s2 = M_ptr;
	int shift = M_size-1;
	for ( i=0 ; i<M_size; i++ ) {
		for ( j = 0 ; j <= i ; j++ )
			*s1++ = (double)*s2++;
		s2 += shift--;
	}
	s1 = s;

	ti = t;
	for ( i=0; i<M_size; i++)
	{
		tj = t;
		for ( j=0; j<i; j++)
		{
			tk = ti; sum = 0.0; k = j;
			while (k--) { sum += *tj++ * *tk++; }
			*tk = (*s++ - sum) / *tj++;
		}
		sum = 0.0; k = i;
		while (k--) { sum += *ti * *ti; ti++; }
		d = *s++ - sum;
		if (d<=0.0)  {
		      printf("Cholesky Factorization failure!\n");
		      exit(1);
		}
		*ti++ = sqrt(d);
	}

	s = s1;

	cvSetZero(F);
	ti = (double*)F->data.ptr;
	s1 = t;
	shift = M_size-1;
	for ( i=0 ; i<M_size; i++ ) {
		for ( j = 0 ; j <= i ; j++ )
			*ti++ = *s1++;
		ti += shift--;
	}

	delete [] M_ptr;
	delete [] s;
	delete [] t;
}

void CCovarianceMatrix::GetNormDiagonalElems(double *norm_diag_elems)
{
	double *tmp_diag = norm_diag_elems;
	double* tmp_SM_ptr = SM_ptr;
	for ( int dy = M_size ; dy > 0 ; dy-- ) {
		if ( *tmp_SM_ptr <= 0 ) {
			*tmp_SM_ptr = MIN_DOUBLE_NUM;
			printf("sqrt error in GetNormDiagonalElems!\n");
			//exit(1);
		}
		*tmp_diag++ = sqrt(*tmp_SM_ptr);
		tmp_SM_ptr += dy;
	}
}

void CCovarianceMatrix::GetDiagonalElems(double *diag_elems)
{
	double *tmp_diag = diag_elems;
	double* tmp_SM_ptr = SM_ptr;
	for ( int dy = M_size ; dy > 0 ; dy-- ) {
		*tmp_diag++ = *tmp_SM_ptr;
		tmp_SM_ptr += dy;
	}
}

void CCovarianceMatrix::NormalizeIllumination(double *parent_win_diag_elems, double *parent_avg)
{
	/* normalize the rows and columns corresponding to position (x,y)-features */
	int dy, dx;
	double* tmp_SM_ptr = SM_ptr;

	double *tmp_dx, *tmp_dy;
	tmp_dy = parent_win_diag_elems;

	if ( parent_avg && this->AVG_ptr ) {
		double* tmp_AVG_ptr = AVG_ptr;
        double* tmp_parent_avg_ptr = parent_avg;

		for ( dy = 0 ; dy < M_size ; dy++ ) {
			tmp_dx = parent_win_diag_elems + dy;
			*tmp_AVG_ptr++ /= *tmp_parent_avg_ptr++;
			for ( dx = dy ; dx < M_size ; dx++ )  {
				*tmp_SM_ptr++ /= (*tmp_dx++) * (*tmp_dy);
			}
			tmp_dy++;
		}
	}
	else {
		for ( dy = 0 ; dy < M_size ; dy++ ) {
			tmp_dx = parent_win_diag_elems + dy;
			for ( dx = dy ; dx < M_size ; dx++ )  {
				*tmp_SM_ptr++ /= (*tmp_dx++) * (*tmp_dy);
			}
			tmp_dy++;
		}
	}
}

void CCovarianceMatrix::NormalizeIllumination(CCovarianceMatrix *parent_win_mat)
{
	double *norm_diag = new double[M_size];
	parent_win_mat->GetNormDiagonalElems(norm_diag);
	this->NormalizeIllumination(norm_diag, parent_win_mat->AVG_ptr);
	delete [] norm_diag;
}

void CCovarianceMatrix::NormalizePositionFeatures(CvRect roi, int x_feature_idx, int y_feature_idx)
{
	double width = (double)roi.width;
	double height = (double)roi.height;

	int dy, dx;

	/* normalize the rows and columns corresponding to position (x,y)-features */
	double* tmp_SM_ptr = SM_ptr;
	for ( dy = 0 ; dy < M_size ; dy++ ) {
		for ( dx = dy ; dx < M_size ; dx++ )  {
			if ( dy == x_feature_idx )
				*tmp_SM_ptr /= width;
			if ( dy == y_feature_idx )
				*tmp_SM_ptr /= height;
			if ( dx == x_feature_idx )
				*tmp_SM_ptr /= width;
			if ( dx == y_feature_idx )
				*tmp_SM_ptr /= height;
			tmp_SM_ptr++;
		}
	}
	/*
	if ( AVG_ptr ) {
        if ( x_feature_idx >= 0 )
            AVG_ptr[x_feature_idx] = 0.5;
        if ( y_feature_idx >= 0 )
            AVG_ptr[y_feature_idx] = 0.5;
	}
	*/
}

void CCovarianceMatrix::ComputeRegionCovFeature(IplImage **feature_imgs, CvRect roi, IplImage *mask_img)
{
	int yx, d, dx, dy;

	/************************************************************************/
	/* covert Opencv IplImage data to one-dimensional double data            */
	/************************************************************************/

	COpencvDataConversion<float,float> ODC1;
	float** feature_data = new float*[M_size];
	for ( d = 0 ; d < M_size ; d++ ) {
		cvSetImageROI(feature_imgs[d], roi);
		feature_data[d] = ODC1.GetImageData(feature_imgs[d]);
		cvResetImageROI(feature_imgs[d]);
	}

	uchar*	mask_data;
	int mask_sum = 0;
	if ( mask_img ) {
		COpencvDataConversion<uchar,uchar> ODC2;
		cvSetImageROI(mask_img, roi);
		CvScalar msum = cvSum(mask_img);
		mask_sum = (int)msum.val[0];
		mask_data = ODC2.GetImageData(mask_img);
		cvResetImageROI(mask_img);
	}
	else {
		mask_data = new uchar[roi.width*roi.height];
		mask_sum = roi.width*roi.height;
		memset(mask_data, 1, sizeof(uchar)*mask_sum);
	}

	/************************************************************************/
	/*  compute the d-dimensional first-order p_vec,			*/
	/*     and d*d-dimensional second-order Q_mat				*/
	/************************************************************************/

	double *roi_p_vec = new double[M_size];
	double *roi_Q_mat = new double[SM_length];
	memset(roi_p_vec, 0, sizeof(double)*M_size);
	memset(roi_Q_mat, 0, sizeof(double)*SM_length);

	double *tmp_p_vec, *tmp_Q_mat;
	double *cur_f_vec_dy, *cur_f_vec_dx, *cur_f_vec, *tmp_f_vec;

	cur_f_vec = new double[M_size];

	float** _feature_data = new float*[M_size];
	for ( d = 0 ; d < M_size ; d++ )
		_feature_data[d] = feature_data[d];

	uchar *cur_mask_data = mask_data;

	for ( yx = 0 ; yx < roi.height*roi.width ; yx++ ) {
		if ( *cur_mask_data++ ) {
			tmp_p_vec = roi_p_vec;
			tmp_Q_mat = roi_Q_mat;
			tmp_f_vec = cur_f_vec;

			for ( d = 0 ; d < M_size ; d++ ) {
				*tmp_f_vec = (double)(*_feature_data[d]++);
				*tmp_p_vec++ += *tmp_f_vec++;
			}

			cur_f_vec_dy = cur_f_vec;
			for ( dy = 0 ; dy < M_size ; dy++ ) {
				cur_f_vec_dx = cur_f_vec + dy;
				for ( dx = dy ; dx < M_size ; dx++ ) {
					*(tmp_Q_mat++) += (*cur_f_vec_dy)*(*cur_f_vec_dx++);
				}
				cur_f_vec_dy++;
			}
		}
	}

	double c1 = (double)(mask_sum)-1.0;
	double c2 = c1*(c1+1);

	double *cov_Q = this->SM_ptr;

	double *p_vec_dx, *p_vec_dy;
	tmp_Q_mat = roi_Q_mat;
	p_vec_dy = roi_p_vec;
	for ( dy = 0 ; dy < M_size ; dy++ ) {
		p_vec_dx = roi_p_vec + dy;
		for ( dx = dy ; dx < M_size ; dx++ ) {
			*cov_Q++ = (double) ( (double)(*tmp_Q_mat++)/c1 - ((double)(*p_vec_dy)*(double)(*p_vec_dx++))/c2 );
		}
		p_vec_dy++;
	}

	CCovarianceMatrix* small_identity_mat = new CCovarianceMatrix(M_size);
	small_identity_mat->SetIdentityMatrix(MIN_IDENTITY_VALUE);
	this->Add(small_identity_mat);

	delete small_identity_mat;

	/************************************************************************/
	/*  release memories                                                    */
	/************************************************************************/
	for ( d = 0 ; d < M_size ; d++ )
		delete [] feature_data[d];
	delete [] _feature_data;
	delete [] feature_data;
	delete [] mask_data;
	delete [] roi_p_vec;
	delete [] roi_Q_mat;
	delete [] cur_f_vec;
}

double CCovarianceMatrix::GeodesicDistanceSum(CCovarianceMatrix **Y_mats, int mat_num, double *weights, int *sel_pts, CCovarianceMatrix *this_mat_inv_sqrt)
{
	double dist_sum = 0;
	int i;
	CCovarianceMatrix *inv_sqrt_X = this_mat_inv_sqrt;
	if ( !inv_sqrt_X ) {
		inv_sqrt_X = new CCovarianceMatrix(this->M_size);
		this->Inverse(inv_sqrt_X);
		inv_sqrt_X->Sqrt();
	}

	if ( sel_pts ) {
		if ( weights ) {
			for ( i = 0 ; i < mat_num ; i++ ) {
				if ( sel_pts[i] )
					dist_sum += weights[i]*this->GeodesicDistance(Y_mats[i], inv_sqrt_X);
			}
		}
		else {
			for ( i = 0 ; i < mat_num ; i++ )
				if ( sel_pts[i] )
					dist_sum += this->GeodesicDistance(Y_mats[i], inv_sqrt_X);
		}
	}
	else {
		if ( weights ) {
			for ( i = 0 ; i < mat_num ; i++ )
				dist_sum += weights[i]*this->GeodesicDistance(Y_mats[i], inv_sqrt_X);
		}
		else {
			for ( i = 0 ; i < mat_num ; i++ )
				dist_sum += this->GeodesicDistance(Y_mats[i], inv_sqrt_X);
		}
	}

	if ( !this_mat_inv_sqrt )
		delete inv_sqrt_X;

	return dist_sum;
}

void CCovarianceMatrix::WeightedMean(CCovarianceMatrix **Y_mats, int mat_num, double *weights, int *sel_pts, int itr_num, double *dist_sums, double itr_stop_value)
{
	int i, j;

	/* initialization for gradient descent iteration */
	double *tmp_weights;
	double tot_weight = 0;
	if ( weights ) {	/* initialization with max weight's covariance matrix */
		int max_weight_idx = 0;
		double max_weight = 0;
		tmp_weights = weights;
		for ( i = 0 ; i < mat_num ; i++ ) {
			if ( sel_pts ) {
				if ( sel_pts[i] ) {
					if ( *tmp_weights > max_weight ) {
						max_weight = *tmp_weights;
						max_weight_idx = i;
					}
					tot_weight += *tmp_weights++;
				}
				else
					tmp_weights++;
			}
			else {
				if ( *tmp_weights > max_weight ) {
					max_weight = *tmp_weights;
					max_weight_idx = i;
				}
				tot_weight += *tmp_weights++;
			}
		}
		Y_mats[max_weight_idx]->Copy(this);
	}
	else			/* initialization with first covariance matrix */
		Y_mats[0]->Copy(this);

	if ( tot_weight > 0 )
		tot_weight = 1.0 / tot_weight;

	CCovarianceMatrix *log_Y = new CCovarianceMatrix(M_size);
	CCovarianceMatrix *sqrt_X = new CCovarianceMatrix(M_size);
	CCovarianceMatrix *inv_sqrt_X = new CCovarianceMatrix(M_size);
	CCovarianceMatrix *sum_log_Y = new CCovarianceMatrix(M_size);
	double one_div_N = 1.0/(double)mat_num;

	CCovarianceMatrix *prev_weighted_mean = new CCovarianceMatrix(M_size);
	this->Copy(prev_weighted_mean);

	/* computing mean avg */
	CCovarianceMatrix **tmp_Y_mats = Y_mats;
	double *tmp_AVG_ptr = this->AVG_ptr;
	double *tmp_Y_AVG_ptr;
	if ( this->AVG_ptr )
        for ( i = 0 ; i < M_size ; i++ )
            *tmp_AVG_ptr++ = 0;
	if ( sel_pts && weights ) {
		tmp_weights = weights;
		double weight_sum = 0;
		for ( i = 0 ; i < mat_num ; i++ ) {
			if ( sel_pts[i] ) {
			    if ( this->AVG_ptr ) {
                    tmp_Y_AVG_ptr = (*tmp_Y_mats)->AVG_ptr;
                    tmp_AVG_ptr = this->AVG_ptr;
                    for ( j = 0 ; j < M_size ; j++ )
                        *tmp_AVG_ptr++ += (*tmp_weights) * (*tmp_Y_AVG_ptr++);
			    }
				weight_sum += *tmp_weights;
			}
			tmp_Y_mats++;
			tmp_weights++;
		}
		if ( this->AVG_ptr ) {
            tmp_AVG_ptr = this->AVG_ptr;
            for ( j = 0 ; j < M_size ; j++ )
                *tmp_AVG_ptr++ /= weight_sum;
		}
	}
	else if ( sel_pts && !weights ) {
		double weight_sum = 0;
		for ( i = 0 ; i < mat_num ; i++ ) {
			if ( sel_pts[i] ) {
                if ( this->AVG_ptr ) {
                    tmp_Y_AVG_ptr = (*tmp_Y_mats)->AVG_ptr;
                    tmp_AVG_ptr = this->AVG_ptr;
                    for ( j = 0 ; j < M_size ; j++ )
                        *tmp_AVG_ptr++ += (*tmp_Y_AVG_ptr++);
                }
				weight_sum += 1.0;
			}
			tmp_Y_mats++;
		}
		if ( this->AVG_ptr ) {
            tmp_AVG_ptr = this->AVG_ptr;
            for ( j = 0 ; j < M_size ; j++ )
                *tmp_AVG_ptr++ /= weight_sum;
		}
	}
	else if ( !sel_pts && weights ) {
		tmp_weights = weights;
		double weight_sum = 0;
		for ( i = 0 ; i < mat_num ; i++ ) {
            if ( this->AVG_ptr ) {
                tmp_Y_AVG_ptr = (*tmp_Y_mats)->AVG_ptr;
                tmp_AVG_ptr = this->AVG_ptr;
                for ( j = 0 ; j < M_size ; j++ )
                    *tmp_AVG_ptr++ += (*tmp_weights) * (*tmp_Y_AVG_ptr++);
            }
			weight_sum += *tmp_weights;
			tmp_Y_mats++;
			tmp_weights++;
		}
		if ( this->AVG_ptr ) {
            tmp_AVG_ptr = this->AVG_ptr;
            for ( j = 0 ; j < M_size ; j++ )
                *tmp_AVG_ptr++ /= weight_sum;
		}
	}
	else {
		double weight_sum = 0;
		for ( i = 0 ; i < mat_num ; i++ ) {
            if ( this->AVG_ptr ) {
                tmp_Y_AVG_ptr = (*tmp_Y_mats)->AVG_ptr;
                tmp_AVG_ptr = this->AVG_ptr;
                for ( j = 0 ; j < M_size ; j++ )
                    *tmp_AVG_ptr++ += (*tmp_Y_AVG_ptr++);
            }
			weight_sum += 1.0;
			tmp_Y_mats++;
		}
		if ( this->AVG_ptr ) {
            tmp_AVG_ptr = this->AVG_ptr;
            for ( j = 0 ; j < M_size ; j++ )
                *tmp_AVG_ptr++ /= weight_sum;
		}
	}

	/* iteration starting */
	for ( i = 0 ; i < itr_num ; i++ ) {
		sum_log_Y->SetZero();

		tmp_Y_mats = Y_mats;

		this->Sqrt(sqrt_X);

		sqrt_X->Inverse(inv_sqrt_X);

		if ( weights ) {	/* different weights for different covariance matrices */
			tmp_weights = weights;
			for ( j = 0 ; j < mat_num ; j++ ) {
				//this->Log_Y(*tmp_Y_mats++, log_Y);

				if ( sel_pts && !sel_pts[j] ) {
					tmp_Y_mats++;
					tmp_weights++;
					continue;
				}

				inv_sqrt_X->Mul2(*tmp_Y_mats++, log_Y);
				log_Y->Log();

				log_Y->Scale(*tmp_weights++);
				sum_log_Y->Add(log_Y);
			}
			sum_log_Y->Scale(tot_weight);
		}
		else {			/* the same weight for all covariance matrices, 1/N */
			for ( j = 0 ; j < mat_num ; j++ ) {
				//this->Log_Y(*tmp_Y_mats++, log_Y);

				inv_sqrt_X->Mul2(*tmp_Y_mats++, log_Y);
				log_Y->Log();

				sum_log_Y->Add(log_Y);
			}

			sum_log_Y->Scale(one_div_N);
		}

		sum_log_Y->Exp();
		sqrt_X->Mul2(sum_log_Y, this);

		//this->Exp_Y(sum_log_Y);

		if ( dist_sums ) {
			dist_sums[i] = this->GeodesicDistanceSum(Y_mats, mat_num, weights, sel_pts, inv_sqrt_X);
			//printf("\t\tWeightedMean itr %d is OK (dist = %.6f)\n", i, dist_sums[i]);
			if ( i > 0 && dist_sums[i-1]-dist_sums[i] < 0 ) {
				prev_weighted_mean->Copy(this);
				break;
			}
			if ( i > 0 && dist_sums[i-1]-dist_sums[i] < itr_stop_value ) {
				break;
			}
		}
		this->Copy(prev_weighted_mean);
	}

	delete log_Y;
	delete sum_log_Y;
	delete sqrt_X;
	delete inv_sqrt_X;
	delete prev_weighted_mean;
}

void CCovarianceMatrix::Add(CCovarianceMatrix *Y_mat, CCovarianceMatrix *sum_mat)
{
	if ( sum_mat ) {
		this->Copy(sum_mat);
		int i;
		double *tmp_sum_SM_ptr = sum_mat->SM_ptr, *tmp_Y_SM_ptr = Y_mat->SM_ptr;
		for ( i = 0 ; i < SM_length ; i++ )
			*tmp_sum_SM_ptr++ += *tmp_Y_SM_ptr++;
	}
	else {
		int i;
		double *tmp_sum_SM_ptr = this->SM_ptr, *tmp_Y_SM_ptr = Y_mat->SM_ptr;
		for ( i = 0 ; i < SM_length ; i++ )
			*tmp_sum_SM_ptr++ += *tmp_Y_SM_ptr++;
	}
}

void CCovarianceMatrix::Scale(double scale)
{
	int i;
	double *tmp_SM_ptr = this->SM_ptr;
	for ( i = 0 ; i < SM_length ; i++ )
		*tmp_SM_ptr++ *= scale;
}

void CCovarianceMatrix::OrthogonalVector(CCovarianceMatrix *y_tangent_vec, CCovarianceMatrix *orth_vec, CCovarianceMatrix *this_mat_inv_sqrt)
{
	CCovarianceMatrix *sqrt_inv_X = this_mat_inv_sqrt;

	if ( !sqrt_inv_X ) {
		sqrt_inv_X = new CCovarianceMatrix(M_size);
		this->Inverse(sqrt_inv_X);
		sqrt_inv_X->Sqrt();
	}

	if ( orth_vec )
		sqrt_inv_X->Mul2(y_tangent_vec, orth_vec);
	else
		sqrt_inv_X->Mul2(y_tangent_vec, this);

	if ( !this_mat_inv_sqrt )
		delete sqrt_inv_X;
}

void CCovarianceMatrix::OrthogonalVectorLog(CCovarianceMatrix *Y_mat, CCovarianceMatrix *vec_log, CCovarianceMatrix *this_mat_sqrt, CCovarianceMatrix *this_mat_inv_sqrt)
{
	CCovarianceMatrix *inv_sqrt_X = this_mat_inv_sqrt;

	if ( !inv_sqrt_X ) {
		inv_sqrt_X = new CCovarianceMatrix(M_size);
		this->Sqrt(inv_sqrt_X);
		inv_sqrt_X->Inverse();
	}

	if ( vec_log ) {
		inv_sqrt_X->Mul2(Y_mat, vec_log);
		vec_log->Log();
	}
	else {
		inv_sqrt_X->Mul2(Y_mat, this);
		this->Log();
	}

	if ( !this_mat_inv_sqrt ) {
		delete inv_sqrt_X;
	}

	/*
	CCovarianceMatrix *log_Y = new CCovarianceMatrix(M_size);
	this->Log_Y(Y_mat, log_Y, this_mat_sqrt, this_mat_inv_sqrt);
	this->OrthogonalVector(log_Y, vec_log, this_mat_inv_sqrt);

	delete log_Y;
	*/
}

void CCovarianceMatrix::SetZero()
{
	memset(SM_ptr, 0, sizeof(double)*SM_length);
	if ( AVG_ptr )
        memset(AVG_ptr, 0, sizeof(double)*M_size);
}

void CCovarianceMatrix::SetIdentityMatrix(double identity_value)
{
	memset(SM_ptr, 0, sizeof(double)*SM_length);

	double *tmp_SM_ptr = SM_ptr;
	for ( int dy = M_size ; dy > 0 ; dy-- ) {
		*tmp_SM_ptr = identity_value;
		tmp_SM_ptr += dy;
	}
	if ( AVG_ptr )
        memset(AVG_ptr, 0, sizeof(double)*M_size);
}

void CCovarianceMatrix::Export(const char* log_msg_fn, bool upper)
{
	if ( log_msg_fn ) {
		ofstream fout(log_msg_fn,ios::app);
		if (fout.fail()) {
			printf("Error opening log output file %s.\n", log_msg_fn);
			fout.close();
			exit(0);
		}

		if ( upper ) {
			double *tmp_SM_ptr = this->SM_ptr;
			for ( int i = 0 ; i < this->SM_length ; i++ ) {
				fout << *tmp_SM_ptr++ << '\t';
			}
			fout << '\n';
		}
		else {
			double *M_ptr = this->GetData(0);
			double *tmp_M_ptr = M_ptr;
			//fout << '\n';
			for ( int i = 0 ; i < M_size ; i++ ) {
				for ( int j = 0 ; j < M_size ; j++ )
					fout << *tmp_M_ptr++ << '\t';
				fout << '\n';
			}
			delete [] M_ptr;
		}

		fout.close();
	}
	else {
		if ( upper ) {
			double *tmp_SM_ptr = this->SM_ptr;
			printf("\n");
			for ( int i = 0 ; i < this->SM_length ; i++ ) {
				printf("%.6f\t", *tmp_SM_ptr++);
			}
			printf("\n");
		}
		else {
			double *M_ptr = this->GetData(0);
			double *tmp_M_ptr = M_ptr;
			printf("\n");
			for ( int i = 0 ; i < M_size ; i++ ) {
				for ( int j = 0 ; j < M_size ; j++ )
					printf("%.6f\t", *tmp_M_ptr++);
				printf("\n");
			}
			delete [] M_ptr;
		}
	}
}

void CCovarianceMatrix::Import(const char* file_name)
{
	// Input the data file of the camera projection matrix
	// and image distortion coefficients
	//--------------------------------------------------------------------
	ifstream fin;
	fin.open(file_name, ios::in);

	if (fin.fail()) {
		printf("Error opening matrix file %s.\n", file_name);
		fin.close();
		exit(0);
	}
	else {
		double *M_ptr = new double[M_size2];
		double *tmp_M_ptr = M_ptr;
		for ( int a = 0 ; a < M_size2 ; a++ )
			fin >> *tmp_M_ptr++;

		fin.close();		// close the file handler

		this->SetData(M_ptr, M_size2);
	}
}

void CCovarianceMatrix::GetSubMatrix(int *rows, int rows_num, CCovarianceMatrix *sub_cov_mat, int *row_start_pos)
{
	if ( sub_cov_mat->M_size != rows_num ) {
		printf("Must be the same dimension!\n");
		exit(1);
	}

	double *sub_SM_ptr = sub_cov_mat->SM_ptr;

	double *row_SM_ptr;
	int row_num_1 = rows_num-1;
	int i, j, row_pos;

	int M2_size_1 = 0;
	if ( !row_start_pos )
		M2_size_1 = this->M_size*2 + 1;

	//int *offset_rows = new int[row_num_1];
	int offset_rows[20];
	for ( i = 0 ; i < row_num_1 ; i++ )
		offset_rows[i] = rows[i+1]-rows[i];

	for ( i = 0 ; i < rows_num ; i++ ) {
		if ( rows[i] >= this->M_size || rows[i] < 0 || ( i>0 && rows[i] <= rows[i-1]) ) {
			printf("The rows of sub-cov-matrix are out of range or not arranged in the increasing order!\n");
			exit(1);
		}

		if ( !row_start_pos )
			row_pos = ((M2_size_1-rows[i])*rows[i])/2;
		else
			row_pos = row_start_pos[rows[i]];


        if ( sub_cov_mat->AVG_ptr && this->AVG_ptr )
            sub_cov_mat->AVG_ptr[i] = this->AVG_ptr[rows[i]];

		row_SM_ptr = this->SM_ptr + row_pos;

		*sub_SM_ptr++ = *row_SM_ptr;
		for ( j = i ; j < row_num_1 ; j++ ) {
			row_SM_ptr += offset_rows[j];
			*sub_SM_ptr++ = *row_SM_ptr;
		}
	}

	//delete [] offset_rows;
}

void CCovarianceMatrix::Resize(int mat_size, bool used_avg)
{
	if ( M_size == mat_size )
		return;

	Delete();

	M_size = mat_size;
	SM_length = (mat_size*(mat_size+1))/2;
	M_size2 = M_size*M_size;

	SM_ptr = new double[SM_length];
	AVG_ptr = NULL;
	if ( used_avg )
        AVG_ptr = new double[M_size];

	//memset(SM_ptr, 0, sizeof(double)*SM_length);
}

int* CCovarianceMatrix::GetRowStartPos()
{
	int* row_start_pos = new int[M_size];
	int M2_size_1 = M_size*2 + 1;

	for ( int a = 0 ; a < M_size ; a++ )
		row_start_pos[a] = ((M2_size_1- a)*a)/2;

	return row_start_pos;
}

void CCovarianceMatrix::Mul(CCovarianceMatrix *mat1, CCovarianceMatrix *mat2, CCovarianceMatrix *mul_mat)
{

}

void CCovarianceMatrix::GeneratedRand(CvRNG rng_state, double range)
{
	double *tmp_SM_ptr = SM_ptr;
	for ( int i = 0 ; i < SM_length ; i++ ) {
		*tmp_SM_ptr++ += (double)(cvRandInt(&rng_state) % 100000)/100000.0 * range;
	}
}

void CCovarianceMatrix::Print(const char* file_name, bool stdprint)
{
    char info[4096];
    double *data = GetData(0);
    sprintf(info, "\nCovariance Matrix:\n");
    double *tmp_data = data;
    for ( int i = 0 ; i < M_size ; i++ ) {
        for ( int j = 0 ; j < M_size ; j++ )
            sprintf(info, "%s  %E", info, *tmp_data++);
        sprintf(info, "%s\n", info);
    }
    if ( AVG_ptr ) {
        sprintf(info, "%sMean:\n", info);
        for ( int i = 0 ; i < M_size ; i++ )
            sprintf(info, "%s  %E", info, AVG_ptr[i]);
        sprintf(info, "%s\n\n", info);
    }

    if ( stdprint || !file_name )
        printf(info);

    delete [] data;

	if ( !file_name ) {
		return;
	}

	ofstream fout(file_name,ios::app);
	if (fout.fail()) {
		printf("Error opening log output file %s.\n", file_name);
		fout.close();
		exit(0);
	}

	fout << info;
	fout.close();
}
