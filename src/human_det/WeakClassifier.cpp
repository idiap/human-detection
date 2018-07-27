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
#include "WeakClassifier.h"
#include "cv.h"

/* Constructor */
WeakClassifier::WeakClassifier(int feature_size_)
{
    feature_size = feature_size_;

    fitting_params      = new double[feature_size+1];
}

/* Destructor */
WeakClassifier::~WeakClassifier()
{
    delete[] fitting_params;
}

/* Test an example sizes are supposed to match with the classifier */
double WeakClassifier::Test(double* example)
{
    double score = 0;
    for(int i=0;i<feature_size;i++){
        score += fitting_params[i]*example[i];
    }
	score += fitting_params[feature_size];

    return score;
}

/* least square regression (weighted or not)
   given a set of examples
*/
void WeakClassifier::WeightedLeastSquareRegression( int n_examples,
                                            double* examples,
                                            double* responses,
                                            double* example_weights)
{
	CvMat* X = cvCreateMat(n_examples, feature_size+1, CV_64F);
	double *tmp_X_ptr = (double*)X->data.ptr;
	double *tmp_examples = examples;
	for ( int dy = 0 ; dy < n_examples ; dy++ ) {
		for ( int dx = 0 ; dx < feature_size ; dx++ )
			*tmp_X_ptr++ = *tmp_examples++;
		*tmp_X_ptr++ = 1.0;
	}

    CvMat y    = cvMat(n_examples,            1, CV_64F, (double*) responses);
    CvMat beta = cvMat(feature_size+1,        1, CV_64F, (double*) this->fitting_params);


    /* Standard least square */
    if(example_weights==NULL){
        cvSolve(X, &y, &beta, CV_SVD);
    }
    /* Weighted least square */
    else {
        /* Weighted least square regression
           beta = (X^tWX)^-1*X^t*W*y
        */
        CvMat* XtW = cvCreateMat(feature_size+1, n_examples, CV_64F);
        cvTranspose(X, XtW);

		/* multiplication by example weights */
		double *tmp_XtW_ptr = (double*)XtW->data.ptr;
		for ( int dy = 0 ; dy < feature_size+1 ; dy++ ) {
			double *tmp_weights = example_weights;
			for ( int dx = 0 ; dx < n_examples ; dx++ )
				*tmp_XtW_ptr++ *= *tmp_weights++;
		}

        CvMat* XtWX = cvCreateMat(feature_size+1, feature_size+1, CV_64F);
        cvMatMul(XtW, X, XtWX);

        CvMat* XtWX_inv = cvCreateMat(feature_size+1, feature_size+1, CV_64F);
        //cvInvert(XtWX, XtWX_inv, CV_SVD);
        cvInv(XtWX, XtWX_inv);

        CvMat* tmp_mat = cvCreateMat(feature_size+1, n_examples, CV_64F);
        cvMatMul(XtWX_inv, XtW, tmp_mat);

        cvMatMul(tmp_mat, &y , &beta);

        cvReleaseMat(&XtW);
        cvReleaseMat(&XtWX_inv);
        cvReleaseMat(&XtWX);
        cvReleaseMat(&tmp_mat);
    }

	cvReleaseMat(&X);
}
