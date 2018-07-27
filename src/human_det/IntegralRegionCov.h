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
// IntegralRegionCov.h: interface for the CIntegralRegionCov class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_INTEGRAL_REGION_COV_H_)
#define _INTEGRAL_REGION_COV_H_

#include "cv.h"
#include "OpencvDataConversion.h"
#include "CovarianceMatrix.h"

#include <ctime>						// clock
#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <fstream>						// file I/O
#include <cmath>			// math includes
#include <iostream>			// I/O streams

using namespace std;					// make std:: accessible

#define MIN_REGION_COV_SIZE	10

#define MAX_SUB_FEATURE_DIM	15
#define MAX_FEATURE_DIM	20


#ifndef CLOCKS_PER_SEC						  /* define clocks-per-second if needed */
#define CLOCKS_PER_SEC 1000000
#endif


/************************************************************************/
/* integral point structure for computing region covariance matrices	*/
/* by using integral images in the first-order and the second-order	 */
/************************************************************************/
typedef struct _IntegralPoint
{
	double*	p;	/* the d-dimensional vector, P(x',y',i) = sum_{x<x',y<y'} F(x,y,i) */
	double*	Q;	/* the d x d dimensional matrix, Q(x',y',i,j) = sum_{x<x',y<y'} F(x,y,i) F(x,y,j),
			   Q is symmetric matrix, here we only keep the elements in the upper triangle,
			   so we only allocate d*(d+1)/2 elements for Q */
	int	mask;	/* the one-dimensional mask value */
}
IntegralPointStruct;

/************************************************************************/
/* feature image structure											  */
/************************************************************************/
typedef struct _FeatureImage
{
	int width;	/* the image width */
	int height;	/* the image height */
	float *data;	/* the float data pointer */
}
FeatureImageStruct;

class CIntegralRegionCov
{
public:
	void GetSubRegionCovNormDiagonalElems(CvRect rect_roi,
            int *sub_feature_idxes, int sub_feature_num, double *norm_diag_elems,
            double *avg_elems=NULL, double min_identity_value=MIN_IDENTITY_VALUE);
	void GetRegionCovNormDiagonalElems(CvRect rect_roi, double *norm_diag_elems,
            double *avg_elems=NULL, double min_identity_value=MIN_IDENTITY_VALUE);
	void GetSubRegionCovFeature(CvRect rect_roi, CCovarianceMatrix* sub_reg_cov,
            int *sub_feature_idxes, double *full_avg_elems=NULL,
            double min_identity_value=MIN_IDENTITY_VALUE, bool compute_avg=true);
	/* for testing */
	void ExportLogMessage(char *msg);

	/* get the region covariance matrix from integral images */
	void GetRegionCovFeature(CvRect rect_roi, CCovarianceMatrix* reg_cov,
            double min_identity_value=MIN_IDENTITY_VALUE, bool compute_avg=true);

	/* compute the integral images */
	void ComputeIntegralImages();

	/* release all the memories */
	void CleanData();

	/* set new multi-dimensional feature images and the mask image */
	void SetNewData(IplImage** new_feature_imgs, IplImage *mask_img = NULL);

	/* save data to a file */
	void Save(char *file_name, CvRect *roi = NULL);

	/* load data from a file */
	void Load(char *file_name);

	/* initialization */
	void Init(int cov_dim, CvSize img_size);

	int		m_nCovDim;			/* the dimension d of the region covariance matrix */

	CvSize		m_szFeatureImage;		/* the size of current input feature images */

	CvSize		m_szMaxFeatureImage;	/* the maximal size of input feature images */

	IplImage**	m_ppFeatureImages;		/* the d-dimensional feature images */

	IplImage*	m_pMaskImage;			/* the mask image (for warped image) */

	IntegralPointStruct*	m_pIntegralImages;	/* the integral images */

	CIntegralRegionCov(int cov_dim, CvSize img_size);
	CIntegralRegionCov();

	virtual ~CIntegralRegionCov();

private:
	int		m_nCovQDim;			/* = d * (d+1) / 2 */
	int*		m_pRowStartPos;
};

#endif // !defined(_INTEGRAL_REGION_COV_H_)

