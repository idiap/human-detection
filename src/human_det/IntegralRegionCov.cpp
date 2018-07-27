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
// IntegralRegionCov.cpp: implementation of the CIntegralRegionCov class.
//
//////////////////////////////////////////////////////////////////////

#include "IntegralRegionCov.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CIntegralRegionCov::CIntegralRegionCov(int cov_dim, CvSize img_size)
{
	m_pMaskImage = NULL;
	m_pIntegralImages = NULL;
	m_pRowStartPos = NULL;

	m_szFeatureImage = cvSize(0,0);

	Init(cov_dim, img_size);
}

CIntegralRegionCov::CIntegralRegionCov()
{
	m_pMaskImage = NULL;
	m_pIntegralImages = NULL;
	m_pRowStartPos = NULL;

	m_szFeatureImage = cvSize(0,0);
}

CIntegralRegionCov::~CIntegralRegionCov()
{
	CleanData();
}

void CIntegralRegionCov::SetNewData(IplImage **new_feature_imgs, IplImage *mask_img)
{
	m_ppFeatureImages = new_feature_imgs;
	if ( mask_img )
		cvCopy(mask_img, m_pMaskImage);
}

void CIntegralRegionCov::Init(int cov_dim, CvSize img_size)
{
	m_nCovDim = cov_dim;
	m_nCovQDim = (cov_dim * (cov_dim+1))/2;

	/* get the size of feature image */
	m_szFeatureImage = img_size;

	/* allocate memory for mask image and set it to ones */
	m_pMaskImage = cvCreateImage(m_szFeatureImage, IPL_DEPTH_8U, 1);
	cvSet(m_pMaskImage, cvScalar(1));

	/* allocate memories for the integral images */
	int length = (m_szFeatureImage.width+1) * (m_szFeatureImage.height+1);
	m_pIntegralImages = new IntegralPointStruct[length];
	int a;
	for ( a = 0 ; a < length ; a++ ) {
		m_pIntegralImages[a].p = new double[m_nCovDim];
		m_pIntegralImages[a].Q = new double[m_nCovQDim];
	}

	m_pRowStartPos = new int[m_nCovDim];
	int M2_size_1 = m_nCovDim*2 + 1;

	for ( a = 0 ; a < m_nCovDim ; a++ )
		m_pRowStartPos[a] = ((M2_size_1- a)*a)/2;

	IntegralPointStruct *cur_integral_pt;

	int p_vec_length = sizeof(double)*m_nCovDim;
	int Q_mat_length = sizeof(double)*m_nCovQDim;

	/* initialize the first-row integral points */
	cur_integral_pt = m_pIntegralImages;
	int x, y;
	for ( x = 0 ; x <= m_szFeatureImage.width ; x++ ) {
		(*cur_integral_pt).mask = 0;
		memset((*cur_integral_pt).p, 0, p_vec_length);
		memset((*cur_integral_pt).Q, 0, Q_mat_length);
		cur_integral_pt++;
	}

	/* initialize the first-column integral points */
	cur_integral_pt = m_pIntegralImages;
	int width_1 = m_szFeatureImage.width+1;
	for ( y = 0 ; y <= m_szFeatureImage.height ; y++ ) {
		(*cur_integral_pt).mask = 0;
		memset((*cur_integral_pt).p, 0, p_vec_length);
		memset((*cur_integral_pt).Q, 0, Q_mat_length);
		cur_integral_pt += width_1;
	}
}

void CIntegralRegionCov::CleanData()
{
	int a;

	/* clean all allocated memories */
	if ( m_pIntegralImages ) {
		int length = (m_szFeatureImage.width+1) * (m_szFeatureImage.height+1);
		for ( a = 0 ; a < length ; a++ ) {
			delete [] m_pIntegralImages[a].p;
			delete [] m_pIntegralImages[a].Q;
		}
		delete [] m_pIntegralImages;
		m_pIntegralImages = NULL;
	}
	if ( m_pMaskImage ) {
		cvReleaseImage(&m_pMaskImage);
		m_pMaskImage = NULL;
	}

	if ( m_pRowStartPos )
		delete [] m_pRowStartPos;
	m_pRowStartPos = NULL;

	m_szFeatureImage = cvSize(0,0);
}

void CIntegralRegionCov::ComputeIntegralImages()
{
	int x, y, d;

	/************************************************************************/
	/* covert Opencv IplImage data to one-dimensional float data			*/
	/************************************************************************/

	COpencvDataConversion<float,float> ODC1;
	float** feature_data = new float*[m_nCovDim];
	for ( d = 0 ; d < m_nCovDim ; d++ ) {
		//cvNormalize(m_ppFeatureImages[d], m_ppFeatureImages[d], 100, 0, CV_MINMAX, m_pMaskImage);
		feature_data[d] = ODC1.GetImageData(m_ppFeatureImages[d]);
	}
	COpencvDataConversion<uchar,uchar> ODC2;
	uchar*	mask_data = ODC2.GetImageData(m_pMaskImage);

	/************************************************************************/
	/* debug																*/
	/************************************************************************/
	/*
	char msg[1024];
	for ( d = 0 ; d < m_nCovDim ; d++ ) {
		sprintf(msg, "\n\nDimension : %d\n\n", d);
		ExportLogMessage(msg);
		yx = 0;
		for ( y = 0 ; y < m_szFeatureImage.height ; y++ ) {
			for ( x = 0 ; x < m_szFeatureImage.width ; x++ ) {
				sprintf(msg, "%f\t", feature_data[d][yx++]);
				ExportLogMessage(msg);
			}
			sprintf(msg, "\n");
			ExportLogMessage(msg);
		}
	}
	*/

	/************************************************************************/
	/*  compute one-dimension integral mask,				*/
	/*	 d-dimensional first-order integral images,			*/
	/*	 and d*d-dimensional second-order integral images			*/
	/************************************************************************/

	IntegralPointStruct *top_integral_pt, *cur_integral_pt;

	float** _feature_data = new float*[m_nCovDim];
	for ( d = 0 ; d < m_nCovDim ; d++ )
		_feature_data[d] = feature_data[d];

	int p_vec_length = sizeof(double)*m_nCovDim;
	int Q_mat_length = sizeof(double)*m_nCovQDim;

	/* compute the p-vectors for the integral points in the rest rows */
	cur_integral_pt = m_pIntegralImages+m_szFeatureImage.width+2;
	top_integral_pt = m_pIntegralImages+1;

	int left_mask;
	double *cur_p_vec, *left_p_vec, *top_p_vec;
	double *cur_Q_mat, *left_Q_mat, *top_Q_mat;
	double *tmp_left_p_vec, *tmp_left_Q_mat;
	uchar* cur_mask_data = mask_data;

	int dx, dy;
	double *cur_f_vec_dy, *cur_f_vec_dx, *cur_f_vec, *tmp_f_vec;

	cur_f_vec = new double[m_nCovDim];

	left_p_vec = new double[m_nCovDim];
	left_Q_mat = new double[m_nCovQDim];

	/* variables related to compute elapsed time for some operation */
	/*
	clock_t stopTime;
	clock_t startTime;
	double elapsedTime;

	startTime = clock();
	*/

	for ( y = 1 ; y <= m_szFeatureImage.height ; y++ ) {
		memset(left_p_vec, 0, p_vec_length);
		memset(left_Q_mat, 0, Q_mat_length);
		left_mask = 0;
		for ( x = 1 ; x <= m_szFeatureImage.width ; x++ ) {
			left_mask += (int)(*cur_mask_data);
			(*cur_integral_pt).mask = left_mask + (*top_integral_pt).mask;
			top_p_vec = (*top_integral_pt).p;
			cur_p_vec = (*cur_integral_pt).p;
			cur_Q_mat = (*cur_integral_pt).Q;
			top_Q_mat = (*top_integral_pt).Q;

			tmp_f_vec = cur_f_vec;
			for ( d = 0 ; d < m_nCovDim ; d++ )
				*tmp_f_vec++ = (double)(*_feature_data[d]++);

			if ( *cur_mask_data++ ) { /* mask = 1 */
				tmp_left_p_vec = left_p_vec;
				tmp_f_vec = cur_f_vec;
				for ( d = 0 ; d < m_nCovDim ; d++ ) {
					(*tmp_left_p_vec) += (*tmp_f_vec++);
					(*cur_p_vec++) = (*top_p_vec++) + (*tmp_left_p_vec++);
				}

				tmp_left_Q_mat = left_Q_mat;
				cur_f_vec_dy = cur_f_vec;
				for ( dy = 0 ; dy < m_nCovDim ; dy++ ) {
					cur_f_vec_dx = cur_f_vec + dy;
					for ( dx = dy ; dx < m_nCovDim ; dx++ ) {
						(*tmp_left_Q_mat) += (*cur_f_vec_dy)*(*cur_f_vec_dx++);
						(*cur_Q_mat++) = (*top_Q_mat++) + (*tmp_left_Q_mat++);
					}
					cur_f_vec_dy++;
				}
			}
			else {	/* mask = 0 */
				tmp_left_p_vec = left_p_vec;
				for ( d = 0 ; d < m_nCovDim ; d++ ) {
					(*cur_p_vec++) = (*top_p_vec++) + (*tmp_left_p_vec++);
				}

				tmp_left_Q_mat = left_Q_mat;
				for ( dy = 0 ; dy < m_nCovDim ; dy++ ) {
					for ( dx = dy ; dx < m_nCovDim ; dx++ ) {
						(*cur_Q_mat++) = (*top_Q_mat++) + (*tmp_left_Q_mat++);
					}
				}
			}

			top_integral_pt++;
			cur_integral_pt++;
		}

		top_integral_pt++;
		cur_integral_pt++;
	}

	/*
	stopTime = clock();

	elapsedTime = (stopTime - startTime) / (double)(CLOCKS_PER_SEC);
	printf("\nElapsed time of computing integral images with size %dx%d for %d dimension region covariances (%.6fs)\n\n",
				m_szFeatureImage.width, m_szFeatureImage.height, m_nCovDim, elapsedTime);
	*/

	/************************************************************************/
	/*  release memories													*/
	/************************************************************************/
	for ( d = 0 ; d < m_nCovDim ; d++ )
		delete [] feature_data[d];
	delete [] feature_data;
	delete [] _feature_data;
	delete [] mask_data;
	delete [] left_p_vec;
	delete [] left_Q_mat;
	delete [] cur_f_vec;
}

void CIntegralRegionCov::GetRegionCovFeature(CvRect rect_roi, CCovarianceMatrix* reg_cov,
        double min_identity_value, bool compute_avg)
{
	if ( reg_cov->M_size != m_nCovDim ) {
		printf("\nMust be the same dimension of region covariance feature!\n");
		printf("reg_cov->M_size = %d, m_nCovDim = %d\n", reg_cov->M_size, m_nCovDim);
		return;
		//exit(1);
	}

	int x1, y1, x2, y2;
	x1 = rect_roi.x;
	y1 = rect_roi.y;
	x2 = x1 + rect_roi.width;
	y2 = y1 + rect_roi.height;

	if ( x2 > m_szFeatureImage.width || y2 > m_szFeatureImage.height ) {
		printf("\nROI must be inside the feature images!\n");
		exit(1);
	}

	int width_1 = m_szFeatureImage.width+1;

	int y1_width_1 = y1*width_1;
	int y2_width_1 = y2*width_1;

	int left_top_pos = y1_width_1 + x1;
	int right_top_pos = y1_width_1 + x2;
	int left_bot_pos = y2_width_1 + x1;
	int right_bot_pos = y2_width_1 + x2;

	IntegralPointStruct *left_top_ipt = m_pIntegralImages + left_top_pos;
	IntegralPointStruct *right_top_ipt = m_pIntegralImages + right_top_pos;
	IntegralPointStruct *left_bot_ipt = m_pIntegralImages + left_bot_pos;
	IntegralPointStruct *right_bot_ipt = m_pIntegralImages + right_bot_pos;

	int n = left_top_ipt->mask + right_bot_ipt->mask - right_top_ipt->mask - left_bot_ipt->mask;

	if ( n < 1 ) {
		printf("\nThe sum of masks in the region of interested = 0!\n");
		return;
		exit(1);
	}

	double c1 = (double)(n-1);
	double c2 = c1 * (double)n;

	double p_vec[MAX_FEATURE_DIM];
	int d;
	double *p00=left_top_ipt->p, *p10=left_bot_ipt->p;
	double *p01=right_top_ipt->p, *p11=right_bot_ipt->p;
	double *p_vec_dy, *p_vec_dx;

	p_vec_dx = p_vec;

	for ( d = 0 ; d < m_nCovDim ; d++ )
		*p_vec_dx++ = (*p00++) + (*p11++) - (*p01++) - (*p10++);

	int dx, dy;
	p_vec_dy = p_vec;

	double *Q00=left_top_ipt->Q, *Q10=left_bot_ipt->Q;
	double *Q01=right_top_ipt->Q, *Q11=right_bot_ipt->Q;
	double *cov_Q = reg_cov->SM_ptr;

	double *identity_cov_Q = reg_cov->SM_ptr;
	double p_vec_dy_c2;
	int i_shift = m_nCovDim;
	double *avg = reg_cov->AVG_ptr;
	for ( dy = 0 ; dy < m_nCovDim ; dy++ ) {
		//p_vec_dx = p_vec + dy;
		p_vec_dx = p_vec_dy;
		p_vec_dy_c2 = (*p_vec_dy)/c2;
		if ( compute_avg )
            *avg++ = (*p_vec_dy)/(double)n + min_identity_value;
		for ( dx = dy ; dx < m_nCovDim ; dx++ ) {
			//*cov_Q++ = ((*Q00++)+(*Q11++)-(*Q01++)-(*Q10++))/c1 - (*p_vec_dy)*(*p_vec_dx++)/c2;
			*cov_Q++ = ((*Q00++)+(*Q11++)-(*Q01++)-(*Q10++))/c1 - p_vec_dy_c2*(*p_vec_dx++);
		}
		*identity_cov_Q += min_identity_value;
		identity_cov_Q += i_shift--;
		p_vec_dy++;
	}

	/*
	cov_Q = reg_cov->SM_ptr;
	for ( dy = m_nCovDim ; dy > 0 ; dy-- ) {
		*cov_Q = MAX(*cov_Q, 0.00001f);
		cov_Q += dy;
	}
	*/

	/*
	CCovarianceMatrix* small_identity_mat = new CCovarianceMatrix(reg_cov->M_size);
	small_identity_mat->SetIdentityMatrix(min_identity_value);
	reg_cov->Add(small_identity_mat);

	delete small_identity_mat;
	*/

	//delete [] p_vec;
}

void CIntegralRegionCov::ExportLogMessage(char *msg)
{
	const char *log_fn = "log_message.txt";
	ofstream fout(log_fn,ios::app);
	if (fout.fail()) {
		printf("Error opening log output file %s.\n", log_fn);
		fout.close();
		exit(0);
	}

	fout << msg;
	fout.close();
}

void CIntegralRegionCov::GetSubRegionCovFeature(CvRect rect_roi, CCovarianceMatrix *sub_reg_cov,
        int *sub_feature_idxes, double *full_avg_elems,
        double min_identity_value, bool compute_avg)
{
	if ( sub_reg_cov->M_size >= m_nCovDim ) {
		printf("sub_reg_cov->M_size >= %d, m_nCovDim = %d\n", sub_reg_cov->M_size, m_nCovDim);
		return;
		//exit(1);
	}

	int x1, y1, x2, y2;
	x1 = rect_roi.x;
	y1 = rect_roi.y;
	x2 = x1 + rect_roi.width;
	y2 = y1 + rect_roi.height;

	if ( x2 > m_szFeatureImage.width || y2 > m_szFeatureImage.height ) {
		printf("\nROI must be inside the feature images!\n");
		exit(1);
	}

	int width_1 = m_szFeatureImage.width+1;

	int y1_width_1 = y1*width_1;
	int y2_width_1 = y2*width_1;

	int left_top_pos = y1_width_1 + x1;
	int right_top_pos = y1_width_1 + x2;
	int left_bot_pos = y2_width_1 + x1;
	int right_bot_pos = y2_width_1 + x2;

	IntegralPointStruct *left_top_ipt = m_pIntegralImages + left_top_pos;
	IntegralPointStruct *right_top_ipt = m_pIntegralImages + right_top_pos;
	IntegralPointStruct *left_bot_ipt = m_pIntegralImages + left_bot_pos;
	IntegralPointStruct *right_bot_ipt = m_pIntegralImages + right_bot_pos;

	double *p00=left_top_ipt->p, *p10=left_bot_ipt->p;
	double *p01=right_top_ipt->p, *p11=right_bot_ipt->p;

	int n = left_top_ipt->mask + right_bot_ipt->mask - right_top_ipt->mask - left_bot_ipt->mask;

	if ( n < 1 ) {
		printf("\nThe sum of masks in the region of interested = 0!\n");
		exit(1);
	}

	double c1 = (double)(n-1);
	double c2 = c1 * (double)n;

	int d;
	int row_pos;
	double *tmp_full_avg_elems;

	if ( full_avg_elems ) {
		tmp_full_avg_elems = full_avg_elems;
		for ( d = 0 ; d < m_nCovDim ; d++ )
			*tmp_full_avg_elems++ = (*p00++) + (*p11++) - (*p01++) - (*p10++);
	}

	if ( sub_reg_cov->M_size == 1 ) {
		int row_idx = sub_feature_idxes[0];
		double p = p00[row_idx] + p11[row_idx] - p01[row_idx] - p10[row_idx];

        if ( compute_avg )
            sub_reg_cov->AVG_ptr[0] = p/(double)n + min_identity_value;

		if ( full_avg_elems ) {
			tmp_full_avg_elems = full_avg_elems;
			for ( d = 0 ; d < m_nCovDim ; d++ ) {
				*tmp_full_avg_elems /= (double)n;
				*tmp_full_avg_elems++ += min_identity_value;
			}
		}

		row_pos = m_pRowStartPos[row_idx];
		sub_reg_cov->SM_ptr[0] = (left_top_ipt->Q[row_pos] + right_bot_ipt->Q[row_pos]
            - left_bot_ipt->Q[row_pos] - right_top_ipt->Q[row_pos])/c1
            - p*p/c2 + min_identity_value;

		return;
	}

	//double *p_vec = new double[sub_reg_cov->M_size];

	double p_vec[MAX_SUB_FEATURE_DIM];

	double *p_vec_dy, *p_vec_dx;

	p_vec_dx = p_vec;

	int *rows_idxes = sub_feature_idxes;

	if ( full_avg_elems ) {
		for ( d = 0 ; d < sub_reg_cov->M_size ; d++ )
			*p_vec_dx++ = full_avg_elems[*rows_idxes++];

		tmp_full_avg_elems = full_avg_elems;
		for ( d = 0 ; d < m_nCovDim ; d++ ){
			*tmp_full_avg_elems /= (double)n;
			*tmp_full_avg_elems++ += min_identity_value;
		}
	}
	else {
		for ( d = 0 ; d < sub_reg_cov->M_size ; d++ ) {
			*p_vec_dx++ = p00[*rows_idxes] + p11[*rows_idxes] - p01[*rows_idxes] - p10[*rows_idxes];
			rows_idxes++;
		}
	}

	int dx, dy;
	p_vec_dy = p_vec;

	double *Q00, *Q10, *Q01, *Q11;

	double *sub_cov_Q = sub_reg_cov->SM_ptr;

	int row_num_1 = sub_reg_cov->M_size-1;

	double *identity_sub_cov_Q = sub_reg_cov->SM_ptr;
	double p_vec_dy_c2;

	int offset_rows[MAX_SUB_FEATURE_DIM];
	int offset_row_dx;
	rows_idxes = offset_rows;
	for ( d = 0 ; d < row_num_1 ; d++ )
		offset_rows[d] = sub_feature_idxes[d+1]-sub_feature_idxes[d];

	rows_idxes = sub_feature_idxes;

	int i_shift = sub_reg_cov->M_size;

	double *avg = sub_reg_cov->AVG_ptr;

	for ( dy = 0 ; dy < sub_reg_cov->M_size ; dy++ ) {
		//p_vec_dx = p_vec + dy;
		p_vec_dx = p_vec_dy;
		p_vec_dy_c2 = (*p_vec_dy)/c2;

        if ( compute_avg )
            *avg++ = (*p_vec_dy)/(double)n + min_identity_value;

		row_pos = m_pRowStartPos[*rows_idxes++];

		Q00=left_top_ipt->Q + row_pos;
		Q10=left_bot_ipt->Q + row_pos;
		Q01=right_top_ipt->Q + row_pos;
		Q11=right_bot_ipt->Q + row_pos;

		*sub_cov_Q++ = ((*Q00)+(*Q11)-(*Q01)-(*Q10))/c1 - p_vec_dy_c2*(*p_vec_dx++);

		for ( dx = dy ; dx < row_num_1 ; dx++ ) {
			//*cov_Q++ = ((*Q00++)+(*Q11++)-(*Q01++)-(*Q10++))/c1 - (*p_vec_dy)*(*p_vec_dx++)/c2;

			offset_row_dx = offset_rows[dx];

			Q00 += offset_row_dx;
			Q01 += offset_row_dx;
			Q10 += offset_row_dx;
			Q11 += offset_row_dx;

			*sub_cov_Q++ = ((*Q00)+(*Q11)-(*Q01)-(*Q10))/c1 - p_vec_dy_c2*(*p_vec_dx++);
		}

		*identity_sub_cov_Q += min_identity_value;
		identity_sub_cov_Q += i_shift--;
		p_vec_dy++;
	}

	/*
	CCovarianceMatrix* small_identity_mat = new CCovarianceMatrix(reg_cov->M_size);
	small_identity_mat->SetIdentityMatrix(min_identity_value);
	reg_cov->Add(small_identity_mat);

	delete small_identity_mat;
	*/

	//delete [] p_vec;
	//delete [] offset_rows;
}

void CIntegralRegionCov::GetRegionCovNormDiagonalElems(CvRect rect_roi, double *norm_diag_elems, double *avg_elems, double min_identity_value)
{
	int x1, y1, x2, y2;
	x1 = rect_roi.x;
	y1 = rect_roi.y;
	x2 = x1 + rect_roi.width;
	y2 = y1 + rect_roi.height;

	if ( x2 > m_szFeatureImage.width || y2 > m_szFeatureImage.height ) {
		printf("\nROI must be inside the feature images!\n");
		exit(1);
	}

	int width_1 = m_szFeatureImage.width+1;

	int y1_width_1 = y1*width_1;
	int y2_width_1 = y2*width_1;

	int left_top_pos = y1_width_1 + x1;
	int right_top_pos = y1_width_1 + x2;
	int left_bot_pos = y2_width_1 + x1;
	int right_bot_pos = y2_width_1 + x2;

	IntegralPointStruct *left_top_ipt = m_pIntegralImages + left_top_pos;
	IntegralPointStruct *right_top_ipt = m_pIntegralImages + right_top_pos;
	IntegralPointStruct *left_bot_ipt = m_pIntegralImages + left_bot_pos;
	IntegralPointStruct *right_bot_ipt = m_pIntegralImages + right_bot_pos;

	int n = left_top_ipt->mask + right_bot_ipt->mask - right_top_ipt->mask - left_bot_ipt->mask;

	if ( n < 1 ) {
		printf("\nThe sum of masks in the region of interested = 0!\n");
		exit(1);
	}

	double c1 = (double)(n-1);
	double c2 = c1 * (double)n;

	double p_vec[MAX_FEATURE_DIM];
	int d;
	double *p00=left_top_ipt->p, *p10=left_bot_ipt->p;
	double *p01=right_top_ipt->p, *p11=right_bot_ipt->p;
	double *p_vec_dy;

	p_vec_dy = p_vec;

	for ( d = 0 ; d < m_nCovDim ; d++ )
		*p_vec_dy++ = (*p00++) + (*p11++) - (*p01++) - (*p10++);

	p_vec_dy = p_vec;

	double *Q00=left_top_ipt->Q, *Q10=left_bot_ipt->Q;
	double *Q01=right_top_ipt->Q, *Q11=right_bot_ipt->Q;

	double *tmp_elems = norm_diag_elems;

	if ( avg_elems ) {
		double *tmp_avg_elems = avg_elems;
		for ( d = m_nCovDim ; d > 0 ; d-- ) {

			*tmp_avg_elems++ = (*p_vec_dy)/(double)n + min_identity_value;

			*tmp_elems = ((*Q00)+(*Q11)-(*Q01)-(*Q10))/c1 - (*p_vec_dy)*(*p_vec_dy)/c2 + min_identity_value;

			*tmp_elems = sqrt(*tmp_elems);
			tmp_elems++;

			Q00 += d;
			Q11 += d;
			Q01 += d;
			Q10 += d;

			p_vec_dy++;
		}
	}
	else {
		for ( d = m_nCovDim ; d > 0 ; d-- ) {

			*tmp_elems = ((*Q00)+(*Q11)-(*Q01)-(*Q10))/c1 - (*p_vec_dy)*(*p_vec_dy)/c2 + min_identity_value;

			*tmp_elems = sqrt(*tmp_elems);
			tmp_elems++;

			Q00 += d;
			Q11 += d;
			Q01 += d;
			Q10 += d;

			p_vec_dy++;
		}
	}

	//delete [] p_vec;
}

void CIntegralRegionCov::GetSubRegionCovNormDiagonalElems(CvRect rect_roi, int *sub_feature_idxes, int sub_feature_num, double *norm_diag_elems, double *avg_elems, double min_identity_value)
{
	int x1, y1, x2, y2;
	x1 = rect_roi.x;
	y1 = rect_roi.y;
	x2 = x1 + rect_roi.width;
	y2 = y1 + rect_roi.height;

	if ( x2 > m_szFeatureImage.width || y2 > m_szFeatureImage.height ) {
		printf("\nROI must be inside the feature images!\n");
		exit(1);
	}

	int width_1 = m_szFeatureImage.width+1;

	int y1_width_1 = y1*width_1;
	int y2_width_1 = y2*width_1;

	int left_top_pos = y1_width_1 + x1;
	int right_top_pos = y1_width_1 + x2;
	int left_bot_pos = y2_width_1 + x1;
	int right_bot_pos = y2_width_1 + x2;

	IntegralPointStruct *left_top_ipt = m_pIntegralImages + left_top_pos;
	IntegralPointStruct *right_top_ipt = m_pIntegralImages + right_top_pos;
	IntegralPointStruct *left_bot_ipt = m_pIntegralImages + left_bot_pos;
	IntegralPointStruct *right_bot_ipt = m_pIntegralImages + right_bot_pos;

	int n = left_top_ipt->mask + right_bot_ipt->mask - right_top_ipt->mask - left_bot_ipt->mask;

	if ( n < 1 ) {
		printf("\nThe sum of masks in the region of interested = 0!\n");
		exit(1);
	}

	double c1 = (double)(n-1);
	double c2 = c1 * (double)n;

	double p_vec[MAX_SUB_FEATURE_DIM];
	int d;
	double *p00=left_top_ipt->p, *p10=left_bot_ipt->p;
	double *p01=right_top_ipt->p, *p11=right_bot_ipt->p;
	double *p_vec_dy;

	p_vec_dy = p_vec;

	int *rows_idxes = sub_feature_idxes;

	for ( d = 0 ; d < sub_feature_num ; d++ ) {
		*p_vec_dy++ = p00[*rows_idxes] + p11[*rows_idxes] - p01[*rows_idxes] - p10[*rows_idxes];
		rows_idxes++;
	}

	p_vec_dy = p_vec;

	double *s_Q00=left_top_ipt->Q, *s_Q10=left_bot_ipt->Q;
	double *s_Q01=right_top_ipt->Q, *s_Q11=right_bot_ipt->Q;

	double *Q00, *Q10, *Q01, *Q11;

	rows_idxes = sub_feature_idxes;
	int row_pos;

	double *tmp_elems = norm_diag_elems;

	if ( avg_elems ) {
		double *tmp_avg_elems = avg_elems;

		for ( d = 0 ; d < sub_feature_num ; d++ ) {

			*tmp_avg_elems++ = (*p_vec_dy)/(double)n + min_identity_value;

			row_pos = m_pRowStartPos[*rows_idxes++];

			Q00 = s_Q00 + row_pos;
			Q01 = s_Q01 + row_pos;
			Q10 = s_Q10 + row_pos;
			Q11 = s_Q11 + row_pos;

			*tmp_elems = ((*Q00)+(*Q11)-(*Q01)-(*Q10))/c1 - (*p_vec_dy)*(*p_vec_dy)/c2 + min_identity_value;

			*tmp_elems = sqrt(*tmp_elems);

			tmp_elems++;

			p_vec_dy++;
		}
	}
	else {
		for ( d = 0 ; d < sub_feature_num ; d++ ) {

			row_pos = m_pRowStartPos[*rows_idxes++];

			Q00 = s_Q00 + row_pos;
			Q01 = s_Q01 + row_pos;
			Q10 = s_Q10 + row_pos;
			Q11 = s_Q11 + row_pos;

			*tmp_elems = ((*Q00)+(*Q11)-(*Q01)-(*Q10))/c1 - (*p_vec_dy)*(*p_vec_dy)/c2 + min_identity_value;

			*tmp_elems = sqrt(*tmp_elems);

			tmp_elems++;

			p_vec_dy++;
		}
	}

	//delete [] p_vec;
}

void CIntegralRegionCov::Save(char *file_name, CvRect *roi)
{
	FILE * pFile = fopen(file_name, "wb");
	if ( !pFile ) {
		printf("Open file error : %s\n", file_name);
		exit(1);
	}

	fwrite(&m_nCovDim, sizeof(int), 1, pFile);

	fwrite(&m_szFeatureImage.width, sizeof(int), 1, pFile);
	fwrite(&m_szFeatureImage.height, sizeof(int), 1, pFile);

	int p_vec_length = sizeof(double)*m_nCovDim;
	int Q_mat_length = sizeof(double)*m_nCovQDim;

	IntegralPointStruct *cur_integral_pt;

	if ( roi ) {
		int width_1 = m_szFeatureImage.width+1;
		int offset = (1+roi->y)*width_1 + (1+roi->x);
		for ( int y = 0 ; y < roi->height ; y++ ) {
			cur_integral_pt = m_pIntegralImages + y*width_1 + offset;
			for ( int x = 0 ; x < roi->width ; x++ ) {
				fwrite((*cur_integral_pt).p, p_vec_length, 1, pFile);
				fwrite((*cur_integral_pt).Q, Q_mat_length, 1, pFile);
				cur_integral_pt++;
			}
		}
	}
	else {
		cur_integral_pt = m_pIntegralImages + m_szFeatureImage.width + 2;
		for ( int y = 0 ; y < m_szFeatureImage.height ; y++ ) {
			for ( int x = 0 ; x < m_szFeatureImage.width ; x++ ) {
				fwrite((*cur_integral_pt).p, p_vec_length, 1, pFile);
				fwrite((*cur_integral_pt).Q, Q_mat_length, 1, pFile);
				cur_integral_pt++;
			}
			cur_integral_pt++;
		}
	}

	fclose(pFile);
}

void CIntegralRegionCov::Load(char *file_name)
{
	FILE * pFile = fopen(file_name, "rb");
	if ( !pFile ) {
		printf("Open file error : %s\n", file_name);
		exit(1);
	}

	//fclose(pFile);
	//return;

	int cov_dim;
	CvSize img_size;

	fread(&cov_dim, sizeof(int), 1, pFile);

	fread(&img_size.width, sizeof(int), 1, pFile);
	fread(&img_size.height, sizeof(int), 1, pFile);


	m_szFeatureImage = img_size;

	//Init(cov_dim, img_size);

	int p_vec_length = sizeof(double)*m_nCovDim;
	int Q_mat_length = sizeof(double)*m_nCovQDim;

	IntegralPointStruct *cur_integral_pt = m_pIntegralImages + m_szFeatureImage.width + 2;
	for ( int y = 0 ; y < m_szFeatureImage.height ; y++ ) {
		for ( int x = 0 ; x < m_szFeatureImage.width ; x++ ) {
			fread((*cur_integral_pt).p, p_vec_length, 1, pFile);
			fread((*cur_integral_pt).Q, Q_mat_length, 1, pFile);
			cur_integral_pt++;
		}
		cur_integral_pt++;
	}

	fclose(pFile);
}
