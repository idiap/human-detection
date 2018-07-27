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
// TestWindows.cpp: implementation of the CTestWindows class.
//
//////////////////////////////////////////////////////////////////////

#include "TestWindows.h"
#include "QuickSort.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

/* constructor method 1 */
CTestWindows::CTestWindows(
        int min_win_width,
		int n_win_width_scales,
		float w_scale_factor,
		float min_hw_ratio,
		int n_hw_scales,
		float h_scale_factor,
		float w_shift_factor,
		float h_shift_factor)
{
	m_nMinWinWidth = min_win_width;
	m_nWinWidthScalesNum = n_win_width_scales;

	m_fMinHWRatio = min_hw_ratio;
	m_nHWScalesNum = n_hw_scales;

	m_fWScaleFactor = w_scale_factor;
	m_fHScaleFactor = h_scale_factor;

	m_fWShiftFactor = w_shift_factor;
	m_fHShiftFactor = h_shift_factor;

	m_pTestWindows = NULL;
	m_nTestWindowsNum = 0;

    m_fMinWidWidthImagePercent = -1.0f;
    m_fMaxWidWidthImagePercent = -1.0f;
};

/* constructor method 2 */
CTestWindows::CTestWindows(
        int min_win_width,
        int max_win_width,
        int n_win_width_scales,
        float min_hw_ratio,
        float max_hw_ratio,
        int n_hw_scales,
        float w_shift_factor,
        float h_shift_factor)
{
	m_nMinWinWidth = min_win_width;
	m_nWinWidthScalesNum = n_win_width_scales;
	m_fWScaleFactor = powf((float)max_win_width/(float)min_win_width,
            1.0f/(float)m_nWinWidthScalesNum) - 1.0f;

	m_fMinHWRatio = min_hw_ratio;
	m_nHWScalesNum = n_hw_scales;
	m_fHScaleFactor = powf(max_hw_ratio/min_hw_ratio,
            1.0f/(float)m_nHWScalesNum) - 1.0f;

	m_fWShiftFactor = w_shift_factor;
	m_fHShiftFactor = h_shift_factor;

	m_pTestWindows = NULL;
	m_nTestWindowsNum = 0;

    m_fMinWidWidthImagePercent = -1.0f;
    m_fMaxWidWidthImagePercent = -1.0f;
}

/* constructor method 3 */
CTestWindows::CTestWindows(
        float min_win_width_percent,
        float max_win_width_percent,
        int n_win_width_scales,
        float min_hw_ratio,
        float max_hw_ratio,
        int n_hw_scales,
        float w_shift_factor,
        float h_shift_factor)
{
    /* for default, we use standard image size */
    CvSize image_size = cvSize(640,480);

    m_fMinWidWidthImagePercent = min_win_width_percent;
    m_fMaxWidWidthImagePercent = max_win_width_percent;

	m_nMinWinWidth = (int)(min_win_width_percent*(float)image_size.width);
	m_nWinWidthScalesNum = n_win_width_scales;
	m_fWScaleFactor = powf(max_win_width_percent/min_win_width_percent,
            1.0f/(float)m_nWinWidthScalesNum) - 1.0f;

	m_fMinHWRatio = min_hw_ratio;
	m_nHWScalesNum = n_hw_scales;
	m_fHScaleFactor = powf(max_hw_ratio/min_hw_ratio,
            1.0f/(float)m_nHWScalesNum) - 1.0f;

	m_fWShiftFactor = w_shift_factor;
	m_fHShiftFactor = h_shift_factor;

	m_pTestWindows = NULL;
	m_nTestWindowsNum = 0;
}


CTestWindows::~CTestWindows()
{
	if ( m_pTestWindows )
		delete [] m_pTestWindows;
}

void CTestWindows::GetTestWindows(int min_win_width, CvRect image_roi)
{
    m_nMinWinWidth = min_win_width;
    GetTestWindows( image_roi );
}

void CTestWindows::GetTestWindows(CvSize image_size)
{
    if ( m_fMinWidWidthImagePercent > 0 && m_fMinWidWidthImagePercent < 1.0f )
        m_nMinWinWidth = (int)(m_fMinWidWidthImagePercent*(float)image_size.width);

    GetTestWindows(cvRect(0,0,image_size.width,image_size.height));
}

void CTestWindows::GetTestWindows(CvRect image_roi)
{
    /* computing the total number of test windows */
    m_nTestWindowsNum = 0;
    for ( int w = 0 ; w < m_nWinWidthScalesNum ; w++ ) {
        int win_width = (int)((float)m_nMinWinWidth*powf(1.0f+m_fWScaleFactor,w));
        if ( win_width > image_roi.width )
            break;
        for ( int h = 0 ; h < m_nHWScalesNum ; h++ ) {
            int win_height = (int)((float)win_width*m_fMinHWRatio*powf(1.0f+m_fHScaleFactor,h));
            if ( win_height > image_roi.height )
                break;
            int x_sample_step = (int)((float)win_width * m_fWShiftFactor);
            int y_sample_step = (int)((float)win_height * m_fHShiftFactor);
            for ( int y = 0 ; y < image_roi.height ; y += y_sample_step ) {
                for ( int x = 0 ; x < image_roi.width ; x += x_sample_step ) {
                    m_nTestWindowsNum++;
                }
            }
        }
    }

    /* get test windows */
    if ( m_pTestWindows )
        delete [] m_pTestWindows;
    m_pTestWindows = new CvRect[m_nTestWindowsNum];
    m_nTestWindowsNum = 0;
    for ( int w = 0 ; w < m_nWinWidthScalesNum ; w++ ) {
        int win_width = (int)((float)m_nMinWinWidth*powf(1.0f+m_fWScaleFactor,w));
        if ( win_width > image_roi.width )
            break;
        for ( int h = 0 ; h < m_nHWScalesNum ; h++ ) {
            int win_height = (int)((float)win_width*m_fMinHWRatio*powf(1.0f+m_fHScaleFactor,h));
            if ( win_height > image_roi.height )
                break;

            int x_sample_step = (int)((float)win_width * m_fWShiftFactor);
            int y_sample_step = (int)((float)win_height * m_fHShiftFactor);

            int x, y;
            int w_left_space, h_left_space;

            for ( y = 0 ; y < image_roi.height - win_height ; y += y_sample_step ) {
                for ( x = 0 ; x < image_roi.width - win_width; x += x_sample_step ) {
                    m_pTestWindows[m_nTestWindowsNum++] = cvRect(x+image_roi.x,
                        y+image_roi.y, win_width, win_height);
                }
                w_left_space = (image_roi.width - win_width) - (x - x_sample_step);
                if ( w_left_space >= 2 ) {
                    m_pTestWindows[m_nTestWindowsNum++] = cvRect(image_roi.x + image_roi.width - win_width - 1 ,
                        y+image_roi.y, win_width, win_height);
                }
            }

            h_left_space = (image_roi.height - win_height) - (y - y_sample_step);
            if ( h_left_space >= 2 ) {
                y = (image_roi.y + image_roi.height - win_height - 1);
                for ( x = 0 ; x < image_roi.width - win_width; x += x_sample_step ) {
                    m_pTestWindows[m_nTestWindowsNum++] = cvRect(x+image_roi.x,
                        y+image_roi.y, win_width, win_height);
                }
                w_left_space = (image_roi.width - win_width) - (x - x_sample_step);
                if ( w_left_space >= 2 ) {
                    m_pTestWindows[m_nTestWindowsNum++] = cvRect(image_roi.x + image_roi.width - win_width - 1,
                        y+image_roi.y, win_width, win_height);
                }
            }


        }
    }
}

double CTestWindows::GetObjTopPointY(CvMat *P, CvMat *H,
        CvPoint2D64f bot_img_pt, double obj_height)
{
	double u_g[2];

	double u_x, u_y, u_z;

    double* H_ptr = (double*)H->data.ptr;

    /* get corresponding ground world point */
	u_x = H_ptr[0]*bot_img_pt.x + H_ptr[1]*bot_img_pt.y + H_ptr[2];
	u_y = H_ptr[3]*bot_img_pt.x + H_ptr[4]*bot_img_pt.y + H_ptr[5];
	u_z = H_ptr[6]*bot_img_pt.x + H_ptr[7]*bot_img_pt.y + H_ptr[8];

	u_g[0] = u_x/u_z;
	u_g[1] = u_y/u_z;

    double* P_ptr = (double*)P->data.ptr;

	u_y = P_ptr[4]*u_g[0] + P_ptr[5]*u_g[1] + P_ptr[6]*obj_height + P_ptr[7];
	u_z = P_ptr[8]*u_g[0] + P_ptr[9]*u_g[1] + P_ptr[10]*obj_height + P_ptr[11];

	double top_img_y = u_y/u_z;

	return top_img_y;
}

CvMat* CTestWindows::LoadCameraParas(char* file_name, CvSize org_img_size,
        double x_scale, double y_scale,
        int left_marg_size, int top_marg_size)
{
    if ( !file_name )
        return NULL;

    CvFileStorage* yaml = cvOpenFileStorage( file_name, 0, CV_STORAGE_READ);
    if(yaml==NULL){
        //printf("File %s could not be read\n", file_name);
        return NULL;
    }
    CvMat* projection_matrix = (CvMat*) cvReadByName(yaml, NULL, "projection_matrix");
    cvReleaseFileStorage(&yaml);

    CvMat ST;
    double scale_t[] = {x_scale, 0.0, left_marg_size ,
                0.0, y_scale, top_marg_size,
                0.0, 0.0, 1.0};
    cvInitMatHeader(&ST, 3, 3, CV_64FC1, scale_t);
    cvMatMul(&ST, projection_matrix, projection_matrix);

    return projection_matrix;
}

void CTestWindows::FilterTestWindows(CvMat* P_double,
        float vertical_bot_factor,
        float min_height, float max_height)
{
    if ( m_nTestWindowsNum == 0 || !P_double )
        return;

    CvMat* wtoi_H = cvCreateMat(3, 3, CV_64F);
    CvMat* itow_H = cvCreateMat(3, 3, CV_64F);

    cvmSet(wtoi_H, 0, 0, cvmGet(P_double, 0, 0));
    cvmSet(wtoi_H, 0, 1, cvmGet(P_double, 0, 1));
    cvmSet(wtoi_H, 0, 2, cvmGet(P_double, 0, 3));
    cvmSet(wtoi_H, 1, 0, cvmGet(P_double, 1, 0));
    cvmSet(wtoi_H, 1, 1, cvmGet(P_double, 1, 1));
    cvmSet(wtoi_H, 1, 2, cvmGet(P_double, 1, 3));
    cvmSet(wtoi_H, 2, 0, cvmGet(P_double, 2, 0));
    cvmSet(wtoi_H, 2, 1, cvmGet(P_double, 2, 1));
    cvmSet(wtoi_H, 2, 2, cvmGet(P_double, 2, 3));

    //cvGetSubRect(P_double, wtoi_H, cvRect(0,0,3,3));
    cvInv(wtoi_H, itow_H);

    bool *labels = new bool[m_nTestWindowsNum];
    int n_filter_test_windows = 0;
    CvPoint2D64f bot_img_pt;
    int min_obj_y, max_obj_y;
    for ( int i = 0 ; i < m_nTestWindowsNum ; i++ ) {
        /* get bottom image point of current test window */
        bot_img_pt.x = (double)m_pTestWindows[i].x +
                (double)(m_pTestWindows[i].width-1)/2.0;
        bot_img_pt.y = (double)m_pTestWindows[i].y +
                (double)m_pTestWindows[i].height * vertical_bot_factor;

        /* get minimal/maximal object Y */
        max_obj_y = (int)GetObjTopPointY(P_double, itow_H, bot_img_pt, min_height);
        min_obj_y = (int)GetObjTopPointY(P_double, itow_H, bot_img_pt, max_height);

        /* check whether the object height inbetween given range */
        int pos_y = m_pTestWindows[i].y+(int)(m_pTestWindows[i].height*(1.0-vertical_bot_factor));
        labels[i] = ( pos_y >= min_obj_y &&
                pos_y <= max_obj_y );

        if ( labels[i] )
            n_filter_test_windows++;
    }

    /* remove invalid test windows */
    if ( n_filter_test_windows > 0 ) {
        CvRect* filter_test_windows = new CvRect[n_filter_test_windows];
        n_filter_test_windows = 0;
        for ( int i = 0 ; i < m_nTestWindowsNum ; i++ ) {
            if ( labels[i] )
                filter_test_windows[n_filter_test_windows++] = m_pTestWindows[i];
        }
        delete [] m_pTestWindows;
        m_pTestWindows = filter_test_windows;
        m_nTestWindowsNum = n_filter_test_windows;
    }
    else {
        delete [] m_pTestWindows;
        m_pTestWindows = NULL;
        m_nTestWindowsNum = 0;
    }

    /* release memories */
    cvReleaseMat(&wtoi_H);
    cvReleaseMat(&itow_H);

    delete [] labels;
}

void CTestWindows::FilterTestWindows(IplImage *mask, float percent_threshold, float top_bot_space_percent)
{
    CIntegralImage* int_img = new CIntegralImage(cvGetSize(mask), true);
    int_img->ComputeIntegralImage(mask);

    bool *labels = new bool[m_nTestWindowsNum];
    int n_filter_test_windows = 0;

    for ( int i = 0 ; i < m_nTestWindowsNum ; i++ ) {
        int x_marg_space = (int)(top_bot_space_percent*(float)m_pTestWindows[i].width);
        int y_marg_space = (int)(top_bot_space_percent*(float)m_pTestWindows[i].height);
        CvRect test_win = cvRect(m_pTestWindows[i].x + x_marg_space,
                            m_pTestWindows[i].y + y_marg_space,
                            m_pTestWindows[i].width - 2*x_marg_space,
                            m_pTestWindows[i].height - 2*y_marg_space);
        float win_mask_sum = int_img->GetRegionIntegralSum(test_win);
        float win_mask_ratio = win_mask_sum/(float)(test_win.width*test_win.height);
        labels[i] = win_mask_ratio > percent_threshold;
        if ( labels[i] )
            n_filter_test_windows++;
    }

    /* remove invalid test windows */
    if ( n_filter_test_windows > 0 ) {
        CvRect* filter_test_windows = new CvRect[n_filter_test_windows];
        n_filter_test_windows = 0;
        for ( int i = 0 ; i < m_nTestWindowsNum ; i++ ) {
            if ( labels[i] )
                filter_test_windows[n_filter_test_windows++] = m_pTestWindows[i];
        }
        delete [] m_pTestWindows;
        m_pTestWindows = filter_test_windows;
        m_nTestWindowsNum = n_filter_test_windows;
    }
    else {
        delete [] m_pTestWindows;
        m_pTestWindows = NULL;
        m_nTestWindowsNum = 0;
    }

    /* release memories */
    delete [] labels;
    delete int_img;
}

void CTestWindows::FilterTestWindows(int max_win_width)
{
    bool *labels = new bool[m_nTestWindowsNum];
    int n_filter_test_windows = 0;

    for ( int i = 0 ; i < m_nTestWindowsNum ; i++ ) {
        labels[i] = m_pTestWindows[i].width < max_win_width;
        if ( labels[i] )
            n_filter_test_windows++;
    }

    /* remove invalid test windows */
    if ( n_filter_test_windows > 0 ) {
        CvRect* filter_test_windows = new CvRect[n_filter_test_windows];
        n_filter_test_windows = 0;
        for ( int i = 0 ; i < m_nTestWindowsNum ; i++ ) {
            if ( labels[i] )
                filter_test_windows[n_filter_test_windows++] = m_pTestWindows[i];
        }
        delete [] m_pTestWindows;
        m_pTestWindows = filter_test_windows;
        m_nTestWindowsNum = n_filter_test_windows;
    }
    else {
        delete [] m_pTestWindows;
        m_pTestWindows = NULL;
        m_nTestWindowsNum = 0;
    }

    /* release memories */
    delete [] labels;
}





//MODIF ALEXANDRE
//ELIMINATE TOO LARGE OR TOO SMALL TEST WINDOWS IN FUNCTION OF LINE IN IMAGE, IF CALIBRATION NOT AVAILABLE


 
 void CTestWindows::FilterTestWindows(int L1, int L2, int h1, int h2, double tolerance)
 {
	 double alpha,beta;
	 alpha = (double)(h2-h1)/(double)(L2-L1);
	 beta = h1 - alpha*L1;
	 bool *labels = new bool[m_nTestWindowsNum]; //false if window needs to be filtered (rejected)
	 int img_line_number;
	 int n_filter_test_windows = 0;
	 
	 for ( int i = 0 ; i < m_nTestWindowsNum ; i++ ) {
		 img_line_number = m_pTestWindows[i].y + m_pTestWindows[i].height; //y=line coordinate of top-left corner, height=height of rectangle
		 if ( m_pTestWindows[i].height > alpha * img_line_number + beta + tolerance*(alpha * img_line_number + beta) || m_pTestWindows[i].height < alpha * img_line_number + beta - tolerance*(alpha * img_line_number + beta))
			 labels[i] = false;
		 else
			 labels[i] = true;
		 if ( labels[i] == true)
			 n_filter_test_windows++;
	 }
	 
	 // remove invalid test windows 
	 if ( n_filter_test_windows > 0 ) {
		 CvRect* filter_test_windows = new CvRect[n_filter_test_windows];
		 n_filter_test_windows = 0;
		 for ( int i = 0 ; i < m_nTestWindowsNum ; i++ ) {
			 if ( labels[i] == true)
				 filter_test_windows[n_filter_test_windows++] = m_pTestWindows[i]; //kept windows
		 }
		 delete [] m_pTestWindows;
		 m_pTestWindows = filter_test_windows;
		 m_nTestWindowsNum = n_filter_test_windows;
	 }
	 else { //case where all the windows have been rejected
		 delete [] m_pTestWindows;
		 m_pTestWindows = NULL;
		 m_nTestWindowsNum = 0;
	 }
	 
	 // release memories 
	 delete [] labels;
 }

 
 // END MODIF ALEXANDRE
 







void CTestWindows::SetTestWindows(CvRect* det_windows, int length)
{
    if ( m_nTestWindowsNum != length ) {
        delete [] m_pTestWindows;
        m_pTestWindows = new CvRect[length];
        m_nTestWindowsNum = length;
    }
    for ( int i = 0 ; i < length ; i++ )
        m_pTestWindows[i] = det_windows[i];
}

void CTestWindows::UpdateTestWindows(int x_offset, int y_offset)
{
    for ( int i = 0 ; i < m_nTestWindowsNum ; i++ ) {
        m_pTestWindows[i].x += x_offset;
        m_pTestWindows[i].y += y_offset;
    }
}
