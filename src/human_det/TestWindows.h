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
// TestWindows.h: interface for the CTestWindows class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_TEST_WINDOWS_H_)
#define _TEST_WINDOWS_H_

#include "cv.h"
#include "IntegralImage.h"

class CTestWindows
{
public:
    /* constructor method 1 */
	CTestWindows(int min_win_width,
		int n_win_width_scales,
		float w_scale_factor,
		float min_hw_ratio,
		int n_hw_scales,
		float h_scale_factor,
		float w_shift_factor,
		float h_shift_factor);

    /* constructor method 2 */
	CTestWindows(int min_win_width,
        int max_win_width,
		int n_win_width_scales,
		float min_hw_ratio,
		float max_hw_ratio,
		int n_hw_scales,
		float w_shift_factor,
		float h_shift_factor);

    /* constructor method 3 */
	CTestWindows(float min_win_width_percent,
        float max_win_width_percent,
		int n_win_width_scales,
		float min_hw_ratio,
		float max_hw_ratio,
		int n_hw_scales,
		float w_shift_factor,
		float h_shift_factor);

	virtual ~CTestWindows();

    /* get test windows */
	void GetTestWindows(int min_win_width, CvRect image_roi);
	void GetTestWindows(CvRect image_roi);
	void GetTestWindows(CvSize image_size);

    /* load the camera parameters */
    CvMat* LoadCameraParas(char* file_name, CvSize org_img_size,
        double x_scale, double y_scale,
        int left_marg_size, int top_marg_size);

    /* filter out test windows
       based on the camera calibration */
	void FilterTestWindows(CvMat *P_double,
        float vertical_bot_factor,
        float min_height, float max_height);

    /* filter out test windows based on
       an mask image and a threshold
       (useful for warped image) */
    void FilterTestWindows(IplImage *mask,
        float percent_threshold,
        float top_bot_space_percent);

    /* remove too large test windows */
    void FilterTestWindows(int max_win_width);
	
	
	
	//MODIF ALEXANDRE
	//ELIMINATE TOO LARGE OR TOO SMALL TEST WINDOWS IN FUNCTION OF LINE IN IMAGE, FOR EXAMPLE WHEN CALIBRATION NOT AVAILABLE
	
	
    void FilterTestWindows(int L1, int L2, int h1, int h2, double tolerance); // simple case where we assume that the relation between the size of the window and 
	//the line in the image plane is linear. 
	//we allow a certain tolerance with regard to the average height obtained by linear regression at a given image line , for example 0.2 (20%). 
	
	//we consider the bottom of bounding boxes
	
	// END MODIF ALEXANDRE
	
	
	
	
	
	

    /* set the test windows */
    void SetTestWindows(CvRect* det_windows, int length);

    /* update the test windows via shift */
    void UpdateTestWindows(int x_offset, int y_offset);

    /* rectagulars of test windows */
	CvRect* m_pTestWindows;

	/* number of test windows */
	int m_nTestWindowsNum;

private:
    /* get image y-value of top poing of
       the object with a given height */
    double GetObjTopPointY(CvMat *P, CvMat *H,
        CvPoint2D64f bot_img_pt, double obj_height);

    /* minimal test window width (in pixels) */
	int m_nMinWinWidth;

    /* number of scales in width */
	int m_nWinWidthScalesNum;

	/* minimal height-width ratio of test windows */
	float m_fMinHWRatio;

	/* number of scales in height-width ratio */
	int m_nHWScalesNum;

	/* minimal/maximal percent of test window width
	   w.r.t. the whole image size */
    float m_fMinWidWidthImagePercent;
    float m_fMaxWidWidthImagePercent;

	/* width/height scale factor (0,1) */
	float m_fWScaleFactor;
	float m_fHScaleFactor;

	/* width/height shift factor of test windows */
	float m_fWShiftFactor;
	float m_fHShiftFactor;
};

#endif // !defined(_TEST_WINDOWS_H_)
