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
// ImageDetections.h: interface for the CImageDetection class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_IMAGE_DETECTION_H_)
#define _IMAGE_DETECTION_H_

#include "BoostClassifierCovariance.h"
#include "NegExamplesGenerator.h"
#include "cv.h"

class CImageDetection
{
public:
    /* human detection on a single image */
    void Detection(CBoostClassifierCovariance* bcc, CTestWindows* tw,
            IplImage* image, IplImage* fp_image=NULL,
            int evualation_start_level=-1,
            int n_cascade_levels=-1, double Xn_prob_alpha=1.0);

    /* draw the merged outputs on the image */
    void Draw(IplImage* image,
            bool show_windows=false,
            int cascade_level=0,
            int radius=2,
            int thinness=1,
            CvScalar center_color=CV_RGB(0,255,0),
            CvScalar win_color=CV_RGB(255,0,255));

    /* export the detection results */
    bool Export(char* file_name);

    /* import the detection results */
    bool Import(char* file_name);

    /* return the found positive detection scores */
    double* GetPosDetScores(int cascade_level=0);

    /* return the found positive detection windows (rectangles) */
    CvRect* GetPosDetWindows(int cascade_level=0);

    /* return the number of found positive detection windows */
    int GetPosDetWindowNum(int cascade_level=0);

    /* update the detection windows */
    void UpdateDetections(int x_marg_size, int y_marg_size,
            double human_left_percent, double human_top_percent,
            int cascade_level=-1);

    /* insert the margin space on the original image for human detection */
    void InsertImageMarginSpace(
            IplImage* image, IplImage** marg_image,
            int left_marg_size, int right_marg_size,
            int top_marg_size, int bottom_marg_size);

    /* constructor */
	CImageDetection();

    /* destructor */
	virtual ~CImageDetection();

private:
    /* the positive scores for the found detection
       candidates at each cascade level */
    double** m_ppDetScores;

    /* the positive windows for the found detection
       candidates at each cascade level */
    CvRect** m_ppDetWindows;

    /* the numbers of the found detection
       candidates at each cascade level  */
    int*     m_pCascadeDetNums;

    /* the maximal (for memory allocation) number of
       the found detection candidates at each cascade level */
    int*     m_pCascadeDetMemNums;

    /* the maximal (for memory allocation) number of
       the cascade levels in the learned model */
    int      m_nMaxMemCascadeLevelsNum;

    /* the number of the cascade levels in the learned model */
    int      m_nCascadeLevelsNum;

    /* the image size */
    CvSize   m_szImage;

    /* update the data */
    void UpdateData(double** scores, CvRect** windows, int cur_len, int &new_len);
};

#endif // !defined(_IMAGE_DETECTION_H_)
