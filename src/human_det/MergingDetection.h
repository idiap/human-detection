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
// MerginDetections.h: interface for the CMergingDetection class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_MERGING_DETECTION_H_)
#define _MERGING_DETECTION_H_

#include "BoostClassifierCovariance.h"
#include "NegExamplesGenerator.h"
#include "cv.h"
#include "Rectangle.h"
#include "QuickSort.h"

class CMergingDetection
{
public:
    /* post-processing the human detection candidates with the windows and scores */
    void Merging(CvRect* det_windows, double* det_scores, int det_windows_num);

    /* return the detection scores of the merged output */
    double* GetMergedDetScores();

    /* return the detection windows of the merged output */
    CvRect* GetMergedDetWindows();

    /* return the number of merged outputs */
    int GetMergedDetWindowNum();

    /* draw the merged outputs on the image */
    void Draw(IplImage* image,
            bool show_windows=true,
            int radius=2,
            int thinness=2,
            CvScalar center_color=CV_RGB(255,0,0),
            CvScalar win_color=CV_RGB(255,0,0));

    /* export the detection results */
    bool Export(char* file_name);

    /* import the detection results */
    bool Import(char* file_name);

    /* import the detection results for multi-person tracking */
	bool Import(char* det_res_fn, int &det_num,
            int neig_det_num_threshold=4, int min_det_win_area=600,
            float *app_probs=NULL, CvRect *obj_rects=NULL);

    void UpdateDetections(int x_marg_size, int y_marg_size,
        double human_left_percent, double human_top_percent);

    /* return the smoothed detection score map (image) */
    IplImage* GetSmoothedScoreImage();

    /* constructor */
	CMergingDetection(CvSize img_size,
            double kernel_size_percent=0.1,
            double sigma=5.0,
            double overlapping_threshold=0.9,
            double front_overlapping_threshold=0.6,
            double front_score_large_scale=1.5);

    CMergingDetection();

    /* destructor */
	virtual ~CMergingDetection();

private:
    /* information associated with the merged outputs */
    double*  m_pMergedDetScores;    /* detection scores */
    int*     m_pMergedNeigDetNums;  /* numbers of neighboring detections */
    CvRect*  m_pMergedDetWindows;   /* detection windows */
    CvPoint* m_pMergedDetCenterPoints;  /* center points of detection windows */
    int      m_nMergedDetNum;       /* the number of the merged outputs */
    int      m_nMergedDetMemNum;    /* the maximal (for memory allocation) number of the merged output */

    /* the smoothing parameters on the detection score map */
    double   m_dbKernelSizePercent;
    double   m_dbSmoothSigma;

    /* the filtering threshold */
    double   m_dbOverlappingThreshold;
    double   m_dbFrontOverlappingThreshold;
    double   m_dbFrontScoreLargeScale;

    /* the detection score map */
    IplImage* m_pDetScoresImage;

    /* the mask image on which the merged outputs are */
    IplImage* m_pDetMaskImage;

    /* update the data */
    void UpdateData();

    /* get the detection windows of the merged outputs */
    void ComputeMergedDetWindows(CvRect* det_windows,
            CvPoint* det_center_pts, double* det_scores, int det_windows_num);

    /* find the local maxima on the smoothed detection score map */
    void FindLocalMaxima(int half_size=2);

    /* remove too closed merged outputs */
    void RemoveTooClosedDetections();

    /* compute the F-measure between two rectangles */
    double GetFMeasure(CvRect rect1, CvRect rect2);

    double GetArea(CvRect rect);
};

#endif // !defined(_MERGING_DETECTION_H_)
