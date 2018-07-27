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
// MerginDetections.cpp: implementation for the CMergingDetection class.
//
//////////////////////////////////////////////////////////////////////

#include "MergingDetection.h"

/* post-processing the human detection candidates with the windows and scores */
void CMergingDetection::Merging(CvRect* det_windows, double* det_scores, int det_windows_num) {
    if ( det_windows_num <= 0 ) {
        m_nMergedDetNum = 0;
        return;
    }

    /* compute the center points of the detection windows */
    CvPoint* det_center_pts = new CvPoint[det_windows_num];
    double avg_width = 0, avg_height = 0;
    for ( int i = 0 ; i < det_windows_num ; i++ ) {
        det_center_pts[i] = cvPoint(det_windows[i].x+det_windows[i].width/2,
                        det_windows[i].y+det_windows[i].height/2);
        avg_width += det_windows[i].width;
        avg_height += det_windows[i].height;
    }

    /* compute the average window */
    avg_width /= (double)det_windows_num;
    avg_height /= (double)det_windows_num;

    /* initialize the data */
    cvSetZero(m_pDetScoresImage);
    cvSetZero(m_pDetMaskImage);

    /* set the scores and mask */
    for ( int i = 0 ; i < det_windows_num ; i++ ) {
        CvPoint pt = det_center_pts[i];
        float* tpt = &(((float*)(m_pDetScoresImage->imageData + m_pDetScoresImage->widthStep*pt.y))[pt.x]);
        *tpt = MAX((float)det_scores[i], *tpt);
        ((uchar*)(m_pDetMaskImage->imageData + m_pDetMaskImage->widthStep*pt.y))[pt.x] = 1;
    }

    /* dilate the mask image which is used for finding local maxima */
    cvDilate(m_pDetMaskImage, m_pDetMaskImage, 0, MAX(m_pDetScoresImage->height/50+1, 3) );

    /* smooth the detection score map */
    cvSmooth( m_pDetScoresImage, m_pDetScoresImage,
           CV_GAUSSIAN,
           (int)(avg_width*m_dbKernelSizePercent)*2+1,
           (int)(avg_height*m_dbKernelSizePercent)*2+1,
           m_dbSmoothSigma );

    /* find the local maxima on the smoothed detection score map */
    FindLocalMaxima();

    /* get the merged detection windows */
    ComputeMergedDetWindows(det_windows, det_center_pts, det_scores, det_windows_num);

    RemoveTooClosedDetections();

    /* release memory */
    delete [] det_center_pts;
}

/* return the detection scores of the merged output */
double* CMergingDetection::GetMergedDetScores() {
    return m_pMergedDetScores;
}

/* return the detection windows of the merged output */
CvRect* CMergingDetection::GetMergedDetWindows() {
    return m_pMergedDetWindows;
}

/* return the number of merged outputs */
int CMergingDetection::GetMergedDetWindowNum() {
    return m_nMergedDetNum;
}

/* draw the merged outputs on the image */
void CMergingDetection::Draw(IplImage* image,
        bool show_windows,
        int radius,
        int thinness,
        CvScalar center_color,
        CvScalar win_color) {
    CvPoint pt1, pt2, pt;
    for ( int i = 0 ; i < m_nMergedDetNum ; i++ ) {
        pt1 = cvPoint(m_pMergedDetWindows[i].x, m_pMergedDetWindows[i].y);
        pt2 = cvPoint(pt1.x+m_pMergedDetWindows[i].width-1,
                pt1.y+m_pMergedDetWindows[i].height-1);
        //cvRectangle(image, pt1, pt2, win_color, thinness);

        pt = m_pMergedDetCenterPoints[i];
        cvCircle( image, pt, radius, center_color, thinness);

        CvSize person_size = cvSize( cvCeil((double)m_pMergedDetWindows[i].width/2.0),
            cvCeil((double)m_pMergedDetWindows[i].height/2.0) );

        if ( show_windows )
            cvRectangle(image, pt1, pt2, win_color, thinness);
        else
            cvEllipse(image, pt, person_size, 0, 0, 360, win_color, thinness);
    }
}

/* export the detection results */
bool CMergingDetection::Export(char* file_name) {
    if ( !file_name ) {
        printf("Please provide file name (XML or YAML format) to export the detection results!\n");
        return false;
    }

    CvFileStorage* fs = cvOpenFileStorage( file_name, 0, CV_STORAGE_WRITE );
    if ( !fs ) {
        printf("Open file error : %s\n", file_name);
        return false;
    }

    char label_name[1024];

    ///////////////////////////////////////////////////////////////////////////////
    cvWriteComment( fs, "Merged human detection outputs on a single image", 0 );
    cvStartWriteStruct( fs, "merged_human_detection", CV_NODE_MAP, NULL, cvAttrList(0,0));

    cvWriteComment( fs, "Size of the input image", 0 );
    cvWriteInt( fs, "image_width", m_pDetScoresImage->width );
    cvWriteInt( fs, "image_height", m_pDetScoresImage->height );

    cvWriteComment( fs, "Number of merged detection outputs", 0 );
    cvWriteInt( fs, "n_det_windows", m_nMergedDetNum );

    for ( int i = 0 ; i < m_nMergedDetNum ; i++ ) {
        sprintf(label_name, "detection_%d", i);
        cvStartWriteStruct( fs, label_name, CV_NODE_MAP, NULL, cvAttrList(0,0));

        cvWriteComment( fs, "Positive detection score", 0 );
        cvWriteReal( fs, "det_score", m_pMergedDetScores[i] );

        cvWriteComment( fs, "Neighboring detection number", 0 );
        cvWriteInt( fs, "neig_det_num", m_pMergedNeigDetNums[i] );

        cvWriteComment( fs, "Positive sub-window", 0 );
        cvWriteInt( fs, "roi_x", m_pMergedDetWindows[i].x );
        cvWriteInt( fs, "roi_y", m_pMergedDetWindows[i].y );
        cvWriteInt( fs, "roi_width", m_pMergedDetWindows[i].width );
        cvWriteInt( fs, "roi_height", m_pMergedDetWindows[i].height );

        cvWriteComment( fs, "Center points", 0 );
        cvWriteInt( fs, "center_pt_x", m_pMergedDetCenterPoints[i].x );
        cvWriteInt( fs, "center_pt_y", m_pMergedDetCenterPoints[i].y );

        cvEndWriteStruct( fs );
    }


    cvEndWriteStruct( fs );

    cvReleaseFileStorage( &fs );

    return true;

}

bool CMergingDetection::Import(char *det_res_fn, int &det_num,
        int neig_det_num_threshold, int min_det_win_area,
        float *app_probs, CvRect *obj_rects)
{
	if ( !det_res_fn || !Import(det_res_fn) ) {
		det_num = 0;
		return false;
	}

	det_num = 0;
	for ( int i = 0 ; i < m_nMergedDetNum ; i++ ) {
	    if ( m_pMergedNeigDetNums[i] >= neig_det_num_threshold &&
                GetArea(m_pMergedDetWindows[i]) >= (double)min_det_win_area ) {
            if ( app_probs )
                app_probs[det_num] = m_pMergedDetScores[i];
            if ( obj_rects )
                obj_rects[det_num] = m_pMergedDetWindows[i];
            det_num++;
        }
	}
	return true;
}

/* import the detection results */
bool CMergingDetection::Import(char* file_name) {
    if ( !file_name ) {
        printf("Please provide file name (XML or YAML format) to import the detection results!\n");
        return false;
    }

    CvFileStorage* fs = cvOpenFileStorage( file_name, 0, CV_STORAGE_READ );
    if ( !fs ) {
        //printf("Open file error : %s\n", file_name);
        return false;
    }

    char label_name[1024];

    ///////////////////////////////////////////////////////////////////////////////
    CvFileNode* merged_human_det_node = cvGetFileNodeByName( fs, NULL, "merged_human_detection" );

    int img_width = cvReadIntByName( fs, merged_human_det_node, "image_width" );
    int img_height = cvReadIntByName( fs, merged_human_det_node, "image_height" );

    //if ( img_width != m_pDetScoresImage->width || img_height != m_pDetScoresImage->height ) {
    //    printf("Not the same image size!\n");
    //    exit(1);
    //}

    m_nMergedDetNum = cvReadIntByName( fs, merged_human_det_node, "n_det_windows" );

    UpdateData();

    for ( int i = 0 ; i < m_nMergedDetNum ; i++ ) {
        sprintf(label_name, "detection_%d", i);
        CvFileNode* det_node = cvGetFileNodeByName( fs, merged_human_det_node, label_name );

        m_pMergedDetScores[i] = cvReadRealByName( fs, det_node, "det_score" );
        m_pMergedNeigDetNums[i] = cvReadIntByName( fs, det_node, "neig_det_num" );

        m_pMergedDetWindows[i].x = cvReadIntByName( fs, det_node, "roi_x" );
        m_pMergedDetWindows[i].y = cvReadIntByName( fs, det_node, "roi_y" );
        m_pMergedDetWindows[i].width = cvReadIntByName( fs, det_node, "roi_width" );
        m_pMergedDetWindows[i].height = cvReadIntByName( fs, det_node, "roi_height" );

        m_pMergedDetCenterPoints[i].x = cvReadIntByName( fs, det_node, "center_pt_x" );
        m_pMergedDetCenterPoints[i].y = cvReadIntByName( fs, det_node, "center_pt_y" );
    }

    cvReleaseFileStorage( &fs );

    return true;
}

void CMergingDetection::UpdateDetections(int x_marg_size, int y_marg_size,
        double human_left_percent, double human_top_percent) {
    if ( x_marg_size == 0 && y_marg_size == 0
        && human_left_percent == 0 && human_top_percent == 0 )
        return;
    for ( int i = 0 ; i < m_nMergedDetNum ; i++ ) {
        m_pMergedDetWindows[i].x -= x_marg_size;
        m_pMergedDetWindows[i].y -= y_marg_size;

        int x_offset = (int)((double)m_pMergedDetWindows[i].width*human_left_percent);
        int y_offset = (int)((double)m_pMergedDetWindows[i].height*human_top_percent);

        m_pMergedDetWindows[i].x += x_offset;
        m_pMergedDetWindows[i].y += y_offset;
        m_pMergedDetWindows[i].width -= 2*x_offset;
        m_pMergedDetWindows[i].height -= 2*y_offset;

        m_pMergedDetCenterPoints[i].x -= x_marg_size;
        m_pMergedDetCenterPoints[i].y -= y_marg_size;
    }
}

/* return the smoothed detection score map (image) */
IplImage* CMergingDetection::GetSmoothedScoreImage() {
    IplImage* image = cvCreateImage(cvGetSize(m_pDetScoresImage), IPL_DEPTH_8U, 1);
    double maxVal = 0;
    double minVal = 0;
    cvMinMaxLoc( m_pDetScoresImage, &minVal, &maxVal, NULL, NULL, m_pDetMaskImage);

    for ( int y = 0 ; y < image->height ; y++ ) {
        uchar* score_uchar = (uchar*)(image->imageData + image->widthStep*y);
        float* score_float = (float*)(m_pDetScoresImage->imageData + m_pDetScoresImage->widthStep*y);
        for ( int x = 0 ; x < image->width ; x++ ) {
            *score_uchar++ = (int)((*score_float++)/maxVal*255.0);
        }
    }
    return image;
}

/* constructor */
CMergingDetection::CMergingDetection() {
    m_pDetScoresImage = NULL;
    m_pDetMaskImage = NULL;

    m_pMergedDetScores = NULL;
    m_pMergedNeigDetNums = NULL;
    m_pMergedDetWindows = NULL;
    m_pMergedDetCenterPoints = NULL;
    m_nMergedDetNum = m_nMergedDetMemNum = 0;
}

CMergingDetection::CMergingDetection(CvSize img_size,
        double kernel_size_percent,
        double sigma,
        double overlapping_threshold,
        double front_overlapping_threshold,
        double front_score_large_scale) {

    m_pDetScoresImage = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
    m_pDetMaskImage = cvCreateImage(img_size, IPL_DEPTH_8U, 1);

    m_pMergedDetScores = NULL;
    m_pMergedNeigDetNums = NULL;
    m_pMergedDetWindows = NULL;
    m_pMergedDetCenterPoints = NULL;
    m_nMergedDetNum = m_nMergedDetMemNum = 0;

    m_dbKernelSizePercent = kernel_size_percent;
    m_dbSmoothSigma = sigma;
    m_dbOverlappingThreshold = overlapping_threshold;

    m_dbFrontOverlappingThreshold = front_overlapping_threshold;
    m_dbFrontScoreLargeScale = front_score_large_scale;
}

/* destructor */
CMergingDetection::~CMergingDetection() {
    if ( m_pMergedDetScores )
        delete [] m_pMergedDetScores;
    if ( m_pMergedNeigDetNums )
        delete [] m_pMergedNeigDetNums;
    if ( m_pMergedDetWindows )
        delete [] m_pMergedDetWindows;
    if ( m_pMergedDetCenterPoints )
        delete [] m_pMergedDetCenterPoints;
    if ( m_pDetScoresImage )
        cvReleaseImage(&m_pDetScoresImage);
    if ( m_pDetMaskImage )
        cvReleaseImage(&m_pDetMaskImage);
}

/* update the data */
void CMergingDetection::UpdateData() {
    if ( m_nMergedDetNum < m_nMergedDetMemNum )
        return;
    if ( m_nMergedDetNum == m_nMergedDetMemNum )
        m_nMergedDetMemNum += 200;
    else
        m_nMergedDetMemNum = m_nMergedDetNum;
    double* new_scores = new double[m_nMergedDetMemNum];
    int*    new_neig_det_nums = new int[m_nMergedDetMemNum];
    CvRect* new_windows = new CvRect[m_nMergedDetMemNum];
    CvPoint* new_center_points = new CvPoint[m_nMergedDetMemNum];
    if ( m_pMergedDetScores ) {
        for ( int i = 0 ; i < m_nMergedDetNum ; i++ ) {
            new_scores[i] = m_pMergedDetScores[i];
            new_neig_det_nums[i] = m_pMergedNeigDetNums[i];
            new_windows[i] = m_pMergedDetWindows[i];
            new_center_points[i] = m_pMergedDetCenterPoints[i];
        }
    }
    if ( m_pMergedDetScores )
        delete [] m_pMergedDetScores;
    if ( m_pMergedNeigDetNums )
        delete [] m_pMergedNeigDetNums;
    if ( m_pMergedDetWindows )
        delete [] m_pMergedDetWindows;
    if ( m_pMergedDetCenterPoints )
        delete [] m_pMergedDetCenterPoints;
    m_pMergedDetScores = new_scores;
    m_pMergedNeigDetNums = new_neig_det_nums;
    m_pMergedDetWindows = new_windows;
    m_pMergedDetCenterPoints = new_center_points;
}

/* find the local maxima on the smoothed detection score map */
void CMergingDetection::FindLocalMaxima(int half_size) {

    int marg_space = half_size*2+1;
    IplImage* img = m_pDetScoresImage;
    IplImage* mask_img = m_pDetMaskImage;
    m_nMergedDetNum = 0;
    float peak_prob, neig_prob;
    uchar* mask_peak;
    for ( int y = marg_space ; y < img->height-2*marg_space ; y++ ) {
        mask_peak = (uchar*)(mask_img->imageData + mask_img->widthStep*y);
        mask_peak += marg_space;
        for ( int x = marg_space ; x < img->width-2*marg_space ; x++ ) {

            if ( *mask_peak == 0 ) {
                mask_peak++;
                continue;
            }

            bool is_peak = true;
            peak_prob = ((float*)(img->imageData + img->widthStep*y))[x];

            if ( peak_prob == 0 ) {
                mask_peak++;
                continue;
            }

            for ( int py = y-half_size ; py <= y+half_size ; py++ ) {
                for ( int px = x-half_size ; px <= x+half_size ; px++ ) {
                    neig_prob = ((float*)(img->imageData + img->widthStep*py))[px];
                    if ( neig_prob > peak_prob ) {
                        is_peak = false;
                        break;
                    }
                }
                if ( !is_peak )
                    break;
            }
            if ( is_peak ) {
                UpdateData();
                m_pMergedDetScores[m_nMergedDetNum] = peak_prob;
                m_pMergedDetCenterPoints[m_nMergedDetNum] = cvPoint(x, y);
                m_nMergedDetNum++;
            }
            mask_peak++;
        }
    }
}

/* get the detection windows of the merged outputs */
void CMergingDetection::ComputeMergedDetWindows(CvRect* det_windows,
        CvPoint* det_center_pts, double* det_scores, int det_windows_num) {
    for ( int i = 0 ; i < m_nMergedDetNum ; i++ ) {
        double width = 0, height = 0;
        double sum_weights = 0;
        double weight;
        m_pMergedNeigDetNums[i] = 0;
        for ( int j = 0 ; j < det_windows_num ; j++ ) {
            int shift_x = abs(det_center_pts[j].x - m_pMergedDetCenterPoints[i].x);
            int shift_y = abs(det_center_pts[j].y - m_pMergedDetCenterPoints[i].y);
            int shift_x_t = MAX(det_windows[j].width/15, 4);
            int shift_y_t = MAX(det_windows[j].height/15, 8);
            if ( shift_x < shift_x_t && shift_y < shift_y_t ) {
                weight = det_scores[j]/(1.0+5.0*(double)(shift_x*shift_y)/(double)(shift_x_t*shift_y_t));
                width += weight*(double)det_windows[j].width;
                height += weight*(double)det_windows[j].height;
                sum_weights += weight;
                m_pMergedNeigDetNums[i]++;
            }
        }
        if ( sum_weights > 0 ) {
            m_pMergedDetWindows[i].width = (int)(width/sum_weights);
            m_pMergedDetWindows[i].height = (int)(height/sum_weights);
            m_pMergedDetWindows[i].x = m_pMergedDetCenterPoints[i].x - m_pMergedDetWindows[i].width/2;
            m_pMergedDetWindows[i].y = m_pMergedDetCenterPoints[i].y - m_pMergedDetWindows[i].height/2;
        }
        else {
            m_pMergedDetWindows[i] = cvRect(0,0,0,0);
        }
    }
}

/* remove too closed merged outputs */
void CMergingDetection::RemoveTooClosedDetections() {
    if ( m_nMergedDetNum <= 0 )
        return;

    double* scores = new double[m_nMergedDetNum];

    memcpy(scores, m_pMergedDetScores, sizeof(double)*m_nMergedDetNum);
    CQuickSort<double, CvRect> QS_rect;
    QS_rect.QuickSort(scores, 0, m_nMergedDetNum-1, false, m_pMergedDetWindows);

    memcpy(scores, m_pMergedDetScores, sizeof(double)*m_nMergedDetNum);
    CQuickSort<double, CvPoint> QS_pt;
    QS_pt.QuickSort(scores, 0, m_nMergedDetNum-1, false, m_pMergedDetCenterPoints);

    memcpy(scores, m_pMergedDetScores, sizeof(double)*m_nMergedDetNum);
    CQuickSort<double, int> QS_int;
    QS_int.QuickSort(scores, 0, m_nMergedDetNum-1, false, m_pMergedNeigDetNums);

    memcpy(m_pMergedDetScores, scores, sizeof(double)*m_nMergedDetNum);

    bool* removed_flags = new bool[m_nMergedDetNum];
    for ( int i = 0 ; i < m_nMergedDetNum ; i++ ) {
        removed_flags[i] = false;
        if ( m_pMergedDetScores[i] == 0 ||
            m_pMergedDetWindows[i].width <= 0 || m_pMergedDetWindows[i].height <= 0 )
            removed_flags[i] = true;
    }
    for ( int i = 0 ; i < m_nMergedDetNum-1 ; i++ ) {
        if ( removed_flags[i] )
            continue;
        for ( int j = i+1 ; j < m_nMergedDetNum ; j++ ) {
            if ( removed_flags[j] )
                continue;
            if ( GetFMeasure(m_pMergedDetWindows[i], m_pMergedDetWindows[j])
                    > m_dbOverlappingThreshold )
                removed_flags[j] = true;
        }
    }

    int* det_front_indices = new int[m_nMergedDetNum];
    int* det_front_y = new int[m_nMergedDetNum];
    for ( int i = 0 ; i < m_nMergedDetNum ; i++ ) {
        det_front_y[i] = m_pMergedDetWindows[i].y;
        det_front_indices[i] = i;
    }

    CQuickSort<int, int> QS_int2;
    QS_int2.QuickSort(det_front_y, 0, m_nMergedDetNum-1, false, det_front_indices);

    for ( int ii = 0 ; ii < m_nMergedDetNum-1 ; ii++ ) {
        int i = det_front_indices[ii];
        if ( removed_flags[i] )
            continue;
        bool unoccluded = true;
        for ( int jj = 0 ; jj < ii ; jj++ ) {
            int j = det_front_indices[jj];
            if ( GetFMeasure(m_pMergedDetWindows[i], m_pMergedDetWindows[j]) > 0.1 ) {
                unoccluded = false;
                break;
            }
        }
        if ( !unoccluded )
            continue;
        for ( int jj = ii+1 ; jj < m_nMergedDetNum ; jj++ ) {
            int j = det_front_indices[jj];
            if ( removed_flags[j] )
                continue;
            if ( m_pMergedDetWindows[i].y - m_pMergedDetWindows[j].y >
                    MAX(m_pMergedDetWindows[i].width/8,5) ) {
                if ( GetFMeasure(m_pMergedDetWindows[i], m_pMergedDetWindows[j])
                    > m_dbFrontOverlappingThreshold &&
                    m_pMergedDetScores[i]*m_dbFrontScoreLargeScale < m_pMergedDetScores[j] ) {
                        removed_flags[i] = true;
                    break;
                }
            }
        }
    }

    int new_merged_det_num = 0;
    for ( int i = 0 ; i < m_nMergedDetNum ; i++ ) {
        if ( removed_flags[i] )
            continue;
        m_pMergedDetScores[new_merged_det_num] = m_pMergedDetScores[i];
        m_pMergedNeigDetNums[new_merged_det_num] = m_pMergedNeigDetNums[i];
        m_pMergedDetWindows[new_merged_det_num] = m_pMergedDetWindows[i];
        m_pMergedDetCenterPoints[new_merged_det_num] = m_pMergedDetCenterPoints[i];
        new_merged_det_num++;
    }
    m_nMergedDetNum = new_merged_det_num;

    delete [] removed_flags;
    delete [] scores;
    delete [] det_front_indices;
    delete [] det_front_y;
}

/* compute the F-measure between two rectangles */
double CMergingDetection::GetFMeasure(CvRect rect1, CvRect rect2) {
    CRectangle RECT;
    CvRect inter_rect = RECT.Intersection(rect1, rect2);
    double measure = 2.0*GetArea(inter_rect)/(GetArea(rect1)+GetArea(rect2));
    if ( GetArea(inter_rect) / GetArea(rect1) > 0.98 ||
        GetArea(inter_rect) / GetArea(rect2) > 0.98 )
        measure = 1.0;
    return measure;
}

double CMergingDetection::GetArea(CvRect rect) {
    return (double)(rect.width*rect.height);
}

