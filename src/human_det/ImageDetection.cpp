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
// ImageDetection.cpp: implementation for the CImageDetection class.
//
//////////////////////////////////////////////////////////////////////

#include "ImageDetection.h"

/* human detection on a single image */
void CImageDetection::Detection(CBoostClassifierCovariance* bcc, CTestWindows* tw,
        IplImage* image, IplImage* fp_image,
        int evualation_start_level,
        int n_cascade_levels, double Xn_prob_alpha) {

    m_nCascadeLevelsNum = bcc->m_nCascadeLevelNum;

    m_szImage = cvGetSize(image);

    int final_pos_level = -1;
    double* cascade_level_probs = NULL;
    if ( evualation_start_level > 0 ) /* for detection evualation */
        cascade_level_probs = new double[bcc->m_nCascadeLevelNum];

    for ( int i = 0 ; i < bcc->m_nCascadeLevelNum ; i++ )
        m_pCascadeDetNums[i] = m_pCascadeDetMemNums[i] = 0;

    /* get the integral region covariance for the given image */
    CNegExamplesGenerator* neg = new CNegExamplesGenerator();
    CIntegralRegionCov* irc = neg->GetIntegralCovariance(
            image,
            bcc->m_pFeatureImgTypes, bcc->m_nFeatureTypeNum,
            fp_image);

    /* get the testing windows */
    //tw->GetTestWindows(cvRect(1,1,irc->m_szFeatureImage.width-2,irc->m_szFeatureImage.height-2));

    /* doing detection */
    for ( int i = 0 ; i < tw->m_nTestWindowsNum ; i++ ) {
        double score;
        if ( evualation_start_level > 0 ) {
            score = bcc->Predictor(irc, tw->m_pTestWindows[i], n_cascade_levels, Xn_prob_alpha,
                            &final_pos_level, cascade_level_probs);
            for ( int j = evualation_start_level ; j <= final_pos_level ; j++ ) {
                UpdateData(&(m_ppDetScores[j]), &(m_ppDetWindows[j]), m_pCascadeDetNums[j], m_pCascadeDetMemNums[j]);
                m_ppDetScores[j][m_pCascadeDetNums[j]] = cascade_level_probs[j];
                m_ppDetWindows[j][m_pCascadeDetNums[j]] = tw->m_pTestWindows[i];
                m_pCascadeDetNums[j]++;
            }
        }
        else {
            score = bcc->Predictor(irc, tw->m_pTestWindows[i], n_cascade_levels, Xn_prob_alpha);
            if ( score > 0 ) {
                UpdateData(&(m_ppDetScores[0]), &(m_ppDetWindows[0]), m_pCascadeDetNums[0], m_pCascadeDetMemNums[0]);
                m_ppDetScores[0][m_pCascadeDetNums[0]] = score;
                m_ppDetWindows[0][m_pCascadeDetNums[0]] = tw->m_pTestWindows[i];
                m_pCascadeDetNums[0]++;
            }
        }
    }

    /* release memories */
    delete neg;
    delete irc;
    if ( cascade_level_probs )
        delete [] cascade_level_probs;
}

/* draw the merged outputs on the image */
void CImageDetection::Draw(IplImage* image,
        bool show_windows,
        int cascade_level,
        int radius,
        int thinness,
        CvScalar center_color,
        CvScalar win_color) {
    CvPoint pt1, pt2, pt;
    for ( int i = 0 ; i < m_pCascadeDetNums[cascade_level] ; i++ ) {
        CvRect rect = m_ppDetWindows[cascade_level][i];
        pt1 = cvPoint(rect.x, rect.y);
        pt2 = cvPoint(pt1.x+rect.width-1, pt1.y+rect.height-1);
        pt = cvPoint((pt1.x+pt2.x)/2, (pt1.y+pt2.y)/2);
        cvCircle( image, pt, radius, center_color, thinness);
        if ( show_windows )
            cvRectangle(image, pt1, pt2, win_color, thinness);
    }
}

/* export the detection results */
bool CImageDetection::Export(char* file_name) {
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
    cvWriteComment( fs, "Human detection candidates on a single image", 0 );
    cvStartWriteStruct( fs, "image_human_detection", CV_NODE_MAP, NULL, cvAttrList(0,0));

    cvWriteComment( fs, "Size of the input image", 0 );
    cvWriteInt( fs, "image_width", m_szImage.width );
    cvWriteInt( fs, "image_height", m_szImage.height );

    cvWriteComment( fs, "Number of detection candidates", 0 );
    cvWriteInt( fs, "n_det_windows", m_pCascadeDetNums[0] );

    for ( int i = 0 ; i < m_pCascadeDetNums[0] ; i++ ) {
        sprintf(label_name, "detection_%d", i);
        cvStartWriteStruct( fs, label_name, CV_NODE_MAP, NULL, cvAttrList(0,0));

        cvWriteComment( fs, "Positive detection score", 0 );
        cvWriteReal( fs, "det_score", m_ppDetScores[0][i] );

        cvWriteComment( fs, "Positive sub-window", 0 );
        cvWriteInt( fs, "roi_x", m_ppDetWindows[0][i].x );
        cvWriteInt( fs, "roi_y", m_ppDetWindows[0][i].y );
        cvWriteInt( fs, "roi_width", m_ppDetWindows[0][i].width );
        cvWriteInt( fs, "roi_height", m_ppDetWindows[0][i].height );

        cvEndWriteStruct( fs );
    }


    cvEndWriteStruct( fs );

    cvReleaseFileStorage( &fs );

    return true;

}

/* import the detection results */
bool CImageDetection::Import(char* file_name) {
    if ( !file_name ) {
        printf("Please provide file name (XML or YAML format) to import the detection results!\n");
        return false;
    }

    CvFileStorage* fs = cvOpenFileStorage( file_name, 0, CV_STORAGE_READ );
    if ( !fs ) {
        printf("Open file error : %s\n", file_name);
        return false;
    }

    char label_name[1024];

    ///////////////////////////////////////////////////////////////////////////////
    CvFileNode* image_human_det_node = cvGetFileNodeByName( fs, NULL, "image_human_detection" );

    m_szImage.width = cvReadIntByName( fs, image_human_det_node, "image_width" );
    m_szImage.height = cvReadIntByName( fs, image_human_det_node, "image_height" );

    m_pCascadeDetNums[0] = cvReadIntByName( fs, image_human_det_node, "n_det_windows" );
    UpdateData(&(m_ppDetScores[0]), &(m_ppDetWindows[0]), m_pCascadeDetNums[0], m_pCascadeDetMemNums[0]);

    for ( int i = 0 ; i < m_pCascadeDetNums[0] ; i++ ) {
        sprintf(label_name, "detection_%d", i);
        CvFileNode* det_node = cvGetFileNodeByName( fs, image_human_det_node, label_name );

        m_ppDetScores[0][i] = cvReadRealByName( fs, det_node, "det_score" );

        m_ppDetWindows[0][i].x = cvReadIntByName( fs, det_node, "roi_x" );
        m_ppDetWindows[0][i].y = cvReadIntByName( fs, det_node, "roi_y" );
        m_ppDetWindows[0][i].width = cvReadIntByName( fs, det_node, "roi_width" );
        m_ppDetWindows[0][i].height = cvReadIntByName( fs, det_node, "roi_height" );
    }

    cvReleaseFileStorage( &fs );

    return true;

}

/* return the found positive detection scores */
double* CImageDetection::GetPosDetScores(int cascade_level) {
    if ( cascade_level < 0 || cascade_level >= m_nCascadeLevelsNum )
        return NULL;
    return m_ppDetScores[cascade_level];
}

/* return the found positive detection windows (rectangles) */
CvRect* CImageDetection::GetPosDetWindows(int cascade_level) {
    if ( cascade_level < 0 || cascade_level >= m_nCascadeLevelsNum )
        return NULL;
    return m_ppDetWindows[cascade_level];
}

/* return the number of found positive detection windows */
int CImageDetection::GetPosDetWindowNum(int cascade_level) {
    if ( cascade_level < 0 || cascade_level >= m_nCascadeLevelsNum )
        return 0;
    return m_pCascadeDetNums[cascade_level];
};

void CImageDetection::UpdateDetections(int x_marg_size, int y_marg_size,
    double human_left_percent, double human_top_percent, int cascade_level) {
    if ( x_marg_size == 0 && y_marg_size == 0
        && human_left_percent == 0 && human_top_percent == 0 )
        return;
    if ( cascade_level >= m_nCascadeLevelsNum )
        return;
    if ( cascade_level < 0 ) {
        for ( int i = 0 ; i < m_nCascadeLevelsNum ; i++ ) {
            for ( int j = 0 ; j < m_pCascadeDetNums[i] ; j++ ) {
                m_ppDetWindows[i][j].x -= x_marg_size;
                m_ppDetWindows[i][j].y -= y_marg_size;

                int x_offset = (int)((double)m_ppDetWindows[i][j].width*human_left_percent);
                int y_offset = (int)((double)m_ppDetWindows[i][j].height*human_top_percent);

                m_ppDetWindows[i][j].x += x_offset;
                m_ppDetWindows[i][j].y += y_offset;
                m_ppDetWindows[i][j].width -= 2*x_offset;
                m_ppDetWindows[i][j].height -= 2*y_offset;
            }
        }
    }
    else {
        int i = cascade_level;
        for ( int j = 0 ; j < m_pCascadeDetNums[i] ; j++ ) {
            m_ppDetWindows[i][j].x -= x_marg_size;
            m_ppDetWindows[i][j].y -= y_marg_size;

            int x_offset = (int)((double)m_ppDetWindows[i][j].width*human_left_percent);
            int y_offset = (int)((double)m_ppDetWindows[i][j].height*human_top_percent);

            m_ppDetWindows[i][j].x += x_offset;
            m_ppDetWindows[i][j].y += y_offset;
            m_ppDetWindows[i][j].width -= 2*x_offset;
            m_ppDetWindows[i][j].height -= 2*y_offset;
        }
    }
}

/* constructor */
CImageDetection::CImageDetection() {
    m_nMaxMemCascadeLevelsNum = 100;
    m_ppDetScores = new double*[m_nMaxMemCascadeLevelsNum];
    m_ppDetWindows = new CvRect*[m_nMaxMemCascadeLevelsNum];
    m_pCascadeDetNums = new int[m_nMaxMemCascadeLevelsNum];
    m_pCascadeDetMemNums = new int[m_nMaxMemCascadeLevelsNum];
    for ( int i = 0 ; i < m_nMaxMemCascadeLevelsNum ; i++ ) {
        m_ppDetScores[i] = NULL;
        m_ppDetWindows[i] = NULL;
        m_pCascadeDetNums[i] = 0;
        m_pCascadeDetMemNums[i] = 0;
    }
    m_nCascadeLevelsNum = 0;
}

/* destructor */
CImageDetection::~CImageDetection() {
    for ( int i = 0 ; i < m_nMaxMemCascadeLevelsNum ; i++ ) {
        if ( m_ppDetScores[i] )
            delete [] m_ppDetScores[i];
        if ( m_ppDetWindows[i] )
            delete [] m_ppDetWindows[i];
    }
    delete [] m_pCascadeDetNums;
    delete [] m_pCascadeDetMemNums;
    delete [] m_ppDetScores;
    delete [] m_ppDetWindows;
}

/* update the data */
void CImageDetection::UpdateData(double** scores, CvRect** windows, int cur_len, int &new_len) {
    if ( cur_len < new_len )
        return;
    if ( cur_len == new_len )
        new_len = cur_len + 10000;
    else
        new_len = cur_len;
    double* new_scores = new double[new_len];
    CvRect* new_windows = new CvRect[new_len];
    if ( (*scores) ) {
        for ( int i = 0 ; i < cur_len ; i++ ) {
            new_scores[i] = (*scores)[i];
            new_windows[i] = (*windows)[i];
        }
    }
    if ( *scores )
        delete [] *scores;
    if ( *windows )
        delete [] *windows;
    *scores = new_scores;
    *windows = new_windows;
}

void CImageDetection::InsertImageMarginSpace(
        IplImage* image, IplImage** marg_image,
        int left_marg_size, int right_marg_size,
        int top_marg_size, int bottom_marg_size)
{
    if ( *marg_image == NULL ) {
        *marg_image = cvCreateImage(cvSize(image->width+left_marg_size+right_marg_size,
                        image->height+top_marg_size+bottom_marg_size),
            image->depth, image->nChannels);
    }

    IplImage* mimg = *marg_image;

    CvRect roi = cvRect(left_marg_size, top_marg_size, image->width, image->height);

    cvSetImageROI(mimg, roi);
    cvCopy(image, mimg);
    cvResetImageROI(mimg);

    /* repeat the left and right edge */
    for ( int y = 0 ; y < image->height ; y++ ) {
        int px = 0;
        int py = y+top_marg_size;
        uchar* marg_data = &((uchar*)(mimg->imageData + mimg->widthStep*py))[px*mimg->nChannels];
        uchar* edge_data = &((uchar*)(mimg->imageData + mimg->widthStep*py))[(px+left_marg_size)*mimg->nChannels];
        for ( int x = 0 ; x < left_marg_size ; x++ ) {
            for ( int i = 0 ; i < mimg->nChannels ; i++ )
                *marg_data++ = edge_data[i];
        }

        px = left_marg_size+image->width;
        marg_data = &((uchar*)(mimg->imageData + mimg->widthStep*py))[px*mimg->nChannels];
        edge_data = &((uchar*)(mimg->imageData + mimg->widthStep*py))[(px-1)*mimg->nChannels];
        for ( int x = 0 ; x < right_marg_size ; x++ ) {
            for ( int i = 0 ; i < mimg->nChannels ; i++ )
                *marg_data++ = edge_data[i];
        }
    }

    /* repeat the top edge */
    for ( int y = 0 ; y < top_marg_size ; y++ ) {
        int px = 0;
        int py = y;
        uchar* marg_data = &((uchar*)(mimg->imageData + mimg->widthStep*py))[px*mimg->nChannels];
        uchar* edge_data = &((uchar*)(mimg->imageData + mimg->widthStep*top_marg_size))[px*mimg->nChannels];
        for ( int x = 0 ; x < mimg->width ; x++ ) {
            for ( int i = 0 ; i < mimg->nChannels ; i++ )
                *marg_data++ = *edge_data++;
        }
    }

    /* repeat the bottom edge */
    for ( int y = 0 ; y < bottom_marg_size ; y++ ) {
        int px = 0;
        int py = y + top_marg_size + image->height;
        uchar* marg_data = &((uchar*)(mimg->imageData + mimg->widthStep*py))[px*mimg->nChannels];
        uchar* edge_data = &((uchar*)(mimg->imageData + mimg->widthStep*(top_marg_size+image->height-1)))[px*mimg->nChannels];
        for ( int x = 0 ; x < mimg->width ; x++ ) {
            for ( int i = 0 ; i < mimg->nChannels ; i++ )
                *marg_data++ = *edge_data++;
        }
    }

}
