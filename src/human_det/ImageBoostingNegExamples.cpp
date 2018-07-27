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
// ImageBoostingNegExamples.cpp: implementation for the CImageBoostingNegExamples class.
//
//////////////////////////////////////////////////////////////////////

#include "ImageBoostingNegExamples.h"

#include "NegExamplesGenerator.h"


/* constructor */
CImageBoostingNegExamples::CImageBoostingNegExamples() {
    m_pImageFileName = NULL;
    m_pFpImageFileName = NULL;


    m_pScaledImagesSizes = NULL;
    m_pImageScalingFactors = NULL;
    m_nScaledImagesNum = 0;

    m_pTestWindowsNums = NULL;
    m_nLastCascadeLevelNo = 0;

    m_ppTestWindows = NULL;

    m_ppScaledImages = NULL;
    m_ppScaledFpImages = NULL;

    m_nTotOriginalTestWindowsNum = 0;

    m_torchRNG = new Torch::Random;
}

/* destructor */
CImageBoostingNegExamples::~CImageBoostingNegExamples() {
    CleanData();
    if ( m_torchRNG )
        delete m_torchRNG;
}

void CImageBoostingNegExamples::CleanData() {
    if ( m_pImageFileName )
        delete [] m_pImageFileName;
    if ( m_pFpImageFileName )
        delete [] m_pFpImageFileName;
    if ( m_pScaledImagesSizes )
        delete [] m_pScaledImagesSizes;
    if ( m_pImageScalingFactors )
        delete [] m_pImageScalingFactors;
    if ( m_ppTestWindows ) {
        for ( int i = 0 ; i < m_nScaledImagesNum ; i++ ) {
            if ( m_ppTestWindows[i] )
                delete [] m_ppTestWindows[i];
        }
        delete [] m_ppTestWindows;
    }
    if ( m_pTestWindowsNums ) {
        delete [] m_pTestWindowsNums;
    }
    if ( m_ppScaledImages ) {
        for ( int i = 0 ; i < m_nScaledImagesNum ; i++ )
            if ( m_ppScaledImages[i] )
                cvReleaseImage(&m_ppScaledImages[i]);
    }
    if ( m_ppScaledFpImages ) {
        for ( int i = 0 ; i < m_nScaledImagesNum ; i++ )
            if ( m_ppScaledFpImages[i] )
                cvReleaseImage(&m_ppScaledFpImages[i]);
    }

    m_pImageFileName = NULL;
    m_pFpImageFileName = NULL;


    m_pScaledImagesSizes = NULL;
    m_pImageScalingFactors = NULL;
    m_nScaledImagesNum = 0;

    m_nTotOriginalTestWindowsNum = 0;

    m_pTestWindowsNums = NULL;
    m_nLastCascadeLevelNo = 0;

    m_ppTestWindows = NULL;

    m_ppScaledImages = NULL;
    m_ppScaledFpImages = NULL;
}

/* export the detection results */
bool CImageBoostingNegExamples::Export(char* file_name) {
    if ( !file_name ) {
        printf("Please provide file name (XML or YAML format) to export the detection results!\n");
        return false;
    }

    CvFileStorage* fs = cvOpenFileStorage( file_name, 0, CV_STORAGE_WRITE );
    if ( !fs ) {
        printf("Open file error : %s\n", file_name);
        return false;
    }

    //char label_name[1024];
    char label_name[1024];

    ///////////////////////////////////////////////////////////////////////////////
    cvWriteComment( fs, "Negative test windows on a single image", 0 );

    cvWriteComment( fs, "Image file name", 0 );
	cvWriteString( fs, "image_file_name", m_pImageFileName );

	if ( m_pFpImageFileName ) {
        cvWriteComment( fs, "Foreground probability image file name", 0 );
        cvWriteString( fs, "fp_image_file_name", m_pFpImageFileName );
	}

    cvWriteComment( fs, "Total number of original test windows", 0 );
    cvWriteInt( fs, "tot_org_test_windows_num", m_nTotOriginalTestWindowsNum );

    cvWriteComment( fs, "Total number of current test windows", 0 );
    cvWriteInt( fs, "tot_cur_test_windows_num", m_nCurTestWindowsNum );

    cvWriteComment( fs, "Last tested cascade level", 0 );
    cvWriteInt( fs, "last_tested_cascade_level", m_nLastCascadeLevelNo );

    cvWriteComment( fs, "Number of scaled images", 0 );
    cvWriteInt( fs, "scaled_image_num", m_nScaledImagesNum );

    for ( int i = 0 ; i < m_nScaledImagesNum ; i++ ) {
        sprintf(label_name, "scaled_image_%d", i);
        cvStartWriteStruct( fs, label_name, CV_NODE_MAP, NULL, cvAttrList(0,0));

        cvWriteComment( fs, "Image scaling factor", 0 );
        cvWriteReal( fs, "image_scaling_factor", m_pImageScalingFactors[i] );

        cvWriteComment( fs, "Scaled image size", 0 );
        cvWriteInt( fs, "scaled_image_width", m_pScaledImagesSizes[i].width );
        cvWriteInt( fs, "scaled_image_height", m_pScaledImagesSizes[i].height );

        cvWriteComment( fs, "Number of negative test windows", 0 );
        cvWriteInt( fs, "neg_test_windows_num", m_pTestWindowsNums[i] );

        if ( m_pTestWindowsNums[i] > 0 ) {
            CvMat* tmp_mat = cvCreateMat(m_pTestWindowsNums[i], 4, CV_16UC1);
            short* ptr = (short*)(tmp_mat->data.s);
            for ( int j = 0 ; j < m_pTestWindowsNums[i] ; j++ ) {
                *ptr++ = (short)(m_ppTestWindows[i][j].x);
                *ptr++ = (short)(m_ppTestWindows[i][j].y);
                *ptr++ = (short)(m_ppTestWindows[i][j].width);
                *ptr++ = (short)(m_ppTestWindows[i][j].height);
            }
            cvWriteComment( fs, "Negative test windows information", 0 );
            cvWriteComment( fs, "[x  y  width height]", 0 );
            cvWrite( fs, "neg_test_windows", tmp_mat );
            cvReleaseMat(&tmp_mat);
        }

        cvEndWriteStruct( fs );
    }

    cvReleaseFileStorage( &fs );

    return true;
}

/* import the detection results */
bool CImageBoostingNegExamples::Import(char* file_name, bool load_all_info) {
    if ( !file_name ) {
        printf("Please provide file name (XML or YAML format) to import the detection results!\n");
        return false;
    }

    CvFileStorage* fs = cvOpenFileStorage( file_name, 0, CV_STORAGE_READ );
    if ( !fs ) {
        //printf("Open file error : %s\n", file_name);
        return false;
    }

    /* clean data */
    CleanData();

    char label_name[1024];

    ///////////////////////////////////////////////////////////////////////////////

    const char* img_file_name = cvReadStringByName( fs, NULL, "image_file_name" );
    const char* fp_img_file_name = cvReadStringByName( fs, NULL, "fp_image_file_name" );

    m_pImageFileName = new char[2048];
    sprintf(m_pImageFileName, "%s", img_file_name);
    if ( fp_img_file_name ) {
        m_pFpImageFileName = new char[2048];
        sprintf(m_pFpImageFileName, "%s", fp_img_file_name);
    }

    m_nTotOriginalTestWindowsNum = cvReadIntByName( fs, NULL, "tot_org_test_windows_num" );
    m_nCurTestWindowsNum = cvReadIntByName( fs, NULL, "tot_cur_test_windows_num" );

    if ( !load_all_info ) {
        cvReleaseFileStorage( &fs );
        return true;
    }

    m_nLastCascadeLevelNo = cvReadIntByName( fs, NULL, "last_tested_cascade_level" );

    m_nScaledImagesNum = cvReadIntByName( fs, NULL, "scaled_image_num" );

    m_pScaledImagesSizes = new CvSize[m_nScaledImagesNum];
    m_pImageScalingFactors = new double[m_nScaledImagesNum];
    m_ppTestWindows = new CvRect*[m_nScaledImagesNum];
    m_pTestWindowsNums = new int[m_nScaledImagesNum];

    for ( int i = 0 ; i < m_nScaledImagesNum ; i++ ) {
        sprintf(label_name, "scaled_image_%d", i);
        CvFileNode* scaled_image_node = cvGetFileNodeByName( fs, NULL, label_name );

        m_pImageScalingFactors[i] = cvReadRealByName( fs, scaled_image_node, "image_scaling_factor" );

        m_pScaledImagesSizes[i].width = cvReadIntByName( fs, scaled_image_node, "scaled_image_width" );
        m_pScaledImagesSizes[i].height = cvReadIntByName( fs, scaled_image_node, "scaled_image_height" );

        m_pTestWindowsNums[i] = cvReadIntByName( fs, scaled_image_node, "neg_test_windows_num" );

        m_ppTestWindows[i] = NULL;
        if ( m_pTestWindowsNums[i] > 0 ) {
            m_ppTestWindows[i] = new CvRect[m_pTestWindowsNums[i]];
            CvMat* tmp_mat = (CvMat*)cvReadByName( fs, scaled_image_node, "neg_test_windows" );
            short* ptr = (short*)(tmp_mat->data.s);
            for ( int j = 0 ; j < m_pTestWindowsNums[i] ; j++ ) {
                m_ppTestWindows[i][j].x = *ptr++;
                m_ppTestWindows[i][j].y = *ptr++;
                m_ppTestWindows[i][j].width = *ptr++;
                m_ppTestWindows[i][j].height = *ptr++;
            }
            cvReleaseMat(&tmp_mat);
        }
    }

    cvReleaseFileStorage( &fs );

    return true;
}

void CImageBoostingNegExamples::FilterTestWindows(int image_marg_space) {
    m_nCurTestWindowsNum = 0;
    for ( int i = 0 ; i < m_nScaledImagesNum ; i++ ) {
        if ( m_pTestWindowsNums[i] > 0 ) {
            bool* valid_flags = new bool[m_pTestWindowsNums[i]];
            for ( int j = 0 ; j < m_pTestWindowsNums[i] ; j++ ) {
                valid_flags[j] = ( m_ppTestWindows[i][j].x >= image_marg_space &&
                    m_ppTestWindows[i][j].y >= image_marg_space &&
                    m_ppTestWindows[i][j].x + m_ppTestWindows[i][j].width <=
                        m_pScaledImagesSizes[i].width-2*image_marg_space &&
                    m_ppTestWindows[i][j].y + m_ppTestWindows[i][j].height <=
                        m_pScaledImagesSizes[i].height-2*image_marg_space );
            }
            int n_valid_wins_num = 0;
            for ( int j = 0 ; j < m_pTestWindowsNums[i] ; j++ )
                if ( valid_flags[j] )
                    m_ppTestWindows[i][n_valid_wins_num++] = m_ppTestWindows[i][j];
            m_pTestWindowsNums[i] = n_valid_wins_num;
            m_nCurTestWindowsNum += m_pTestWindowsNums[i];
        }
    }
}

void CImageBoostingNegExamples::Test(CBoostClassifierCovariance *BCC,
        int last_cascade_level, double Xn_prob_alpha,
        char* image_fn, char* fp_image_fn)
{
    if ( BCC->m_nCascadeLevelNum == 0 )
        return;

    char* _image_fn = image_fn;
    char* _fp_image_fn = fp_image_fn;
    if ( !_image_fn )
        _image_fn = m_pImageFileName;
    if ( !_fp_image_fn )
        _fp_image_fn = m_pFpImageFileName;

    /* load images */
    IplImage* org_img = cvLoadImage(_image_fn);
    if ( !org_img ) {
        printf("Load image error : %s\n", _image_fn);
        exit(1);
    }
    IplImage* org_fp_img = NULL;
    if ( _fp_image_fn ) {
        org_fp_img = cvLoadImage(_fp_image_fn);
        if ( !org_fp_img ) {
            printf("Load image error : %s\n", _fp_image_fn);
            exit(1);
        }
    }

    CNegExamplesGenerator* NEG = new CNegExamplesGenerator();

    m_ppScaledImages = new IplImage*[m_nScaledImagesNum];
    m_ppScaledFpImages = new IplImage*[m_nScaledImagesNum];
    for ( int i = 0 ; i < m_nScaledImagesNum ; i++ )
        m_ppScaledImages[i] = m_ppScaledFpImages[i] = NULL;

    int n_used_cascade_levels = last_cascade_level;
    if ( n_used_cascade_levels < 0 )
        n_used_cascade_levels = BCC->m_nCascadeLevelNum;
    n_used_cascade_levels = MIN(BCC->m_nCascadeLevelNum, n_used_cascade_levels);

    for ( int i = 0 ; i < m_nScaledImagesNum ; i++ ) {
        if ( m_pTestWindowsNums[i] > 0 ) {
            /* get the scaled images */
            m_ppScaledImages[i] = cvCreateImage(m_pScaledImagesSizes[i], org_img->depth, org_img->nChannels);
            cvResize(org_img, m_ppScaledImages[i]);
            if ( org_fp_img ) {
                m_ppScaledFpImages[i] = cvCreateImage(m_pScaledImagesSizes[i], org_fp_img->depth, org_fp_img->nChannels);
                cvResize(org_fp_img, m_ppScaledFpImages[i]);
            }

            /* get the integral region covariance */
            CIntegralRegionCov* IRC = NEG->GetIntegralCovariance(m_ppScaledImages[i],
                    BCC->m_pFeatureImgTypes, BCC->m_nFeatureTypeNum, m_ppScaledFpImages[i]);

            bool* test_flags = new bool[m_pTestWindowsNums[i]];
            int pos_test_num = 0;
            for ( int j = 0 ; j < m_pTestWindowsNums[i] ; j++ ) {
                test_flags[j] = BCC->CascadeLevelPredictor(IRC, m_ppTestWindows[i][j],
                        m_nLastCascadeLevelNo, n_used_cascade_levels, Xn_prob_alpha);
                if (test_flags[j] )
                    pos_test_num++;
            }

            delete IRC;

            if ( pos_test_num == 0 ) {
                delete [] m_ppTestWindows[i];
                m_ppTestWindows[i] = NULL;
            }
            else {
                CvRect* new_test_windows = new CvRect[pos_test_num];
                pos_test_num = 0;
                for ( int j = 0 ; j < m_pTestWindowsNums[i] ; j++ ) {
                    if ( test_flags[j] ) {
                        new_test_windows[pos_test_num++] = m_ppTestWindows[i][j];
                    }
                }
                delete [] m_ppTestWindows[i];
                m_ppTestWindows[i] = new_test_windows;
            }

            m_pTestWindowsNums[i] = pos_test_num;

            delete [] test_flags;

        }
    }

    m_nLastCascadeLevelNo = n_used_cascade_levels;

    /* release memories */
    delete NEG;
    cvReleaseImage(&org_img);
    if ( org_fp_img )
        cvReleaseImage(&org_fp_img);

    m_nCurTestWindowsNum = 0;
    for ( int i = 0 ; i < m_nScaledImagesNum ; i++ )
        m_nCurTestWindowsNum += m_pTestWindowsNums[i];

    return;
}

void CImageBoostingNegExamples::GenerateInitTestWindows(
        char* image_fn, char* fp_image_fn,
        double* img_scale_factors, int n_scale_factors,
        int min_test_win_width, int max_test_win_width,
        int test_win_marg_space, CTestWindows* TW) {

    CleanData();

    m_pImageFileName = new char[2048];

    sprintf(m_pImageFileName, "%s", image_fn);

    if ( fp_image_fn ) {
        m_pFpImageFileName = new char[2048];
        sprintf(m_pFpImageFileName, "%s", fp_image_fn);
    }

    m_pScaledImagesSizes = new CvSize[n_scale_factors];
    m_pImageScalingFactors = new double[n_scale_factors];
    m_nScaledImagesNum = n_scale_factors;

    IplImage* org_img = cvLoadImage(image_fn);
    if ( !org_img ) {
        printf("Load image error : %s!\n", image_fn);
        exit(1);
    }
    CvSize img_size = cvGetSize(org_img);
    for ( int i = 0 ; i < n_scale_factors ; i++ ) {
        m_pImageScalingFactors[i] = img_scale_factors[i];
        m_pScaledImagesSizes[i].width = cvCeil(img_scale_factors[i]*(double)img_size.width);
        m_pScaledImagesSizes[i].height = cvCeil(img_scale_factors[i]*(double)img_size.height);
    }
    cvReleaseImage(&org_img);

    m_ppTestWindows = new CvRect*[n_scale_factors];
    m_pTestWindowsNums = new int[n_scale_factors];

    m_nLastCascadeLevelNo = 0;

    int min_win_width = min_test_win_width;
    float w_scale_factor = 0.1f;
    float min_hw_ratio = 1.5f;
    float h_scale_factor = 0.1f;
    int n_hw_scales = 8;
    float w_shift_factor = 0.08f;
    float h_shift_factor = 0.08f;
    int n_win_width_scales = 100;

    bool input_TW = (TW != NULL);
    if ( !TW )
        TW = new CTestWindows(min_win_width,
            n_win_width_scales,
            w_scale_factor,
            min_hw_ratio,
            n_hw_scales,
            h_scale_factor,
            w_shift_factor,
            h_shift_factor);

    m_nTotOriginalTestWindowsNum = 0;
    for ( int i = 0 ; i < n_scale_factors ; i++ ) {
        CvRect image_roi = cvRect( test_win_marg_space+1, test_win_marg_space+1,
            m_pScaledImagesSizes[i].width-2*test_win_marg_space-2,
            m_pScaledImagesSizes[i].height-2*test_win_marg_space-2 );

        TW->GetTestWindows(image_roi);
        TW->FilterTestWindows(max_test_win_width);

        m_ppTestWindows[i] = NULL;

        m_pTestWindowsNums[i] = TW->m_nTestWindowsNum;

        m_nTotOriginalTestWindowsNum += TW->m_nTestWindowsNum;
        if ( m_pTestWindowsNums[i] > 0 ) {
            m_ppTestWindows[i] = new CvRect[m_pTestWindowsNums[i]];

            for ( int j = 0 ; j < m_pTestWindowsNums[i] ; j++ ) {
                m_ppTestWindows[i][j] = TW->m_pTestWindows[j];
            }
        }
    }

    if ( !input_TW )
        delete TW;

    m_nCurTestWindowsNum = m_nTotOriginalTestWindowsNum;

    /* save the test windows */
    //Export(test_wins_fn);
}

int CImageBoostingNegExamples::GetTestWindows(CvRect* test_windows, int n_test_windows, int *img_scale_indices)
{
    int* start_locations = new int[m_nScaledImagesNum];
    start_locations[0] = 0;
    int tot_n_test_windows = m_pTestWindowsNums[0];
    for ( int i = 1 ; i < m_nScaledImagesNum ; i++ ) {
        start_locations[i] = start_locations[i-1] + m_pTestWindowsNums[i-1];
        tot_n_test_windows += m_pTestWindowsNums[i];
    }

    int n_valid_test_windows = 0;
    if ( tot_n_test_windows <= n_test_windows ) {
        for ( int i = 0 ; i < m_nScaledImagesNum ; i++ ) {
            for ( int j = 0 ; j < m_pTestWindowsNums[i] ; j++ ) {
                img_scale_indices[n_valid_test_windows] = i;
                test_windows[n_valid_test_windows++] = m_ppTestWindows[i][j];
            }
        }
    }
    else {
        int* shuffledIndices = new int[tot_n_test_windows];
        m_torchRNG->getShuffledIndices(shuffledIndices, tot_n_test_windows);
        for ( int i = 0 ; i < n_test_windows ; i++ ) {
            for ( int j = m_nScaledImagesNum-1 ; j >= 0 ; j-- ) {
                int shift = shuffledIndices[i] - start_locations[j];
                if ( shift >= 0 ) {
                    img_scale_indices[n_valid_test_windows] = j;
                    test_windows[n_valid_test_windows++] = m_ppTestWindows[j][shift];
                    break;
                }
            }
        }
        delete [] shuffledIndices;
    }

    delete [] start_locations;

    return n_valid_test_windows;
}
