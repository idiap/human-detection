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
// ImageBoostingNegExampless.h: interface for the CImageBoostingNegExamples class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_IMAGE_BOOSTING_NEG_EXAMPLES_H_)
#define _IMAGE_BOOSTING_NEG_EXAMPLES_H_

#include "BoostClassifierCovariance.h"
#include "TestWindows.h"
#include "FileList.h"
#include "cv.h"


class CImageBoostingNegExamples
{
public:
    /* export the detection results */
    bool Export(char* file_name);

    /* import the detection results */
    bool Import(char* file_name, bool load_all_info=true);

    /* generate initial test windows */
    void GenerateInitTestWindows(char* image_fn, char* fp_image_fn,
        double* img_scale_factors, int n_scale_factors,
        int min_test_win_width, int max_test_win_width,
        int test_win_marg_space,
        CTestWindows* TW=NULL);

    /* filtering the test windows using current learned model */
    void Test(CBoostClassifierCovariance *BCC,
        int last_cascade_level=-1,
        double Xn_prob_alpha=1.0,
        char* image_fn=NULL, char* fp_image_fn=NULL);

    /* randomly get the test windows */
    int GetTestWindows(CvRect* test_windows, int n_test_windows, int *img_scale_indices);

    /* filter the test windows */
    void FilterTestWindows(int image_marg_space);

    /* constructor */
	CImageBoostingNegExamples();

    /* destructor */
	virtual ~CImageBoostingNegExamples();

    /* clean the data */
    void CleanData();

    /* scaled images */
    IplImage**  m_ppScaledImages;
    IplImage**  m_ppScaledFpImages;

    /* numbers of test windows at the generation step and current test step */
    int  m_nTotOriginalTestWindowsNum;
    int  m_nCurTestWindowsNum;

    /* the last cascade level numbers of the detection windows */
    int     m_nLastCascadeLevelNo;

private:
    /* Torch random generator */
    Torch::Random *m_torchRNG;

    /* the image file name */
    char*    m_pImageFileName;
    char*    m_pFpImageFileName;

    /* the used image scale factors */
    CvSize*  m_pScaledImagesSizes;
    double*  m_pImageScalingFactors;
    int      m_nScaledImagesNum;

    /* the test windows */
    CvRect** m_ppTestWindows;

    int*     m_pTestWindowsNums;
};

///////////////////////////////////////////////////////////////////

class CNegTestWindows
{
public:
    CNegTestWindows(char* test_wins_fn_list) {
        m_pTestWindowsFileList = new CFileList(test_wins_fn_list);
        if ( m_pTestWindowsFileList->GetListLength() == 0 )
            m_nLastCascadeLevelNo = 0;
        else {
            CImageBoostingNegExamples* IBNE = new CImageBoostingNegExamples();
            IBNE->Import(m_pTestWindowsFileList->GetFileName(0));
            m_nLastCascadeLevelNo = IBNE->m_nLastCascadeLevelNo;
            delete IBNE;

            m_pOrgTestWindowsNums = new int[m_pTestWindowsFileList->GetListLength()];
            m_pCurTestWindowsNums = new int[m_pTestWindowsFileList->GetListLength()];
            for ( int i = 0 ; i < m_pTestWindowsFileList->GetListLength() ; i++ ) {
                CImageBoostingNegExamples* IBNE = new CImageBoostingNegExamples();
                IBNE->Import(m_pTestWindowsFileList->GetFileName(i), false);
                m_pOrgTestWindowsNums[i] = IBNE->m_nTotOriginalTestWindowsNum;
                m_pCurTestWindowsNums[i] = IBNE->m_nCurTestWindowsNum;
                delete IBNE;
            }
        }
    };

	virtual ~CNegTestWindows() {
	    if ( m_pTestWindowsFileList );
            delete m_pTestWindowsFileList;
        if ( m_pOrgTestWindowsNums )
            delete [] m_pOrgTestWindowsNums;
        if ( m_pCurTestWindowsNums )
            delete [] m_pCurTestWindowsNums;
	};

    void SetNeededNegExamples(int n_needed_neg_examples) {
        m_dbSelectPercent = 0;
        for ( int i = 0 ; i < m_pTestWindowsFileList->GetListLength() ; i++ )
            m_dbSelectPercent += (double)(m_pCurTestWindowsNums[i]);
        if ( m_dbSelectPercent > 0 )
            m_dbSelectPercent = (double)n_needed_neg_examples/m_dbSelectPercent;
    };

    double GetFalsePositiveRate() {
        long tot_org_test_windows_num = 0;
        long tot_cur_test_windows_num = 0;
        for ( int i = 0 ; i < m_pTestWindowsFileList->GetListLength() ; i++ ) {
            tot_cur_test_windows_num += m_pCurTestWindowsNums[i];
            tot_org_test_windows_num += m_pOrgTestWindowsNums[i];
        }
        return (double)tot_cur_test_windows_num/(double)tot_org_test_windows_num;
    };

    int GetMaxSelectExamplesNum(int frame_index) {
        if ( m_pCurTestWindowsNums[frame_index] == 0 )
            return 0;
        return cvCeil((double)m_pCurTestWindowsNums[frame_index]*m_dbSelectPercent);
    };

    CFileList* m_pTestWindowsFileList;
    int    m_nLastCascadeLevelNo;
    int*   m_pOrgTestWindowsNums;
    int*   m_pCurTestWindowsNums;
    double m_dbSelectPercent;
};


#endif // !defined(_IMAGE_BOOSTING_NEG_EXAMPLES_H_)
