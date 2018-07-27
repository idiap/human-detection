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
// NegExamplesGenerator.h: interface for the CNegExamplesGenerator class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_NEG_EXAMPLES_GENERATOR_H_)
#define _NEG_EXAMPLES_GENERATOR_H_

#include "stdio.h"
#include <stdlib.h>
#include "BoostClassifierCovariance.h"
#include "FileList.h"
#include "Timer.h"
#include "highgui.h"
#include "TestWindows.h"
#include "Rectangle.h"
#include "ImageBoostingNegExamples.h"

class CNegExamplesGenerator
{
public:
	CNegExamplesGenerator(char *neg_sample_sg_images_list_fn, char *neg_sample_fp_sg_images_list_fn,
			   char *neg_sample_cp_images_list_fn, char *neg_sample_fp_cp_images_list_fn, char *neg_test_windows_list_fn,
			   char *neg_sample_noncp_images_list_fn, char *neg_sample_fp_noncp_images_list_fn,
               char *pos_sample_images_list_fn, char *pos_sample_fp_images_list_fn,
			   double *height_width_ratios, int min_det_win_width, int max_det_win_width,
			   int max_neg_examples_num, int max_2_neg_examples_num, bool neg_examples_from_pos=false,
			   int flipping_prev_cropped_neg_start_level=0,
			   int flipping_neg_img_start_level=1000,
			   int neg_boosting_cascade_step=0,
			   double neg_det_Xn_prob_alpha=1.0,
			   int end_noncp_cascade_level=1000);

    CNegExamplesGenerator();

	virtual ~CNegExamplesGenerator();

    /* set the input/output information for cropped negative examples */
    void SetNegCroppedExampleList(char* output_dir);

    /* generate a list of cropped negative examples */
	void GenerateNegExamples(char* boost_classifier_fn, char* log_msg_fn=NULL);
	void GenerateNegExamples(CBoostClassifierCovariance* boost_classifier, char* log_msg_fn=NULL);

	/* compute integral region covariances for inputed image(s) */
	CIntegralRegionCov* GetIntegralCovariance(IplImage* image,
            int* feature_types, int n_feature_types,
            IplImage* fp_image=NULL);
	CIntegralRegionCov* GetIntegralCovariance(char* image_fn, IplImage** image,
            int* feature_types, int n_feature_types,
            char* fp_image_fn=NULL, IplImage** fp_image=NULL);
	CIntegralRegionCov* GetIntegralCovariance(char* image_fn,
            int* feature_types, int n_feature_types,
            char* fp_image_fn=NULL);

    /* export log message to a text file */
    void ExportLogMessage(char *log_file_name, char *msg);

    /* positive image lists */
	CFileList *POS_OLIST;
	CFileList *POS_FLIST;

    /* file lists for cropped negative examples: color and foreground probability */
	CFileList *m_pNegCroppedExColorList;
	CFileList *m_pNegCroppedExFgList;

    /* output directory for negative cropped examples */
    char* m_pOutputDir;

    /* save the image with roi */
    void SaveImageROI(IplImage* img, CvRect roi, char* file_name);

    /* generate file name */
	char* GenerateFileName(const char* data_dir, const char* file_prefix, const char* file_ext, int idx);

    /* maximal negative example number */
	int m_nMaxNegExamplesNum;
    int m_nMax2NegExamplesNum;

    /* the balancing alpha using negative rejection threshold for boosting generation */
    double m_dbNegDetXnProbAlpha;

private:
    int m_nEndNoncpCascadeLevel;

    /* the class instant of negative test windows */
    CNegTestWindows* NTW;

    CRectangle m_cRect;

    /* store currently found negative cropped image lists */
    void SaveNegCroppedImageLists(int n_examples, bool finish=false);

    /* get test windows */
    CTestWindows* GetTestWindows(CvSize img_size);

	/* clean data */
	void CleanData();

    /* generate negative cropped examples */
	void GenerateNegExamples_Prev(CBoostClassifierCovariance* boost_classifier,
            int &start_example_loc, int end_example_loc,
            char* log_msg_fn=NULL);
	void GenerateNegExamples_POS(CBoostClassifierCovariance* boost_classifier,
            int &start_example_loc, int end_example_loc,
            char* log_msg_fn=NULL);
	void GenerateNegExamples_NEG_SG(CBoostClassifierCovariance* boost_classifier,
            int &start_example_loc, int end_example_loc,
            char* log_msg_fn=NULL);
	void GenerateNegExamples_NEG_CP(CBoostClassifierCovariance* boost_classifier,
            int &start_example_loc, int end_example_loc,
            char* log_msg_fn=NULL);
	void GenerateNegExamples_NEG_NONCP(CBoostClassifierCovariance* boost_classifier,
            int &start_example_loc, int end_example_loc,
            char* log_msg_fn=NULL);

    void TestNegExamples_NEG_CP(CBoostClassifierCovariance* boost_classifier,
            int &start_example_loc, int end_example_loc, char* log_msg_fn);

	/* save negative cropped example */
	void SaveNegCroppedExample(IplImage* image, CvRect roi, int example_idx,
            int edge_margin_size, int down_scale_factor, IplImage* fp_image=NULL);

    /* update the rectangle (sub-window) of interest */
    CvRect UpdateRectangle(CvSize img_size, CvRect roi, int &edge_margin_size);

    /* generate random negative example window */
    floatRect GetRandNegativeSgExample(CvSize img_size);
    floatRect GetRandNegativeExample(CvSize img_size);
    floatRect GetRandNegativePosExample(CvSize img_size);

    /* Torch random generator */
    Torch::Random *m_torchRNG;

    /* timer class instant */
    CTimer* m_pTimer;

	/* file lists */
	CFileList *NEG_SG_OLIST;
	CFileList *NEG_CP_OLIST;
	CFileList *NEG_NONCP_OLIST;

	CFileList *NEG_SG_FLIST;
	CFileList *NEG_CP_FLIST;
	CFileList *NEG_NONCP_FLIST;

    /* minimal and maximal height/width ratios to be used to generate negative window */
	double m_pHeightWidthRatios[2];

	/* minimal and maximal detection window width */
	int m_nMinDetWinWidth;
	int m_nMaxDetWinWidth;

    /* do we use positive examples to generate negative examples */
    bool m_bNegExamplesFromPos;

    /* starting cascade level for flipping negative images */
    int m_nFlippingPrevCroppedNegStartLevel;
    int m_nFlippingNegImgStartLevel;

    /* starting cascade level where not enough negative are generated */
    int m_nNotEnoughNegCascadeStartLevel;

    /* the negative boosting cascade step
       while no enough negative examples are generated at last cascade levels */
    int m_nNegBoostingCascadeLevelStep;
};

#endif // !defined(_NEG_EXAMPLES_GENERATOR_H_)
