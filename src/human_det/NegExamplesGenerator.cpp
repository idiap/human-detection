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
// NegExamplesGenerator.cpp: implementation of the CNegExamplesGenerator class.
//
//////////////////////////////////////////////////////////////////////

#include "NegExamplesGenerator.h"
#include "QuickSort.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CNegExamplesGenerator::CNegExamplesGenerator(
            char *neg_example_sg_images_list_fn, char *neg_example_fp_sg_images_list_fn,
            char *neg_example_cp_images_list_fn, char *neg_example_fp_cp_images_list_fn, char *neg_test_windows_list_fn,
            char *neg_example_noncp_images_list_fn, char *neg_example_fp_noncp_images_list_fn,
            char *pos_sample_images_list_fn, char *pos_sample_fp_images_list_fn,
            double *height_width_ratios, int min_det_win_width, int max_det_win_width,
            int max_neg_examples_num, int max_2_neg_examples_num, bool neg_examples_from_pos,
            int flipping_prev_cropped_neg_start_level,
            int flipping_neg_img_start_level,
            int neg_boosting_cascade_step,
            double neg_det_Xn_prob_alpha,
            int end_noncp_cascade_level)
{
	/* load training file lists */
	NEG_SG_OLIST = new CFileList( neg_example_sg_images_list_fn );
	NEG_CP_OLIST = new CFileList( neg_example_cp_images_list_fn );
	NEG_NONCP_OLIST = new CFileList( neg_example_noncp_images_list_fn );

	NEG_SG_FLIST = new CFileList( neg_example_fp_sg_images_list_fn );
	NEG_CP_FLIST = new CFileList( neg_example_fp_cp_images_list_fn );
	NEG_NONCP_FLIST = new CFileList( neg_example_fp_noncp_images_list_fn );

	POS_OLIST = new CFileList(pos_sample_images_list_fn);
	POS_FLIST = new CFileList(pos_sample_fp_images_list_fn);

	m_nEndNoncpCascadeLevel = end_noncp_cascade_level;

	m_pHeightWidthRatios[0] = height_width_ratios[0];
	m_pHeightWidthRatios[1] = height_width_ratios[1];

	m_nMinDetWinWidth = min_det_win_width;
	m_nMaxDetWinWidth = max_det_win_width;

	m_nMaxNegExamplesNum = max_neg_examples_num;

	m_nMax2NegExamplesNum = MAX(max_neg_examples_num, max_2_neg_examples_num);

	m_bNegExamplesFromPos = neg_examples_from_pos;

	m_torchRNG = new Torch::Random;

    m_nFlippingPrevCroppedNegStartLevel = flipping_prev_cropped_neg_start_level;
    m_nFlippingNegImgStartLevel = flipping_neg_img_start_level;

    m_nNotEnoughNegCascadeStartLevel = -1;

    m_nNegBoostingCascadeLevelStep = neg_boosting_cascade_step;

    m_dbNegDetXnProbAlpha = neg_det_Xn_prob_alpha;

    m_pTimer = new CTimer();

    NTW = new CNegTestWindows(neg_test_windows_list_fn);
}

CNegExamplesGenerator::CNegExamplesGenerator()
{
	NEG_SG_OLIST = NULL;
	NEG_CP_OLIST = NULL;
	NEG_NONCP_OLIST = NULL;

	NEG_SG_FLIST = NULL;
	NEG_CP_FLIST = NULL;
	NEG_NONCP_FLIST = NULL;

	POS_OLIST = NULL;
	POS_FLIST = NULL;

	m_torchRNG = NULL;

    m_pTimer = NULL;

    m_bNegExamplesFromPos = false;

    NTW = NULL;
}

CNegExamplesGenerator::~CNegExamplesGenerator()
{
	CleanData();
}

void CNegExamplesGenerator::CleanData()
{
    if ( NEG_SG_OLIST )
        delete NEG_SG_OLIST;
    if ( NEG_CP_OLIST )
        delete NEG_CP_OLIST;
    if ( NEG_NONCP_OLIST )
        delete NEG_NONCP_OLIST;

    if ( NEG_SG_FLIST )
        delete NEG_SG_FLIST;
    if ( NEG_CP_FLIST )
        delete NEG_CP_FLIST;
    if ( NEG_NONCP_FLIST )
        delete NEG_NONCP_FLIST;

    if ( POS_OLIST )
        delete POS_OLIST;
    if ( POS_FLIST )
        delete POS_FLIST;

    if ( m_torchRNG )
        delete m_torchRNG;
    if ( m_pTimer )
        delete m_pTimer;
    if ( NTW )
        delete NTW;
}

char* CNegExamplesGenerator::GenerateFileName(const char* data_dir, const char* file_prefix, const char* file_ext, int idx)
{
	char *file_name = new char[1024];

    if ( data_dir[strlen(data_dir)-1] == '/' ||
        data_dir[strlen(data_dir)-1] == '\\' )
        sprintf(file_name, "%s%s_%06d.%s", data_dir, file_prefix, idx, file_ext);
    else
#ifndef _WIN32  /* under Linux system */
        sprintf(file_name, "%s/%s_%06d.%s", data_dir, file_prefix, idx, file_ext);
#else           /* under Window system */
        sprintf(file_name, "%s\\%s_%06d.%s", data_dir, file_prefix, idx, file_ext);
#endif

	return file_name;
}

void CNegExamplesGenerator::SetNegCroppedExampleList(char* output_dir)
{
    m_pOutputDir = output_dir;

    if ( !m_pOutputDir || strlen(m_pOutputDir) == 0 ) {
        printf("Please provide a temporary directory\n");
        printf("       for storing intermediate negative cropped example images and others!\n");
        exit(1);
    }

    char neg_cropped_ex_color_list_fn[2048];
    char neg_cropped_ex_fg_list_fn[2048];

    if ( output_dir[strlen(output_dir)-1] == '/' ||
        output_dir[strlen(output_dir)-1] == '\\' ) {
        sprintf(neg_cropped_ex_color_list_fn, "%sneg_cropped_color_images.lst", output_dir);
        sprintf(neg_cropped_ex_fg_list_fn, "%sneg_cropped_fp_images.lst", output_dir);
    }
    else {
#ifndef _WIN32  /* under Linux system */
        sprintf(neg_cropped_ex_color_list_fn, "%s/neg_cropped_color_images.lst", output_dir);
        sprintf(neg_cropped_ex_fg_list_fn, "%s/neg_cropped_fp_images.lst", output_dir);
#else           /* under Window system */
        sprintf(neg_cropped_ex_color_list_fn, "%s\\neg_cropped_color_images.lst", output_dir);
        sprintf(neg_cropped_ex_fg_list_fn, "%s\\neg_cropped_fp_images.lst", output_dir);
#endif
    }

    m_pNegCroppedExColorList = new CFileList(neg_cropped_ex_color_list_fn);
    m_pNegCroppedExFgList = new CFileList(neg_cropped_ex_fg_list_fn);

    m_pNegCroppedExColorList->Init(m_nMax2NegExamplesNum);
    if ( POS_FLIST && POS_FLIST->GetListLength() )
        m_pNegCroppedExFgList->Init(m_nMax2NegExamplesNum);
}

void CNegExamplesGenerator::ExportLogMessage(char *log_file_name, char *msg)
{
	if ( !log_file_name )
		return;

	ofstream fout(log_file_name,ios::app);
	if (fout.fail()) {
		printf("Error opening log output file %s.\n", log_file_name);
		fout.close();
		exit(0);
	}

	fout << msg;
	fout.close();
}

void CNegExamplesGenerator::GenerateNegExamples(char* boost_classifier_fn, char* log_msg_fn)
{
    CBoostClassifierCovariance* boost_classifier = new CBoostClassifierCovariance();
    boost_classifier->Load(boost_classifier_fn);

    GenerateNegExamples(boost_classifier, log_msg_fn);

    delete boost_classifier;
}

CIntegralRegionCov* CNegExamplesGenerator::GetIntegralCovariance(IplImage* image,
        int* feature_types, int n_feature_types, IplImage* fp_image)
{
    /* Related to the computation of integral feature images */
    CBuildFeatureImages* BFI = new CBuildFeatureImages(feature_types, n_feature_types);
    int feature_img_num = BFI->GetFeatureImageNumber();

    CvSize img_size = cvGetSize(image);
    IplImage** feature_imgs = new IplImage*[feature_img_num];
    for ( int a = 0 ; a < feature_img_num ; a++ ) {
        feature_imgs[a] = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
    }
    BFI->BuildFeatureImages(feature_imgs, image, fp_image);
    delete BFI;

    /* Compute the integral covariance images */
    CIntegralRegionCov *IRC = new CIntegralRegionCov(feature_img_num, img_size);
    IRC->SetNewData(feature_imgs);
    IRC->ComputeIntegralImages();

    /* Release memory */
    for ( int a = 0 ; a < feature_img_num ; a++ ) {
        cvReleaseImage(&feature_imgs[a]);
    }
    delete[] feature_imgs;

    return IRC;
}

CIntegralRegionCov* CNegExamplesGenerator::GetIntegralCovariance(
            char* image_fn, IplImage** image,
            int* feature_types, int n_feature_types,
            char* fp_image_fn, IplImage** fp_image)
{
    /* load images */
    *image = cvLoadImage(image_fn);
    if ( !(*image) ) {
        printf("Load image error: %s\n", image_fn);
        return NULL;
    }
    *fp_image = NULL;
    if ( fp_image_fn ) {
        *fp_image = cvLoadImage(fp_image_fn);
        if ( !(*fp_image) ) {
            printf("Load image error: %s\n", fp_image_fn);
            if ( (*image) )
                delete (*image);
            return NULL;
        }
    }

    CIntegralRegionCov* IRC = GetIntegralCovariance(*image, feature_types, n_feature_types, *fp_image);

    return IRC;
}

CIntegralRegionCov* CNegExamplesGenerator::GetIntegralCovariance(
            char* image_fn,
            int* feature_types, int n_feature_types,
            char* fp_image_fn)
{
    /* load images */
    IplImage* image = cvLoadImage(image_fn);
    if ( !image ) {
        printf("Load image error: %s\n", image_fn);
        return NULL;
    }
    IplImage* fp_image = NULL;
    if ( fp_image_fn ) {
        fp_image = cvLoadImage(fp_image_fn);
        if ( !fp_image ) {
            printf("Load image error: %s\n", fp_image_fn);
            if ( image )
                delete image;
            return NULL;
        }
    }

    CIntegralRegionCov* IRC = GetIntegralCovariance(image, feature_types, n_feature_types, fp_image);

    cvReleaseImage(&image);
    if ( fp_image )
        cvReleaseImage(&fp_image);

    return IRC;
}

void CNegExamplesGenerator::SaveNegCroppedExample(IplImage* image, CvRect roi, int example_idx,
        int edge_margin_size, int down_scale_factor, IplImage* fp_image)
{
	const char* neg_color_img_file_prefix = "neg_color";
	const char* neg_fp_img_file_prefix = "neg_fp";
	const char* neg_cropped_ext = "bmp";

    CvSize sub_img_size;

    /* save the negative cropped images to files */
    char* file_name = GenerateFileName(m_pOutputDir, neg_color_img_file_prefix, neg_cropped_ext, example_idx);

    CvRect updated_roi = roi;

    updated_roi.x -= edge_margin_size;
    updated_roi.y -= edge_margin_size;
    updated_roi.width += 2*edge_margin_size;
    updated_roi.height += 2*edge_margin_size;

    //printf("\tROI = [%d, %d, %d, %d]\n", updated_roi.x, updated_roi.y, updated_roi.width, updated_roi.height);

    if ( down_scale_factor == 1 ) {
        SaveImageROI(image, updated_roi, file_name);
    }
    else {
        sub_img_size.width = updated_roi.width/down_scale_factor;
        sub_img_size.height = updated_roi.height/down_scale_factor;
        IplImage* sub_image = cvCreateImage(sub_img_size, image->depth, image->nChannels);
        cvSetImageROI(image, updated_roi);
        cvResize(image, sub_image);
        cvSaveImage(file_name, sub_image);
        cvResetImageROI(image);
        cvReleaseImage(&sub_image);
    }

    m_pNegCroppedExColorList->SetFileName(example_idx, file_name);
    delete [] file_name;

    if ( fp_image ) {
        file_name = GenerateFileName(m_pOutputDir, neg_fp_img_file_prefix, neg_cropped_ext, example_idx);
        if ( down_scale_factor == 1 ) {
            SaveImageROI(fp_image, updated_roi, file_name);
        }
        else {
            IplImage* sub_image = cvCreateImage(sub_img_size, fp_image->depth, fp_image->nChannels);
            cvSetImageROI(fp_image, updated_roi);
            cvResize(fp_image, sub_image);
            cvSaveImage(file_name, sub_image);
            cvResetImageROI(fp_image);
            cvReleaseImage(&sub_image);
        }
        m_pNegCroppedExFgList->SetFileName(example_idx, file_name);
        delete [] file_name;
    }
}

void CNegExamplesGenerator::GenerateNegExamples_Prev(CBoostClassifierCovariance* boost_classifier,
        int &start_example_loc, int end_example_loc, char* log_msg_fn)
{
    if ( !m_pNegCroppedExColorList->GetListLength() )
        return;

    if ( start_example_loc >= end_example_loc )
        return;

    if ( NTW->m_nLastCascadeLevelNo > 0 || boost_classifier->m_nCascadeLevelNum > NTW->m_nLastCascadeLevelNo+2 )
        return;

    CIntegralRegionCov *IRC = NULL;
    IplImage *image=NULL, *fp_image=NULL;
    int* feature_img_types = boost_classifier->m_pFeatureImgTypes;
    int feature_type_num = boost_classifier->m_nFeatureTypeNum;
    char elapsed_time_msg[1024];
    char log_msg[2048];

    int used_prev_neg_examples_num = 0;
    int used_flip_prev_neg_examples_num = 0;
    int pos_flip_prev_neg_examples_num = 0;

	/* find some negative examples exampled from negative images obtained at the (k-1)-th cascade level */
    m_pTimer->Start();

    int prev_neg_cropped_images_num = m_pNegCroppedExColorList->GetListLength();

    int init_start_example_loc = start_example_loc;

    for ( int i = 0 ; i <  prev_neg_cropped_images_num ; i++ ) {
        if ( start_example_loc >= end_example_loc )
            break;

        IRC = GetIntegralCovariance(m_pNegCroppedExColorList->GetFileName(i), &image,
                feature_img_types, feature_type_num,
                m_pNegCroppedExFgList->GetFileName(i), &fp_image );
        if ( !IRC )
            continue;

        used_prev_neg_examples_num++;

        bool used_cur_example = false;
        CvRect neg_int_example_rect = cvRect(boost_classifier->m_nNegCroppedEdgeMarginSize,
            boost_classifier->m_nNegCroppedEdgeMarginSize,
            IRC->m_szFeatureImage.width-boost_classifier->m_nNegCroppedEdgeMarginSize*2,
            IRC->m_szFeatureImage.height-boost_classifier->m_nNegCroppedEdgeMarginSize*2);
//        if ( boost_classifier->Predictor(IRC, neg_int_example_rect, -1, 1.0, NULL, NULL, i) > 0  ) {
        if ( boost_classifier->Predictor(IRC, neg_int_example_rect, -1, m_dbNegDetXnProbAlpha) > 0  ) {
            /* save files */
            SaveNegCroppedExample(image, neg_int_example_rect, start_example_loc,
                    boost_classifier->m_nNegCroppedEdgeMarginSize, 1, fp_image);

            start_example_loc++;

            used_cur_example = true;

            SaveNegCroppedImageLists(start_example_loc);
        }

        /* flip the image for testing */
        if ( !used_cur_example &&
            boost_classifier->m_nCascadeLevelNum > m_nFlippingPrevCroppedNegStartLevel ) {
            if ( start_example_loc >= end_example_loc )
                break;

            delete IRC;
            cvFlip(image, image, 1);
            if ( fp_image )
                cvFlip(fp_image, fp_image, 1);

            IRC = GetIntegralCovariance(image,
                feature_img_types, feature_type_num,
                fp_image );
            if ( !IRC )
                continue;

            used_prev_neg_examples_num++;
            used_flip_prev_neg_examples_num++;

            neg_int_example_rect = cvRect(boost_classifier->m_nNegCroppedEdgeMarginSize,
                boost_classifier->m_nNegCroppedEdgeMarginSize,
                IRC->m_szFeatureImage.width-boost_classifier->m_nNegCroppedEdgeMarginSize*2,
                IRC->m_szFeatureImage.height-boost_classifier->m_nNegCroppedEdgeMarginSize*2);
    //        if ( boost_classifier->Predictor(IRC, neg_int_example_rect, -1, 1.0, NULL, NULL, i) > 0  ) {
            if ( boost_classifier->Predictor(IRC, neg_int_example_rect, -1, m_dbNegDetXnProbAlpha) > 0  ) {
                /* save files */
                SaveNegCroppedExample(image, neg_int_example_rect, start_example_loc,
                        boost_classifier->m_nNegCroppedEdgeMarginSize, 1, fp_image);

                start_example_loc++;
                pos_flip_prev_neg_examples_num++;

                SaveNegCroppedImageLists(start_example_loc);
            }
        }

        /* release memories */
        cvReleaseImage(&image);
        if ( fp_image )
            cvReleaseImage(&fp_image);
        delete IRC;
    }

    m_pTimer->Stop(true);
    m_pTimer->PrintElapsedTimeMsg( elapsed_time_msg );

    /* output message */
    sprintf(log_msg, "\n\tselecting %d examples from %d false postive examples at %d-th cascade level - (%s)\n",
        start_example_loc, used_prev_neg_examples_num,
        boost_classifier->m_nCascadeLevelNum, elapsed_time_msg);
    sprintf(log_msg, "%s\t%.4f(%%) false positive rate \n", log_msg,
        (double)(start_example_loc-init_start_example_loc)/(double)MAX(used_prev_neg_examples_num,1)*100.0);
    if ( used_flip_prev_neg_examples_num > 0 ) {
        int n_pos_neg_examples = start_example_loc-init_start_example_loc-pos_flip_prev_neg_examples_num;
        int n_used_neg_examples = MAX(used_prev_neg_examples_num-used_flip_prev_neg_examples_num,1);
        sprintf(log_msg, "%s\t%.4f(%%) false positive rate without flipping, %d of %d\n", log_msg,
            (double)n_pos_neg_examples/(double)n_used_neg_examples*100.0,
            n_pos_neg_examples, n_used_neg_examples);
        sprintf(log_msg, "%s\t%.4f(%%) false positive rate with flipping, %d of %d \n", log_msg,
            (double)(pos_flip_prev_neg_examples_num)/
            (double)MAX(used_flip_prev_neg_examples_num,1)*100.0,
            pos_flip_prev_neg_examples_num,
            used_flip_prev_neg_examples_num);
    }
    printf(log_msg);
    ExportLogMessage(log_msg_fn, log_msg);
}

void CNegExamplesGenerator::GenerateNegExamples_POS(CBoostClassifierCovariance* boost_classifier,
        int &start_example_loc, int end_example_loc, char* log_msg_fn)
{
    if ( start_example_loc >= end_example_loc )
        return;

    /* generate shuffled indices */
    int* shuffledIndices = NULL;
    if ( POS_OLIST->GetListLength() ) {
        shuffledIndices = new int[POS_OLIST->GetListLength()];
        m_torchRNG->getShuffledIndices(shuffledIndices, POS_OLIST->GetListLength());
    }
    else
        return;

    CIntegralRegionCov *IRC = NULL;
    IplImage *image=NULL, *fp_image=NULL;
    int* feature_img_types = boost_classifier->m_pFeatureImgTypes;
    int feature_type_num = boost_classifier->m_nFeatureTypeNum;
    char elapsed_time_msg[1024];
    char log_msg[2048];

    int level_check_neg_pos_examples_num = cvCeil(MIN(
        (float)boost_classifier->m_nCascadeLevelNum/(float)boost_classifier->m_nMaxCascadeLevelNum,1.0f)*
        (float)m_nMaxNegExamplesNum);

    if ( level_check_neg_pos_examples_num <= 0 ) {
        delete [] shuffledIndices;
        return;
    }

    int eff_neg_pos_examples_num = 0;
    long checked_neg_pos_det_windows_num = 0;
	int possible_neg_examples_per_pos_image = 1000;
    int neg_extracted_examples_num_per_sg_pos_image = 1;
    int used_pos_examples_num = 0;

    /* find some negative examples exampled from positive example images */
    m_pTimer->Start();

    for ( int i = 0 ; i < POS_OLIST->GetListLength() ; i++ ) {
        if ( start_example_loc >= end_example_loc || eff_neg_pos_examples_num >= level_check_neg_pos_examples_num )
            break;

        IRC = GetIntegralCovariance(POS_OLIST->GetFileName(shuffledIndices[i]), &image,
                feature_img_types, feature_type_num,
                POS_FLIST->GetFileName(shuffledIndices[i]), &fp_image );

        if ( !IRC )
            continue;

        used_pos_examples_num++;

        /* randomly generate negative double detection window examples */
        CvRect neg_int_example_rect;
        int tot_neg_examples_num_per_image = 0;
        for ( int j = 0 ; j < possible_neg_examples_per_pos_image ; j++ ) {
            /* compute the integer window example */
            neg_int_example_rect = m_cRect.GetIntRect(GetRandNegativePosExample(IRC->m_szFeatureImage),
                    cvRect(0,0,IRC->m_szFeatureImage.width,IRC->m_szFeatureImage.height));

            int margin_size = boost_classifier->m_nNegCroppedEdgeMarginSize;
            neg_int_example_rect = UpdateRectangle(IRC->m_szFeatureImage,
                    neg_int_example_rect,
                    margin_size);

            checked_neg_pos_det_windows_num++;

            /* classify it using previous k-1 cascade boost classifiers */
            if ( boost_classifier->Predictor(IRC, neg_int_example_rect, -1, m_dbNegDetXnProbAlpha) == 0 )
                continue;

            tot_neg_examples_num_per_image++;
            eff_neg_pos_examples_num++;

            /* save the negative cropped images to files */
            SaveNegCroppedExample(image, neg_int_example_rect, start_example_loc,
                    margin_size, margin_size/boost_classifier->m_nNegCroppedEdgeMarginSize, fp_image);

            start_example_loc++;

            SaveNegCroppedImageLists(start_example_loc);

            if ( tot_neg_examples_num_per_image >= neg_extracted_examples_num_per_sg_pos_image ||
                start_example_loc >= end_example_loc || eff_neg_pos_examples_num >= level_check_neg_pos_examples_num )
                break;
        }

        /* release memories */
        cvReleaseImage(&image);
        if ( fp_image )
            cvReleaseImage(&fp_image);
        delete IRC;
    }
    m_pTimer->Stop(true);
    m_pTimer->PrintElapsedTimeMsg( elapsed_time_msg );

    delete [] shuffledIndices;

    /* output message */
    sprintf(log_msg, "\n\tselecting %d examples from %d positive images - (%s)\n",
        eff_neg_pos_examples_num, used_pos_examples_num, elapsed_time_msg);
    sprintf(log_msg, "%s\t%.4f(%%) false positive rate\n", log_msg,
        (double)eff_neg_pos_examples_num/(double)MAX(1,checked_neg_pos_det_windows_num)*100.0);
    printf(log_msg);
    ExportLogMessage(log_msg_fn, log_msg);
}
void CNegExamplesGenerator::GenerateNegExamples_NEG_SG(CBoostClassifierCovariance* boost_classifier,
        int &start_example_loc, int end_example_loc, char* log_msg_fn)
{
    if ( start_example_loc >= end_example_loc )
        return;

    /* generate shuffled indices */
    int* shuffledIndices = NULL;
    if ( NEG_SG_OLIST->GetListLength() ) {
        shuffledIndices = new int[NEG_SG_OLIST->GetListLength()];
        m_torchRNG->getShuffledIndices(shuffledIndices, NEG_SG_OLIST->GetListLength());
    }
    else
        return;

    int level_check_neg_sg_examples_num = cvCeil(MIN(1.5f*(float)boost_classifier->m_nCascadeLevelNum/
            (float)boost_classifier->m_nMaxCascadeLevelNum,1.0f)
            *(float)NEG_SG_OLIST->GetListLength());

    if ( level_check_neg_sg_examples_num <= 0 ) {
        delete [] shuffledIndices;
        return;
    }

    CIntegralRegionCov *IRC = NULL;
    IplImage *image=NULL, *fp_image=NULL;
    int* feature_img_types = boost_classifier->m_pFeatureImgTypes;
    int feature_type_num = boost_classifier->m_nFeatureTypeNum;
    char elapsed_time_msg[1024];
    char log_msg[2048];

    long checked_neg_sg_det_windows_num = 0;
    int eff_neg_sg_examples_num = 0;
    int possible_neg_examples_per_sg_image = 100;
    int neg_extracted_examples_num_per_sg_neg_image = 1;
    int used_image_examples_num = 0;

    /* find some negative examples exampled from negative single example images */
    m_pTimer->Start();

    for ( int i = 0 ; i < NEG_SG_OLIST->GetListLength() ; i++ ) {
        if ( start_example_loc >= end_example_loc || eff_neg_sg_examples_num >= level_check_neg_sg_examples_num )
            break;

        IRC = GetIntegralCovariance(NEG_SG_OLIST->GetFileName(shuffledIndices[i]), &image,
                feature_img_types, feature_type_num,
                NEG_SG_FLIST->GetFileName(shuffledIndices[i]), &fp_image );

        if ( !IRC )
            continue;

        used_image_examples_num++;

        /* randomly generate negative double detection window examples */
        CvRect neg_int_example_rect;
        int tot_neg_examples_num_per_image = 0;
        for ( int j = 0 ; j < possible_neg_examples_per_sg_image ; j++ ) {
            /* compute the integer window example */
            neg_int_example_rect = m_cRect.GetIntRect(GetRandNegativeSgExample(IRC->m_szFeatureImage),
                    cvRect(0,0,IRC->m_szFeatureImage.width,IRC->m_szFeatureImage.height));

            int margin_size = boost_classifier->m_nNegCroppedEdgeMarginSize;
            neg_int_example_rect = UpdateRectangle(IRC->m_szFeatureImage,
                    neg_int_example_rect,
                    margin_size);

            checked_neg_sg_det_windows_num++;

            /* classify it using previous k-1 cascade boost classifiers */
            if ( boost_classifier->Predictor(IRC, neg_int_example_rect, -1, m_dbNegDetXnProbAlpha) == 0 )
                continue;

            tot_neg_examples_num_per_image++;
            eff_neg_sg_examples_num++;

            /* save the negative cropped images to files */
            SaveNegCroppedExample(image, neg_int_example_rect, start_example_loc,
                    margin_size, margin_size/boost_classifier->m_nNegCroppedEdgeMarginSize, fp_image);

            start_example_loc++;

            SaveNegCroppedImageLists(start_example_loc);

            if ( tot_neg_examples_num_per_image >= neg_extracted_examples_num_per_sg_neg_image ||
                start_example_loc >= end_example_loc || eff_neg_sg_examples_num >= level_check_neg_sg_examples_num )
                break;
        }

        /* release memories */
        cvReleaseImage(&image);
        if ( fp_image )
            cvReleaseImage(&fp_image);
        delete IRC;
    }

    m_pTimer->Stop(true);
    m_pTimer->PrintElapsedTimeMsg( elapsed_time_msg );

    delete [] shuffledIndices;

    /* output message */
    sprintf(log_msg, "\n\tselecting %d examples from %d negative single example images - (%s)\n",
        eff_neg_sg_examples_num, used_image_examples_num, elapsed_time_msg);
    sprintf(log_msg, "%s\t%.4f(%%) false positive rate \n", log_msg,
        (double)eff_neg_sg_examples_num/(double)MAX(checked_neg_sg_det_windows_num,1)*100.0);
    printf(log_msg);
    ExportLogMessage(log_msg_fn, log_msg);
}

void CNegExamplesGenerator::TestNegExamples_NEG_CP(CBoostClassifierCovariance* boost_classifier,
        int &start_example_loc, int end_example_loc, char* log_msg_fn)
{
    if ( start_example_loc >= end_example_loc )
        return;

    int* shuffledIndices = NULL;
    if ( NEG_CP_OLIST->GetListLength() ) {
        shuffledIndices = new int[NEG_CP_OLIST->GetListLength()];
        m_torchRNG->getShuffledIndices(shuffledIndices, NEG_CP_OLIST->GetListLength());
    }
    else
        return;

    char elapsed_time_msg[1024];
    char log_msg[2048];

   /* for negative cp_images */
    int possible_neg_examples_per_image = MAX(cvCeil(1.5*
            (double)(end_example_loc-start_example_loc)/
            (double)NEG_CP_OLIST->GetListLength()), 1);

    long tot_org_test_windows_num = 0;
    long tot_cur_test_windows_num = 0;

    printf("\nTotal false positive rate = %f\n\n", NTW->GetFalsePositiveRate());

    NTW->SetNeededNegExamples(end_example_loc-start_example_loc);

    printf("\nSelection percent = %f\n\n", NTW->m_dbSelectPercent);

    /* find some negative examples exampled from negative cp example images */
    m_pTimer->Start();

    int used_image_examples_num = 0;
    for ( int i = 0 ; i < NEG_CP_OLIST->GetListLength() ; i++ ) {
        if ( start_example_loc >= end_example_loc )
            break;

	    /* load images */
		char* image_fn = NEG_CP_OLIST->GetFileName(shuffledIndices[i]);
		char* fp_image_fn = NEG_CP_FLIST->GetFileName(shuffledIndices[i]);

        /* save the negative cropped images to files */
        char* test_wins_fn = GenerateFileName(m_pOutputDir, "neg_cp_test_windows", "yml", shuffledIndices[i]);

	    CImageBoostingNegExamples* IBNE = new CImageBoostingNegExamples();

        if ( !IBNE->Import(test_wins_fn) ) {
            if ( !IBNE->Import(NTW->m_pTestWindowsFileList->GetFileName(shuffledIndices[i])) ) {
                double scale_factors[] = {1.0, 0.5};
                int n_scale = 2;
                IBNE->GenerateInitTestWindows(image_fn, fp_image_fn,
                    scale_factors, n_scale,
                    m_nMinDetWinWidth,
                    m_nMaxDetWinWidth,
                    boost_classifier->m_nNegCroppedEdgeMarginSize);
            }
            else {
                IBNE->FilterTestWindows(boost_classifier->m_nNegCroppedEdgeMarginSize);
                NTW->m_pTestWindowsFileList->SetFileName(shuffledIndices[i], test_wins_fn);
            }
        }

        IBNE->Test(boost_classifier, -1, m_dbNegDetXnProbAlpha,
            NEG_CP_OLIST->GetFileName(shuffledIndices[i]),
            NEG_CP_FLIST->GetFileName(shuffledIndices[i]) );


        IBNE->Export(test_wins_fn);

        delete [] test_wins_fn;

        used_image_examples_num++;

        tot_org_test_windows_num += IBNE->m_nTotOriginalTestWindowsNum;
        tot_cur_test_windows_num += IBNE->m_nCurTestWindowsNum;

        if ( NTW->m_pTestWindowsFileList->GetListLength() > 0 )
            possible_neg_examples_per_image = NTW->GetMaxSelectExamplesNum(shuffledIndices[i]);

        if ( possible_neg_examples_per_image ) {
            CvRect* test_windows = new CvRect[possible_neg_examples_per_image];
            int* scale_indices = new int[possible_neg_examples_per_image];

            int n_valid_test_windows = IBNE->GetTestWindows(test_windows, possible_neg_examples_per_image, scale_indices);

            printf("Frame %06d    -   %d / %d \n", i,
                possible_neg_examples_per_image, IBNE->m_nCurTestWindowsNum);

            for ( int j = 0 ; j < n_valid_test_windows ; j++ ) {
                if ( start_example_loc >= end_example_loc )
                    break;

                /* save the negative cropped images to files */
                SaveNegCroppedExample(IBNE->m_ppScaledImages[scale_indices[j]],
                        test_windows[j], start_example_loc,
                        boost_classifier->m_nNegCroppedEdgeMarginSize,
                        1,
                        IBNE->m_ppScaledFpImages[scale_indices[j]]);

                start_example_loc++;

                SaveNegCroppedImageLists(start_example_loc);
            }

            delete [] test_windows;
            delete [] scale_indices;
        }

        delete IBNE;

        if ( start_example_loc >= end_example_loc )
            break;
    }

    m_pTimer->Stop(true);
    m_pTimer->PrintElapsedTimeMsg( elapsed_time_msg );

    delete [] shuffledIndices;

	if ( boost_classifier->m_pFalsePosRates ) {
		boost_classifier->m_pFalsePosRates[boost_classifier->m_nCascadeLevelNum] =
			(double)tot_cur_test_windows_num/(double)MAX(tot_org_test_windows_num,1);
	}

    /* output message */
    sprintf(log_msg, "\n\tselecting %d examples from %d negative cp example images - (%s)\n",
        end_example_loc-start_example_loc, used_image_examples_num, elapsed_time_msg);
    sprintf(log_msg, "%s\t%.4f(%%) false positive rate \n", log_msg,
        (double)tot_cur_test_windows_num/(double)MAX(tot_org_test_windows_num,1)*100.0);
    printf(log_msg);
    ExportLogMessage(log_msg_fn, log_msg);
}

void CNegExamplesGenerator::GenerateNegExamples_NEG_CP(CBoostClassifierCovariance* boost_classifier,
        int &start_example_loc, int end_example_loc, char* log_msg_fn)
{
    if ( start_example_loc >= end_example_loc )
        return;

    int* shuffledIndices = NULL;
    if ( NEG_CP_OLIST->GetListLength() ) {
        shuffledIndices = new int[NEG_CP_OLIST->GetListLength()];
        m_torchRNG->getShuffledIndices(shuffledIndices, NEG_CP_OLIST->GetListLength());
    }
    else
        return;

   /* for negative cp_images */
    int possible_neg_examples_per_image = MAX(cvCeil(1.5*
            (double)(end_example_loc-start_example_loc)/
            (double)NEG_CP_OLIST->GetListLength()), 1);

    CIntegralRegionCov *IRC = NULL;
	IplImage **scale_images=NULL, **scale_fp_images=NULL;
    int* feature_img_types = boost_classifier->m_pFeatureImgTypes;
    int feature_type_num = boost_classifier->m_nFeatureTypeNum;
    char elapsed_time_msg[1024];
    char log_msg[2048];

    long checked_det_windows_num = 0;
    int eff_neg_cp_examples_num = 0;

    int used_image_examples_num = 0;

	double scale_factors[] = {1.0, 0.5};
	int n_scale = 2;

	CIntegralRegionCov **IRCs = new CIntegralRegionCov*[n_scale];
	scale_images = new IplImage*[n_scale];
	scale_fp_images = new IplImage*[n_scale];

    /* find some negative examples exampled from negative cp example images */
    m_pTimer->Start();

    for ( int i = 0 ; i < NEG_CP_OLIST->GetListLength() ; i++ ) {
        if ( start_example_loc >= end_example_loc )
            break;

	    /* load images */
		char* image_fn = NEG_CP_OLIST->GetFileName(shuffledIndices[i]);
		char* fp_image_fn = NEG_CP_FLIST->GetFileName(shuffledIndices[i]);
	    scale_images[0] = cvLoadImage(image_fn);
	    if ( !scale_images[0] ) {
	        printf("Load image error: %s\n", image_fn);
			continue;
	    }
	    scale_fp_images[0] = NULL;
	    if ( fp_image_fn ) {
	        scale_fp_images[0] = cvLoadImage(fp_image_fn);
	        if ( !scale_fp_images[0] ) {
	            printf("Load image error: %s\n", fp_image_fn);
	            if ( scale_images[0] )
	                delete scale_images[0];
				continue;
	        }
	    }

	    /* flipping the images */
	    if ( boost_classifier->m_nCascadeLevelNum > m_nFlippingNegImgStartLevel ) {
	        cvFlip(scale_images[0], scale_images[0], 1);
	        if ( scale_fp_images[0] )
                cvFlip(scale_fp_images[0], scale_fp_images[0], 1);
	    }

		for ( int s = 1 ; s < n_scale ; s++ ) {
			double image_scale = scale_factors[s];
			CvSize scale_img_size = cvSize((int)((double)scale_images[0]->width*image_scale),
					(int)((double)scale_images[0]->height*image_scale));
			scale_images[s] = cvCreateImage(scale_img_size, scale_images[0]->depth, scale_images[0]->nChannels);
			cvResize(scale_images[0], scale_images[s]);
			scale_fp_images[s] = NULL;
			if ( scale_fp_images[0] ) {
				scale_fp_images[s] = cvCreateImage(scale_img_size, scale_fp_images[0]->depth, scale_fp_images[0]->nChannels);
				cvResize(scale_fp_images[0], scale_fp_images[s]);
			}
		}
		for ( int s = 0 ; s < n_scale ; s++ )
			IRCs[s] = GetIntegralCovariance(scale_images[s], feature_img_types, feature_type_num, scale_fp_images[s]);

        used_image_examples_num++;

        for ( int itr = 0 ; itr < 6 ; itr++ ) {
			int rnd_scale = m_torchRNG->random()%n_scale;
			IRC = IRCs[rnd_scale];

            /* get test windows */
            CTestWindows* TW = GetTestWindows(cvGetSize(scale_images[rnd_scale]));
            int* indices = new int[TW->m_nTestWindowsNum];
            m_torchRNG->getShuffledIndices(indices, TW->m_nTestWindowsNum);

            bool loop_break = false;
            /* randomly generate negative double detection window examples */
            CvRect neg_int_example_rect;
            int tot_neg_examples_num_per_image = 0;
            for ( int j = 0 ; j < TW->m_nTestWindowsNum ; j++ ) {
                neg_int_example_rect = TW->m_pTestWindows[indices[j]];
				if ( neg_int_example_rect.width > m_nMaxDetWinWidth )
					continue;

                int margin_size = boost_classifier->m_nNegCroppedEdgeMarginSize;
                neg_int_example_rect = UpdateRectangle(IRC->m_szFeatureImage,
                        neg_int_example_rect,
                        margin_size);

                checked_det_windows_num++;

                /* classify it using previous k-1 cascade boost classifiers */
                if ( boost_classifier->Predictor(IRC, neg_int_example_rect, -1, m_dbNegDetXnProbAlpha) == 0 )
                    continue;

                tot_neg_examples_num_per_image++;
                eff_neg_cp_examples_num++;

				if ( margin_size/boost_classifier->m_nNegCroppedEdgeMarginSize > 1 ) {
					printf("too large negative example image !\n");
				}

                /* save the negative cropped images to files */
                SaveNegCroppedExample(scale_images[rnd_scale], neg_int_example_rect, start_example_loc,
                        margin_size, margin_size/boost_classifier->m_nNegCroppedEdgeMarginSize, scale_fp_images[rnd_scale]);

                start_example_loc++;

                SaveNegCroppedImageLists(start_example_loc);

                if ( tot_neg_examples_num_per_image >= possible_neg_examples_per_image ||
                    start_example_loc >= end_example_loc ) {
                    loop_break = true;
                    break;
                }
            }

            /* release memories */
            delete TW;
            delete [] indices;

            if ( loop_break || tot_neg_examples_num_per_image >
                MAX(possible_neg_examples_per_image/3,1) )
                break;
        }

        /* release memories */
		for ( int s = 0 ; s < n_scale ; s++ ) {
			cvReleaseImage(&scale_images[s]);
			if ( scale_fp_images[s] )
				cvReleaseImage(&scale_fp_images[s]);
			delete IRCs[s];
		}
    }
    m_pTimer->Stop(true);
    m_pTimer->PrintElapsedTimeMsg( elapsed_time_msg );

    delete [] shuffledIndices;
	delete [] scale_images;
	delete [] scale_fp_images;
	delete [] IRCs;

	if ( boost_classifier->m_pFalsePosRates ) {
		boost_classifier->m_pFalsePosRates[boost_classifier->m_nCascadeLevelNum] =
			(double)eff_neg_cp_examples_num/(double)MAX(checked_det_windows_num,1);
	}

    /* output message */
    sprintf(log_msg, "\n\tselecting %d examples from %d negative cp example images - (%s)\n",
        eff_neg_cp_examples_num, used_image_examples_num, elapsed_time_msg);
    sprintf(log_msg, "%s\t%.4f(%%) false positive rate \n", log_msg,
        (double)eff_neg_cp_examples_num/(double)MAX(checked_det_windows_num,1)*100.0);
    printf(log_msg);
    ExportLogMessage(log_msg_fn, log_msg);
}

void CNegExamplesGenerator::GenerateNegExamples_NEG_NONCP(CBoostClassifierCovariance* boost_classifier,
        int &start_example_loc, int end_example_loc, char* log_msg_fn)
{
    if ( start_example_loc >= end_example_loc )
        return;

    if ( boost_classifier->m_nCascadeLevelNum >= m_nEndNoncpCascadeLevel )
        return;

    int* shuffledIndices = NULL;
    if ( NEG_NONCP_OLIST->GetListLength() ) {
        shuffledIndices = new int[NEG_NONCP_OLIST->GetListLength()];
        m_torchRNG->getShuffledIndices(shuffledIndices, NEG_NONCP_OLIST->GetListLength());
    }
    else
        return;

    CIntegralRegionCov *IRC = NULL;
	IplImage **scale_images=NULL, **scale_fp_images=NULL;
    int* feature_img_types = boost_classifier->m_pFeatureImgTypes;
    int feature_type_num = boost_classifier->m_nFeatureTypeNum;
    char elapsed_time_msg[1024];
    char log_msg[2048];

    int possible_neg_examples_per_image = MAX(cvCeil(1.5*
            (double)(end_example_loc-start_example_loc)/
            (double)NEG_CP_OLIST->GetListLength()), 1);

    long checked_det_windows_num = 0;
    int eff_neg_noncp_examples_num = 0;

	double scale_factors[] = {1.0, 0.5, 0.25};
	int n_scale = 3;

	CIntegralRegionCov **IRCs = new CIntegralRegionCov*[n_scale];
	scale_images = new IplImage*[n_scale];
	scale_fp_images = new IplImage*[n_scale];

    int used_image_examples_num = 0;

    m_pTimer->Start();


    /* for negative noncp_images */
    for ( int i = 0 ; i < NEG_NONCP_OLIST->GetListLength() ; i++ ) {
        if ( start_example_loc >= end_example_loc )
            break;

	    /* load images */
		char* image_fn = NEG_NONCP_OLIST->GetFileName(shuffledIndices[i]);
        int rand_fp_idx = cvFloor(m_torchRNG->uniform() * (double)NEG_NONCP_OLIST->GetListLength());
		char* fp_image_fn = NEG_NONCP_FLIST->GetFileName(rand_fp_idx);
	    scale_images[0] = cvLoadImage(image_fn);
	    if ( !scale_images[0] ) {
	        printf("Load image error: %s\n", image_fn);
			continue;
	    }
	    scale_fp_images[0] = NULL;
	    if ( fp_image_fn ) {
	        scale_fp_images[0] = cvLoadImage(fp_image_fn);
	        if ( !scale_fp_images[0] ) {
	            printf("Load image error: %s\n", fp_image_fn);
	            if ( scale_images[0] )
	                delete scale_images[0];
				continue;
	        }
	    }

	    /* flipping the images */
	    if ( boost_classifier->m_nCascadeLevelNum > m_nFlippingNegImgStartLevel ) {
	        cvFlip(scale_images[0], scale_images[0], 1);
	        if ( scale_fp_images[0] )
                cvFlip(scale_fp_images[0], scale_fp_images[0], 1);
	    }

		for ( int s = 1 ; s < n_scale ; s++ ) {
			double image_scale = scale_factors[s];
			CvSize scale_img_size = cvSize((int)((double)scale_images[0]->width*image_scale),
					(int)((double)scale_images[0]->height*image_scale));
			scale_images[s] = cvCreateImage(scale_img_size, scale_images[0]->depth, scale_images[0]->nChannels);
			cvResize(scale_images[0], scale_images[s]);
			scale_fp_images[s] = NULL;
			if ( scale_fp_images[0] ) {
				scale_fp_images[s] = cvCreateImage(scale_img_size, scale_fp_images[0]->depth, scale_fp_images[0]->nChannels);
				cvResize(scale_fp_images[0], scale_fp_images[s]);
			}
		}
		for ( int s = 0 ; s < n_scale ; s++ )
			IRCs[s] = GetIntegralCovariance(scale_images[s], feature_img_types, feature_type_num, scale_fp_images[s]);

        used_image_examples_num++;

        for ( int itr = 0 ; itr < 2 ; itr++ ) {
			int rnd_scale = m_torchRNG->random()%n_scale;
			IRC = IRCs[rnd_scale];

            /* get test windows */
            CTestWindows* TW = GetTestWindows(cvGetSize(scale_images[rnd_scale]));
            int* indices = new int[TW->m_nTestWindowsNum];
            m_torchRNG->getShuffledIndices(indices, TW->m_nTestWindowsNum);

            bool loop_break = false;

            /* randomly generate negative double detection window examples */
            CvRect neg_int_example_rect;
            int tot_neg_examples_num_per_image = 0;
            for ( int j = 0 ; j < TW->m_nTestWindowsNum ; j++ ) {
                neg_int_example_rect = TW->m_pTestWindows[indices[j]];
				if ( neg_int_example_rect.width > m_nMaxDetWinWidth )
					continue;

                int margin_size = boost_classifier->m_nNegCroppedEdgeMarginSize;
                neg_int_example_rect = UpdateRectangle(IRC->m_szFeatureImage,
                        neg_int_example_rect,
                        margin_size);

                checked_det_windows_num++;

                /* classify it using previous k-1 cascade boost classifiers */
                if ( boost_classifier->Predictor(IRC, neg_int_example_rect, -1, m_dbNegDetXnProbAlpha) == 0 )
                    continue;

                tot_neg_examples_num_per_image++;
                eff_neg_noncp_examples_num++;

				if ( margin_size/boost_classifier->m_nNegCroppedEdgeMarginSize > 1 ) {
					printf("too large negative example image !\n");
				}

                /* save the negative cropped images to files */
                SaveNegCroppedExample(scale_images[rnd_scale], neg_int_example_rect, start_example_loc,
                        margin_size, margin_size/boost_classifier->m_nNegCroppedEdgeMarginSize, scale_fp_images[rnd_scale]);

                start_example_loc++;

                SaveNegCroppedImageLists(start_example_loc);

                if ( tot_neg_examples_num_per_image >= possible_neg_examples_per_image ||
                    start_example_loc >= end_example_loc ) {
                    loop_break = true;
                    break;
                }
            }

            /* release memories */
            delete TW;
            delete [] indices;

            if ( loop_break || tot_neg_examples_num_per_image >
                MAX(possible_neg_examples_per_image/3,1) )
                break;
        }

        /* release memories */
		for ( int s = 0 ; s < n_scale ; s++ ) {
			cvReleaseImage(&scale_images[s]);
			if ( scale_fp_images[s] )
				cvReleaseImage(&scale_fp_images[s]);
			delete IRCs[s];
		}
    }

    m_pTimer->Stop(true);
    m_pTimer->PrintElapsedTimeMsg( elapsed_time_msg );

    delete [] shuffledIndices;
	delete [] scale_images;
	delete [] scale_fp_images;
	delete [] IRCs;

    /* output message */
    sprintf(log_msg, "\n\tselecting %d examples from %d negative noncp example images - (%s)\n",
        eff_neg_noncp_examples_num, used_image_examples_num, elapsed_time_msg);
    sprintf(log_msg, "%s\t%.4f(%%) false positive rate \n", log_msg,
        (double)eff_neg_noncp_examples_num/(double)MAX(checked_det_windows_num,1)*100.0);
    printf(log_msg);
    ExportLogMessage(log_msg_fn, log_msg);
}

void CNegExamplesGenerator::GenerateNegExamples(CBoostClassifierCovariance* boost_classifier, char* log_msg_fn)
{
    int start_example_loc = 0;
    int max_neg_examples_num = boost_classifier->m_nCascadeLevelNum == 0 ?
        m_nMaxNegExamplesNum: m_nMax2NegExamplesNum;
    int end_example_loc = max_neg_examples_num;

    char log_msg[2048];

    sprintf(log_msg, "\n\nTraining Logit-Boost Level %d ...\n", boost_classifier->m_nCascadeLevelNum);
    printf(log_msg);
    ExportLogMessage(log_msg_fn, log_msg);

    int n_prev_neg_examples = m_pNegCroppedExColorList->GetListLength();

    GenerateNegExamples_Prev(boost_classifier, start_example_loc, end_example_loc, log_msg_fn);

    if ( boost_classifier->m_nCascadeLevelNum > 15 &&
        n_prev_neg_examples < max_neg_examples_num &&
        m_nNotEnoughNegCascadeStartLevel < 0 ) {
        m_nNotEnoughNegCascadeStartLevel = boost_classifier->m_nCascadeLevelNum;
    }

    if ( m_nNegBoostingCascadeLevelStep == 0 ||
        m_nNotEnoughNegCascadeStartLevel < 0 ||
        ( m_nNotEnoughNegCascadeStartLevel > 0 &&
        (m_nNotEnoughNegCascadeStartLevel - boost_classifier->m_nCascadeLevelNum)%
        (m_nNegBoostingCascadeLevelStep+1) == 0 ) ) {
        GenerateNegExamples_NEG_SG(boost_classifier, start_example_loc, end_example_loc, log_msg_fn);

        if ( m_bNegExamplesFromPos )
            GenerateNegExamples_POS(boost_classifier, start_example_loc, end_example_loc, log_msg_fn);

        if ( NEG_NONCP_OLIST->GetListLength() )
            end_example_loc = start_example_loc + cvFloor((float)(max_neg_examples_num-start_example_loc)*0.95f);
        if ( NTW->m_nLastCascadeLevelNo <= 0 ||
            boost_classifier->m_nCascadeLevelNum <= NTW->m_nLastCascadeLevelNo+2 )
            GenerateNegExamples_NEG_CP(boost_classifier, start_example_loc, end_example_loc, log_msg_fn);
        else
            TestNegExamples_NEG_CP(boost_classifier, start_example_loc, end_example_loc, log_msg_fn);

        end_example_loc = max_neg_examples_num;
        GenerateNegExamples_NEG_NONCP(boost_classifier, start_example_loc, end_example_loc, log_msg_fn);
    }

    /* save the covariance data for negative examples */
    SaveNegCroppedImageLists(start_example_loc, true);

    sprintf(log_msg, "\nStore the images for %d negative cropped examples \n", start_example_loc);
    sprintf(log_msg, "%s\t:%s\n", log_msg, m_pNegCroppedExColorList->GetFileListName());
    printf(log_msg);
    ExportLogMessage(log_msg_fn, log_msg);
}

void CNegExamplesGenerator::SaveNegCroppedImageLists(int n_examples, bool finish)
{
    if ( !finish && n_examples%5 == 0 ) {
        printf("=");
        if ( n_examples%250 == 0 ) {
			char buffer [80];
            m_pTimer->PrintLocalTime(buffer, 80, "Now is %Y-%m-%d %H:%M:%S");
            printf("%6d  ->  %s\n", n_examples, buffer);
		}
    }

    if ( !finish && n_examples%100 != 0 )
        return;

    /* save the covariance data for negative examples */
    m_pNegCroppedExColorList->SetListLength(n_examples);
    m_pNegCroppedExColorList->Write();
    if ( m_pNegCroppedExFgList->GetListLength() ) {
        m_pNegCroppedExFgList->SetListLength(n_examples);
        m_pNegCroppedExFgList->Write();
    }
}

void CNegExamplesGenerator::SaveImageROI(IplImage* img, CvRect roi, char* file_name)
{
    if ( roi.x < 0 || roi.y < 0 ||
        roi.x+roi.width > img->width ||
        roi.y+roi.height > img->height ) {
        printf("Saving negative cropped image error - beyound image boundary!\n");
        exit(1);
    }

	cvSetImageROI(img, roi);
	cvSaveImage(file_name, img);
	cvResetImageROI(img);
}

floatRect CNegExamplesGenerator::GetRandNegativeSgExample(CvSize img_size)
{
    floatRect neg_win_example;

	double max_win_percent = 1.0;
	double min_win_percent = 0.4;

	double hw_ratio = (double)img_size.height / (double)img_size.width;

    for ( ; ; ) {
        neg_win_example.width = m_torchRNG->boundedUniform(0.6, 1.0);
        neg_win_example.height = neg_win_example.width *
                m_torchRNG->boundedUniform(m_pHeightWidthRatios[0], m_pHeightWidthRatios[1])/hw_ratio;

        if ( neg_win_example.width * neg_win_example.height >= max_win_percent ||
            neg_win_example.width * neg_win_example.height <= min_win_percent )
            continue;

        neg_win_example.x = m_torchRNG->boundedUniform(0, 1.0-neg_win_example.width);
        neg_win_example.y = m_torchRNG->boundedUniform(0, 1.0-neg_win_example.height);

        break;
    }

    return neg_win_example;
}

floatRect CNegExamplesGenerator::GetRandNegativeExample(CvSize img_size)
{
    floatRect neg_win_example;

	double min_hw_ratio = m_pHeightWidthRatios[0];
	double max_hw_ratio = m_pHeightWidthRatios[1];

	double min_float_height, max_float_height;
	min_float_height = (double)(m_nMinDetWinWidth*min_hw_ratio)/(double)img_size.height;
	min_float_height = MIN(min_float_height, 0.9);

	max_float_height = (double)(m_nMaxDetWinWidth*max_hw_ratio)/(double)img_size.height;
	max_float_height = MIN(max_float_height, 1.0);

	double max_min_float_height = max_float_height - min_float_height;
	double max_min_hw_ratio = max_hw_ratio-min_hw_ratio;

    neg_win_example.height = min_float_height + max_min_float_height * m_torchRNG->uniform();
    double win_neg_width = neg_win_example.height * (double)img_size.height / ( min_hw_ratio + max_min_hw_ratio * m_torchRNG->uniform() );
    if ( win_neg_width < m_nMinDetWinWidth ) {
        double scale = (double)m_nMinDetWinWidth/(double)win_neg_width;
        win_neg_width = m_nMinDetWinWidth;
        neg_win_example.height *= scale;
        neg_win_example.height = MIN(neg_win_example.height, 1.0);
    }
    win_neg_width = MIN(win_neg_width, img_size.width);
    neg_win_example.width = (double)win_neg_width/(double)img_size.width;

    neg_win_example.x = (1.0-neg_win_example.width) * m_torchRNG->uniform();
    neg_win_example.y = (1.0-neg_win_example.height) * m_torchRNG->uniform();

    return neg_win_example;
}


floatRect CNegExamplesGenerator::GetRandNegativePosExample(CvSize img_size)
{
    floatRect neg_win_example;

	double max_win_percent = 0.4;
	double min_win_percent = 0.15;
	double max_edge_marg_percent = 0.1;

	double hw_ratio = (double)img_size.height / (double)img_size.width;

    for ( ; ; ) {
        neg_win_example.width = m_torchRNG->boundedUniform(0.35, 0.6);
        neg_win_example.height = neg_win_example.width *
                m_torchRNG->boundedUniform(m_pHeightWidthRatios[0], m_pHeightWidthRatios[1])/hw_ratio;

        if ( neg_win_example.width * neg_win_example.height >= max_win_percent ||
            neg_win_example.width * neg_win_example.height <= min_win_percent )
            continue;

        switch(m_torchRNG->random()%4) {
        case 0:   /* left-top */
            neg_win_example.x = m_torchRNG->boundedUniform(0, max_edge_marg_percent);
            neg_win_example.y = m_torchRNG->boundedUniform(0, max_edge_marg_percent);
            break;
        case 1:   /* left-bottom */
            neg_win_example.x = m_torchRNG->boundedUniform(0, max_edge_marg_percent);
            neg_win_example.y = 1.0-neg_win_example.height-m_torchRNG->boundedUniform(0, max_edge_marg_percent);
            break;
        case 2:   /* right-top */
            neg_win_example.x = 1.0-neg_win_example.width-m_torchRNG->boundedUniform(0, max_edge_marg_percent);
            neg_win_example.y = m_torchRNG->boundedUniform(0, max_edge_marg_percent);
            break;
        case 3:   /* right-bottom */
            neg_win_example.x = 1.0-neg_win_example.width-m_torchRNG->boundedUniform(0, max_edge_marg_percent);
            neg_win_example.y = 1.0-neg_win_example.height-m_torchRNG->boundedUniform(0, max_edge_marg_percent);
            break;
        }

        if ( neg_win_example.x < 0 )
            neg_win_example.x = 0.0;
        if ( neg_win_example.y < 0 )
            neg_win_example.y = 0.0;

        /*
        neg_win_example.x = m_torchRNG->boundedUniform(0, 1.0-neg_win_example.width);
        neg_win_example.y = m_torchRNG->boundedUniform(0, 1.0-neg_win_example.height);
        */

        break;
    }

    return neg_win_example;
}

CvRect CNegExamplesGenerator::UpdateRectangle(CvSize img_size, CvRect roi, int &edge_margin_size)
{

    if ( roi.width > m_nMaxDetWinWidth )
        edge_margin_size *= (roi.width/m_nMaxDetWinWidth + 1);


    if ( img_size.width < edge_margin_size*2 ||
        img_size.height < edge_margin_size*2 ) {
        printf("Edge margin size is too larger than image size!\n");
        exit(1);
    }

    CvRect updated_roi = roi;

    CvPoint ipt1 = cvPoint(edge_margin_size,edge_margin_size);
    CvPoint ipt2 = cvPoint(img_size.width-edge_margin_size*2,
            img_size.height-edge_margin_size*2);

    CvPoint rpt1 = cvPoint(updated_roi.x,updated_roi.y);
    CvPoint rpt2 = cvPoint(updated_roi.x+updated_roi.width, updated_roi.y+updated_roi.height);

    updated_roi.x = MAX(ipt1.x, rpt1.x);
    updated_roi.y = MAX(ipt1.y, rpt1.y);

    updated_roi.width = MIN(ipt2.x, rpt2.x) - updated_roi.x;
    updated_roi.height = MIN(ipt2.y, rpt2.y) - updated_roi.y;

    return updated_roi;
}

CTestWindows* CNegExamplesGenerator::GetTestWindows(CvSize img_size)
{
    int min_win_width = m_nMinDetWinWidth;
    float w_scale_factor = 0.1f + m_torchRNG->boundedUniform(0.0, 0.1);
    float min_hw_ratio = m_pHeightWidthRatios[0];
    float h_scale_factor = 0.1f + m_torchRNG->boundedUniform(0.0, 0.1);
    int n_hw_scales = cvRound( log(m_pHeightWidthRatios[1]/m_pHeightWidthRatios[0])/
                    log(1.0f+h_scale_factor) );
    float w_shift_factor = 0.1f + m_torchRNG->boundedUniform(0.0, 0.1);
    float h_shift_factor = 0.05f + m_torchRNG->boundedUniform(0.0, 0.1);
    int n_win_width_scales = 15;

    CTestWindows* TW = new CTestWindows(min_win_width,
            n_win_width_scales,
            w_scale_factor,
            min_hw_ratio,
            n_hw_scales,
            h_scale_factor,
            w_shift_factor,
            h_shift_factor);

    CvRect image_roi;
    image_roi.x = (int)(m_torchRNG->boundedUniform(0.0, 1.0)*4.0);
    image_roi.y = (int)(m_torchRNG->boundedUniform(0.0, 1.0)*4.0);
    image_roi.width = img_size.width - image_roi.x;
    image_roi.height = img_size.height - image_roi.y;

    TW->GetTestWindows(image_roi);

    return TW;
}
