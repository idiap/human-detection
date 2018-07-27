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
#include "BoostClassifierCovariance.h"

CCascadeClassifierCovariance::CCascadeClassifierCovariance(int features_num, int tot_features_num,
			int weak_classifiers_num, bool map, bool used_avg)
{
	m_nMaxWeakClassifiersNum = weak_classifiers_num;
	m_nWeakClassifiersNum = 0;

	m_nFeaturesNum = features_num;
	m_nTotFeaturesNum = tot_features_num;
	m_bMap = map;
	m_bUsedFeatureMean = used_avg;

	m_ppWeakClassifers = new WeakClassifierCovariance*[m_nMaxWeakClassifiersNum];

	for ( int i = 0 ; i < m_nMaxWeakClassifiersNum ; i++ ) {
		//m_ppWeakClassifers[i] = new WeakClassifierCovariance(-1,-1,-1,-1,features_num,tot_features_num,NULL,map,used_avg);
		m_ppWeakClassifers[i] = NULL;
	}

}

CCascadeClassifierCovariance::~CCascadeClassifierCovariance()
{
	for ( int i = 0 ; i < m_nMaxWeakClassifiersNum ; i++ ) {
	    if ( m_ppWeakClassifers[i] )
            delete m_ppWeakClassifers[i];
	}
	delete [] m_ppWeakClassifers;
}

void CCascadeClassifierCovariance::GetClassificationRates(double* classifier_signs, int* labels, int n_examples)
{
    m_dbPositiveDetectionRate = m_dbFalsePositiveRate = 0;

    int n_pos_examples = 0, n_neg_examples = 0;
    for ( int i = 0 ; i < n_examples ; i++ ) {
        if ( labels[i] ) {
            n_pos_examples++;
            if ( classifier_signs[i] > 0 )
                m_dbPositiveDetectionRate += 1.0;
        }
        else {
            n_neg_examples++;
            if ( classifier_signs[i] > 0 )
                m_dbFalsePositiveRate += 1.0;
        }
    }

    m_dbPositiveDetectionRate /= (double)n_pos_examples;
    m_dbFalsePositiveRate /= (double)n_neg_examples;
}

bool CCascadeClassifierCovariance::AddWeakClassifier(WeakClassifierCovariance* weak_classifier)
{
    if ( weak_classifier == NULL ) {
        if ( m_nWeakClassifiersNum >= m_nMaxWeakClassifiersNum )
            return false;
        if ( m_ppWeakClassifers[m_nWeakClassifiersNum] )
            delete m_ppWeakClassifers[m_nWeakClassifiersNum];
		m_ppWeakClassifers[m_nWeakClassifiersNum] = new
            WeakClassifierCovariance(-1,-1,-1,-1,m_nFeaturesNum,m_nTotFeaturesNum,NULL,m_bMap,m_bUsedFeatureMean);
        m_nWeakClassifiersNum++;
    }
    else {
        if ( weak_classifier->negative_log_likelihood <
            m_ppWeakClassifers[m_nWeakClassifiersNum-1]->negative_log_likelihood )
            weak_classifier->Copy(m_ppWeakClassifers[m_nWeakClassifiersNum-1]);
    }

    return true;
}

/************************************************************************************/
/************************************************************************************/
/************************************************************************************/

CBoostClassifierCovariance::CBoostClassifierCovariance(int *feature_img_types, int feature_type_num,
		int *cascade_subset_sizes, int max_cascade_level_num, int max_weak_classifiers_num,
		bool map, bool identity_map, bool used_avg, bool mean_cov_fall_ex, bool normalization,
		double level_reject_percent,
		double pos_example_edge_margin_x_percent,
		double pos_example_edge_margin_y_percent,
        double pos_example_shift_percent,
        int pos_example_shift_itr_num,
        int neg_cropped_edge_margin_size)
{
	m_nMaxCascadeLevelNum = max_cascade_level_num;
	m_nMaxWeakClassifiersNum = max_weak_classifiers_num;
	m_nCascadeLevelNum = 0;

	m_dbLevelRejectPercent = level_reject_percent;

	m_dbPosExEdgeMarginXPercent = pos_example_edge_margin_x_percent;
	m_dbPosExEdgeMarginYPercent = pos_example_edge_margin_y_percent;
    m_dbPosExShiftPercent = pos_example_shift_percent;
    m_nPosExShiftItrNum = pos_example_shift_itr_num;
	m_nNegCroppedEdgeMarginSize = neg_cropped_edge_margin_size;

	m_nFeatureTypeNum = feature_type_num;
	m_pFeatureImgTypes = NULL;
    if ( feature_img_types ) {
        m_pFeatureImgTypes = new int[m_nFeatureTypeNum];
        memcpy(m_pFeatureImgTypes, feature_img_types, sizeof(int)*m_nFeatureTypeNum);
    }

	m_bMap = map;
	m_bIdentityMap = identity_map;
	if ( m_bIdentityMap )
        m_bMap = true;
	m_bUsedFeatureMean = used_avg;
    m_bMeanCovFromAllExamples = mean_cov_fall_ex;
    m_bNormalization = normalization;

	/* get the actual number of features */
	CBuildFeatureImages* BFI = new CBuildFeatureImages(m_pFeatureImgTypes, m_nFeatureTypeNum);
	m_nTotFeatureNum = BFI->GetFeatureImageNumber();
	delete BFI;

    m_pCascadeSubsetSizes = new int[m_nMaxCascadeLevelNum];
	if ( cascade_subset_sizes ) {
        memcpy(m_pCascadeSubsetSizes, cascade_subset_sizes, sizeof(int)*m_nMaxCascadeLevelNum);
	}
	else {
	    for ( int i = 0 ; i < m_nMaxCascadeLevelNum ; i++ )
            m_pCascadeSubsetSizes[i] = m_nTotFeatureNum;
	}

	m_ppCascadeClassifers = new CCascadeClassifierCovariance*[m_nMaxCascadeLevelNum];

	for ( int i = 0 ; i < m_nMaxCascadeLevelNum ; i++ ) {
		m_ppCascadeClassifers[i] = new CCascadeClassifierCovariance(
						m_pCascadeSubsetSizes?m_pCascadeSubsetSizes[i]:m_nTotFeatureNum,
                        m_nTotFeatureNum, max_weak_classifiers_num, map, used_avg);
	}

	m_pCascadeLevelWeights = new double[m_nMaxCascadeLevelNum];
	ComputeCascadeLevelWeights(level_reject_percent);

	m_bUpdatePositiveWeightScale = false;
	m_dbPositiveExampleWeightScale = 1.0;
	m_dbAdapativePositiveExampleWeightScale = 1.0;

    m_dbBetaFilterPercent = -1.0;

    m_dbMaxFalsePosRate = 0.65;
    m_dbMinPosDetRate = 0.9985;
    m_dbMarginThreshold = 0.2;

	m_pPosDetRates = new double[m_nMaxWeakClassifiersNum];
	m_pFalsePosRates = new double[m_nMaxWeakClassifiersNum];
	for ( int i = 0 ; i < m_nMaxCascadeLevelNum ; i++ ) {
		m_pPosDetRates[i] = 0;
		m_pFalsePosRates[i] = 0;
	}
}

CBoostClassifierCovariance::CBoostClassifierCovariance()
{
    m_ppCascadeClassifers = NULL;
    m_pCascadeLevelWeights = NULL;
	m_pCascadeSubsetSizes = NULL;
    m_pFeatureImgTypes = NULL;

	m_pPosDetRates = NULL;
	m_pFalsePosRates = NULL;

	m_nCascadeLevelNum = 0;
	m_bUpdatePositiveWeightScale = false;
	m_dbPositiveExampleWeightScale = 1.0;
	m_dbAdapativePositiveExampleWeightScale = 1.0;

    m_dbBetaFilterPercent = -1.0;
}

CBoostClassifierCovariance::~CBoostClassifierCovariance()
{
	CleanData();
}

void CBoostClassifierCovariance::SetPositiveExampleWeightScale(double pos_example_weight_scale)
{
    m_dbPositiveExampleWeightScale = pos_example_weight_scale;
    m_dbAdapativePositiveExampleWeightScale = pos_example_weight_scale;
}

void CBoostClassifierCovariance::SetUpdatePositiveWeightScale(bool update_pos_weight_scale)
{
    m_bUpdatePositiveWeightScale = update_pos_weight_scale;
}

void CBoostClassifierCovariance::SetMaxTrainingExampleNumbers(
    int max_pos_examples_num, int max_neg_examples_num,
    int max_sel_pos_examples_num, int max_sel_neg_examples_num)
{
    m_nMaxPosExamplesNum = max_pos_examples_num;
    m_nMaxNegExamplesNum = max_neg_examples_num;

    m_nMaxSelNegExamplesNum = (max_sel_neg_examples_num<=0) ?
                            m_nMaxNegExamplesNum:
                            MIN(m_nMaxNegExamplesNum, max_sel_neg_examples_num);
    m_nMaxSelPosExamplesNum = (max_sel_pos_examples_num<=0) ?
                            m_nMaxPosExamplesNum:
                            MIN(m_nMaxPosExamplesNum, max_sel_pos_examples_num);
}

void CBoostClassifierCovariance::SetTrainingStopConditions(double max_false_pos_rate,
        double min_pos_det_rate, double margin_threshold )
{
    m_dbMinPosDetRate = min_pos_det_rate;
    m_dbMaxFalsePosRate = max_false_pos_rate;
    m_dbMarginThreshold = margin_threshold;
}

void CBoostClassifierCovariance::ComputeAdapativePosExampleWeightScale(
        int n_pos_examples, int n_neg_examples)
{
    m_dbAdapativePositiveExampleWeightScale = m_dbPositiveExampleWeightScale *
            (double)n_neg_examples/(double)n_pos_examples + 1.0;
}

bool CBoostClassifierCovariance::AddCascadeLevel()
{
    if ( m_nCascadeLevelNum >= m_nMaxCascadeLevelNum ) {
        return false;
    }
    else
        m_nCascadeLevelNum++;

    return true;
}

bool CBoostClassifierCovariance::AddWeakClassifier(WeakClassifierCovariance* weak_classifier)
{
    if ( m_nCascadeLevelNum < 1 || m_nCascadeLevelNum > m_nMaxCascadeLevelNum )
        return false;

    return m_ppCascadeClassifers[m_nCascadeLevelNum-1]->AddWeakClassifier(weak_classifier);
}

WeakClassifierCovariance* CBoostClassifierCovariance::
    GetCurrentWeakClassifier()
{
    if ( m_nCascadeLevelNum < 1 || m_nCascadeLevelNum > m_nMaxCascadeLevelNum )
        return NULL;

    return m_ppCascadeClassifers[m_nCascadeLevelNum-1]->m_ppWeakClassifers
        [m_ppCascadeClassifers[m_nCascadeLevelNum-1]->m_nWeakClassifiersNum];
}

bool CBoostClassifierCovariance::Save(const char* file_name)
{
	if ( !file_name ) {
		printf("Please provide file name to save the trained model!\n");
		return false;
	}

	CvFileStorage* fs = cvOpenFileStorage( file_name, 0, CV_STORAGE_WRITE );
	if ( !fs ) {
		printf("Open file error : %s\n", file_name);
		return false;
	}

	CvMat tmp_mat;
	char msg[1024], label_name[1024];

	const char* s_true = "true";
	const char* s_false = "false";

    ///////////////////////////////////////////////////////////////////////////////
	cvWriteComment( fs, "Image features to be used for computing covariances", 0 );
    cvStartWriteStruct( fs, "image_cov_features", CV_NODE_MAP, NULL, cvAttrList(0,0));

	sprintf(msg, "Image features: ");
	for ( int i = 0 ; i < m_nFeatureTypeNum ; i++ ) {
		switch( m_pFeatureImgTypes[i] ) {
		case F_POSITION_XY: sprintf(msg, "%s F_POSITION_XY", msg); break;
		case F_RGB_COLOR: sprintf(msg, "%s F_RGB_COLOR", msg); break;
		case F_HSV_COLOR: sprintf(msg, "%s F_HSV_COLOR", msg); break;
		case F_GRAY_LEVEL: sprintf(msg, "%s F_GRAY_LEVEL", msg); break;
		case F_IMG_GRAD_XY: sprintf(msg, "%s F_IMG_GRAD_XY", msg); break;
		case F_IMG_GRAD_XY_SQRT: sprintf(msg, "%s F_IMG_GRAD_XY_SQRT", msg); break;
		case F_IMG_GRAD_XY2: sprintf(msg, "%s F_IMG_GRAD_XY2", msg); break;
		case F_IMG_GRAD_XY2_SQRT: sprintf(msg, "%s F_IMG_GRAD_XY2_SQRT", msg); break;
		case F_IMG_EDGE_XY_ORI: sprintf(msg, "%s F_IMG_EDGE_XY_ORI", msg); break;
		case F_IMG_EDGE_XY2_ORI: sprintf(msg, "%s F_IMG_EDGE_XY2_ORI", msg); break;
		case F_IMG_LBP: sprintf(msg, "%s F_IMG_LBP", msg); break;
		case F_IMG_CANNY_EDGE: sprintf(msg, "%s F_IMG_CANNY_EDGE", msg); break;
		case F_FG_PROB_MAP: sprintf(msg, "%s F_FG_PROB_MAP", msg); break;
		case F_FG_BINARY_MAP: sprintf(msg, "%s F_FG_BINARY_MAP", msg); break;
		case F_MAP_GRAD_XY: sprintf(msg, "%s F_MAP_GRAD_XY", msg); break;
		case F_MAP_GRAD_XY_SQRT: sprintf(msg, "%s F_MAP_GRAD_XY_SQRT", msg); break;
		case F_MAP_GRAD_XY2: sprintf(msg, "%s F_MAP_GRAD_XY2", msg); break;
		case F_MAP_GRAD_XY2_SQRT: sprintf(msg, "%s F_MAP_GRAD_XY2_SQRT", msg); break;
		case F_MAP_EDGE_XY_ORI: sprintf(msg, "%s F_MAP_EDGE_XY_ORI", msg); break;
		case F_MAP_EDGE_XY2_ORI: sprintf(msg, "%s F_MAP_EDGE_XY2_ORI", msg); break;
		case F_BIN_GRAD_XY: sprintf(msg, "%s F_BIN_GRAD_XY", msg); break;
		case F_BIN_GRAD_XY_SQRT: sprintf(msg, "%s F_BIN_GRAD_XY_SQRT", msg); break;
		case F_BIN_EDGE_XY_ORI: sprintf(msg, "%s F_BIN_EDGE_XY_ORI", msg); break;
		}
		if ( i < m_nFeatureTypeNum-1 )
            sprintf(msg, "%s, ", msg);
	}
	cvWriteComment( fs, msg, 0 );

	cvWriteInt( fs, "feature_types_num", m_nFeatureTypeNum );

	tmp_mat = cvMat( 1, m_nFeatureTypeNum, CV_32S, m_pFeatureImgTypes );
	cvWrite( fs, "feature_types", &tmp_mat );

	cvWriteComment( fs, "Total number of features to be used for training and detection", 0 );
	cvWriteInt( fs, "total_feature_num", m_nTotFeatureNum );

    cvEndWriteStruct( fs );

    ///////////////////////////////////////////////////////////////////////////////
    cvWriteComment( fs, "Logit-Boost training parameters", 0 );
    cvStartWriteStruct( fs, "logit_boost_paras", CV_NODE_MAP, NULL, cvAttrList(0,0));

	cvWriteComment( fs, "Maximal cascade level number", 0 );
	cvWriteInt( fs, "max_cascade_level_num", m_nMaxCascadeLevelNum );

	cvWriteComment( fs, "Maximal weak classifiers at each cascade level", 0 );
	cvWriteInt( fs, "max_weak_classifiers_num", m_nMaxWeakClassifiersNum );

	cvWriteComment( fs, "Maximal number of positive examples to be used at each cascade level", 0 );
	cvWriteInt( fs, "max_pos_examples_num", m_nMaxPosExamplesNum );

	cvWriteComment( fs, "Maximal number of negative examples to be used at each cascade level", 0 );
	cvWriteInt( fs, "max_neg_examples_num", m_nMaxNegExamplesNum );

	cvWriteComment( fs, "Maximal number of selected positive examples to be used at each cascade level", 0 );
	cvWriteInt( fs, "max_sel_pos_examples_num", m_nMaxSelPosExamplesNum );

	cvWriteComment( fs, "Maximal number of selected negative examples to be used at each cascade level", 0 );
	cvWriteInt( fs, "max_sel_neg_examples_num", m_nMaxSelNegExamplesNum );

    bool full_features = true;
    for ( int i = 0 ; i < m_nMaxCascadeLevelNum ; i++ )
        full_features &= m_pCascadeSubsetSizes?m_pCascadeSubsetSizes[i]==m_nTotFeatureNum:true;

	cvWriteComment( fs, "Do the subset features used?", 0 );
	cvWriteString( fs, "used_subset_features", full_features ? s_false : s_true );
	if ( !full_features ) {
		cvWriteComment( fs, "Subset feature size at all the cascade levels", 0 );
		tmp_mat = cvMat( 1, m_nMaxCascadeLevelNum, CV_32S, m_pCascadeSubsetSizes );
		cvWrite( fs, "cascade_subset_sizes", &tmp_mat );
	}

	cvWriteComment( fs, "Positive example weight for unbalanced training", 0 );
	cvWriteReal( fs, "pos_example_weight_scale", m_dbPositiveExampleWeightScale );

	cvWriteComment( fs, "Do we update positive example weight based on the same weights sum of both positive/negative examples?", 0 );
	cvWriteString( fs, "update_pos_weight_scale", m_bUpdatePositiveWeightScale ? s_true : s_false );

	cvWriteComment( fs, "Is it mapped to tangent space?", 0 );
	cvWriteString( fs, "map_tangent_space", m_bMap ? s_true : s_false );

	cvWriteComment( fs, "Is it mapped to tangent space based on identity matrix or mean covariance matrix?", 0 );
	cvWriteString( fs, "cov_identity_map", m_bIdentityMap ? s_true : s_false );

	cvWriteComment( fs, "Does it use the mean features?", 0 );
	cvWriteString( fs, "used_feature_mean", m_bUsedFeatureMean ? s_true : s_false );

	cvWriteComment( fs, "Will it be normalized based on the parent window?", 0 );
	cvWriteString( fs, "normalization", m_bNormalization ? s_true : s_false );

	cvWriteComment( fs, "Do we compute mean covariance matrix for mapping from all the examples or only positive examples", 0 );
	cvWriteString( fs, "mean_cov_from_all_examples", m_bMeanCovFromAllExamples ? s_true : s_false );

	cvWriteComment( fs, "False example reject percent at each cascade level", 0 );
	cvWriteReal( fs, "level_reject_percent", m_dbLevelRejectPercent );

	cvWriteComment( fs, "To avoid boundary effect", 0 );
	cvWriteComment( fs, "Margin size percent for positive image examples w.r.t. their image sizes", 0 );
	cvWriteReal( fs, "pos_ex_edge_margin_x_percent", m_dbPosExEdgeMarginXPercent );
	cvWriteReal( fs, "pos_ex_edge_margin_y_percent", m_dbPosExEdgeMarginYPercent );

	cvWriteComment( fs, "the edge shift percent of positive image examples w.r.t. their image sizes", 0 );
	cvWriteReal( fs, "pos_ex_shift_percent", m_dbPosExShiftPercent );

	cvWriteComment( fs, "the iteration of randomly selecting positive examples from positive cropped images", 0 );
	cvWriteInt( fs, "pos_ex_shift_itr_num", m_nPosExShiftItrNum );

	cvWriteComment( fs, "Margin size (in pixels) of saved negative cropped image examples", 0 );
	cvWriteInt( fs, "neg_cropped_edge_margin_size", m_nNegCroppedEdgeMarginSize );

	cvWriteComment( fs, "Filtering percent of examples for next weak classifier training", 0 );
	cvWriteReal( fs, "beta_filter_percent", m_dbBetaFilterPercent );

    cvEndWriteStruct( fs );

    ///////////////////////////////////////////////////////////////////////////////
	cvWriteComment( fs, "Logit-boost training stop constraints", 0 );
    cvStartWriteStruct( fs, "logit_boost_stop_constraints", CV_NODE_MAP, NULL, cvAttrList(0,0));

	cvWriteComment( fs, "Minimal positive detection rate for each cascade level", 0 );
	cvWriteReal( fs, "min_pos_det_rate", m_dbMinPosDetRate );

	cvWriteComment( fs, "Maximal false positive rate for each cascade level", 0 );
	cvWriteReal( fs, "max_false_pos_rate", m_dbMaxFalsePosRate );

	cvWriteComment( fs, "Margin threshold for each cascade level", 0 );
	cvWriteReal( fs, "margin_threshold", m_dbMarginThreshold );

    cvEndWriteStruct( fs );

    ///////////////////////////////////////////////////////////////////////////////
	cvWriteComment( fs, "Learned cascade level number", 0 );
	cvWriteInt( fs, "cascade_level_num",  m_nCascadeLevelNum );

    ///////////////////////////////////////////////////////////////////////////////
	cvStartWriteStruct( fs, "stat_infos", CV_NODE_MAP, NULL, cvAttrList(0,0));

	int *weak_classifier_nums = new int[m_nCascadeLevelNum];
	for ( int i = 0 ; i < m_nCascadeLevelNum ; i++ )
		weak_classifier_nums[i] = m_ppCascadeClassifers[i]->m_nWeakClassifiersNum;
	tmp_mat = cvMat( 1, m_nCascadeLevelNum, CV_32S, weak_classifier_nums );
	cvWriteComment( fs, "Numbers of weak classifiers at all trained cascade levels", 0 );
	cvWrite( fs, "cascade_weak_classifier_nums", &tmp_mat );
	delete [] weak_classifier_nums;

	tmp_mat = cvMat( 1, m_nMaxCascadeLevelNum, CV_64F, m_pPosDetRates );
	cvWriteComment( fs, "Detection rates of positive examples at all cascade levels", 0 );
	cvWrite( fs, "cascade_pos_det_rates", &tmp_mat );

	tmp_mat = cvMat( 1, m_nMaxCascadeLevelNum, CV_64F, m_pFalsePosRates );
	cvWriteComment( fs, "False positive rates while generating negative cropped image for training at all cascade levels", 0 );
	cvWrite( fs, "cascade_false_pos_rates", &tmp_mat );

    cvEndWriteStruct( fs );

	for ( int i = 0 ; i < m_nCascadeLevelNum ; i++ ) {
        ///////////////////////////////////////////////////////////////////////////////
		sprintf(label_name, "cascade_%d_classifiers", i);
		cvStartWriteStruct( fs, label_name, CV_NODE_MAP, NULL, cvAttrList(0,0));

		sprintf(msg, "The %d-th cascade level : %d weak classifiers", i, m_ppCascadeClassifers[i]->m_nWeakClassifiersNum );
		cvWriteComment( fs, msg, 0 );

		cvWriteInt( fs, "weak_classifiers_num", m_ppCascadeClassifers[i]->m_nWeakClassifiersNum );

		cvWriteComment( fs, "Thresholds for current logit-boost cascade level", 0 );

		cvWriteComment( fs, "positive probability threshold of Xn", 0 );
		cvWriteReal( fs, "Xn_pos_prob_threshold", m_ppCascadeClassifers[i]->m_dbXnPosProbThreshold );

		cvWriteComment( fs, "classifier function threshold of Xn", 0 );
		cvWriteReal( fs, "Xn_classifier_fun_threshold", m_ppCascadeClassifers[i]->m_dbXnClassifierFunThreshold );

		cvWriteComment( fs, "positive probability threshold of Xp", 0 );
		cvWriteReal( fs, "Xp_pos_prob_threshold", m_ppCascadeClassifers[i]->m_dbXpPosProbThreshold );

		cvWriteComment( fs, "classifier function threshold of Xp", 0 );
		cvWriteReal( fs, "Xp_classifier_fun_threshold", m_ppCascadeClassifers[i]->m_dbXpClassifierFunThreshold );

		for ( int j = 0 ; j < m_ppCascadeClassifers[i]->m_nWeakClassifiersNum ; j++ ) {

			WeakClassifierCovariance* wcc = m_ppCascadeClassifers[i]->m_ppWeakClassifers[j];

            ///////////////////////////////////////////////////////////////////////////////
			sprintf(label_name, "weak_classifier_%d", j);
			cvStartWriteStruct( fs, label_name, CV_NODE_MAP, NULL, cvAttrList(0,0));

			sprintf(msg, "Information in the %d-th weak classifier at the %d-th cascade level", j, i );
			cvWriteComment( fs, msg, 0 );

            cvWriteComment( fs, "dimension number of covariance features", 0 );
			cvWriteInt( fs, "n_features", wcc->n_features );

            cvWriteComment( fs, "indices of covariance features", 0 );
			tmp_mat = cvMat( 1, wcc->n_features, CV_32S, wcc->feature_indices );
			cvWrite( fs, "feature_indices", &tmp_mat );

			cvWriteComment( fs, "normalized float rectangular sub-window [x, y, width, height]", 0 );
			cvWriteReal( fs, "roi_x", wcc->roi.x );
			cvWriteReal( fs, "roi_y", wcc->roi.y );
			cvWriteReal( fs, "roi_width", wcc->roi.width );
			cvWriteReal( fs, "roi_height", wcc->roi.height );

			cvWriteComment( fs, "learned fitting function parameters", 0 );
			tmp_mat = cvMat( 1, wcc->feature_size+1, CV_64F, wcc->fitting_params );
			cvWrite( fs, "fitting_params", &tmp_mat );

			if ( m_bMap && !m_bIdentityMap ) {
			    cvWriteComment( fs, "weighted mean covariance matrix used for mapping", 0 );

			    cvWriteComment( fs, "upper triangular data of weighted mean covariance", 0 );
				tmp_mat = cvMat( 1, wcc->weighted_mean->SM_length, CV_64F, wcc->weighted_mean->SM_ptr );
				cvWrite( fs, "weighted_mean_SM_ptr", &tmp_mat );
			}

			cvEndWriteStruct( fs );
		}
		cvEndWriteStruct( fs );
	}

	cvReleaseFileStorage( &fs );

	printf("\nExport %d logit-boost classifiers to %s\n", m_nCascadeLevelNum, file_name);

	return true;
}

bool CBoostClassifierCovariance::Load(const char* file_name)
{
	if ( !file_name ) {
		printf("Please provide file name to load the trained model!\n");
		return false;
	}

	CvFileStorage* fs = cvOpenFileStorage( file_name, 0, CV_STORAGE_READ );
	if ( !fs ) {
		//printf("Open file error : %s\n", file_name);
		return false;
	}

	// Clean the data
	CleanData();

	CvMat* tmp_mat;
	char label_name[1024];

    //////////////////////////////////////////////////////////////////////////////
	CvFileNode* cov_features_node = cvGetFileNodeByName( fs, NULL, "image_cov_features" );

	m_nFeatureTypeNum = cvReadIntByName( fs, cov_features_node, "feature_types_num" );
	tmp_mat = (CvMat*)cvReadByName( fs, cov_features_node, "feature_types" );
	m_pFeatureImgTypes = new int[m_nFeatureTypeNum];
	memcpy(m_pFeatureImgTypes, tmp_mat->data.ptr, sizeof(int)*m_nFeatureTypeNum);
	cvReleaseMat(&tmp_mat);

	m_nTotFeatureNum = cvReadIntByName( fs, cov_features_node, "total_feature_num" );

    //////////////////////////////////////////////////////////////////////////////
	CvFileNode* logit_boost_paras_node = cvGetFileNodeByName( fs, NULL, "logit_boost_paras" );

	m_nMaxCascadeLevelNum = cvReadIntByName( fs, logit_boost_paras_node, "max_cascade_level_num" );
	m_nMaxWeakClassifiersNum = cvReadIntByName( fs, logit_boost_paras_node, "max_weak_classifiers_num" );

	m_nMaxPosExamplesNum = cvReadIntByName( fs, logit_boost_paras_node, "max_pos_examples_num" );
	m_nMaxNegExamplesNum = cvReadIntByName( fs, logit_boost_paras_node, "max_neg_examples_num" );

	m_nMaxSelPosExamplesNum = cvReadIntByName( fs, logit_boost_paras_node, "max_sel_pos_examples_num" );
	m_nMaxSelNegExamplesNum = cvReadIntByName( fs, logit_boost_paras_node, "max_sel_neg_examples_num" );

    bool used_subset_features = bool_check( cvReadStringByName( fs, logit_boost_paras_node, "used_subset_features" ) );

    m_pCascadeSubsetSizes = NULL;
	if ( used_subset_features ) {
        tmp_mat = (CvMat*)cvReadByName( fs, logit_boost_paras_node, "cascade_subset_sizes" );
        if ( tmp_mat ) {
            m_pCascadeSubsetSizes = new int[m_nMaxCascadeLevelNum];
            memcpy(m_pCascadeSubsetSizes, tmp_mat->data.ptr, sizeof(int)*m_nMaxCascadeLevelNum);
            cvReleaseMat(&tmp_mat);
        }
        else {
            printf("Loading subset sizes error!\n");
            exit(1);
        }
	}
	else {
        m_pCascadeSubsetSizes = new int[m_nMaxCascadeLevelNum];
        for ( int i = 0 ; i < m_nMaxCascadeLevelNum ; i++ )
            m_pCascadeSubsetSizes[i] = m_nTotFeatureNum;
	}

	m_dbPositiveExampleWeightScale = cvReadRealByName( fs, logit_boost_paras_node, "pos_example_weight_scale" );

	m_bUpdatePositiveWeightScale = bool_check( cvReadStringByName( fs, logit_boost_paras_node, "update_pos_weight_scale" ) );

    m_bMap = bool_check( cvReadStringByName( fs, logit_boost_paras_node, "map_tangent_space" ) );
    m_bIdentityMap = bool_check( cvReadStringByName( fs, logit_boost_paras_node, "cov_identity_map" ) );
    if ( m_bIdentityMap )
        m_bMap = true;

    m_bUsedFeatureMean = bool_check( cvReadStringByName( fs, logit_boost_paras_node, "used_feature_mean" ) );
    m_bNormalization = bool_check( cvReadStringByName( fs, logit_boost_paras_node, "normalization" ) );
	m_bMeanCovFromAllExamples = bool_check( cvReadStringByName( fs, logit_boost_paras_node, "mean_cov_from_all_examples" ) );

    m_dbLevelRejectPercent = cvReadRealByName( fs, logit_boost_paras_node, "level_reject_percent" );

    m_dbPosExEdgeMarginXPercent = cvReadRealByName( fs, logit_boost_paras_node, "pos_ex_edge_margin_x_percent" );
    m_dbPosExEdgeMarginYPercent = cvReadRealByName( fs, logit_boost_paras_node, "pos_ex_edge_margin_y_percent" );

    m_dbPosExShiftPercent = cvReadRealByName( fs, logit_boost_paras_node, "pos_ex_shift_percent" );
    m_nPosExShiftItrNum = cvReadIntByName( fs, logit_boost_paras_node, "pos_ex_shift_itr_num" );

    m_nNegCroppedEdgeMarginSize = cvReadIntByName( fs, logit_boost_paras_node, "neg_cropped_edge_margin_size" );

    m_dbBetaFilterPercent = cvReadRealByName( fs, logit_boost_paras_node, "beta_filter_percent" );

    //////////////////////////////////////////////////////////////////////////////
	CvFileNode* logit_boost_stop_node = cvGetFileNodeByName( fs, NULL, "logit_boost_stop_constraints" );
    m_dbMinPosDetRate = cvReadRealByName( fs, logit_boost_stop_node, "min_pos_det_rate" );
    m_dbMaxFalsePosRate = cvReadRealByName( fs, logit_boost_stop_node, "max_false_pos_rate" );
    m_dbMarginThreshold = cvReadRealByName( fs, logit_boost_stop_node, "margin_threshold" );

    //////////////////////////////////////////////////////////////////////////////
	m_nCascadeLevelNum = cvReadIntByName( fs, NULL, "cascade_level_num" );

	m_pCascadeLevelWeights = new double[m_nMaxCascadeLevelNum];
	ComputeCascadeLevelWeights(m_dbLevelRejectPercent);

    ///////////////////////////////////////////////////////////////////////////////
	CvFileNode* stat_infos_node = cvGetFileNodeByName( fs, NULL, "stat_infos" );
	if ( stat_infos_node ) {
        m_pPosDetRates = new double[m_nMaxCascadeLevelNum];
        tmp_mat = (CvMat*)cvReadByName( fs, stat_infos_node, "cascade_pos_det_rates" );
        if ( tmp_mat ) {
            memcpy(m_pPosDetRates, tmp_mat->data.ptr, sizeof(double)*m_nMaxCascadeLevelNum);
            cvReleaseMat(&tmp_mat);
        }

        m_pFalsePosRates = new double[m_nMaxCascadeLevelNum];
        tmp_mat = (CvMat*)cvReadByName( fs, stat_infos_node, "cascade_false_pos_rates" );
        if ( tmp_mat ) {
            memcpy(m_pFalsePosRates, tmp_mat->data.ptr, sizeof(double)*m_nMaxCascadeLevelNum);
            cvReleaseMat(&tmp_mat);
        }
	}

    // allocate memories
	m_ppCascadeClassifers = new CCascadeClassifierCovariance*[m_nMaxCascadeLevelNum];
	for ( int i = 0 ; i < m_nMaxCascadeLevelNum ; i++ ) {
		m_ppCascadeClassifers[i] = new CCascadeClassifierCovariance(m_pCascadeSubsetSizes?m_pCascadeSubsetSizes[i]:m_nTotFeatureNum,
                        m_nTotFeatureNum, m_nMaxWeakClassifiersNum, m_bMap, m_bUsedFeatureMean);
	}

	for ( int i = 0 ; i < m_nCascadeLevelNum ; i++ ) {
        //////////////////////////////////////////////////////////////////////////////
		sprintf(label_name, "cascade_%d_classifiers", i);
	    CvFileNode* cascade_node = cvGetFileNodeByName( fs, NULL, label_name );

	    m_ppCascadeClassifers[i]->m_nWeakClassifiersNum = cvReadIntByName( fs, cascade_node, "weak_classifiers_num" );
	    m_ppCascadeClassifers[i]->m_dbXnPosProbThreshold = cvReadRealByName( fs, cascade_node, "Xn_pos_prob_threshold" );
	    m_ppCascadeClassifers[i]->m_dbXnClassifierFunThreshold = cvReadRealByName( fs, cascade_node, "Xn_classifier_fun_threshold" );
	    m_ppCascadeClassifers[i]->m_dbXpPosProbThreshold = cvReadRealByName( fs, cascade_node, "Xp_pos_prob_threshold" );
	    m_ppCascadeClassifers[i]->m_dbXpClassifierFunThreshold = cvReadRealByName( fs, cascade_node, "Xp_classifier_fun_threshold" );

		for ( int j = 0 ; j < m_ppCascadeClassifers[i]->m_nWeakClassifiersNum ; j++ ) {
            int feature_indices[100];
            int n_features;
            floatRect roi;

            //////////////////////////////////////////////////////////////////////////////
			sprintf(label_name, "weak_classifier_%d", j);
			CvFileNode* weak_node = cvGetFileNodeByName( fs, cascade_node, label_name );

			n_features = cvReadIntByName( fs, weak_node, "n_features" );
			tmp_mat = (CvMat*)cvReadByName( fs, weak_node, "feature_indices" );
			memcpy(feature_indices, tmp_mat->data.ptr, sizeof(int)*n_features);
			cvReleaseMat(&tmp_mat);

			roi.x = (float)cvReadRealByName( fs, weak_node, "roi_x" );
			roi.y = (float)cvReadRealByName( fs, weak_node, "roi_y" );
			roi.width = (float)cvReadRealByName( fs, weak_node, "roi_width" );
			roi.height = (float)cvReadRealByName( fs, weak_node, "roi_height" );

            m_ppCascadeClassifers[i]->m_ppWeakClassifers[j] = new WeakClassifierCovariance(
                roi.x, roi.y, roi.width, roi.height,
                m_pCascadeSubsetSizes?m_pCascadeSubsetSizes[i]:m_nTotFeatureNum,
                m_nTotFeatureNum,
                feature_indices,
                m_bMap, m_bUsedFeatureMean);
			WeakClassifierCovariance* wcc = m_ppCascadeClassifers[i]->m_ppWeakClassifers[j];
            wcc->setPositiveExampleWeightScale(m_dbAdapativePositiveExampleWeightScale);

			tmp_mat = (CvMat*)cvReadByName( fs, weak_node, "fitting_params" );
			memcpy(wcc->fitting_params, tmp_mat->data.ptr, sizeof(double)*(wcc->feature_size+1));
			cvReleaseMat(&tmp_mat);

			if ( m_bMap && !m_bIdentityMap ) {
				tmp_mat = (CvMat*)cvReadByName( fs, weak_node, "weighted_mean_SM_ptr" );
				memcpy(wcc->weighted_mean->SM_ptr, tmp_mat->data.ptr, sizeof(double)*wcc->weighted_mean->SM_length);
				cvReleaseMat(&tmp_mat);

				wcc->weighted_mean->Sqrt(wcc->weighted_mean_sqrt);
				wcc->weighted_mean_sqrt->Inverse(wcc->weighted_mean_inv_sqrt);
			}
			if ( m_bIdentityMap ) {
				wcc->weighted_mean->SetIdentityMatrix(1.0);
				wcc->weighted_mean_sqrt->SetIdentityMatrix(1.0);
				wcc->weighted_mean_inv_sqrt->SetIdentityMatrix(1.0);
			}
		}
	}

	cvReleaseFileStorage( &fs );

	printf("\nImport %d logit-boost classifiers from %s\n", m_nCascadeLevelNum, file_name);

	return true;
}

void CBoostClassifierCovariance::CleanData()
{
    if ( m_ppCascadeClassifers ) {
        for ( int i = 0 ; i < m_nMaxCascadeLevelNum ; i++ ) {
            delete m_ppCascadeClassifers[i];
        }
        delete [] m_ppCascadeClassifers;
        m_ppCascadeClassifers = NULL;
    }

    if ( m_pCascadeLevelWeights )
        delete [] m_pCascadeLevelWeights;
	m_pCascadeLevelWeights = NULL;

	if ( m_pCascadeSubsetSizes )
        delete [] m_pCascadeSubsetSizes;
    m_pCascadeSubsetSizes = NULL;

    if ( m_pFeatureImgTypes )
        delete [] m_pFeatureImgTypes;
    m_pFeatureImgTypes = NULL;

	if ( m_pPosDetRates )
		delete [] m_pPosDetRates;
	if ( m_pFalsePosRates )
		delete [] m_pFalsePosRates;
}

void CBoostClassifierCovariance::ComputeCascadeLevelWeights(double level_reject_percent)
{
	double tot_weight = 0;
	int i;
	for ( i = 0 ; i < m_nMaxCascadeLevelNum ; i++ ) {
		m_pCascadeLevelWeights[i] = pow(level_reject_percent, m_nMaxCascadeLevelNum-i+1);
		tot_weight += m_pCascadeLevelWeights[i];
	}
	for ( i = 0 ; i < m_nMaxCascadeLevelNum ; i++ )
		m_pCascadeLevelWeights[i] /= tot_weight;
}

WeakClassifierCovariance* CBoostClassifierCovariance::
    CreateWeakClassifier(floatRect roi, int n_features, int *feature_indices)
{
    if ( m_nCascadeLevelNum < 1 || m_nCascadeLevelNum > m_nMaxCascadeLevelNum )
        return NULL;

    WeakClassifierCovariance* weak_classifier = new WeakClassifierCovariance(
        roi.x, roi.y, roi.width, roi.height,
        n_features,
        m_nTotFeatureNum,
        feature_indices,
        m_bMap, m_bUsedFeatureMean);

    weak_classifier->setPositiveExampleWeightScale(m_dbAdapativePositiveExampleWeightScale);

    return weak_classifier;
}

void CBoostClassifierCovariance::ExportLogMessage(const char* log_file_name, char *msg)
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

double CBoostClassifierCovariance::Predictor(CIntegralRegionCov *irc, CvRect roi,
				int boost_classifiers_num,
				double Xn_prob_alpha,
				int *final_pos_level, double *cascade_level_probs, int neg_obj_id)
{
	if ( boost_classifiers_num == 0 || m_nCascadeLevelNum == 0 )
		return 1.0;

    if ( boost_classifiers_num < 0 )
        boost_classifiers_num = m_nCascadeLevelNum;
    boost_classifiers_num = MIN(boost_classifiers_num, m_nCascadeLevelNum);

	double F_X = 0;
	int mat_size = irc->m_nCovDim;
	CvRect sub_int_win;
	struct floatRect sub_float_win;
	CCovarianceMatrix* reg_cov = new CCovarianceMatrix(mat_size);

	double parent_diag_elems[20];
	double sub_parent_diag_elems[20];

	double parent_avg_elems[20];
	double sub_parent_avg_elems[20];

	if ( m_bNormalization )
		irc->GetRegionCovNormDiagonalElems(roi, parent_diag_elems, parent_avg_elems, PARENT_MIN_IDENTITY_VALUE);

	double tot_p_X = 0.0;
	double e_F_X, p_X = 0;
	int i;

	WeakClassifierCovariance* wcc;

	double Xp_prob_alpha = 1.0 - Xn_prob_alpha;

	for ( int k = 0 ; k < boost_classifiers_num ; k++ ) {
		F_X = 0;
		for ( int l = 0 ; l < m_ppCascadeClassifers[k]->m_nWeakClassifiersNum ; l++ ) {
			wcc = m_ppCascadeClassifers[k]->m_ppWeakClassifers[l];
			sub_float_win = wcc->roi;
			sub_int_win = m_cRect.GetIntRect(sub_float_win, roi);

			double f_l;
			reg_cov->Resize(wcc->n_features);

			if ( wcc->n_features == mat_size )
                irc->GetRegionCovFeature(sub_int_win, reg_cov, MIN_IDENTITY_VALUE, m_bUsedFeatureMean);
                //irc->GetRegionCovFeature(sub_int_win, reg_cov, MIN_IDENTITY_VALUE);
            else
                irc->GetSubRegionCovFeature(sub_int_win, reg_cov, wcc->feature_indices,
                    NULL, MIN_IDENTITY_VALUE, m_bUsedFeatureMean);

			if ( m_bNormalization ) {
				if ( wcc->n_features < mat_size ) {
					for ( i = 0 ; i < wcc->n_features ; i++ ) {
						sub_parent_diag_elems[i] = parent_diag_elems[wcc->feature_indices[i]];
						if ( m_bUsedFeatureMean )
							sub_parent_avg_elems[i] = parent_avg_elems[wcc->feature_indices[i]];
					}
					reg_cov->NormalizeIllumination(sub_parent_diag_elems, m_bUsedFeatureMean?sub_parent_avg_elems:NULL);
				}
				else
					reg_cov->NormalizeIllumination(parent_diag_elems, m_bUsedFeatureMean?parent_avg_elems:NULL);
			}
			else {
				int x_feature_idx = wcc->feature_indices[0]==0?0:-1;
				int y_feature_idx = wcc->feature_indices[1]==1?1:-1;
				reg_cov->NormalizePositionFeatures(roi, x_feature_idx, y_feature_idx);
			}
			f_l = wcc->TestCov(reg_cov);
			F_X += 0.5 * f_l;
/*
			if ( neg_obj_id == 0 ) {
			    char log_msg[1024];
                sprintf(log_msg, "best-sub-win : [%.3f, %.3f, %.3f, %.3f]\n",
                    wcc->roi.x,
                    wcc->roi.y,
                    wcc->roi.width,
                    wcc->roi.height );
                ExportLogMessage("debug2.log", log_msg);
			    reg_cov->Print("debug2.log");
			}
*/
		}

        e_F_X = exp(F_X);
        p_X = e_F_X / ( e_F_X + 1.0/e_F_X );

        if ( boost_classifiers_num < m_nMaxCascadeLevelNum )
            tot_p_X += m_pCascadeLevelWeights[k]/m_pCascadeLevelWeights[boost_classifiers_num]*p_X;
        else
            tot_p_X += m_pCascadeLevelWeights[k]*p_X;

        if ( cascade_level_probs )
            cascade_level_probs[k] = tot_p_X;
        if ( p_X < Xn_prob_alpha*m_ppCascadeClassifers[k]->m_dbXnPosProbThreshold +
            Xp_prob_alpha*m_ppCascadeClassifers[k]->m_dbXpPosProbThreshold ) {
            if ( final_pos_level )
                *final_pos_level = k - 1;
            tot_p_X = 0;
            if ( cascade_level_probs )
                cascade_level_probs[k] = tot_p_X;
            delete reg_cov;

            /* debug output */
            if ( neg_obj_id >= 0 ) {
            char log_msg[1024];
            sprintf(log_msg, "%05d  :  F_X  %E\t P_X  %E\n",
                    neg_obj_id, F_X, p_X);
            ExportLogMessage("debug2.log", log_msg);
            }

            return tot_p_X;
        }
	}

	delete reg_cov;

	if ( final_pos_level )
		*final_pos_level = boost_classifiers_num-1;

    /* debug output */

	if ( neg_obj_id >= 0 ) {
    char log_msg[1024];
    sprintf(log_msg, "%05d  :  F_X  %E\t P_X  %E\n",
            neg_obj_id, F_X, p_X);
    ExportLogMessage("debug2.log", log_msg);
	}

	return tot_p_X;
}

bool CBoostClassifierCovariance::CascadeLevelPredictor(CIntegralRegionCov *irc, CvRect roi,
                int start_cascade_level,
				int end_cascade_level,
				double Xn_prob_alpha)
{
    start_cascade_level = MAX(0, start_cascade_level);
    end_cascade_level = MIN(m_nCascadeLevelNum, end_cascade_level);

	if ( start_cascade_level < 0 || end_cascade_level <= 0 ||
        end_cascade_level <= start_cascade_level )
		return true;

	double F_X = 0;
	int mat_size = irc->m_nCovDim;
	CvRect sub_int_win;
	struct floatRect sub_float_win;
	CCovarianceMatrix* reg_cov = new CCovarianceMatrix(mat_size);

	double parent_diag_elems[20];
	double sub_parent_diag_elems[20];

	double parent_avg_elems[20];
	double sub_parent_avg_elems[20];

	if ( m_bNormalization )
		irc->GetRegionCovNormDiagonalElems(roi, parent_diag_elems, parent_avg_elems, PARENT_MIN_IDENTITY_VALUE);

	double e_F_X, p_X = 0;
	int i;

	WeakClassifierCovariance* wcc;

	double Xp_prob_alpha = 1.0 - Xn_prob_alpha;

	for ( int k = start_cascade_level ; k < end_cascade_level ; k++ ) {
		F_X = 0;
		for ( int l = 0 ; l < m_ppCascadeClassifers[k]->m_nWeakClassifiersNum ; l++ ) {
			wcc = m_ppCascadeClassifers[k]->m_ppWeakClassifers[l];
			sub_float_win = wcc->roi;
			sub_int_win = m_cRect.GetIntRect(sub_float_win, roi);

			double f_l;
			reg_cov->Resize(wcc->n_features);

			if ( wcc->n_features == mat_size )
                irc->GetRegionCovFeature(sub_int_win, reg_cov, MIN_IDENTITY_VALUE, m_bUsedFeatureMean);
                //irc->GetRegionCovFeature(sub_int_win, reg_cov, MIN_IDENTITY_VALUE);
            else
                irc->GetSubRegionCovFeature(sub_int_win, reg_cov, wcc->feature_indices,
                    NULL, MIN_IDENTITY_VALUE, m_bUsedFeatureMean);

			if ( m_bNormalization ) {
				if ( wcc->n_features < mat_size ) {
					for ( i = 0 ; i < wcc->n_features ; i++ ) {
						sub_parent_diag_elems[i] = parent_diag_elems[wcc->feature_indices[i]];
						if ( m_bUsedFeatureMean )
							sub_parent_avg_elems[i] = parent_avg_elems[wcc->feature_indices[i]];
					}
					reg_cov->NormalizeIllumination(sub_parent_diag_elems, m_bUsedFeatureMean?sub_parent_avg_elems:NULL);
				}
				else
					reg_cov->NormalizeIllumination(parent_diag_elems, m_bUsedFeatureMean?parent_avg_elems:NULL);
			}
			else {
				int x_feature_idx = wcc->feature_indices[0]==0?0:-1;
				int y_feature_idx = wcc->feature_indices[1]==1?1:-1;
				reg_cov->NormalizePositionFeatures(roi, x_feature_idx, y_feature_idx);
			}
			f_l = wcc->TestCov(reg_cov);
			F_X += 0.5 * f_l;
		}

        e_F_X = exp(F_X);
        p_X = e_F_X / ( e_F_X + 1.0/e_F_X );

        if ( p_X < Xn_prob_alpha*m_ppCascadeClassifers[k]->m_dbXnPosProbThreshold +
            Xp_prob_alpha*m_ppCascadeClassifers[k]->m_dbXpPosProbThreshold ) {
            delete reg_cov;
            return false;
        }
	}

	delete reg_cov;

	return true;
}


double CBoostClassifierCovariance::GetCascadeClassifierFun(
                CIntegralRegionCov *irc, CvRect roi,
				int boost_classifiers_idx)
{
	if ( boost_classifiers_idx < 0 || m_nCascadeLevelNum == 0 )
		return 1.0;

    boost_classifiers_idx = MIN(boost_classifiers_idx, m_nCascadeLevelNum-1);

	double F_X = 0;
	int mat_size = irc->m_nCovDim;
	CvRect sub_int_win;
	struct floatRect sub_float_win;
	CCovarianceMatrix* reg_cov = new CCovarianceMatrix(mat_size);

	double parent_diag_elems[20];
	double sub_parent_diag_elems[20];

	double parent_avg_elems[20];
	double sub_parent_avg_elems[20];

	if ( m_bNormalization )
		irc->GetRegionCovNormDiagonalElems(roi, parent_diag_elems, parent_avg_elems, PARENT_MIN_IDENTITY_VALUE);

	int i;

	WeakClassifierCovariance* wcc;

    int k = boost_classifiers_idx;
    F_X = 0;
    for ( int l = 0 ; l < m_ppCascadeClassifers[k]->m_nWeakClassifiersNum ; l++ ) {
        wcc = m_ppCascadeClassifers[k]->m_ppWeakClassifers[l];
        sub_float_win = wcc->roi;
        sub_int_win = m_cRect.GetIntRect(sub_float_win, roi);

        double f_l;
        reg_cov->Resize(wcc->n_features);

        if ( wcc->n_features == mat_size )
            irc->GetRegionCovFeature(sub_int_win, reg_cov, MIN_IDENTITY_VALUE, m_bUsedFeatureMean);
            //irc->GetRegionCovFeature(sub_int_win, reg_cov, MIN_IDENTITY_VALUE);
        else
            irc->GetSubRegionCovFeature(sub_int_win, reg_cov, wcc->feature_indices,
                NULL, MIN_IDENTITY_VALUE, m_bUsedFeatureMean);

        if ( m_bNormalization ) {
            if ( wcc->n_features < mat_size ) {
                for ( i = 0 ; i < wcc->n_features ; i++ ) {
                    sub_parent_diag_elems[i] = parent_diag_elems[wcc->feature_indices[i]];
                    if ( m_bUsedFeatureMean )
                        sub_parent_avg_elems[i] = parent_avg_elems[wcc->feature_indices[i]];
                }
                reg_cov->NormalizeIllumination(sub_parent_diag_elems, m_bUsedFeatureMean?sub_parent_avg_elems:NULL);
            }
            else
                reg_cov->NormalizeIllumination(parent_diag_elems, m_bUsedFeatureMean?parent_avg_elems:NULL);
        }
        else {
            int x_feature_idx = wcc->feature_indices[0]==0?0:-1;
            int y_feature_idx = wcc->feature_indices[1]==1?1:-1;
            reg_cov->NormalizePositionFeatures(roi, x_feature_idx, y_feature_idx);
        }
        f_l = wcc->TestCov(reg_cov);
        F_X += 0.5 * f_l;
    }

	delete reg_cov;

	return F_X;
}
