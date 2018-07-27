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
// BoostClassifierCovariance.h: interface for the CCascadeClassifierCovariance and CBoostClassifierCovariance classes.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_BOOST_CLASSIFIER_COVARIANCE_H_)
#define _BOOST_CLASSIFIER_COVARIANCE_H_

#include "WeakClassifierCovariance.h"
#include "BuildFeatureImages.h"
#include "IntegralRegionCov.h"
#include "cv.h"
#include "Rectangle.h"

inline bool bool_check(const char* x) {
	return !strcmp(x, "true");
};

class CCascadeClassifierCovariance
{
public:
	CCascadeClassifierCovariance(int features_num,   /* number of subset features */
			int tot_features_num,                    /* total number of features */
			int weak_classifiers_num=500,   /* maximal number of weak classifiers */
			bool map=true,          /* do we map the covariance features onto the tangent space */
			bool used_avg=false);   /* do we use the feature mean */

	~CCascadeClassifierCovariance();

	/* add a new weak classifier (weak_classifier=NULL),
	   or update current weak classifier */
	bool AddWeakClassifier(WeakClassifierCovariance* weak_classifier=NULL);

	/* get the classification rates:
	   positive detection rate and false positive rate */
	void GetClassificationRates(double* classifier_signs, int* labels, int n_examples);

    /* weak classifiers in one cascade level */
	WeakClassifierCovariance** m_ppWeakClassifers;

    /* feature number to consider in one cascade level */
	int m_nFeaturesNum;

	/* total feature number to be used */
	int m_nTotFeaturesNum;

    /* map the covariance to tangent space or not */
	bool m_bMap;

	/* integrate the feature mean with covariance or not */
	bool m_bUsedFeatureMean;

    /* current learned weak classifier number */
	int m_nWeakClassifiersNum;

    /* maximal weak classifier number */
	int m_nMaxWeakClassifiersNum;

    /* classification rates */
    double m_dbPositiveDetectionRate;
    double m_dbFalsePositiveRate;

    /* some thresholds for classification */
	double m_dbXnPosProbThreshold;
	double m_dbXnClassifierFunThreshold;
	double m_dbXpPosProbThreshold;
	double m_dbXpClassifierFunThreshold;
private:
};

class CBoostClassifierCovariance
{
public:
    /* constructor for training */
	CBoostClassifierCovariance(int *feature_img_types,
				int feature_type_num,
				int *cascade_subset_sizes=NULL,
				int max_cascade_level_num=35,
				int max_weak_classifiers_num=500,
				bool map=true,
				bool identity_map=false,
				bool used_avg=false,
				bool mean_cov_fall_ex=false,
				bool normalization=true,
				double level_reject_percent=0.35,
				double pos_example_edge_margin_x_percent=0.15,
				double pos_example_edge_margin_y_percent=0.1,
                double pos_example_shift_percent=0.0,
                int pos_example_shift_itr_num=1,
				int neg_cropped_edge_margin_size=3);

    /* constructor for detection */
    CBoostClassifierCovariance();

	~CBoostClassifierCovariance();

	/* set positive example weight scale (default = 1.0) */
	void SetPositiveExampleWeightScale(double pos_example_weight_scale=1.0);

	/* set whether do we update the positive weight scale */
	void SetUpdatePositiveWeightScale(bool update_pos_weight_scale=true);

	void ComputeAdapativePosExampleWeightScale(int n_pos_examples, int n_neg_examples);

    /* set the maximal negative/positive example numbers */
    void SetMaxTrainingExampleNumbers(int max_pos_examples_num, int max_neg_examples_num,
        int max_sel_pos_examples_num=0, int max_sel_neg_examples_num=0);

	/* add a new cascade level for training */
	bool AddCascadeLevel();

	/* add a new weak classifier (weak_classifier=NULL),
	   or update current weak classifier */
	bool AddWeakClassifier(WeakClassifierCovariance* weak_classifier=NULL);

	/* get the currently trained weak classifier */
	WeakClassifierCovariance* GetCurrentWeakClassifier();

    /* test the given sample */
	double Predictor(CIntegralRegionCov *irc,
				CvRect roi,
				int *final_pos_level,
				double *cascade_level_probs,
                double Xn_prob_alpha=0.5) {
		return Predictor(irc, roi, m_nCascadeLevelNum, Xn_prob_alpha,
                    final_pos_level, cascade_level_probs);
	};

	double Predictor(CIntegralRegionCov *irc,
				CvRect roi,
                double Xn_prob_alpha=0.5) {
		return Predictor(irc, roi, m_nCascadeLevelNum, Xn_prob_alpha,
                    NULL, NULL);
	};

    /* test the given sample */
	double Predictor(CIntegralRegionCov *irc,
				CvRect roi,
				int boost_classifiers_num,
                double Xn_prob_alpha=1.0,
                int *final_pos_level=NULL,
                double *cascade_level_probs=NULL,
                int neg_obj_id=-1);

    /* test the given sample using some cascade levels */
    bool CascadeLevelPredictor(CIntegralRegionCov *irc,
                CvRect roi,
                int start_cascade_level,
				int end_cascade_level,
				double Xn_prob_alpha);

    /* get classifier function at some cascade level */
    double GetCascadeClassifierFun(
                CIntegralRegionCov *irc, CvRect roi,
				int boost_classifiers_idx);

    void ExportLogMessage(const char* log_file_name, char *msg);

    /* save the learned model (.xml or .yml) */
	bool Save(const char* file_name);

    /* load the learned model (.xml or .yml) */
	bool Load(const char* file_name);

    /* create a weak classifier for current cascade level */
    WeakClassifierCovariance* CreateWeakClassifier(floatRect roi,
            int n_features, int *feature_indices);

    /* set the training stop conditions */
    void SetTrainingStopConditions(double max_false_pos_rate,
            double min_pos_det_rate,
            double margin_threshold );


    /* model data */
	CCascadeClassifierCovariance** m_ppCascadeClassifers;

    /* positive example weight for unbalanced training */
    double m_dbPositiveExampleWeightScale;
    double m_dbAdapativePositiveExampleWeightScale;

    /* whether to update positive example weight based on
       the same weights sum of both positive/negative examples */
    bool m_bUpdatePositiveWeightScale;

    /* the beta filter percent, (1-beta)*100% examples are kept
       in next weak classifier training process */
    double m_dbBetaFilterPercent;

	/* detection rates of positive examples at each cascade level */
	double* m_pPosDetRates;

	/* false positive rates while generating negative cropped image for training
	   at each cascade level */
	double* m_pFalsePosRates;

    /* feature image types */
	int* m_pFeatureImgTypes;

	/* number of feature image types */
	int m_nFeatureTypeNum;

    /* total feature number to be used */
	int m_nTotFeatureNum;

    /* map the covariance to tangent space or not */
	bool m_bMap;

	/* map the covairance to tangent space based on
	   identity matrix or covariance mean matrix */
	bool m_bIdentityMap;

	/* computing mean covariance matrix for mapping
	   from all the examples or only positive examples */
	bool m_bMeanCovFromAllExamples;

	/* integrate the feature mean with covariance or not */
	bool m_bUsedFeatureMean;

    /* normalize the covariance and mean features based on the parent window */
	bool m_bNormalization;

    /* subset feature sizes to be used in each cascade level */
	int* m_pCascadeSubsetSizes;

    /* computed weights for each cascade level */
	double* m_pCascadeLevelWeights;

    /* currently learned cascade level number */
	int m_nCascadeLevelNum;

    /* maximal cascade level number */
	int m_nMaxCascadeLevelNum;

	/* maximal weak classifier number in each cascade level */
	int m_nMaxWeakClassifiersNum;

	/* maximal number of positive examples to be used for training at each cascade level */
    int m_nMaxPosExamplesNum;

	/* maximal number of negative examples to be used for training at each cascade level */
    int m_nMaxNegExamplesNum;

	/* maximal number of selected positive examples to be used for training at each cascade level */
    int m_nMaxSelPosExamplesNum;

	/* maximal number of selected negative examples to be used for training at each cascade level */
    int m_nMaxSelNegExamplesNum;

    /* to avoid to boundary effect */
    /* margin percent for positive image examples */
	double m_dbPosExEdgeMarginXPercent;
	double m_dbPosExEdgeMarginYPercent;

    double m_dbPosExShiftPercent;
    int m_nPosExShiftItrNum;

	/* margin size of saved negative cropped image examples */
	int m_nNegCroppedEdgeMarginSize;

    /* training stop conditions */
    double m_dbMaxFalsePosRate;
    double m_dbMinPosDetRate;
    double m_dbMarginThreshold;

private:
    CRectangle m_cRect;

    /* clean the data */
	void CleanData();

    /* computing the weights for each cascade level */
	void ComputeCascadeLevelWeights(double level_reject_percent=0.35);

    /* false example reject percent at each cascade level */
    double m_dbLevelRejectPercent;
};

#endif // !defined(_BOOST_CLASSIFIER_COVARIANCE_H_)
