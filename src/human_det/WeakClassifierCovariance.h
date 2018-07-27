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
#if !defined(_WEAK_CLASSIFIER_COVARIANCE_H_)
#define _WEAK_CLASSIFIER_COVARIANCE_H_

#include "WeakClassifier.h"
#include "CovarianceMatrix.h"
#include "IntegralRegionCov.h"
#include "Random.h"
#include "Rectangle.h"


// Weak classifier for covariance features
class WeakClassifierCovariance : public WeakClassifier {

public:
	/* Torch random generator */
    Torch::Random *random;

	/* do we map the covariance features
	   on the tangent space */
    bool map_features;

	/* do we need feature mean */
    bool used_feature_mean;

	/* number of features to be considered */
    int n_features;

	/* indices of features of interest */
    int *feature_indices;

	/* total number of features */
    int total_n_features;

    /* total number of valid symmetric covariance matrix elements */
    int total_valid_SM_size;
    int *valid_SM_indices;

    /* total number of valid mean feature elements */
    int total_valid_AVG_size;
    int *valid_AVG_indices;

	/* float region of interest (0 <= roi.x, roi.y <= 1) */
    struct floatRect roi;

    /* positive example weight for unbalanced training */
    double pos_ex_wscale;
    double updated_pos_ex_wscale;

    /* negative log likelihood */
    double negative_log_likelihood;

    /* Needed to map the covariance matrix to the tangent space */
    CCovarianceMatrix *weighted_mean;
    CCovarianceMatrix *weighted_mean_sqrt;
    CCovarianceMatrix *weighted_mean_inv_sqrt;
    CCovarianceMatrix *vec_log_x;

    /* If no feature indices are given, they are selected as been
	   the first n_features indices */
    WeakClassifierCovariance(float x, float y,
                           float width, float height,
                           int n_features,
                           int total_features,
                           int *feature_indices=NULL,
                           bool map=true,
                           bool used_avg=false);

    ~WeakClassifierCovariance();

	/* set positive example weight scale,
	   useful for unbalanced training */
    void setPositiveExampleWeightScale(double pos_ex_wscale_);

    /* Test a given example */
    double TestCov(CCovarianceMatrix* reg_cov);

	/* copy current data to another class instance */
    void Copy(WeakClassifierCovariance* dst);

	/* computing weights and responses of examples */
    void computeWeights(int n_examples, int* labels, double* cur_positive_probs, double* w, double* z );

	/* computing the mean covariance matrix for mapping */
    void computeMappingParameters(CCovarianceMatrix **X_mats, int n_matrices, double* weights, int* labels);

	/* get the negative binormial log-likelihood */
    double getNegativeBinormialLogLikelihood(double* cur_positive_probs, int n_examples, int* labels);

	/* update the positive weight scale, assuming that
	   sum of positive examples' weights equals to
	   sum of negative examples' weights */
    void updatePositiveWeightScale(double* w, int* labels, int n_examples);

	/* update the positive weights */
	void updatePositiveWeights(double* w, int* labels, int n_examples);

	/* update the posibitive probabilities and classifier functions */
    void updatePositiveProbFuncs(double* examples_data, int n_examples,
        double* positive_probs, double* classifier_funs);

    /* extract the feature vector for a given image. if
	   map_features is set to true, features are projected
	   onto the tangent space. Features are extracted based
	   on the current feature_indices values. Data is stored in
	   a table of n_features*n_examples  */
    double* extractFeatures(CCovarianceMatrix **cov_mats, int n_matrices);

	/* get covariance (full- or sub-) matrix from integral covariance data */
    CCovarianceMatrix* getCovMatrix(CIntegralRegionCov* cov_int, bool is_full_mat=false);

	/* randomly generate a new sub-window sample */
    void sampleWindow();
};

#endif // !defined(_WEAK_CLASSIFIER_COVARIANCE_H_)
