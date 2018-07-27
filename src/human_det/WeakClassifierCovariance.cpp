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
#include "WeakClassifierCovariance.h"

int compare(const void * a, const void * b)
{
    if (*(int*)a == *(int*)b)
        return 0;
    else
        if (*(int*)a < *(int*)b)
            return -1;
        else
            return 1;
}

/* Constructor */
WeakClassifierCovariance::WeakClassifierCovariance(float x, float y,
                         float width, float height,
                         int n_features_,
                         int total_features,
                         int *feature_indices_,
                         bool map,
                         bool used_avg):
                         WeakClassifier( used_avg ?
                                (n_features_*(n_features_+3))/2 :
                                (n_features_*(n_features_+1))/2 )
{

    random = new Torch::Random;

    map_features = map;
    n_features = n_features_;
    used_feature_mean = used_avg;

    pos_ex_wscale = 1.0;
	updated_pos_ex_wscale = 1.0;

    total_n_features = total_features;

    /* no indices are given -> take n_features at random */
    if(feature_indices_ == NULL) {
        int* all_feature_indices = new int[total_features];
        feature_indices    = new int[n_features];
        random->getShuffledIndices(all_feature_indices, total_features);
        /* Warning, the feature indices must be sorted !!!*/
        qsort(all_feature_indices, n_features, sizeof(int), compare);
        for (int i=0;i<n_features;i++){
            feature_indices[i] = all_feature_indices[i];
        }
        delete[] all_feature_indices;
    }
    else{
        feature_indices = new int[n_features];
        for (int i=0;i<n_features;i++){
            feature_indices[i] = feature_indices_[i];
        }
    }

    int SM_indices[1000];
    int AVG_indices[1000];

    total_valid_AVG_size = 0;
    total_valid_SM_size = 0;
    int SM_index = 0;
    int AVG_index = 0;
    for ( int a = 0 ; a < n_features ; a++ ) {
        if ( used_feature_mean ) {
            if ( ( feature_indices[a] == 0 || feature_indices[a] == 1 ) )
                AVG_index++;
            else
                AVG_indices[total_valid_AVG_size++] = AVG_index++;
        }

        for ( int b = a ; b < n_features ; b++ ) {
            if ( ( feature_indices[a] == 0 || feature_indices[a] == 1 ) &&
                ( feature_indices[b] == 0 || feature_indices[b] == 1 ) )
                SM_index++;
            else
                SM_indices[total_valid_SM_size++] = SM_index++;
        }
    }

    feature_size = total_valid_SM_size + total_valid_AVG_size;

    valid_SM_indices = new int[total_valid_SM_size];
    if ( total_valid_AVG_size )
        valid_AVG_indices = new int[total_valid_AVG_size];
    else
        valid_AVG_indices = NULL;
    for ( int i = 0 ; i < total_valid_SM_size ; i++ )
        valid_SM_indices[i] = SM_indices[i];
    for ( int i = 0 ; i < total_valid_AVG_size ; i++ )
        valid_AVG_indices[i] = AVG_indices[i];

    /* if no region of interest is given, sample one */
    if(x==-1 | y==-1 | width==-1 | height==-1){
        sampleWindow();
    }
    else{
        this->roi.x = x;
        this->roi.y = y;
        this->roi.width  = width;
        this->roi.height = height;
    }

    /* weighted mean related information is allocated only if mapping */
    if(map) {
        weighted_mean          = new CCovarianceMatrix(this->n_features);
        weighted_mean_sqrt     = new CCovarianceMatrix(this->n_features);
        weighted_mean_inv_sqrt = new CCovarianceMatrix(this->n_features);

		vec_log_x              = new CCovarianceMatrix(this->n_features);

		weighted_mean->SetIdentityMatrix(1.0);
		weighted_mean_sqrt->SetIdentityMatrix(1.0);
		weighted_mean_inv_sqrt->SetIdentityMatrix(1.0);
    }
    else{
        weighted_mean          = NULL;
        weighted_mean_sqrt     = NULL;
        weighted_mean_inv_sqrt = NULL;
        vec_log_x              = NULL;
    }

    negative_log_likelihood = 1.0E+100;
}

void WeakClassifierCovariance::
    setPositiveExampleWeightScale(double pos_ex_wscale_)
{
    pos_ex_wscale = pos_ex_wscale_;
    updated_pos_ex_wscale = pos_ex_wscale_;
}

/* Sample a window */
void WeakClassifierCovariance::sampleWindow()
{
	double min_subwin_area = 0.02;

    this->roi.width  = random->boundedUniform(0.1, 1.0);
    this->roi.height = random->boundedUniform(0.1, 1.0);

    this->roi.height = MAX( this->roi.height, min_subwin_area/this->roi.width );

    this->roi.x      = random->boundedUniform(0.0, 1.0-this->roi.width);
    this->roi.y      = random->boundedUniform(0.0, 1.0-this->roi.height);
}

/* Destructor */
WeakClassifierCovariance::~WeakClassifierCovariance()
{
    delete random;

    if(feature_indices != NULL){
        delete[] feature_indices;
    }

    if(weighted_mean != NULL){
        delete weighted_mean;
    }

    if(weighted_mean_sqrt != NULL){
        delete weighted_mean_sqrt;
    }

    if(weighted_mean_inv_sqrt != NULL){
        delete weighted_mean_inv_sqrt;
    }

    if(vec_log_x != NULL){
        delete vec_log_x;
    }

    if(valid_SM_indices != NULL){
        delete[] valid_SM_indices;
    }

    if(valid_AVG_indices != NULL){
        delete[] valid_AVG_indices;
    }
}

void WeakClassifierCovariance
     ::computeWeights(int n_examples, int* labels, double* cur_positive_probs, double* w, double* z )
{
    double* pX = cur_positive_probs;

    for (int i=0;i<n_examples;i++){
        w[i] = MAX(pX[i] * (1.0 - pX[i]), MIN_DOUBLE_NUM);

        if(labels[i]){
            z[i] = 1.0 / MAX(pX[i], MIN_DOUBLE_NUM);
        }
        else {
            z[i] = -1.0 / MAX(1.0 - pX[i], MIN_DOUBLE_NUM);
        }
    }
}

void WeakClassifierCovariance::updatePositiveWeightScale(double* w, int* labels, int n_examples)
{
    double neg_w_sum = 0, pos_w_sum = 0;
    for ( int i = 0 ; i < n_examples ; i++ ) {
        if ( labels[i] ) {
            pos_w_sum += w[i];
        }
        else {
            neg_w_sum += w[i];
        }
    }

    updated_pos_ex_wscale = pos_ex_wscale * neg_w_sum / pos_w_sum;
}

void WeakClassifierCovariance::updatePositiveWeights(double* w, int* labels, int n_examples)
{
	if ( updated_pos_ex_wscale == 1.0 )
		return;

    for (int i=0;i<n_examples;i++){
        if(labels[i])
	        w[i] *= updated_pos_ex_wscale;
    }
}

void WeakClassifierCovariance
    ::updatePositiveProbFuncs(double* examples_data, int n_examples,
        double* positive_probs, double* classifier_funs)
{
    double* tmp_p_X = positive_probs;
    double* tmp_F_X = classifier_funs;
    double e_F_X = 0;
	double* tmp_data = examples_data;
    for ( int i = 0 ; i < n_examples ; i++ ) {
        *tmp_F_X += 0.5*this->Test(tmp_data);

        e_F_X = exp(*tmp_F_X++);
        *tmp_p_X++ = e_F_X / ( e_F_X + 1.0/e_F_X );

		tmp_data += feature_size;
    }
}

/* Test a given example */
double WeakClassifierCovariance::TestCov(CCovarianceMatrix* reg_cov)
{
	double score = 0;

	if ( map_features ) {
		CCovarianceMatrix *w_inv_sqrt = this->weighted_mean_inv_sqrt;
		w_inv_sqrt->Mul2(reg_cov, vec_log_x);
		vec_log_x->Log();
	}
	else
		vec_log_x = reg_cov;

	double *SM_ptr = vec_log_x->SM_ptr;
	int i;
	double *tmp_fitting_params = fitting_params;
	for ( i = 0 ; i < total_valid_SM_size ; i++ )
       	score += SM_ptr[valid_SM_indices[i]] * (*tmp_fitting_params++);

    double *AVG_ptr = reg_cov->AVG_ptr;
    for ( i = 0 ; i < total_valid_AVG_size ; i++ ) {
        score += AVG_ptr[valid_AVG_indices[i]] * (*tmp_fitting_params++);
    }

	score += *tmp_fitting_params;

    if ( !map_features )
        vec_log_x = NULL;

	return score;
}

/* Computing the mapping parameters on the tangent space */
void WeakClassifierCovariance::
     computeMappingParameters(CCovarianceMatrix **X_mats, int n_matrices, double* weights, int* labels)
{
    /* Test first covariance matrix to verify its size */
    if(X_mats[0]->M_size < n_features){
        fprintf(stderr, "Problem with feature size\n");
    }

    CCovarianceMatrix **sub_cov_mats = NULL;
    if ( X_mats[0]->M_size > n_features ) {
        sub_cov_mats = new CCovarianceMatrix*[n_matrices];
        for ( int i = 0 ; i < n_matrices ; i++ ) {
            sub_cov_mats[i] = new CCovarianceMatrix(n_features);
            X_mats[i]->GetSubMatrix(feature_indices, n_features, sub_cov_mats[i]);
        }
    }
    else
        sub_cov_mats = X_mats;

    double dist_sums[MAX_WEIGHTED_MEAN_ITR_NUM];

    weighted_mean->WeightedMean(sub_cov_mats, n_matrices,
                              weights, labels,
                              MAX_WEIGHTED_MEAN_ITR_NUM, dist_sums);

    weighted_mean->Sqrt(weighted_mean_sqrt);
    weighted_mean_sqrt->Inverse(weighted_mean_inv_sqrt);

    if ( X_mats[0]->M_size > n_features ) {
        for ( int i = 0 ; i < n_matrices ; i++ )
            delete sub_cov_mats[i];
        delete [] sub_cov_mats;
    }
}

/* compute negative binormial log likelihood */
double WeakClassifierCovariance::
    getNegativeBinormialLogLikelihood(double* cur_positive_probs, int n_examples, int* labels)
{
	negative_log_likelihood = 0;
	int *tmp_y = labels;
	double *tmp_p_X = cur_positive_probs;

	for ( int i = 0 ; i < n_examples ; i++ ) {
		if ( *tmp_y++ ) {
			if ( *tmp_p_X <= 0 ) {
				//printf("p_X = %.20f - log error in getNegativeBinormialLogLikelihood!\n", *tmp_p_X);

				*tmp_p_X = MIN_DOUBLE_NUM;
                negative_log_likelihood += pos_ex_wscale*log(*tmp_p_X++) - 10000.0;
			}
			else
                negative_log_likelihood += pos_ex_wscale*log(*tmp_p_X++);
		}
		else {
			if ( *tmp_p_X >= 1.0 ) {
				//printf("p_X = %.20f - log error in getNegativeBinormialLogLikelihood!\n", *tmp_p_X);
				*tmp_p_X = 1.0-MIN_DOUBLE_NUM;

                *tmp_p_X = MIN(*tmp_p_X, 1.0-MIN_DOUBLE_NUM);
                negative_log_likelihood += log(1.0 - *tmp_p_X++) - 10000.0;
			}
			else {
                *tmp_p_X = MIN(*tmp_p_X, 1.0-MIN_DOUBLE_NUM);
                negative_log_likelihood += log(1.0 - *tmp_p_X++);
			}
		}
	}

	negative_log_likelihood = -negative_log_likelihood;

	return negative_log_likelihood;
}

/* Extract the vector of double given a set of covariance matrix
   returns feature vectors stored in a single dimensional double array */
double* WeakClassifierCovariance::
        extractFeatures(CCovarianceMatrix **cov_mats, int n_matrices)
{
    if( cov_mats[0]->M_size != total_n_features ){
        fprintf(stderr, "Error, the size of the covariance features do not correspond to the classifier size\n");
    }

    CCovarianceMatrix* sub_cov_mat = NULL;
    if ( cov_mats[0]->M_size > n_features )
        sub_cov_mat = new CCovarianceMatrix(n_features);

    double* data = new double[n_matrices*this->feature_size];

    double* tmp_data = data;
    for(int i=0;i<n_matrices;i++){
        if ( cov_mats[0]->M_size > n_features )
            cov_mats[i]->GetSubMatrix(feature_indices, n_features, sub_cov_mat);
        else
            sub_cov_mat = cov_mats[i];

        /* Project the covariance matrix data onto tangent space */
        if(map_features){
            weighted_mean->OrthogonalVectorLog(sub_cov_mat, vec_log_x,
                                             this->weighted_mean_sqrt,
                                             this->weighted_mean_inv_sqrt);
            for(int j=0;j<total_valid_SM_size;j++){
                *tmp_data++ = vec_log_x->SM_ptr[valid_SM_indices[j]];
            }
        }
        /* No projection here, just copying data */
        else{
            for(int j=0;j<total_valid_SM_size;j++){
                *tmp_data++ = sub_cov_mat->SM_ptr[valid_SM_indices[j]];
            }
        }

        for ( int j = 0 ; j < total_valid_AVG_size ; j++ )
            *tmp_data++ = sub_cov_mat->AVG_ptr[valid_SM_indices[j]];
    }

    if ( cov_mats[0]->M_size > n_features )
        delete sub_cov_mat;

    return data;
}

void WeakClassifierCovariance::
    Copy(WeakClassifierCovariance* dst)
{
	/* from WeakClassifierCovariance */
    dst->map_features = this->map_features;
    dst->used_feature_mean = this->used_feature_mean;
    dst->n_features = this->n_features;

    if ( dst->feature_indices )
        delete[] dst->feature_indices;
    dst->feature_indices = new int[dst->n_features];
    for ( int i = 0 ; i < this->n_features; i++ )
        dst->feature_indices[i] = this->feature_indices[i];

    dst->total_n_features = this->total_n_features;
    dst->roi = this->roi;

    dst->total_valid_SM_size = this->total_valid_SM_size;
    dst->total_valid_AVG_size = this->total_valid_AVG_size;

    if ( dst->valid_SM_indices )
        delete [] dst->valid_SM_indices;
    if ( dst->valid_AVG_indices )
        delete [] dst->valid_AVG_indices;
    dst->valid_SM_indices = NULL;
    dst->valid_AVG_indices = NULL;

    if ( dst->total_valid_SM_size ) {
        dst->valid_SM_indices = new int[dst->total_valid_SM_size];
        for ( int i = 0 ; i < dst->total_valid_SM_size ; i++ )
            dst->valid_SM_indices[i] = this->valid_SM_indices[i];
    }
    if ( dst->total_valid_AVG_size ) {
        dst->valid_AVG_indices = new int[dst->total_valid_AVG_size];
        for ( int i = 0 ; i < dst->total_valid_AVG_size ; i++ )
            dst->valid_AVG_indices[i] = this->valid_AVG_indices[i];
    }

    dst->pos_ex_wscale = this->pos_ex_wscale;
    dst->negative_log_likelihood = this->negative_log_likelihood;

    if ( this->map_features ) {
        this->weighted_mean->Copy(dst->weighted_mean);
        this->weighted_mean_sqrt->Copy(dst->weighted_mean_sqrt);
        this->weighted_mean_inv_sqrt->Copy(dst->weighted_mean_inv_sqrt);
    }

	/* from WeakClassifier */
    dst->feature_size = this->feature_size;
    for ( int i = 0 ; i < this->feature_size+1; i++ )
        dst->fitting_params[i] = this->fitting_params[i];

	dst->threshold = this->threshold;
}

/* Extract covariance matrix based on roi of the classifier */
CCovarianceMatrix* WeakClassifierCovariance::
                   getCovMatrix(CIntegralRegionCov* cov_int, bool is_full_mat)
{
    CvRect rect_roi;
    rect_roi.x = cvFloor(roi.x*cov_int->m_szFeatureImage.width);
    rect_roi.y = cvFloor(roi.y*cov_int->m_szFeatureImage.height);
    rect_roi.width = cvFloor(roi.width*cov_int->m_szFeatureImage.width);
    rect_roi.height = cvFloor(roi.height*cov_int->m_szFeatureImage.height);

    if ( is_full_mat ) {
        CCovarianceMatrix* reg_cov = new CCovarianceMatrix(total_n_features);
        cov_int->GetRegionCovFeature(rect_roi, reg_cov);
        return reg_cov;
    }
    else {
        CCovarianceMatrix* sub_reg_cov = new CCovarianceMatrix(n_features);
        cov_int->GetSubRegionCovFeature(rect_roi, sub_reg_cov, this->feature_indices);
        return sub_reg_cov;
    }

    return NULL;
}
