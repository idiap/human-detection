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
// BuildFeatureImages.cpp: implementation of the CBuildFeatureImages class.
//
//////////////////////////////////////////////////////////////////////

#include "BuildFeatureImages.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CBuildFeatureImages::CBuildFeatureImages(int *feature_img_types, int feature_type_num)
{
	m_nFeatureTypeNum = feature_type_num;
	m_pFeatureImageTypes = feature_img_types;

	int tot_feature_img_num = 0;
	int a;
	for ( a = 0 ; a < m_nFeatureTypeNum ; a++ ) {
		switch(m_pFeatureImageTypes[a]) {
		case F_POSITION_XY:		/* the 2D position (x,y) */
			tot_feature_img_num += 2;
			break;
		case F_RGB_COLOR:		/* the RGB color values (R, G, B) */
			tot_feature_img_num += 3;
			break;
		case F_HSV_COLOR:		/* the HSV color values (H, S, V) */
			tot_feature_img_num += 3;
			break;
		case F_GRAY_LEVEL:		/* the grayscale values */
			tot_feature_img_num++;
			break;
		case F_IMG_GRAD_XY:		/* the first-order gradients at x- and y-directions (I_x, I_y) */
			tot_feature_img_num += 2;
			break;
		case F_IMG_GRAD_XY_SQRT:	/* the first-order gradients, sqrt(I_x^2+I_y^2) */
			tot_feature_img_num++;
			break;
		case F_IMG_GRAD_XY2:		/* the second-order gradients at x- and y-directions (I_xx, I_yy) */
			tot_feature_img_num += 2;
			break;
		case F_IMG_GRAD_XY2_SQRT:	/* the second-order gradients, sqrt(I_xx^2+I_yy^2) */
			tot_feature_img_num++;
			break;
		case F_IMG_EDGE_XY_ORI:		/* the first-order edge orientation, arctan ( |I_x| / |I_y|) */
			tot_feature_img_num++;
			break;
		case F_IMG_EDGE_XY2_ORI:	/* the second-order edge orientation, arctan ( |I_xx| / |I_yy| ) */
			tot_feature_img_num++;
			break;
		case F_IMG_LBP:			/* the local binary pattern feature of the gray-level image */
			tot_feature_img_num++;
			break;
		case F_IMG_CANNY_EDGE:		/* the canny edge image */
			tot_feature_img_num++;
			break;
		case F_FG_PROB_MAP:		/* the foreground probability map */
			tot_feature_img_num++;
			break;
		case F_FG_BINARY_MAP:		/* the binary foreground mask map */
			tot_feature_img_num++;
			break;
		case F_MAP_GRAD_XY:		/* the first-order gradients on the foreground probability map, (F_x, F_y) */
			tot_feature_img_num += 2;
			break;
		case F_MAP_GRAD_XY_SQRT:	/* the first-order gradients, sqrt(F_x^2+F_y^2) */
			tot_feature_img_num++;
			break;
		case F_MAP_GRAD_XY2:		/* the second-order gradients at x- and y-directions (F_xx, F_yy) */
			tot_feature_img_num += 2;
			break;
		case F_MAP_GRAD_XY2_SQRT:	/* the second-order gradients, sqrt(F_xx^2+F_yy^2) */
			tot_feature_img_num++;
			break;
		case F_MAP_EDGE_XY_ORI:		/* the first-order edge orientation, arctan ( |F_x| / |F_y|) */
			tot_feature_img_num++;
			break;
		case F_MAP_EDGE_XY2_ORI:	/* the second-order edge orientation, arctan ( |F_xx| / |F_yy| ) */
			tot_feature_img_num++;
			break;
		case F_BIN_GRAD_XY:		/* the first-order gradients on the foreground binary mask map, (B_x, B_y) */
			tot_feature_img_num += 2;
			break;
		case F_BIN_GRAD_XY_SQRT:	/* the first-order gradients, sqrt(B_x^2+B_y^2) */
			tot_feature_img_num++;
			break;
		case F_BIN_EDGE_XY_ORI:		/* the first-order edge orientation, arctan ( |B_x| / |B_y|) */
			tot_feature_img_num++;
			break;
		}
	}

	m_nFeatureImageNum = tot_feature_img_num;
}

CBuildFeatureImages::~CBuildFeatureImages()
{
}

int CBuildFeatureImages::GetFeatureImageNumber()
{
	return m_nFeatureImageNum;
}

void CBuildFeatureImages::BuildFeatureImages(IplImage **feature_imgs, IplImage *org_img, IplImage *fg_prob_img, float fg_threshold, int one_existing_pos_feature)
{
	CvSize img_size = cvGetSize(org_img);
	bool bColor = (org_img->nChannels == 3);

	IplImage *org_float_img = cvCreateImage(img_size, IPL_DEPTH_32F, org_img->nChannels);
	IplImage *gray_float_img = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
	IplImage *gray_img = cvCreateImage(img_size, IPL_DEPTH_8U, 1);
	IplImage *fg_prob_float_img = 0;
	IplImage *fg_mask_img = 0;

	IplImage *canny_edge_img;

	cvConvert(org_img, org_float_img);

	if ( bColor ) {
		cvCvtColor(org_img, gray_img, CV_BGR2GRAY);
		cvConvert(gray_img, gray_float_img);
	}
	else
		cvCopy(org_float_img, gray_float_img);

	if ( fg_prob_img ) {
		fg_prob_float_img = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
		if ( fg_prob_img->nChannels == 3 && fg_prob_img->depth != IPL_DEPTH_32F ) {
			IplImage* tmp_fg_prob_img = cvCreateImage(img_size, IPL_DEPTH_32F, 3);
			cvConvert(fg_prob_img, tmp_fg_prob_img);
			cvCvtColor(tmp_fg_prob_img, fg_prob_float_img, CV_BGR2GRAY);
			cvReleaseImage(&tmp_fg_prob_img);
			fg_threshold *= 255.0f;
		}
		else if ( fg_prob_img->depth != IPL_DEPTH_32F ) { 	/* uchar type */
			cvConvert(fg_prob_img, fg_prob_float_img);
			fg_threshold *= 255.0f;
		}
		else
			cvCopy(fg_prob_img, fg_prob_float_img);

        bool used_fg_mask_img = false;
        for ( int i = 0 ; i < m_nFeatureTypeNum ; i++ ) {
            if ( m_pFeatureImageTypes[i] == F_BIN_GRAD_XY ||
                m_pFeatureImageTypes[i] == F_BIN_GRAD_XY_SQRT ||
                m_pFeatureImageTypes[i] == F_BIN_EDGE_XY_ORI ||
                m_pFeatureImageTypes[i] == F_FG_BINARY_MAP ) {
                used_fg_mask_img = true;
                break;
            }
        }
        if ( used_fg_mask_img ) {
            fg_mask_img = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
			cvThreshold(fg_prob_float_img, fg_mask_img, fg_threshold, 1, CV_THRESH_BINARY);
        }
	}

	bool bImgGrad_XY = false;
	bool bImgGrad_XY2 = false;
	bool bMapGrad_XY = false;
	bool bMapGrad_XY2 = false;
	bool bMaskGrad_XY = false;

	IplImage *img_grad_x, *img_grad_y, *img_grad_xx, *img_grad_yy;
	IplImage *map_grad_x, *map_grad_y, *map_grad_xx, *map_grad_yy;
	IplImage *mask_grad_x, *mask_grad_y;

	int a;
	int fimg_idx = 0;
	for ( a = 0 ; a < m_nFeatureTypeNum ; a++ ) {
		switch(m_pFeatureImageTypes[a]) {
		case F_POSITION_XY:		/* the 2D position (x,y) */
			if ( !one_existing_pos_feature )
				ComputePositionMaps(img_size, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
			fimg_idx += 2;
			break;
		case F_RGB_COLOR:		/* the RGB color values (R, G, B) */
			if ( bColor )
				cvSplit(org_float_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1], feature_imgs[fimg_idx+2], NULL);
			else {
				printf("Must input color images for building feature images!\n");
				exit(1);
			}
			fimg_idx += 3;
			break;
		case F_HSV_COLOR:		/* the HSV color values (H, S, V) */
			if ( bColor ) {
				IplImage *hsv_img = cvCreateImage(img_size, org_float_img->depth, 3);
				cvCvtColor( org_float_img, hsv_img, CV_BGR2HSV);
				cvSplit(org_float_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1], feature_imgs[fimg_idx+2], NULL);
				cvReleaseImage(&hsv_img);
			}
			else {
				printf("Must input color images for building feature images!\n");
				exit(1);
			}
			fimg_idx += 3;
			break;
		case F_GRAY_LEVEL:		/* the grayscale values */
			cvCopy(gray_float_img, feature_imgs[fimg_idx]);
			fimg_idx++;
			break;
		case F_IMG_GRAD_XY:		/* the first-order gradients at x- and y-directions (I_x, I_y) */
			ComputeFirstOrderGradients(gray_float_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
			img_grad_x = feature_imgs[fimg_idx];
			img_grad_y = feature_imgs[fimg_idx+1];
			bImgGrad_XY = true;
			fimg_idx += 2;
			break;
		case F_IMG_GRAD_XY_SQRT:	/* the first-order gradients, sqrt(I_x^2+I_y^2) */
			if ( bImgGrad_XY )
				ComputeGradients(img_grad_x, img_grad_y, feature_imgs[fimg_idx]);
			else {
				img_grad_x = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
				img_grad_y = cvCreateImage(img_size, IPL_DEPTH_32F, 1);

				ComputeFirstOrderGradients(gray_float_img, img_grad_x, img_grad_y);
				ComputeGradients(img_grad_x, img_grad_y, feature_imgs[fimg_idx]);

				cvReleaseImage(&img_grad_x);
				cvReleaseImage(&img_grad_y);
			}
			fimg_idx++;
			break;
		case F_IMG_GRAD_XY2:		/* the second-order gradients at x- and y-directions (I_xx, I_yy) */
			if ( bImgGrad_XY )
				ComputeSecondOrderGradients(img_grad_x, img_grad_y, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
			else
				ComputeSecondOrderGradients(gray_float_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
			img_grad_xx = feature_imgs[fimg_idx];
			img_grad_yy = feature_imgs[fimg_idx+1];
			bImgGrad_XY2 = true;
			fimg_idx += 2;
			break;
		case F_IMG_GRAD_XY2_SQRT:	/* the second-order gradients, sqrt(I_xx^2+I_yy^2) */
			if ( bImgGrad_XY2 )
				ComputeGradients(img_grad_xx, img_grad_yy, feature_imgs[fimg_idx]);
			else {
				img_grad_xx = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
				img_grad_yy = cvCreateImage(img_size, IPL_DEPTH_32F, 1);

				if ( bImgGrad_XY )
					ComputeSecondOrderGradients(img_grad_x, img_grad_y, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
				else
					ComputeSecondOrderGradients(gray_float_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);

				ComputeGradients(img_grad_xx, img_grad_yy, feature_imgs[fimg_idx]);

				cvReleaseImage(&img_grad_xx);
				cvReleaseImage(&img_grad_yy);
			}
			fimg_idx++;
			break;
		case F_IMG_EDGE_XY_ORI:		/* the first-order edge orientation, arctan ( |I_x| / |I_y|) */
			if ( bImgGrad_XY )
				ComputeEdgeOrientation(img_grad_x, img_grad_y, feature_imgs[fimg_idx]);
			else {
				img_grad_x = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
				img_grad_y = cvCreateImage(img_size, IPL_DEPTH_32F, 1);

				ComputeFirstOrderGradients(gray_float_img, img_grad_x, img_grad_y);
				ComputeEdgeOrientation(img_grad_x, img_grad_y, feature_imgs[fimg_idx]);

				cvReleaseImage(&img_grad_x);
				cvReleaseImage(&img_grad_y);
			}
			fimg_idx++;
			break;
		case F_IMG_EDGE_XY2_ORI:	/* the second-order edge orientation, arctan ( |I_xx| / |I_yy| ) */
			if ( bImgGrad_XY2 )
				ComputeEdgeOrientation(img_grad_xx, img_grad_yy, feature_imgs[fimg_idx]);
			else {
				img_grad_xx = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
				img_grad_yy = cvCreateImage(img_size, IPL_DEPTH_32F, 1);

				if ( bImgGrad_XY )
					ComputeSecondOrderGradients(img_grad_x, img_grad_y, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
				else
					ComputeSecondOrderGradients(gray_float_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);

				ComputeEdgeOrientation(img_grad_xx, img_grad_yy, feature_imgs[fimg_idx]);

				cvReleaseImage(&img_grad_xx);
				cvReleaseImage(&img_grad_yy);
			}
			fimg_idx++;
			break;
		case F_IMG_LBP:			/* the local binary pattern feature of the gray-level image */
			ComputeLBP(gray_float_img, feature_imgs[fimg_idx]);
			fimg_idx++;
			break;
		case F_IMG_CANNY_EDGE:		/* the canny edge image */
			canny_edge_img = cvCreateImage(img_size, IPL_DEPTH_8U, 1);
			cvCanny(gray_img, canny_edge_img, 20, 100);
			cvConvert(canny_edge_img, feature_imgs[fimg_idx]);
			cvReleaseImage(&canny_edge_img);
			fimg_idx++;
			break;
		case F_FG_PROB_MAP:		/* the foreground probability map */
			cvCopy(fg_prob_float_img, feature_imgs[fimg_idx]);
			fimg_idx++;
			break;
		case F_FG_BINARY_MAP:		/* the binary foreground mask map */
			/* get the foreground mask by thresholding */
			//cvThreshold(fg_prob_float_img, feature_imgs[fimg_idx], fg_threshold, 1, CV_THRESH_BINARY);
			cvCopy(fg_mask_img, feature_imgs[fimg_idx]);
			fimg_idx++;
			break;
		case F_MAP_GRAD_XY:		/* the first-order gradients on the foreground probability map, (F_x, F_y) */
			ComputeFirstOrderGradients(fg_prob_float_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
			map_grad_x = feature_imgs[fimg_idx];
			map_grad_y = feature_imgs[fimg_idx+1];
			bMapGrad_XY = true;
			fimg_idx += 2;
			break;
		case F_MAP_GRAD_XY_SQRT:	/* the first-order gradients, sqrt(F_x^2+F_y^2) */
			if ( bMapGrad_XY )
				ComputeGradients(map_grad_x, map_grad_y, feature_imgs[fimg_idx]);
			else {
				map_grad_x = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
				map_grad_y = cvCreateImage(img_size, IPL_DEPTH_32F, 1);

				ComputeFirstOrderGradients(fg_prob_float_img, map_grad_x, map_grad_y);
				ComputeGradients(map_grad_x, map_grad_y, feature_imgs[fimg_idx]);

				cvReleaseImage(&map_grad_x);
				cvReleaseImage(&map_grad_y);
			}
			fimg_idx++;
			break;
		case F_MAP_GRAD_XY2:		/* the second-order gradients at x- and y-directions (F_xx, F_yy) */
			if ( bMapGrad_XY )
				ComputeSecondOrderGradients(map_grad_x, map_grad_y, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
			else
				ComputeSecondOrderGradients(fg_prob_float_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
			map_grad_xx = feature_imgs[fimg_idx];
			map_grad_yy = feature_imgs[fimg_idx+1];
			bMapGrad_XY2 = true;
			fimg_idx += 2;
			break;
		case F_MAP_GRAD_XY2_SQRT:	/* the second-order gradients, sqrt(F_xx^2+F_yy^2) */
			if ( bMapGrad_XY2 )
				ComputeGradients(map_grad_xx, map_grad_yy, feature_imgs[fimg_idx]);
			else {
				map_grad_xx = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
				map_grad_yy = cvCreateImage(img_size, IPL_DEPTH_32F, 1);

				if ( bMapGrad_XY )
					ComputeSecondOrderGradients(map_grad_x, map_grad_y, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
				else
					ComputeSecondOrderGradients(fg_prob_float_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);

				ComputeGradients(map_grad_xx, map_grad_yy, feature_imgs[fimg_idx]);

				cvReleaseImage(&map_grad_xx);
				cvReleaseImage(&map_grad_yy);
			}
			fimg_idx++;
			break;
		case F_MAP_EDGE_XY_ORI:		/* the first-order edge orientation, arctan ( |F_x| / |F_y|) */
			if ( bMapGrad_XY )
				ComputeEdgeOrientation(map_grad_x, map_grad_y, feature_imgs[fimg_idx]);
			else {
				map_grad_x = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
				map_grad_y = cvCreateImage(img_size, IPL_DEPTH_32F, 1);

				ComputeFirstOrderGradients(fg_prob_float_img, map_grad_x, map_grad_y);
				ComputeEdgeOrientation(map_grad_x, map_grad_y, feature_imgs[fimg_idx]);

				cvReleaseImage(&map_grad_x);
				cvReleaseImage(&map_grad_y);
			}
			fimg_idx++;
			break;
		case F_MAP_EDGE_XY2_ORI:	/* the second-order edge orientation, arctan ( |F_xx| / |F_yy| ) */
			if ( bMapGrad_XY2 )
				ComputeEdgeOrientation(map_grad_xx, map_grad_yy, feature_imgs[fimg_idx]);
			else  {
				map_grad_xx = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
				map_grad_yy = cvCreateImage(img_size, IPL_DEPTH_32F, 1);

				if ( bMapGrad_XY )
					ComputeSecondOrderGradients(map_grad_x, map_grad_y, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
				else
					ComputeSecondOrderGradients(fg_prob_float_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);

				ComputeEdgeOrientation(map_grad_xx, map_grad_yy, feature_imgs[fimg_idx]);

				cvReleaseImage(&map_grad_xx);
				cvReleaseImage(&map_grad_yy);
			}
			fimg_idx++;
			break;


		case F_BIN_GRAD_XY:		/* the first-order gradients on the foreground binary mask map, (B_x, B_y) */
			ComputeFirstOrderGradients(fg_mask_img, feature_imgs[fimg_idx], feature_imgs[fimg_idx+1]);
			mask_grad_x = feature_imgs[fimg_idx];
			mask_grad_y = feature_imgs[fimg_idx+1];
			bMaskGrad_XY = true;
			fimg_idx += 2;
			break;
		case F_BIN_GRAD_XY_SQRT:	/* the first-order gradients, sqrt(B_x^2+B_y^2) */
			if ( bMaskGrad_XY )
				ComputeGradients(mask_grad_x, mask_grad_y, feature_imgs[fimg_idx]);
			else {
				mask_grad_x = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
				mask_grad_y = cvCreateImage(img_size, IPL_DEPTH_32F, 1);

				ComputeFirstOrderGradients(fg_mask_img, mask_grad_x, mask_grad_y);
				ComputeGradients(mask_grad_x, mask_grad_y, feature_imgs[fimg_idx]);

				cvReleaseImage(&mask_grad_x);
				cvReleaseImage(&mask_grad_y);
			}
			fimg_idx++;
			break;
		case F_BIN_EDGE_XY_ORI:		/* the first-order edge orientation, arctan ( |B_x| / |B_y|) */
			if ( bMaskGrad_XY )
				ComputeEdgeOrientation(mask_grad_x, mask_grad_y, feature_imgs[fimg_idx]);
			else {
				mask_grad_x = cvCreateImage(img_size, IPL_DEPTH_32F, 1);
				mask_grad_y = cvCreateImage(img_size, IPL_DEPTH_32F, 1);

				ComputeFirstOrderGradients(fg_mask_img, mask_grad_x, mask_grad_y);
				ComputeEdgeOrientation(mask_grad_x, mask_grad_y, feature_imgs[fimg_idx]);

				cvReleaseImage(&mask_grad_x);
				cvReleaseImage(&mask_grad_y);
			}
			fimg_idx++;
			break;

		}
	}

	for ( a = 0 ; a < m_nFeatureImageNum ; a++ )
		cvAbs(feature_imgs[a], feature_imgs[a]);

	cvReleaseImage(&org_float_img);
	cvReleaseImage(&gray_float_img);
	cvReleaseImage(&gray_img);
	if ( fg_prob_float_img )
		cvReleaseImage(&fg_prob_float_img);
    if ( fg_mask_img )
		cvReleaseImage(&fg_mask_img);
}

void CBuildFeatureImages::ComputeFirstOrderGradients(IplImage *gray_float_img, IplImage *grad_x, IplImage *grad_y)
{
	int aperture_size = 3;

	cvSobel(gray_float_img, grad_x, 1, 0, aperture_size);
	cvSobel(gray_float_img, grad_y, 0, 1, aperture_size);

//	cvAbs(grad_x, grad_x);
//	cvAbs(grad_y, grad_y);
}

void CBuildFeatureImages::ComputeSecondOrderGradients(IplImage *grad_x, IplImage *grad_y, IplImage *grad_xx, IplImage *grad_yy)
{
	int aperture_size = 3;
	cvSobel(grad_x, grad_xx, 1, 0, aperture_size);
	cvSobel(grad_y, grad_yy, 0, 1, aperture_size);

//	cvAbs(grad_xx, grad_xx);
//	cvAbs(grad_yy, grad_yy);
}


void CBuildFeatureImages::ComputeSecondOrderGradients(IplImage *gray_float_img, IplImage *grad_xx, IplImage *grad_yy)
{
	IplImage *grad_x = cvCreateImage(cvGetSize(gray_float_img), IPL_DEPTH_32F, 1);
	IplImage *grad_y = cvCreateImage(cvGetSize(gray_float_img), IPL_DEPTH_32F, 1);

	ComputeFirstOrderGradients(gray_float_img, grad_x, grad_y);
	ComputeSecondOrderGradients(grad_x, grad_y, grad_xx, grad_yy);

//	cvAbs(grad_xx, grad_xx);
//	cvAbs(grad_yy, grad_yy);

	cvReleaseImage(&grad_x);
	cvReleaseImage(&grad_y);
}

void CBuildFeatureImages::ComputeGradients(IplImage *grad_x, IplImage *grad_y, IplImage *grad_xy)
{
	int y, x;
	float* x_grad_x_data;
	float* x_grad_y_data;
	float* x_grad_xy_data;

	for ( y = 0 ; y < grad_x->height ; y++ ) {
		x_grad_x_data = (float*)(grad_x->imageData + grad_x->widthStep*y);
		x_grad_y_data = (float*)(grad_y->imageData + grad_y->widthStep*y);
		x_grad_xy_data = (float*)(grad_xy->imageData + grad_xy->widthStep*y);

		for ( x = 0 ; x < grad_x->width ; x++ )  {
			*x_grad_xy_data++ = sqrtf(((*x_grad_x_data) * (*x_grad_x_data) + (*x_grad_y_data) * (*x_grad_y_data)));
			x_grad_x_data++;
			x_grad_y_data++;
		}
	}
}

void CBuildFeatureImages::ComputeEdgeOrientation(IplImage *grad_x, IplImage *grad_y, IplImage *edge_ori)
{
	int y, x;
	float* x_grad_x_data;
	float* x_grad_y_data;
	float* x_edge_ori_data;

	for ( y = 0 ; y < grad_x->height ; y++ ) {
		x_grad_x_data = (float*)(grad_x->imageData + grad_x->widthStep*y);
		x_grad_y_data = (float*)(grad_y->imageData + grad_y->widthStep*y);
		x_edge_ori_data = (float*)(edge_ori->imageData + edge_ori->widthStep*y);

		for ( x = 0 ; x < grad_x->width ; x++ )  {
			*x_edge_ori_data++ = cvFastArctan(fabsf(*x_grad_y_data++), fabsf(*x_grad_x_data++));
			//*x_edge_ori_data++ = atan2f(fabsf(*x_grad_y_data++), fabsf(*x_grad_x_data++));
		}
	}
}

void CBuildFeatureImages::ComputePositionMaps(CvSize img_size, IplImage *map_x, IplImage *map_y)
{
	float *map_x_data = new float[img_size.width*img_size.height];
	float *map_y_data = new float[img_size.width*img_size.height];

	float *_x_data = map_x_data, *_y_data = map_y_data;

	if ( !map_y_data )
		return;

	int y, x;
	for ( y = 0 ; y < img_size.height ; y++ )
		for ( x = 0 ; x < img_size.width ; x++ ) {
			*_x_data++ = (float)x;
			*_y_data++ = (float)y;
		}

	COpencvDataConversion<float,float> ODC;
	ODC.SetImageData(map_x, map_x_data);
	ODC.SetImageData(map_y, map_y_data);

	delete [] map_x_data;
	delete [] map_y_data;
}

//void CBuildFeatureImages::ConvertImage(IplImage *src_uchar_Img, IplImage *dst_float_Img)
//{
//	if ( src_uchar_Img->width != dst_float_Img->width || src_uchar_Img->height != dst_float_Img->height ) {
//		printf("ConvertImage size is not correct!\n");
//	}
//	int y, x;
//	uchar* x_src_data;
//	float* x_dst_data;
//	for ( y = 0 ; y < src_uchar_Img->height ; y++ ) {
//		x_src_data = (uchar*)(src_uchar_Img->imageData + src_uchar_Img->widthStep*y);
//		x_dst_data = (float*)(dst_float_Img->imageData + dst_float_Img->widthStep*y);

//		for ( x = 0 ; x < src_uchar_Img->width*src_uchar_Img->nChannels ; x++ )
//			*x_dst_data++ = (float)(*x_src_data++);
//	}
//}

void CBuildFeatureImages::ComputeLBP(IplImage *gray_img, IplImage *lbp_img, int radius, int neig_pt_num, float robust_noise)
{
	IplImage* new_img = cvCreateImage(cvSize(gray_img->width+radius*2, gray_img->height+radius*2), gray_img->depth, 1);
	cvSetZero(new_img);

	CvRect img_roi = cvRect(radius,radius,gray_img->width,gray_img->height);
	cvSetImageROI(new_img, img_roi);
	cvCopy(gray_img, new_img);
	cvResetImageROI(new_img);

	COpencvDataConversion<float,float> ODC;
	float* img_data = ODC.GetImageData(new_img);

	int *offset_x = new int[neig_pt_num];
	int *offset_y = new int[neig_pt_num];

	int a;
	float angle;
	for ( a = 0 ; a < neig_pt_num ; a++ ) {
		angle = (float)a/(float)neig_pt_num*2.0f*PI;
		offset_x[a] = cvRound((float)radius*(cosf(angle)+1.0f));
		offset_y[a] = cvRound((float)radius*(-sinf(angle)+1.0f));
	}

	float **neig_pts_data = new float*[neig_pt_num];

	for ( a = 0 ; a < neig_pt_num ; a++ )
		neig_pts_data[a] = img_data+offset_y[a]*new_img->width+offset_x[a];

	float *cur_pt_data = img_data+radius*new_img->width+radius;

	float* lbp_data = new float[gray_img->width*gray_img->height];
	memset(lbp_data,0,sizeof(float)*gray_img->width*gray_img->height);
	float* _lbp_data = lbp_data;

	int* bit_OR = new int[neig_pt_num];
	bit_OR[0] = 1;
	for ( a = 1 ; a < neig_pt_num ; a++ )
		bit_OR[a] = bit_OR[a-1]<<1;
	int* _bit_OR = bit_OR;

	int bit_lbp;

	int y, x;
	for ( y = radius ; y < gray_img->height+radius ; y++ ) {
		for ( x = radius ; x < gray_img->width+radius ; x++ ) {
			_bit_OR = bit_OR;
			bit_lbp = 0;
			for ( a = 0 ; a < neig_pt_num ; a++ ) {
				if ( BIN_LBP_ELEM(*neig_pts_data[a], *cur_pt_data, robust_noise) )
					bit_lbp |= *_bit_OR;
				neig_pts_data[a]++;
				_bit_OR++;
			}
			*_lbp_data++ = (float)bit_lbp;
			cur_pt_data++;
		}
		for ( a = 0 ; a < neig_pt_num ; a++ )
			neig_pts_data[a] += radius;
		cur_pt_data += radius;
	}

	ODC.SetImageData(lbp_img, lbp_data);

	/************************************************************************/
	/* release memories                                                     */
	/************************************************************************/
	delete [] lbp_data;
	delete [] neig_pts_data;
	delete [] img_data;
	delete [] bit_OR;
	delete [] offset_x;
	delete [] offset_y;
	cvReleaseImage(&new_img);
}

void CBuildFeatureImages::SetBinFeature(IplImage *feature_img, int bins_num, float *range)
{
	int y, x;
	float* x_data;
	float low, high_low;
	low = range[0];
	high_low = range[1]-range[0];
	float bins_num_1 = (float)bins_num - 1.0f;
	float factor = bins_num_1/high_low;

	for ( y = 0 ; y < feature_img->height ; y++ ) {
		x_data = (float*)(feature_img->imageData + feature_img->widthStep*y);
		for ( x = 0 ; x < feature_img->width ; x++ ) {
			*x_data -= low;
			if ( *x_data <= 0 )
				*x_data++ = 0;
			else if ( *x_data >= high_low )
				*x_data++ = bins_num_1;
			else {
				*x_data = (float)((int)(*x_data*factor));
				x_data++;
			}
		}
	}
}
