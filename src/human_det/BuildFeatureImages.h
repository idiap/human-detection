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
// BuildFeatureImages.h: interface for the CBuildFeatureImages class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_BUILD_FEATURE_IMAGES_H_)
#define _BUILD_FEATURE_IMAGES_H_

#include "cv.h"
#include <cstdio>
#include "OpencvDataConversion.h"

/* define a set of feature mapping types */
#define F_POSITION_XY		0		/* the 2D position (x,y) */
#define F_RGB_COLOR		1		/* the RGB color values (R, G, B) */
#define F_HSV_COLOR		2		/* the HSV color values (H, S, V) */
#define F_GRAY_LEVEL		3		/* the grayscale values */
#define F_IMG_GRAD_XY		4		/* the first-order gradients at x- and y-directions (I_x, I_y) */
#define F_IMG_GRAD_XY_SQRT	6		/* the first-order gradients, sqrt(I_x^2+I_y^2) */
#define F_IMG_GRAD_XY2		7		/* the second-order gradients at x- and y-directions (I_xx, I_yy) */
#define F_IMG_GRAD_XY2_SQRT	9		/* the second-order gradients, sqrt(I_xx^2+I_yy^2) */
#define F_IMG_EDGE_XY_ORI	10		/* the first-order edge orientation, arctan ( |I_x| / |I_y|) */
#define F_IMG_EDGE_XY2_ORI	11		/* the second-order edge orientation, arctan ( |I_xx| / |I_yy| ) */
#define F_IMG_LBP		12		/* the local binary pattern feature of the gray-level image */
#define F_IMG_CANNY_EDGE	13		/* the canny edge image of gray-scale image */
#define F_FG_PROB_MAP		20		/* the foreground probability map */
#define F_FG_BINARY_MAP		21		/* the binary foreground mask map */
#define F_MAP_GRAD_XY		22		/* the first-order gradients on the foreground probability map, (F_x, F_y) */
#define F_MAP_GRAD_XY_SQRT	23		/* the first-order gradients, sqrt(F_x^2+F_y^2) */
#define F_MAP_GRAD_XY2		24		/* the second-order gradients at x- and y-directions (F_xx, F_yy) */
#define F_MAP_GRAD_XY2_SQRT	25		/* the second-order gradients, sqrt(F_xx^2+F_yy^2) */
#define F_MAP_EDGE_XY_ORI	26		/* the first-order edge orientation, arctan ( |F_x| / |F_y|) */
#define F_MAP_EDGE_XY2_ORI	27		/* the second-order edge orientation, arctan ( |F_xx| / |F_yy| ) */
#define F_BIN_GRAD_XY		28		/* the first-order gradients on the foreground mask (binary) image, (B_x, B_y) */
#define F_BIN_GRAD_XY_SQRT	29		/* the first-order gradients, sqrt(B_x^2+B_y^2) */
#define F_BIN_EDGE_XY_ORI	30		/* the first-order edge orientation, arctan ( |B_x| / |B_y|) */

#ifndef PI
#define PI	3.1415926f
#endif

/* computing Local Binary Pattern element */
#define BIN_LBP_ELEM(c1, c2, offset)  (c2 > c1 ? c2 - c1 >= offset : c1 - c2 < offset )
//#define BIN_LBP_ELEM(c1, c2, offset)  (c2 > c1)

class CBuildFeatureImages
{
public:
	/* sampling feature image with a set of bins and ranges */
	void SetBinFeature(IplImage* feature_img, int bins_num, float *range);

	/* computing local binary pattern image */
	void ComputeLBP(IplImage *gray_img, IplImage *lbp_img, int radius=3, int neig_pt_num=6, float robust_noise=3.0f);

	/* convert image from uchar to float */
	//void ConvertImage(IplImage *src_uchar_Img, IplImage *dst_float_Img);

	/* computing position (x,y) maps */
	void ComputePositionMaps(CvSize img_size, IplImage *map_x, IplImage *map_y);

	/* computing edge orientation, arctan(|I_y|/|I_x|) */
	void ComputeEdgeOrientation(IplImage *grad_x, IplImage *grad_y, IplImage *edge_ori);

	/* computing first-order intensity gradient, sqrt(I_x^2+I_y^2) */
	void ComputeGradients(IplImage *grad_x, IplImage *grad_y, IplImage *grad_xy);

	/* computing second-order intensity gradients, I_xx, I_yy */
	void ComputeSecondOrderGradients(IplImage *grad_x, IplImage *grad_y, IplImage *grad_xx, IplImage *grad_yy);

	/* computing second-order intensity gradients, I_xx, I_yy */
	void ComputeSecondOrderGradients(IplImage *gray_img, IplImage *grad_xx, IplImage *grad_yy);

	/* computing first-order intensity gradients, I_x, I_y */
	void ComputeFirstOrderGradients(IplImage* gray_img, IplImage *grad_x, IplImage *grad_y);

	int m_nFeatureTypeNum;		/* the number of feature types */
	int m_nFeatureImageNum;		/* the dimension number of feature images */
	int* m_pFeatureImageTypes;	/* feature type vector */

	/* building feature images */
	void BuildFeatureImages(IplImage **feature_imgs, IplImage *org_img, IplImage *fg_prob_img = NULL, float fg_threshold = 0.15f, int one_existing_pos_feature = 0);

	/* get feature image number */
	int GetFeatureImageNumber();

	CBuildFeatureImages(int *feature_img_types, int feature_type_num);
	virtual ~CBuildFeatureImages();
};

#endif // !defined(_BUILD_FEATURE_IMAGES_H_)

