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
#ifndef _CROSS_BILATERAL_FILTER_H_
#define _CROSS_BILATERAL_FILTER_H_

//#define CHRONO
#include "geom.h"
#include "fast_lbf.h"
#include "cv.h"

typedef Array_2D<float> image_type;

class CCrossBilateralFilter
{
public:
	image_type m_iEdgeImage;		/* the edge image */
	image_type m_iInputImage;		/* the input image */
	image_type m_iFilteredImage;		/* the filtered image */
	
	float	m_fSigmaS;			/* the gaussian sigma in the spatial domain */
	float	m_fSigmaR;			/* the gaussian sigma in the intensity domain */
	
	CCrossBilateralFilter() 
	{
	}
	
	~CCrossBilateralFilter()
	{
		ClearData();	/* clear memories */
	}
    
	void ClearData()
	{
	}
 
	/* both edge and input images must be float type */
	void Initialization(IplImage *first_input_img, IplImage *first_edge_img, float sigma_s = 5.0f, float sigma_r = 0.1f) 
	{
		ClearData();
		
		if ( first_edge_img->depth != IPL_DEPTH_32F || first_input_img->depth != IPL_DEPTH_32F )
		{
		    printf("Error: both opencv edge and input images must be IPL_DEPTH_32F (float) type!\n");
		    exit(1);
		}
		
		if ( first_edge_img->nChannels != 1 || first_input_img->nChannels != 1 )
		{
		    printf("Error: both opencv edge and input images must be single channel images!\n");
		    exit(1);
		}
		
		m_fSigmaS = sigma_s;
		m_fSigmaR = sigma_r;
		
		CvSize img_size = cvGetSize(first_edge_img);
		
		// allocate memories
		m_iEdgeImage = image_type(img_size.width, img_size.height);
		m_iInputImage = image_type(img_size.width, img_size.height);
		m_iFilteredImage = image_type(img_size.width, img_size.height);
	}
    
	/* both edge and input images must be float type */
	void SetNewImages(IplImage *new_input_img, IplImage *new_edge_img)
	{
		if ( new_edge_img->depth != IPL_DEPTH_32F || new_input_img->depth != IPL_DEPTH_32F )
		{
		    printf("Error: both opencv edge and input images must be IPL_DEPTH_32F (float) type!\n");
		    exit(1);
		}
		
		if ( new_edge_img->nChannels != 1 || new_input_img->nChannels != 1 )
		{
		    printf("Error: both opencv edge and input images must be single channel images!\n");
		    exit(1);
		}
		
		SetImageData(new_edge_img, &m_iEdgeImage);
		SetImageData(new_input_img, &m_iInputImage);
	}
    
	/* do fast cross bilateral filter */
	void FastCrossBF()
	{
		Image_filter::fast_LBF( m_iInputImage, m_iEdgeImage,
				 m_fSigmaS, m_fSigmaR,
				 false,
				 &m_iFilteredImage, &m_iFilteredImage);
	}
    
	/* return the filtered image via cross bilateral filter */
	void GetFilteredImage(IplImage *filtered_img)
	{
		GetImageData(&m_iFilteredImage, filtered_img);
	}
    
    
private:
	/* get image data of IplImage from image_type */
	void GetImageData(image_type *src, IplImage *dst)
	{
		float *x_data;
		int x, y;
		for ( y = 0 ; y < dst->height ; y++ ) {
		    x_data = (float*)(dst->imageData + dst->widthStep*y);
		    for ( x = 0 ; x < dst->width ; x++ ) {
			*x_data++ = (*src)(x,y);
		    }
		}
	}
    
	/* set image data of IplImage to image_type */
	void SetImageData(IplImage *src, image_type *dst)
	{
		int x, y;
		float *x_data;
		for ( y = 0 ; y < src->height ; y++ ) {
		    x_data = (float*)(src->imageData + src->widthStep*y);
		    for ( x = 0 ; x < src->width ; x++ ) {
			(*dst)(x,y) = *x_data++;
		    }
		}
	}
};

#endif
