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
// IntegralMaskImage.cpp: implementation of the CIntegralImage class.
//
//////////////////////////////////////////////////////////////////////

#include "IntegralImage.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CIntegralImage::CIntegralImage(CvSize img_size, bool is_mask_img)
{
	m_szImage = img_size;
	int length = (m_szImage.width+1)*(m_szImage.height+1);
	m_pIntegralData = new float[length];
	memset(m_pIntegralData, 0, sizeof(float)*length);
	m_bMaskImg = is_mask_img;

	m_cvIntegralROI = cvRect(0,0,m_szImage.width,m_szImage.height);

	m_pImage = NULL;
}

CIntegralImage::~CIntegralImage()
{
	delete [] m_pIntegralData;
}

void CIntegralImage::ComputeIntegralImage(IplImage *img)
{
	m_pImage = img;

	if ( m_cvIntegralROI.width == 0 || m_cvIntegralROI.height == 0 )
		return;

	float left_data;
	float *cur_integral_pt, *top_integral_pt;
	int x, y, a;

	cur_integral_pt = m_pIntegralData+m_szImage.width+2;
	top_integral_pt = m_pIntegralData+1;

	cvSetImageROI(img, m_cvIntegralROI);

	int xoffset = m_szImage.width - m_cvIntegralROI.width + 1;

	if ( img->depth == IPL_DEPTH_8U ) {
		COpencvDataConversion<uchar,uchar> ODC;
		uchar* img_data;
		uchar* cur_img_data;
		if ( img )
			img_data = ODC.GetImageData(img);
		else {
			img_data = new uchar[m_cvIntegralROI.width*m_cvIntegralROI.height];
			cur_img_data = img_data;
			for ( a = 0 ; a < m_cvIntegralROI.width*m_cvIntegralROI.height ; a++ )
				*cur_img_data++ = 1;
		}
		
		cur_img_data = img_data;
		
		if ( m_bMaskImg ) {
			for ( y = 1 ; y <= m_cvIntegralROI.height ; y++ ) {
				left_data = 0;
				for ( x = 1 ; x <= m_cvIntegralROI.width ; x++ ) {
					left_data += (*cur_img_data++?1:0);
					*cur_integral_pt++ = left_data + *top_integral_pt++;
				}
				top_integral_pt += xoffset;
				cur_integral_pt += xoffset;
			}
		}
		else {
			for ( y = 1 ; y <= m_cvIntegralROI.height ; y++ ) {
				left_data = 0;
				for ( x = 1 ; x <= m_cvIntegralROI.width ; x++ ) {
					left_data += *cur_img_data++;
					*cur_integral_pt++ = left_data + *top_integral_pt++;
				}
				top_integral_pt += xoffset;
				cur_integral_pt += xoffset;
			}
		}
		
		delete [] img_data;
	}
	else {
		COpencvDataConversion<float,float> ODC;
		float* img_data;
		float* cur_img_data;
		if ( img )
			img_data = ODC.GetImageData(img);
		else {
			img_data = new float[m_cvIntegralROI.width*m_cvIntegralROI.height];
			cur_img_data = img_data;
			for ( a = 0 ; a < m_cvIntegralROI.width*m_cvIntegralROI.height ; a++ )
				*cur_img_data++ = 1;
		}
		
		cur_img_data = img_data;
		
		if ( m_bMaskImg ) {
			for ( y = 1 ; y <= m_cvIntegralROI.height ; y++ ) {
				left_data = 0;
				for ( x = 1 ; x <= m_cvIntegralROI.width ; x++ ) {
					left_data += (*cur_img_data++?1:0);
					*cur_integral_pt++ = left_data + *top_integral_pt++;
				}
				top_integral_pt += xoffset;
				cur_integral_pt += xoffset;
			}
		}
		else {
			for ( y = 1 ; y <= m_cvIntegralROI.height ; y++ ) {
				left_data = 0;
				for ( x = 1 ; x <= m_cvIntegralROI.width ; x++ ) {
					left_data += *cur_img_data++;
					*cur_integral_pt++ = left_data + *top_integral_pt++;
				}
				top_integral_pt += xoffset;
				cur_integral_pt += xoffset;
			}
		}
		
		delete [] img_data;
	}

	cvResetImageROI(img);
}

float CIntegralImage::GetRegionIntegralSum(CvRect rect_roi)
{
	if ( m_cvIntegralROI.width == 0 || m_cvIntegralROI.height == 0 )
		return 0;

	int x1, y1, x2, y2;
	x1 = rect_roi.x - m_cvIntegralROI.x;
	y1 = rect_roi.y - m_cvIntegralROI.y;
	x2 = x1 + rect_roi.width;
	y2 = y1 + rect_roi.height;	
	
	if ( x1 < 0 || y1 < 0 || x2 > m_szImage.width || y2 > m_szImage.height ) {
		printf("\nROI [%d,%d,%d,%d] must be inside the images!\n", rect_roi.x, rect_roi.y, rect_roi.width, rect_roi.height);
		printf("m_cvIntegralROI = [%d, %d, %d, %d]\n", m_cvIntegralROI.x, m_cvIntegralROI.y, m_cvIntegralROI.width, m_cvIntegralROI.height);
		printf("m_szImage = [%d, %d]\n", m_szImage.width, m_szImage.height);
		return -1.0f;
	}

	int width_1 = m_szImage.width+1;
	
	int left_top_pos = y1*width_1 + x1;
	int right_top_pos = y1*width_1 + x2;
	int left_bot_pos = y2*width_1 + x1;
	int right_bot_pos = y2*width_1 + x2;
	
	float mask_area = m_pIntegralData[right_bot_pos] + m_pIntegralData[left_top_pos] 
		- m_pIntegralData[right_top_pos] - m_pIntegralData[left_bot_pos];

	return mask_area;
}

void CIntegralImage::SetIntegralROI(CvRect roi)
{
	m_cvIntegralROI = roi;

}
