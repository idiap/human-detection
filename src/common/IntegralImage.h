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
// IntegralMaskImage.h: interface for the CIntegralMaskImage class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_INTEGRAL_MASK_IMAGE_H_)
#define _INTEGRAL_MASK_IMAGE_H_

#include <cstdio>
#include "stdlib.h"
#include <memory>
#include <cmath>

#include <ctime>						// clock
#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <fstream>						// file I/O
#include <cmath>						// math includes
#include <iostream>						// I/O streams

using namespace std;

#include "cv.h"

#include "OpencvDataConversion.h"


class CIntegralImage
{
public:
	void SetIntegralROI(CvRect roi);
	float GetRegionIntegralSum(CvRect rect_roi);
	void ComputeIntegralImage(IplImage *img);
	CIntegralImage(CvSize img_size, bool is_mask_img=true);
	virtual ~CIntegralImage();

	float *m_pIntegralData;
	CvSize m_szImage;
	bool m_bMaskImg;
	IplImage *m_pImage;

	CvRect m_cvIntegralROI;
};

#endif // !defined(_INTEGRAL_MASK_IMAGE_H_)
