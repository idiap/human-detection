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
// Rectangle.h: interface for the CRectangle class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_RECTANGLE_H_)
#define _RECTANGLE_H_

#include "cv.h"

// CvRect equivalent with floating values
struct floatRect {
    float x;
    float y;
    float width;
    float height;
};

class CRectangle
{
public:
    /* return the intersection of two OpenCV rectangles */
	CvRect Intersection(CvRect rect1, CvRect rect2)
	{
		CvPoint r1pt1 = cvPoint(rect1.x, rect1.y);
		CvPoint r1pt2 = cvPoint(rect1.x + rect1.width, rect1.y + rect1.height);

		CvPoint r2pt1 = cvPoint(rect2.x, rect2.y);
		CvPoint r2pt2 = cvPoint(rect2.x + rect2.width, rect2.y + rect2.height);

		CvPoint pt1, pt2;

		/* Test if rectangle overlaps, if not return empty Rectangle */
		pt1 = cvPoint(MAX(r1pt1.x, r2pt1.x), MAX(r1pt1.y, r2pt1.y));
		pt2 = cvPoint(MIN(r1pt2.x, r2pt2.x), MIN(r1pt2.y, r2pt2.y));

		CvRect r = cvRect(pt1.x, pt1.y, pt2.x - pt1.x, pt2.y - pt1.y);

		if ( r.width <= 0 || r.height <= 0 )
			r = cvRect(0,0,0,0);   // I expect to be empty

		return r;
	};

    /* return the intersection of one OpenCV rectangle
       and another whole image boundary */
	CvRect Intersection(CvRect rect, CvSize img_size) {
		return Intersection(rect, cvRect(0,0,img_size.width,img_size.height));
	};

    /* return the union of two OpenCV rectangles */
	CvRect Union(CvRect rect1, CvRect rect2) {
		CvPoint r1pt1 = cvPoint(rect1.x, rect1.y);
		CvPoint r1pt2 = cvPoint(rect1.x + rect1.width, rect1.y + rect1.height);

		CvPoint r2pt1 = cvPoint(rect2.x, rect2.y);
		CvPoint r2pt2 = cvPoint(rect2.x + rect2.width, rect2.y + rect2.height);

		CvPoint pt1, pt2;
		pt1 = cvPoint(MIN(r1pt1.x, r2pt1.x), MIN(r1pt1.y, r2pt1.y));
		pt2 = cvPoint(MAX(r1pt2.x, r2pt2.x), MAX(r1pt2.y, r2pt2.y));

		CvRect r = cvRect(pt1.x, pt1.y, pt2.x - pt1.x, pt2.y - pt1.y);

		return r;
	};

    /* return the intersection of one OpenCV rectangle
       and another whole image boundary */
	CvRect Union(CvRect rect, CvSize img_size) {
		return Union(rect, cvRect(0,0,img_size.width,img_size.height));
	};

    /* return the integer rectangle from a normalized float rectangle,
       whose parent window is the whole image boundary */
	CvRect GetIntRect(floatRect float_rect, CvSize img_size) {
	    return GetIntRect(float_rect, cvRect(0,0,img_size.width,img_size.height));
	};

    /* return the integer rectangle from a normalized float rectangle,
       whose parent window is another integer rectangle */
	CvRect GetIntRect(floatRect float_rect, CvRect rect_roi) {
        CvRect int_rect;

        int_rect.x = rect_roi.x + cvFloor(float_rect.x*(float)rect_roi.width);
        int_rect.y = rect_roi.y + cvFloor(float_rect.y*(float)rect_roi.height);

        int_rect.x = MAX(int_rect.x, rect_roi.x);
        int_rect.y = MAX(int_rect.y, rect_roi.y);

        int_rect.x = MIN(int_rect.x, rect_roi.x+rect_roi.width-1);
        int_rect.y = MIN(int_rect.y, rect_roi.y+rect_roi.height-1);

        int_rect.width = cvCeil(float_rect.width*(float)rect_roi.width);
        int_rect.height = cvCeil(float_rect.height*(float)rect_roi.height);

        int_rect.width = MAX(int_rect.width, 1);
        int_rect.height = MAX(int_rect.height, 1);

        int_rect.width = MIN(int_rect.width, rect_roi.x+rect_roi.width-int_rect.x);
        int_rect.height = MIN(int_rect.height, rect_roi.y+rect_roi.height-int_rect.y);

        return int_rect;
	};

	/* get the rectangle area */
	int Area(CvRect rect) {
	    return (rect.width*rect.height);
	};

	float FArea(CvRect rect) {
	    return (float)(rect.width*rect.height);
	};

	CRectangle() {};
	virtual ~CRectangle() {};
};

#endif // !defined(_RECTANGLE_H_)


