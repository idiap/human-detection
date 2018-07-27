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
#include "cv.h"

class Camera {
  private:
    void Scale(double x, double y);
    void VanishingPointHomography();
    void VanishingPointUndistortedImageROI();
    void UndistortedImageROI();
    CvRect GetROI(CvPoint2D64f *pts, int n_points);

  public:
	  void SaveUpdatedCalibration(char *filename);
    IplImage* current_image;     // pointer to current image, no image data copy!!
    IplImage* undistorted_image;
    IplImage* mask;

    CvMat* projection_matrix;
    CvMat* homography;
    CvMat* off_plane_homography; // to project points in the image on a plane parallel
                                 // to the ground plane, with z=off_plane_height

    CvMat* vanishing_point_homography;
    CvMat* inv_vanishing_point_homography;

    IplImage* mapx;              // used to store the precomputed undistortion map for x
    IplImage* mapy;              // used to store the precomputed undistortion map for y

    double k;                 // first order distortion corefficient
    double cx;                // position y of center = cx*image_width
    double cy;                // position y of center = cy*image_height
    double s;                 // set factor between x dimension and y dimension

    double off_plane_height;  // store the current height of the off plane homography

    bool vertical_vp_mapping_flag;
    bool distortion_removal_flag;

    Camera(IplImage* image);
    ~Camera();

    CvRect undist_roi;
    CvRect wrapped_roi;

    /*!Load the existing camera ibration file that consists of a 3x4 camera
       projection matrix and 4 projection parameters
    */
    void loadCalibration(char* filename);
    /*!Based on the estimated homography and the loaded k, cx, cy and s
       distortion parameters, makes the mapx and mapy LUT pixel positions
    */
    void WarpingMapXY();

    void WrappedToUndistorted(CvPoint2D64f* src, CvPoint2D64f* dst);
    void UndistortedToOriginal(CvPoint2D64f* src, CvPoint2D64f* dst);
    void OriginalToUndistorted(CvPoint2D64f* src, CvPoint2D64f* dst);
    void UndistortedToWrapped(CvPoint2D64f* src, CvPoint2D64f* dst);

    void WrappedToOriginal(CvPoint2D64f* src, CvPoint2D64f* dst);
    void OriginalToWrapped(CvPoint2D64f* src, CvPoint2D64f* dst);

    void getMask();

    /*!Undistort the current image based on mapx and mapy LUT
    */
    void undistort(IplImage* image, CvRect *roi=NULL);
    void undistort(IplImage* image, IplImage* undis_image, CvRect *roi=NULL);
    void computeTransforms(double scale_x=1.0, double scale_y=1.0, bool vp_flag=true);
    void DrawGrid(bool to_original=false);
    void DrawAxes(bool to_original=false);

    void Proj3DTo2D(CvPoint3D64f* src, CvPoint2D64f* dst, bool to_original=false);
    void Proj2DTo3D(CvPoint2D64f* src, CvPoint3D64f* dst, bool from_original=false, double height=0.0);

    void ComputeHomographyHeight(CvMat* mat, double height);
};
