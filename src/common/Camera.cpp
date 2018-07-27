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
#include "Camera.h"
#include "stdio.h"

Camera::Camera(IplImage* image = NULL){
  if(image == NULL){
    fprintf(stderr, "No image given for constructing the camera object, exiting...");
    exit(0);
  }

  // pointer to the image to be processed
  current_image = image;

  // this is crucial !!!
  undist_roi  = cvRect(0, 0, current_image->width, current_image->height);
  wrapped_roi = cvRect(0, 0, current_image->width, current_image->height);

  // will be allocated depending on the undistortion parameters
  // their size is unknown yet
  undistorted_image = NULL;
  mapx = NULL;
  mapy = NULL;
  mask = NULL;

  // Camera projection matrix, will be allocated
  // when loading the calibration file
  // the calibration is related to the *undistorted*
  // image 3D -> 2D
  projection_matrix = NULL;
  homography           = cvCreateMat(3, 3, CV_64F);
  off_plane_homography = cvCreateMat(3, 3, CV_64F);
  off_plane_height = 0;

  // homography to make vertical lines appear vertical
  // forward and backward transformation
  vanishing_point_homography     = cvCreateMat(3, 3, CV_64F);
  inv_vanishing_point_homography = cvCreateMat(3, 3, CV_64F);
}

Camera::~Camera(){
  if(undistorted_image)
    cvReleaseImage(&undistorted_image);
  if(projection_matrix)
    cvReleaseMat(&projection_matrix);
  if(homography)
    cvReleaseMat(&homography);
  if(inv_vanishing_point_homography)
    cvReleaseMat(&inv_vanishing_point_homography);
  if(vanishing_point_homography)
    cvReleaseMat(&vanishing_point_homography);
  if(mapx)
    cvReleaseImage(&mapx);
  if(mapy)
    cvReleaseImage(&mapy);
  if(mask)
    cvReleaseImage(&mask);
}

// load projection matrix and distortion coefficents from yml file
void Camera::loadCalibration(char* filename){
  CvFileStorage* yaml = cvOpenFileStorage( filename, 0, CV_STORAGE_READ);
  CvMat* distortion_params;
  if(yaml==NULL){
    printf("File %s could not be read\n", filename);
    exit(1);
  }
  else{
    projection_matrix   = (CvMat*) cvReadByName(yaml, NULL, "projection_matrix");
    distortion_params   = (CvMat*) cvReadByName(yaml, NULL, "distortion_coeffs");
  }

  if ( !distortion_params ) {
    distortion_removal_flag = false;
    k = 0;  s = 1.0;
    cx = cy = 0.5;
  }
  else {
    distortion_removal_flag = true;
    // copy parameters based on Jian's setup
    k  = cvmGet(distortion_params,0,0);
    cx = cvmGet(distortion_params,0,1);
    cy = cvmGet(distortion_params,0,2);
    s  = cvmGet(distortion_params,0,3);
    cvReleaseMat(&distortion_params);
  }

  cvReleaseFileStorage(&yaml);
}

void Camera::computeTransforms(double scale_x, double scale_y, bool vp_flag){
  vertical_vp_mapping_flag = vp_flag;
  UndistortedImageROI();
  VanishingPointHomography();
  if((scale_x != 1.0) || (scale_y != 1.0)){
    Scale(scale_x, scale_y);
  }
  VanishingPointUndistortedImageROI();
  WarpingMapXY();
}

void Camera::VanishingPointHomography(){
  if ( !vertical_vp_mapping_flag ) {
    cvSetZero(vanishing_point_homography);
    cvmSet(vanishing_point_homography, 0, 0, 1.0);
    cvmSet(vanishing_point_homography, 1, 1, 1.0);
    cvmSet(vanishing_point_homography, 2, 2, 1.0);

    cvSetZero(inv_vanishing_point_homography);
    cvmSet(inv_vanishing_point_homography, 0, 0, 1.0);
    cvmSet(inv_vanishing_point_homography, 1, 1, 1.0);
    cvmSet(inv_vanishing_point_homography, 2, 2, 1.0);

    return;
  }

  // T matrix moves the origin point to the image center
  // (x_o, y_o) = (undistorted_image->width/2, undistorted_image->height/2)
  CvMat T, R, G;
  double t[] = {1.0, 0.0,  -undist_roi.width*cx,
                0.0, 1.0, -undist_roi.height*cy,
                0.0, 0.0,                   1.0};
  cvInitMatHeader(&T, 3, 3, CV_64FC1, t);

  // get the vanishing point from the third column of the projection matrix
  // and translate it w.r.t the new origin
  CvMat P;
  double p[] = { cvmGet(projection_matrix, 0, 2),
                 cvmGet(projection_matrix, 1, 2),
                 cvmGet(projection_matrix, 2, 2)};
  cvInitMatHeader(&P, 3, 1, CV_64F, p);
  cvMatMul(&T, &P, &P);

  p[0] /= p[2];
  p[1] /= p[2];

  double norm = cvSqrt(p[0]*p[0] + p[1]*p[1]);
  if(p[1]<0)
    norm = -norm;

  double r[] = {p[1]/norm, -p[0]/norm, 0.0,
                p[0]/norm,  p[1]/norm, 0.0,
                      0.0,        0.0, 1.0};
  cvInitMatHeader(&R, 3, 3, CV_64FC1, r);

  double g[] = {1.0,       0.0, 0.0,
                0.0,       1.0, 0.0,
                0.0, -1.0/norm, 1.0};
  cvInitMatHeader(&G, 3, 3, CV_64FC1, g);

  CvMat* GR = cvCreateMat(3, 3, CV_64F);
  cvMatMul(&G, &R, GR);
  cvMatMul(GR, &T, vanishing_point_homography);
  cvInvert(vanishing_point_homography, inv_vanishing_point_homography);
  cvReleaseMat(&GR);
}

// Compute the size of the undistorted image
void Camera::UndistortedImageROI() {
  if ( !distortion_removal_flag ) {
    undist_roi = cvRect(0,0,current_image->width,current_image->height);
    return;
  }

  double cols = (double) current_image->width;
  double rows = (double) current_image->height;

  CvPoint2D64f src[8];
  src[0].x = 0.0;       src[0].y =      0.0;
  src[1].x = cols/2.0f; src[1].y =      0.0;
  src[2].x =      cols; src[2].y =      0.0;
  src[3].x =      cols; src[3].y = rows/2.0;
  src[4].x =      cols; src[4].y =     rows;
  src[5].x =  cols/2.0; src[5].y =     rows;
  src[6].x =       0.0; src[6].y =     rows;
  src[7].x =       0.0; src[7].y = rows/2.0;
  CvPoint2D64f dst[8];

  for(int i=0;i<8;i++){
    OriginalToUndistorted(&src[i], &dst[i]);
  }
  undist_roi = GetROI(dst, 8);
}

void Camera::VanishingPointUndistortedImageROI(){
  CvPoint2D64f src[4];
  src[0].x =              0.0; src[0].y =               0.0;
  src[1].x = undist_roi.width; src[1].y = undist_roi.height;
  src[2].x =              0.0; src[2].y = undist_roi.height;
  src[3].x = undist_roi.width; src[3].y =               0.0;
  CvPoint2D64f dst[4];
  for(int i=0;i<4;i++){
    UndistortedToWrapped(&src[i], &dst[i]);
  }

  wrapped_roi = GetROI(dst, 4);
  wrapped_roi.x += cvRound(undist_roi.width*cx);
  wrapped_roi.y += cvRound(undist_roi.height*cy);

  // allocate image and LUTs memory
  undistorted_image = cvCreateImage(cvSize(wrapped_roi.width, wrapped_roi.height), IPL_DEPTH_8U, 3);
  mapx = cvCreateImage(cvSize(wrapped_roi.width, wrapped_roi.height), IPL_DEPTH_32F, 1);
  mapy = cvCreateImage(cvSize(wrapped_roi.width, wrapped_roi.height), IPL_DEPTH_32F, 1);
  mask = cvCreateImage(cvSize(wrapped_roi.width, wrapped_roi.height), IPL_DEPTH_8U, 1);

  // update vanishing point matrix and camera projection matrix
  CvMat T;
  double t[] = {1.0, 0.0, -undist_roi.width*cx  + wrapped_roi.x,
                0.0, 1.0, -undist_roi.height*cy + wrapped_roi.y,
                0.0, 0.0,                                 1.0};
  cvInitMatHeader(&T, 3, 3, CV_64FC1, t);
  cvMatMul( inv_vanishing_point_homography, &T, inv_vanishing_point_homography);

  cvInvert(inv_vanishing_point_homography, vanishing_point_homography);
  cvMatMul(vanishing_point_homography, projection_matrix, projection_matrix);

  ComputeHomographyHeight(homography, 0);

  //cvmSet(homography, 0, 0, cvmGet(projection_matrix, 0, 0));
  //cvmSet(homography, 0, 1, cvmGet(projection_matrix, 0, 1));
  //cvmSet(homography, 0, 2, cvmGet(projection_matrix, 0, 3));
  //cvmSet(homography, 1, 0, cvmGet(projection_matrix, 1, 0));
  //cvmSet(homography, 1, 1, cvmGet(projection_matrix, 1, 1));
  //cvmSet(homography, 1, 2, cvmGet(projection_matrix, 1, 3));
  //cvmSet(homography, 2, 0, cvmGet(projection_matrix, 2, 0));
  //cvmSet(homography, 2, 1, cvmGet(projection_matrix, 2, 1));
  //cvmSet(homography, 2, 2, cvmGet(projection_matrix, 2, 3));
  //cvInvert(homography, homography);
}

// Computes the mapping between the final undistorted +
// vanishing point corrected image and the original image
void Camera::WarpingMapXY(){
  CvPoint2D64f src, dst;
  for ( int y = 0 ; y < undistorted_image->height; y++ ) {
    for ( int x = 0 ; x < undistorted_image->width; x++ ) {
      src = cvPoint2D64f(x,y);
      WrappedToUndistorted(&src, &dst);
      UndistortedToOriginal(&dst, &dst);

      cvSetReal2D(mapx, y, x, dst.x);
      cvSetReal2D(mapy, y, x, dst.y);
    }
  }
  getMask();
}

// Draw a grid from the 3D world
// on the current view (to verify the camera calibration)
void Camera::DrawGrid(bool to_original){
  CvPoint2D64f p_2D;
  CvPoint3D64f p_3D;
  p_3D.z = 180;
  for(int i=-850;i<850;i=i+50){
    for(int j=-200;j<550;j=j+50){
      p_3D.x = i;
      p_3D.y = j;
      Proj3DTo2D(&p_3D, &p_2D, to_original);
      if(!to_original){
        cvCircle(undistorted_image, cvPoint((int) p_2D.x, (int) p_2D.y), 1, cvScalar(60,60,60), 1);
      }
      else{
        cvCircle(current_image, cvPoint((int) p_2D.x, (int) p_2D.y), 1, cvScalar(60,60,60), 1);
      }
    }
  }
}

void Camera::DrawAxes(bool to_original){
  CvPoint3D64f p_3D[4];
  p_3D[0] = cvPoint3D64f( 0,  0,  0);
  p_3D[1] = cvPoint3D64f(80,  0,  0);
  p_3D[2] = cvPoint3D64f( 0, 80,  0);
  p_3D[3]= cvPoint3D64f( 0,  0, 80);

  CvPoint2D64f p_2D[5];          // fifth is the origin projected onto the original view
  Proj3DTo2D(&p_3D[0], &p_2D[0]);

  WrappedToUndistorted(&p_2D[0], &p_2D[4]);
  UndistortedToOriginal(&p_2D[4], &p_2D[4]);
  for(int i=1;i<4;i++){
    Proj3DTo2D(&p_3D[i], &p_2D[i], to_original);
    if(!to_original){
      cvLine(undistorted_image, cvPoint(cvRound(p_2D[0].x), cvRound(p_2D[0].y)),
                                cvPoint(cvRound(p_2D[i].x), cvRound(p_2D[i].y)),
                                CV_RGB(0, 0, 0), 1);
    }
    else{
      cvLine(current_image, cvPoint(cvRound(p_2D[4].x), cvRound(p_2D[4].y)),
                            cvPoint(cvRound(p_2D[i].x), cvRound(p_2D[i].y)),
                            CV_RGB(0, 0, 0), 1);
    }
   }
}

void Camera::undistort(IplImage* image, CvRect *roi){
  current_image = image;

  if ( roi && (roi->width<=0 || roi->height<=0) )
	  return;

  if ( roi ) {
    cvSetImageROI(undistorted_image, *roi);
    cvSetImageROI(mapx, *roi);
    cvSetImageROI(mapy, *roi);
  }

  if ( distortion_removal_flag || vertical_vp_mapping_flag )
	  cvRemap(current_image, undistorted_image,
		mapx, mapy,
		CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS,
		cvScalarAll(0));
  else
	  cvResize(image, undistorted_image);

  if ( roi ) {
    cvResetImageROI(undistorted_image);
    cvResetImageROI(mapx);
    cvResetImageROI(mapy);
  }
}

void Camera::undistort(IplImage* image, IplImage* undis_image, CvRect *roi) {
  current_image = image;

  if ( roi && (roi->width<=0 || roi->height<=0) )
	  return;

  if ( roi ) {
    cvSetImageROI(undis_image, *roi);
    cvSetImageROI(mapx, *roi);
    cvSetImageROI(mapy, *roi);
  }

  if ( distortion_removal_flag || vertical_vp_mapping_flag )
	  cvRemap(current_image, undis_image,
		mapx, mapy,
		CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS,
		cvScalarAll(0));
  else
	  cvResize(image, undis_image);

  if ( roi ) {
    cvResetImageROI(undis_image);
    cvResetImageROI(mapx);
    cvResetImageROI(mapy);
  }
}

void Camera::WrappedToUndistorted(CvPoint2D64f* src, CvPoint2D64f* dst){
  CvMat point;
  double p[] = {src->x, src->y, 1.0};
  cvInitMatHeader(&point, 3, 1, CV_64F, p);
  cvMatMul(inv_vanishing_point_homography, &point, &point);
  dst->x = p[0]/p[2];
  dst->y = p[1]/p[2];
}

void Camera::UndistortedToOriginal(CvPoint2D64f* src, CvPoint2D64f* dst){
  if ( !distortion_removal_flag ) {
    dst->x = src->x;
    dst->y = src->y;
    return;
  }

  CvMat coeffs;
  double c[] = {0.0, 0.0, 1.0, 0.0};
  cvInitMatHeader(&coeffs, 1, 4, CV_64F, c);

  CvMat roots;
  double r[] = {0.0, 0.0, 0.0};
  cvInitMatHeader(&roots, 1, 3, CV_64F, r);

  //adding the x and y position of the computed rectangle
  double xu = src->x + (double) undist_roi.x;
  double yu = src->y + (double) undist_roi.y;

  //double ru = sqrt( pow((xu-cx*current_image->width)/s, 2.0) + pow((yu-cy*current_image->height), 2.0) ); // Undistorted radius
  double xc = (xu-cx*current_image->width)/s;
  double yc = (yu-cy*current_image->height);
  double ru = sqrt( xc*xc + yc*yc ); // Undistorted radius

  double rmax = sqrt(current_image->width*current_image->width/4.0 + current_image->height*current_image->height/4.0);

  // Normalize radius so that it varies from 0 at the centre to a maximum
  // of 1 at the corners.
  ru /= rmax;
  double rd = ru;

  cvmSet(&coeffs, 0, 0, k);
  cvmSet(&coeffs, 0, 3, -rd);
  double eig_val = 1000000;
  if(cvSolveCubic( &coeffs, &roots ) > 0){
    for(int i=0;i<3;i++){
      if(cvmGet(&roots, 0, i) < eig_val && cvmGet(&roots, 0, i) > 0){
        eig_val = cvmGet(&roots, 0, i);
      }
    }
  }
  if ( eig_val >= 0 )
    rd = eig_val;

  ru += FLT_EPSILON;

  dst->x = cx*undist_roi.width  + (src->x - cx*undist_roi.width)*rd/ru;
  dst->y = cy*undist_roi.height + (src->y - cy*undist_roi.height)*rd/ru;

  dst->x += undist_roi.x;
  dst->y += undist_roi.y;
}

void Camera::WrappedToOriginal(CvPoint2D64f* src, CvPoint2D64f* dst)
{
  CvPoint2D64f undist_pt;
  WrappedToUndistorted(src, &undist_pt);
  UndistortedToOriginal(&undist_pt, dst);
}

void Camera::OriginalToWrapped(CvPoint2D64f* src, CvPoint2D64f* dst)
{
  CvPoint2D64f undist_pt;
  OriginalToUndistorted(src, &undist_pt);
  UndistortedToWrapped(&undist_pt, dst);
}

void Camera::OriginalToUndistorted(CvPoint2D64f* src, CvPoint2D64f* dst){
  if ( !distortion_removal_flag ) {
    dst->x = src->x;
    dst->y = src->y;
    return;
  }

  double rmax2 = current_image->width*current_image->width/4.0 +
                current_image->height*current_image->height/4.0;

  /*
  double rd2 = pow((src->x-cx*current_image->width)/s, 2.0) +
               pow((src->y-cy*current_image->height), 2.0);
  */
  double xc = (src->x-cx*current_image->width)/s;
  double yc = (src->y-cy*current_image->height);
  double rd2 = xc*xc + yc*yc;

  rd2 *= k;
  rd2 /= rmax2;

  dst->x = src->x + (src->x - cx*current_image->width)*rd2 - undist_roi.x;
  dst->y = src->y + (src->y - cy*current_image->height)*rd2 - undist_roi.y;
}

void Camera::UndistortedToWrapped(CvPoint2D64f* src, CvPoint2D64f* dst){
  CvMat point;
  double p[] = {src->x, src->y, 1.0};
  cvInitMatHeader(&point, 3, 1, CV_64F, p);
  cvMatMul(vanishing_point_homography, &point, &point);
  dst->x = p[0]/p[2];
  dst->y = p[1]/p[2];
}

void Camera::ComputeHomographyHeight(CvMat* homography_mat, double height){
  if(height != 0)
    off_plane_height= height;

  cvmSet(homography_mat, 0, 0, cvmGet(projection_matrix, 0, 0));
  cvmSet(homography_mat, 0, 1, cvmGet(projection_matrix, 0, 1));
  cvmSet(homography_mat, 0, 2, cvmGet(projection_matrix, 0, 3)
         + height*cvmGet(projection_matrix, 0, 2));
  cvmSet(homography_mat, 1, 0, cvmGet(projection_matrix, 1, 0));
  cvmSet(homography_mat, 1, 1, cvmGet(projection_matrix, 1, 1));
  cvmSet(homography_mat, 1, 2, cvmGet(projection_matrix, 1, 3)
         + height*cvmGet(projection_matrix, 1, 2));
  cvmSet(homography_mat, 2, 0, cvmGet(projection_matrix, 2, 0));
  cvmSet(homography_mat, 2, 1, cvmGet(projection_matrix, 2, 1));
  cvmSet(homography_mat, 2, 2, cvmGet(projection_matrix, 2, 3)
         + height*cvmGet(projection_matrix, 2, 2));

  cvInvert(homography_mat, homography_mat);
}

void Camera::Proj2DTo3D(CvPoint2D64f* src, CvPoint3D64f* dst, bool from_original, double height){
  CvPoint2D64f tmp;
  double p1[3];
  if(from_original){
    OriginalToUndistorted(src, &tmp);
    UndistortedToWrapped(&tmp, &tmp);
    p1[0] = tmp.x; p1[1] = tmp.y; p1[2] = 1.0;
  }
  else{
    p1[0] = src->x; p1[1] = src->y; p1[2] = 1.0;
  }

  if((off_plane_height != height) & (height!= 0)){
    ComputeHomographyHeight(off_plane_homography, height);
  }

  double p2[] = {0.0, 0.0, 0.0};
  CvMat point2D, point3D;
  cvInitMatHeader(&point2D, 3, 1, CV_64F, p1);
  cvInitMatHeader(&point3D, 3, 1, CV_64F, p2);
  if(height != 0){
    cvMatMul(off_plane_homography, &point2D, &point3D);
  }
  else{
    cvMatMul(homography, &point2D, &point3D);
  }

  dst->x = p2[0]/p2[2];
  dst->y = p2[1]/p2[2];
  dst->z = height;
}

void Camera::Proj3DTo2D(CvPoint3D64f* src, CvPoint2D64f* dst, bool to_original){
  double p1[] = {src->x, src->y, src->z, 1.0};
  double p2[] = {0.0, 0.0, 0.0};
  CvMat point2D, point3D;
  cvInitMatHeader(&point3D, 4, 1, CV_64F, p1);
  cvInitMatHeader(&point2D, 3, 1, CV_64F, p2);
  cvMatMul(projection_matrix, &point3D, &point2D);

  dst->x = p2[0]/p2[2];
  dst->y = p2[1]/p2[2];
  if(to_original){
    WrappedToUndistorted(dst, dst);
    UndistortedToOriginal(dst, dst);
  }
}

CvRect Camera::GetROI(CvPoint2D64f* pts, int n_points){
  double colmin = pts[0].x, colmax = pts[0].x;
  double rowmin = pts[0].y, rowmax = pts[0].y;
  for (int i = 1 ; i < n_points; i++ ) {
    colmin = MIN(colmin, pts[i].x);
    colmax = MAX(colmax, pts[i].x);
    rowmin = MIN(rowmin, pts[i].y);
    rowmax = MAX(rowmax, pts[i].y);
  }

  int icolmin = cvFloor(colmin);
  int icolmax = cvCeil(colmax) ;
  int irowmin = cvFloor(rowmin);
  int irowmax = cvCeil(rowmax) ;

  return cvRect(icolmin, irowmin, icolmax - icolmin, irowmax - irowmin);
}

void Camera::Scale(double x, double y){
  CvMat scale;
  double s[] = {  x, 0.0, 0.0,
                0.0,   y, 0.0,
                0.0, 0.0, 1.0};
  cvInitMatHeader(&scale, 3, 3, CV_64F, s);
  cvMatMul(&scale, vanishing_point_homography, vanishing_point_homography);
  cvInvert(vanishing_point_homography, inv_vanishing_point_homography);
}
void Camera::getMask(){
  cvSet(mask, cvScalar(0));
  for(int i=0;i<mapx->width;i++){
    for(int j=0;j<mapx->height;j++){
      if((cvGetReal2D(mapx,j,i) <  current_image->width) &
         (cvGetReal2D(mapx,j,i) >=                    0) &
         (cvGetReal2D(mapy,j,i) < current_image->height) &
         (cvGetReal2D(mapy,j,i) >=                    0)){
        cvSetReal2D(mask, j, i, 255);
      }
    }
  }
}

void Camera::SaveUpdatedCalibration(char *filename)
{
  CvFileStorage* yaml = cvOpenFileStorage( filename, 0, CV_STORAGE_WRITE);
  if(yaml==NULL){
    printf("File %s could not be written\n", filename);
    exit(1);
  }

  cvWriteComment( yaml, "image size", 0);

  cvWriteInt( yaml, "image_width", mask->width);
  cvWriteInt( yaml, "image_height", mask->height);

  cvWriteComment( yaml, "the 3x4 projection matrix of warped image after distortion removal and infinite vp mapping", 0);
  cvWrite( yaml, "projection_matrix", projection_matrix, cvAttrList(0,0) );

  cvReleaseFileStorage(&yaml);
}
