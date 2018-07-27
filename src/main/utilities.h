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

/* get next image frame data from images name list or a video file */
static bool get_next_image(CvCapture* video, CFileList* list, Camera* cam,
                           IplImage* image, double noise_sigma,
                           int cur_frame_idx, int n_skipped_frames) {
    IplImage* org_image = NULL;
    if ( video ) {
        org_image = cvQueryFrame(video);
        if ( !org_image )
            return false;
#ifdef _WIN32
        // commented by remi emonet, because it does not seem to be a real need
        //        if ( org_image ) {
        //            cvFlip(org_image, org_image, 0);
        //            org_image->origin = 0;
        //        }
#endif
        for ( int i = 0 ; i < n_skipped_frames ; i++ )
            cvGrabFrame(video);
    }
    else {
        if ( cur_frame_idx >= list->GetListLength() )
            return false;
        org_image = cvLoadImage(list->GetFileName(cur_frame_idx));
    }
    
    if ( cam ) {
        cam->undistort(org_image);
        cvCopy(cam->undistorted_image, image);
    }
    else
        cvResize(org_image, image);
    
    if ( noise_sigma > 0 )
        cvSmooth(image, image, CV_GAUSSIAN, 0, 0, noise_sigma);
    
    if ( !video )
        cvReleaseImage(&org_image);
    return true;
}


static char* get_file_name(const char* data_dir, const char* file_prefix, const char* file_ext, int idx, bool no_underscore = false) {
    char *file_name = new char[1024];
    
    if ( data_dir[strlen(data_dir)-1] == '/' ||
         data_dir[strlen(data_dir)-1] == '\\' )
    {
    	if (no_underscore)
    		sprintf(file_name, "%s%s%06d.%s", data_dir, file_prefix, idx, file_ext);
    	else
    		sprintf(file_name, "%s%s_%06d.%s", data_dir, file_prefix, idx, file_ext);
    }
    else
    {
#ifndef _WIN32  /* under Linux system */
    	if (no_underscore)
    		sprintf(file_name, "%s/%s%06d.%s", data_dir, file_prefix, idx, file_ext);
    	else
    		sprintf(file_name, "%s/%s_%06d.%s", data_dir, file_prefix, idx, file_ext);
#else           /* under Window system */
    if (no_underscore)
    	sprintf(file_name, "%s\\%s%06d.%s", data_dir, file_prefix, idx, file_ext);
    else
    	sprintf(file_name, "%s\\%s_%06d.%s", data_dir, file_prefix, idx, file_ext);
#endif
    }
    return file_name;
}


void skip_first_frames(CvCapture* VIDEO, int start_frame_idx) {
#ifdef _win32
    cvSetCaptureProperty(VIDEO, CV_CAP_PROP_POS_FRAMES, start_frame_idx);
#else
    for (int i = 0; i < start_frame_idx; i++)
        cvGrabFrame(VIDEO);
#endif
}
