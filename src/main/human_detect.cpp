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

#include "MultiLayerBGS.h"
#include "CmdLine.h"
#include "cv.h"
#include "highgui.h"
#include "FileList.h"
#include "BoostClassifierCovariance.h"
#include "Timer.h"
#include "TestWindows.h"
#include "NegExamplesGenerator.h"
#include "ImageDetection.h"
#include "MergingDetection.h"
#include "Camera.h"

#include "utilities.h"

#ifndef CHECK_STRING
#define CHECK_STRING(x)  {if(x && strlen(x)==0) {x=NULL;} }
#endif

/* get next image frame data from images name list or a video file */
bool get_next_image(CvCapture* video, CFileList* img_list, CFileList* fp_img_list,
        Camera* cam, CMultiLayerBGS* bgs,
        IplImage* image, IplImage* fp_image, double noise_sigma,
        int cur_frame_idx, int n_skipped_frames);

int main (int argc, char* argv[])
{
    Torch::CmdLine cmd;

    cmd.addText("\n******************************************************************************");
    cmd.addText("**   PROGRAM: Human Detection on a Set of Images or a Video Stream          **");
    cmd.addText("**                                                                          **");
    cmd.addText("**                                                                          **");
    cmd.addText("**             Algorithm:            Dr Jian Yao & Jean-Marc Odobez         **");
    cmd.addText("**             Main program author:  Dr. Jian Yao                           **");
    cmd.addText("**                       IDIAP Research Institute                           **");
    cmd.addText("**                       CH-1920, Martigny, Switzerland                     **");
    cmd.addText("**             Emails:   odobez@idiap.ch                                    **");
    cmd.addText("**                       jianyao75@gmail.com                                **");
    cmd.addText("**                                                                          **");
    cmd.addText("**             Date:     November 02, 2008                                  **");
    cmd.addText("******************************************************************************");

    cmd.addText("\nArguments:");
    char*   logit_boost_fn;
    char*   images_list_or_video_fn;
    cmd.addSCmdArg("logit boost model", &logit_boost_fn,
			"the file name of logit boost model (XML or YAML format) ");
    cmd.addSCmdArg("images list or video", &images_list_or_video_fn,
			"the list of images or an video to be tested ");

    cmd.addText("\nOptions:");

    cmd.addText("\n *Frames selection:");
    int     detection_start_frame_idx = 0;
    int     detection_end_frame_idx = -1;
    int     skipped_frames_num = 0;
    cmd.addICmdOption("-sf", &detection_start_frame_idx, 0,
            "the start frame index for human detection");
    cmd.addICmdOption("-ef", &detection_end_frame_idx, -1,
            "the end frame index for human detection");
    cmd.addICmdOption("-sfn", &skipped_frames_num, 0,
            "the number of skipped frames between two used frames for detection");

    cmd.addText("\n *Foreground map:");

    char*   fp_images_list_fn = NULL;
    cmd.addSCmdOption("-fplist", &fp_images_list_fn, "",
                        "the images file list of the foreground probability images corresponding to the tested images");

    cmd.addText("\n *Pre-processing options:");

    real    width_scale = 1.0;
    real    height_scale = 1.0;
    real    noise_sigma = 0;
    cmd.addRCmdOption("-xs", &width_scale, 1.0,
            "the image scaling factor at the x-direction (width)");
    cmd.addRCmdOption("-ys", &height_scale, 1.0,
            "the image scaling factor at the y-direction (height)");
    cmd.addRCmdOption("-nsig", &noise_sigma, 0.0,
                        "the Gaussian noise sigma used to remove the noise on the input image");

    cmd.addText("");
    int     left_marg_space = 0;
    int     right_marg_space = 0;
    int     top_marg_space = 0;
    int     bottom_marg_space = 0;
    cmd.addICmdOption("-lms", &left_marg_space, 0,
                        "the expanded left-margin space");
    cmd.addICmdOption("-rms", &right_marg_space, 0,
                        "the expanded right-margin space");
    cmd.addICmdOption("-tms", &top_marg_space, 0,
                        "the expanded top-margin space");
    cmd.addICmdOption("-bms", &bottom_marg_space, 0,
                        "the expanded bottom-margin space");

    cmd.addText("\n *Foreground detection options:");
    char*   bg_model_fn = NULL;
    real    mode_learn_rate_per_second = 0.01;
    real    weight_learn_rate_per_second = 0.01;
    real    init_mode_weight = 0.001;
    real    frame_duration = 1.0/25.0;
    bool    nolearn = false;
    bool    nobilateral = false;
    cmd.addSCmdOption("-bgm", &bg_model_fn, "",
                    "the learned background model file");
    cmd.addRCmdOption("-mlr", &mode_learn_rate_per_second, 0.01,
            "the learning rate per second for the background modes in the detection processing");
    cmd.addRCmdOption("-wlr", &weight_learn_rate_per_second, 0.01,
            "the learning rate per second for the background mode weights in the detection processing");
    cmd.addRCmdOption("-initw", &init_mode_weight, 0.001,
            "the initial mode weight for a newly added background mode in the detection processing");
    cmd.addRCmdOption("-fdur", &frame_duration, 1.0/25.0,
            "the duration between two used frames");
    cmd.addBCmdOption("-nolearn", &nolearn, false,
                      "disable the learning in background subtraction completely [false]");
    cmd.addBCmdOption("-nobilateral", &nobilateral, false,
                      "disable the bilateral filtering on the distance map (speed up the process) [false]");

    cmd.addText("\n *Detector model options:");

    int     n_cascade_levels = -1;
    real    Xn_prob_alpha = 1.0;
    cmd.addICmdOption("-clevel", &n_cascade_levels, -1,
            "the total number of used cascade levels for detection");
    cmd.addRCmdOption("-xnpa", &Xn_prob_alpha, 1.0,
            "the balancing alpha using negative rejection threshold");

    cmd.addText("\n *Testing windows scanning parameters:");

    int     min_det_win_width = 60;
    real    w_scale_factor = 0.1;
    int     n_win_width_scales = 25;
    cmd.addICmdOption("-wmin", &min_det_win_width, 60,
                        "minimal test window width (in pixels)");
    cmd.addRCmdOption("-wsf", &w_scale_factor, 0.1,
                        "scaling factor in width");
    cmd.addICmdOption("-wns", &n_win_width_scales, 25,
                        "number of scales in width");

    real    min_hw_ratio = 2.0;
    int     n_hw_scales = 1;
    real    h_scale_factor = 0.1;
    cmd.addRCmdOption("-hwmin", &min_hw_ratio, 2.0,
                        "minimal height/width ratio");
    cmd.addRCmdOption("-hwsf", &h_scale_factor, 0.1,
                        "height/width ratio scaling factor");
    cmd.addICmdOption("-hwns", &n_hw_scales, 1,
                        "number of scales in height-width ratio");

    real    w_shift_factor = 0.1;
    real    h_shift_factor = 0.1;
    cmd.addRCmdOption("-wshift", &w_shift_factor, 0.1,
                        "width shift factor of test windows");
    cmd.addRCmdOption("-hshift", &h_shift_factor, 0.1,
                        "height shift factor of test windows");

    cmd.addText("\n *Testing windows pre-filtering options:");
    char*   input_mask_fn = NULL;
    real    fg_percent_threshold = 0;
    char*   cam_paras_fn = NULL;
    char*   updated_cam_paras_fn = NULL;
    bool    warping_vp_flag = false;
    real    min_human_height = 150.0;
    real    max_human_height = 230.0;
    real    vertical_bot_factor = 0.885;
    real    horizon_left_factor = 0.25;
	
	
	
    cmd.addSCmdOption("-mask", &input_mask_fn, "",
                        "the input mask image");
    cmd.addRCmdOption("-fpt", &fg_percent_threshold, 0,
                        "the foreground percent threshold used to remove sub-windows");
    cmd.addSCmdOption("-cam", &cam_paras_fn, "",
                        "the camera parameters file (3x4 projection matrix + distortion coefficients) (YAML format)");
    cmd.addSCmdOption("-ucam", &updated_cam_paras_fn, "",
                        "the updated camera parameters file (YAML format)");
    cmd.addBCmdOption("-vp", &warping_vp_flag, false,
                        "warping the vertical vanishing point to infinity [false]");
    cmd.addRCmdOption("-hhmin", &min_human_height, 150.0,
                        "the minimal human height");
    cmd.addRCmdOption("-hhmax", &max_human_height, 230.0,
                        "the maximal human height");
    cmd.addRCmdOption("-vbf", &vertical_bot_factor, 0.885,
                        "the relative position of the bottom point w.r.t. test window height");
    cmd.addRCmdOption("-hlf", &horizon_left_factor, 0.25,
                        "the relative position of the left point w.r.t. test window width");
	
	
	
	//MODIF1 ALEXANDRE
	
	int     L1;
	int     L2;
	int     h1;
	int     h2;
	real    tolerance;
	
	cmd.addICmdOption("-L1", &L1, -1,
					  "the 1st reference line in the image (y-coordinate of the bottom of a reference test window, given in pixels from the top)");
	cmd.addICmdOption("-h1", &h1, -1,
					  "the test window reference height (in pixels) at line L1");
	cmd.addICmdOption("-L2", &L2, -1,
					  "the 2nd reference line in the image (y-coordinate of the bottom of a reference test window, given in pixels from the top)");
	cmd.addICmdOption("-h2", &h2, -1,
					  "the test window reference height (in pixels) at line L2");
	
	cmd.addRCmdOption("-tol", &tolerance, -1,
					  "the fraction of reference window height allowed to be kept under and above the reference size");

	
	
	//END MODIF1 ALEXANDRE
	
	
	
	
	

    cmd.addText("\n *Post-processing options:");
    real    kernel_size_percent = 0.1;
    real    kernel_sigma = 4.0;
    cmd.addRCmdOption("-ksp", &kernel_size_percent, 0.1,
                        "the kernel size percent w.r.t. the detection window");
    cmd.addRCmdOption("-ksig", &kernel_sigma, 4.0,
                        "the kernel Gaussian sigma used to smooth the detection score map for finding the local maxima");

    cmd.addText("");
    real    overlapping_threshold = 0.8;
    real    front_overlapping_threshold = 0.5;
    real    front_score_large_scale = 2.0;
    cmd.addRCmdOption("-ovt", &overlapping_threshold, 0.8,
                        "the overlapping F-measure threshold used to remove too closed merged outputs");
    cmd.addRCmdOption("-fovt", &front_overlapping_threshold, 0.5,
                        "the overlapping F-measure threshold used to remove closed front merged outputs");
    cmd.addRCmdOption("-fsls", &front_score_large_scale, 2.0,
                        "the detection score scale factor threshold used to remove closed front merged outputs");

    cmd.addText("\n *Outputs:");
    char*   output_dir = NULL;
    bool    output_image_flag = false;
    bool    output_score_image_flag = false;
    bool    output_det_res_flag = false;
    bool    output_merged_det_res_flag = false;
    cmd.addSCmdOption("-od", &output_dir, "",
                        "the output directory for the detection results");
    cmd.addBCmdOption("-oimg", &output_image_flag, false,
                        "export the image on which the detection results display [false]");
    cmd.addBCmdOption("-oscore", &output_score_image_flag, false,
                        "export the smoothed detection score image [false]");
    cmd.addBCmdOption("-odet", &output_det_res_flag, false,
                        "export the detection results (XML or YAML format) [false]");
    cmd.addBCmdOption("-omdet", &output_merged_det_res_flag, false,
                        "export the merged detection outputs (XML or YAML format) [false]");

    cmd.addText("\n *Display results:");
    bool    show_det_results = false;
    bool    show_bg_sub_results = false;
    cmd.addBCmdOption("-ddet", &show_det_results, false,
                        "display the human detection results [false]");
    cmd.addBCmdOption("-dbgs", &show_bg_sub_results, false,
                        "display the background subtraction results (foreground probability map) [false]");

    cmd.addText(" \n");
    cmd.addText("WARNING: if only a video is given, a background model MUST be given.");
    cmd.addText("");
    cmd.addText("");

    /* Read the command line */
    cmd.read(argc, argv);

    //////////////////////////////////////////////////////////////////////////////////

    /* test */
    /*
    logit_boost_fn = "/home/holoyaja/Data/IDIAPPerson/TRAIN/models/color_u_234.yml";
    images_list_or_video_fn = "/home/holoyaja/Sources/tracker/projects/codeblocks/boost_classisifer/big00.olist";
    fp_images_list_fn = "/home/holoyaja/Sources/tracker/projects/codeblocks/boost_classisifer/big00.flist";
    cam_paras_fn = "/home/holoyaja/CARETAKER/Torino_0/Bigliettatrice_00_updated_cam_paras.yml";
    Xn_prob_alpha = 0.75;
    n_cascade_levels = 30;
    min_det_win_width = 30;
    n_win_width_scales = 15;
    w_scale_factor = 0.1;
    min_human_height = 160;
    max_human_height = 220;
    shift_factor = 0.1;
    //kernel_sigma = 5.0;
    //noise_sigma = 0.5;
    left_marg_space = 0;
    right_marg_space = 0;
    bottom_marg_space = 0;
    h_scale_factor = 0.2;
    n_hw_scales = 2;
    kernel_sigma = 5;
    kernel_size_percent = 0.1;
    overlapping_threshold = 0.8;
    front_overlapping_threshold = 0.6;
    fg_percent_threshold = 0.2;
    */


    CHECK_STRING(fp_images_list_fn);
    CHECK_STRING(cam_paras_fn);
    CHECK_STRING(updated_cam_paras_fn);
    CHECK_STRING(output_dir);
    CHECK_STRING(bg_model_fn);
    CHECK_STRING(input_mask_fn);

    //////////////////////////////////////////////////////////////////////////////////

    /* load logit-boost model */
	CBoostClassifierCovariance* BCC = new CBoostClassifierCovariance();
	BCC->Load(logit_boost_fn);

	CTimer timer;

    /* define the class instant of test windows */
    CTestWindows* CAM_TW = new CTestWindows(min_det_win_width,
            n_win_width_scales,
            w_scale_factor,
            min_hw_ratio,
            n_hw_scales,
            h_scale_factor,
            w_shift_factor,
            h_shift_factor);
    CTestWindows* DET_TW = new CTestWindows(min_det_win_width,
            n_win_width_scales,
            w_scale_factor,
            min_hw_ratio,
            n_hw_scales,
            h_scale_factor,
            w_shift_factor,
            h_shift_factor);

    /* define the class instants for CImageDetection  */
	CImageDetection* DET = new CImageDetection();

    /* check the input data: image list or avi video */
    CvCapture*  VIDEO = NULL;     /* the opencv video capture handle */
    VIDEO = cvCaptureFromFile( images_list_or_video_fn );
    CFileList* IMG_LIST = NULL;
    if ( VIDEO == NULL ) {
        IMG_LIST = new CFileList(images_list_or_video_fn);
    }
    else {
#ifdef _win32
        cvSetCaptureProperty(VIDEO, CV_CAP_PROP_POS_FRAMES, detection_start_frame_idx);
#else
        for ( int i = 0 ; i < detection_start_frame_idx ; i++ )
                cvGrabFrame(VIDEO);
#endif
    }

    CFileList* FP_IMG_LIST = new CFileList(fp_images_list_fn);

	/* load the images */
	IplImage* org_img = VIDEO ? cvQueryFrame(VIDEO) : cvLoadImage(IMG_LIST->GetFileName(0));

    Camera* CAM = NULL;
    CvSize img_size;
    if ( cam_paras_fn ) {
        CAM = new Camera(org_img);
        CAM->loadCalibration(cam_paras_fn);
        CAM->computeTransforms(width_scale, height_scale, warping_vp_flag);
        img_size = cvGetSize(CAM->undistorted_image);
        if ( updated_cam_paras_fn )
            CAM->SaveUpdatedCalibration(updated_cam_paras_fn);
    }
    else
        img_size = cvSize( cvCeil((double)org_img->width*width_scale),
            cvCeil((double)org_img->height*height_scale) );

    CvSize test_img_size = cvSize(img_size.width + left_marg_space + right_marg_space,
            img_size.height + top_marg_space + bottom_marg_space);

    /* define the class instants for CMerginDetection */
    CMergingDetection* MERG = new CMergingDetection(img_size,
            kernel_size_percent, kernel_sigma, overlapping_threshold,
            front_overlapping_threshold, front_score_large_scale);

	/* load the camera parameters */
	CvMat* proj_mat = NULL;
	if ( CAM ) {
        proj_mat = cvCreateMat(3, 4, CV_64F);
        cvCopy(CAM->projection_matrix, proj_mat);
        CvMat T;
        double t[] = {1.0, 0.0, -left_marg_space ,
                    0.0, 1.0, -top_marg_space,
                    0.0, 0.0, 1.0};
        cvInitMatHeader(&T, 3, 3, CV_64FC1, t);
        cvMatMul(&T, proj_mat, proj_mat);
	}




    /* get the testing windows */
    CAM_TW->GetTestWindows(cvRect(1,1,test_img_size.width-2,test_img_size.height-2));

    /* pre-filtering the testing windows using limited heights */
    CAM_TW->FilterTestWindows(proj_mat, vertical_bot_factor, min_human_height, max_human_height);

	
	
	
	
	//MODIF2 ALEXANDRE
    
    //prototype : FilterTestWindows(int L1, int L2, int h1, int h2, double tolerance)
	
	if (L1!=-1 && L2!=-1 && h1!=-1 && h2!=-1 && tolerance!=-1) {
		
		CAM_TW->FilterTestWindows(L1, L2, h1, h2, tolerance);
	}	
    else {
		
		if (L1!=-1 || L2!=-1 || h1!=-1 || h2!=-1 || tolerance!=-1) { //at least one of these arguments is initialized
			printf("You need to provide 5 arguments (4 integers: L1,L2,h1,h2  and  1 real number: tolerance (-tol))\n");
			exit(EXIT_FAILURE);
		}
	
	}
		
    //END MODIF2 ALEXANDRE	
	
	


	/* allocate memories */
    IplImage* det_image = cvCreateImage(img_size, org_img->depth, org_img->nChannels);
    IplImage* det_fp_image = NULL;
    IplImage* det_fg_mask_image = NULL;

    /* load background model and set detection parameters */
    IplImage* fg_prob_image = NULL;
    IplImage* fg_mask_image = NULL;
    CMultiLayerBGS* BGS = NULL;
    if ( bg_model_fn ) {
        fg_prob_image = cvCreateImage(img_size, det_image->depth, 1);
        fg_mask_image = cvCreateImage(img_size, det_image->depth, 1);

        BGS = new CMultiLayerBGS();
        BGS->Init(img_size.width, img_size.height);
        BGS->SetForegroundMaskImage(fg_mask_image);
        BGS->SetForegroundProbImage(fg_prob_image);
        BGS->Load(bg_model_fn);
        BGS->SetFrameDuration(frame_duration);
        BGS->SetBackgroundModeLearningRate(mode_learn_rate_per_second);
        BGS->SetModeWeightLearningRate(weight_learn_rate_per_second);
        BGS->SetLowInitialModWeight(init_mode_weight);
        
        BGS->m_disableLearning = nolearn;
        BGS->m_disableBilateral = nobilateral;
    }

    if ( fg_percent_threshold > 0 )
            det_fg_mask_image = cvCreateImage(test_img_size, det_image->depth, 1);

    if ( BGS || (FP_IMG_LIST && FP_IMG_LIST->GetListLength() > 0) )
        det_fp_image = cvCreateImage(img_size, IPL_DEPTH_8U, 1);

    if ( fg_percent_threshold > 0 && det_fp_image && !fg_mask_image )
        fg_mask_image = cvCreateImage(img_size, det_image->depth, 1);

    if ( input_mask_fn ) {
            IplImage* input_mask_img = cvLoadImage(input_mask_fn, 0);
            CAM_TW->FilterTestWindows(input_mask_img, 0.7, (1.0-vertical_bot_factor));
            BGS->SetBkMaskImage(input_mask_img);
            cvReleaseImage(&input_mask_img);
        }

    /* human detection on input stream */
    const char* disp_det_win_name = "IDIAP Human Detection";
    const char* disp_bgs_win_name = "Foreground Probability Map";
	char lost_time[1024];
    if ( detection_end_frame_idx <= 0 ) {
        detection_end_frame_idx = 999999;
    }
    if ( show_det_results )
        cvNamedWindow(disp_det_win_name);
    if ( show_bg_sub_results && det_fp_image )
        cvNamedWindow(disp_bgs_win_name);

    for ( int frame = detection_start_frame_idx ; frame < detection_end_frame_idx ; frame += skipped_frames_num+1 ) {
        /* load images */

            timer.Start();
	    if ( !get_next_image(VIDEO, IMG_LIST, FP_IMG_LIST, CAM, BGS,
                det_image, det_fp_image, noise_sigma, frame, skipped_frames_num) )
            break;

		if ( IMG_LIST && IMG_LIST->GetListLength() > 0 )
            printf("\nDetecting people on the %d-th image : %s\n", frame, IMG_LIST->GetFileName(frame) );
        else
            printf("\nDetecting people on the %d-th video frame image \n", frame );

                timer.Stop();
                timer.PrintElapsedTimeMsg(lost_time);
                printf("Bg subtraction time : %s \n", lost_time);
                
                /* record the start time of human detection */
		timer.Start();

        /* get the expanded images */
        if ( left_marg_space > 0 || right_marg_space > 0 || top_marg_space > 0 || bottom_marg_space > 0 ) {
            IplImage* ex_image = NULL;
            DET->InsertImageMarginSpace(det_image, &ex_image,
                    left_marg_space, right_marg_space, top_marg_space, bottom_marg_space);
            cvReleaseImage(&det_image);
            det_image = ex_image;
            IplImage* ex_fp_image = NULL;
            if ( det_fp_image ) {
                DET->InsertImageMarginSpace(det_fp_image, &ex_fp_image,
                        left_marg_space, right_marg_space, top_marg_space, bottom_marg_space);
                cvReleaseImage(&det_fp_image);
                det_fp_image = ex_fp_image;
            }
        }


	// MODIFICATION JMO: DON'T RECOPY THE TEST WINDOWS IF THERE ARE NO FILTERING
	// OF THESE TEST WINDOWS

	{
	  int  filter_test_windows=0;
	  int  NtestWindows;

	  if ( fg_percent_threshold > 0 && det_fp_image )
	    filter_test_windows=1;

	  if(filter_test_windows){
	    DET_TW->SetTestWindows(CAM_TW->m_pTestWindows, CAM_TW->m_nTestWindowsNum);
	    /* pre-filtering the test windows using foreground information */
	    if ( fg_percent_threshold > 0 && det_fp_image ) {
	      cvThreshold(det_fp_image, det_fg_mask_image, 51, 255, CV_THRESH_BINARY);
	      DET_TW->FilterTestWindows(det_fg_mask_image, fg_percent_threshold, (1.0-vertical_bot_factor));
	    }
	    /* human detection */
	    DET->Detection(BCC, DET_TW, det_image, det_fp_image, -1, n_cascade_levels, Xn_prob_alpha);
	    
	    NtestWindows=DET_TW->m_nTestWindowsNum;
	  }
	  else {
	    // direclty provide the CAM_TW windows
	    /* human detection */
	    DET->Detection(BCC, CAM_TW, det_image, det_fp_image, -1, n_cascade_levels, Xn_prob_alpha);
	  
	    NtestWindows=CAM_TW->m_nTestWindowsNum;
	  }
	  
	  
	  /* merging detection candidates */
	  MERG->Merging(DET->GetPosDetWindows(), DET->GetPosDetScores(), DET->GetPosDetWindowNum());

	  timer.Stop();
	  timer.PrintElapsedTimeMsg(lost_time);
	  printf("Detection time : %s - %d test windows - %d detection candidates - %d merged detection outputs \n",
		 lost_time,NtestWindows, DET->GetPosDetWindowNum(), MERG->GetMergedDetWindowNum());
	}


        DET->UpdateDetections(0, 0, horizon_left_factor, 1.0-vertical_bot_factor);
        MERG->UpdateDetections(0, 0, horizon_left_factor, 1.0-vertical_bot_factor);

        /* draw the results */
        DET->Draw(det_image);
        MERG->Draw(det_image);

        /* display the detection results */
        if ( show_det_results ) {
            cvShowImage(disp_det_win_name, det_image);
            cvWaitKey(10);
        }
        if ( show_bg_sub_results && det_fp_image ) {
            cvShowImage(disp_bgs_win_name, det_fp_image);
            cvWaitKey(10);
        }

        /* save the detection results */
        if ( output_dir ) {
            if ( output_image_flag ) {
		        char* file_name = get_file_name(output_dir, "det_image", "jpg", frame);
                cvSaveImage(file_name, det_image);
                delete [] file_name;
            }
            if ( output_score_image_flag ) {
		        char* file_name = get_file_name(output_dir, "det_score_image", "jpg", frame);
                IplImage* score_image = MERG->GetSmoothedScoreImage();
                cvSaveImage(file_name, score_image);
                cvReleaseImage(&score_image);
                delete [] file_name;
            }
            if ( output_det_res_flag ) {
		        char* file_name = get_file_name(output_dir, "det_res", "yml", frame);
                DET->UpdateDetections(left_marg_space, top_marg_space, 0, 0);
                DET->Export(file_name);
                delete [] file_name;
            }
            if ( output_merged_det_res_flag ) {
		        char* file_name = get_file_name(output_dir, "det_merged_res", "yml", frame);
                MERG->UpdateDetections(left_marg_space, top_marg_space, 0, 0);
                MERG->Export(file_name);
                delete [] file_name;
            }
        }
    }

    /* release memories */
    delete CAM_TW;
    delete DET_TW;
    delete BCC;
    cvReleaseImage(&det_image);
    if ( det_fp_image )
        cvReleaseImage(&det_fp_image);
    if ( det_fg_mask_image )
        cvReleaseImage(&det_fg_mask_image);
    delete DET;
    delete MERG;
    if ( BGS )
        delete BGS;
    if ( IMG_LIST )
        delete IMG_LIST;
    if ( FP_IMG_LIST )
        delete FP_IMG_LIST;
    if ( VIDEO )
        cvReleaseCapture(&VIDEO);
    if ( CAM )
        delete CAM;
    if ( fg_prob_image )
        cvReleaseImage(&fg_prob_image);
    if ( fg_mask_image )
        cvReleaseImage(&fg_mask_image);

	return 1;
}

/* get next image frame data from images name list or a video file */
bool get_next_image(CvCapture* video, CFileList* img_list, CFileList* fp_img_list,
        Camera* cam, CMultiLayerBGS* bgs,
        IplImage* image, IplImage* fp_image, double noise_sigma,
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
        if ( cur_frame_idx >= img_list->GetListLength() )
            return false;
        org_image = cvLoadImage(img_list->GetFileName(cur_frame_idx));
    }

    if ( cam ) {
        cam->undistort(org_image, image);
    }
    else
        cvResize(org_image, image);

    if ( fp_img_list && fp_img_list->GetListLength() > 0 ) {
        IplImage* org_fp_image = cvLoadImage(fp_img_list->GetFileName(cur_frame_idx), 0);
        if ( cam ) {
            cam->undistort(org_fp_image, fp_image);
        }
        else
            cvResize(org_fp_image, fp_image);
        cvReleaseImage(&org_fp_image);
    }
    else if ( bgs && fp_image ) {
		/* set the new image data */
		bgs->SetRGBInputImage(image);

		/* do background subtraction process */
		bgs->Process();

		cvCopy(bgs->m_pFgProbImg, fp_image);
    }

    if ( noise_sigma > 0 ) {
        cvSmooth(image, image, CV_GAUSSIAN, 0, 0, noise_sigma);
        if ( fp_image )
            cvSmooth(fp_image, fp_image, CV_GAUSSIAN, 0, 0, noise_sigma);
    }

    if ( !video )
        cvReleaseImage(&org_image);
    return true;
}
