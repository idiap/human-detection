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
// bgsub_detect.cpp : Defines the entry point for the console application.
//


//#include "BackgroundSubtractionAPI.h"		/* API header file of background subtraction algorithm */
//#include "BackgroundSubtraction.h"		/* header file of background subtraction algorithm */
#include "MultiLayerBGS.h"		/* header file of background subtraction algorithm */

#include "CmdLine.h"
#include "Timer.h"
#include "FileList.h"
#include "Camera.h"

#include "cv.h"
#include "highgui.h"

#ifndef CHECK_STRING
#define CHECK_STRING(x)  {if(x && strlen(x)==0) {x=NULL;} }
#endif

#include "utilities.h"

#if !(defined(BGSUB_LEARN) || defined(BGSUB_DETECT))
#error "BGSUB_LEARN or BGSUB_DETECT should be defined."
#endif
#if (defined(BGSUB_LEARN) && defined(BGSUB_DETECT))
#error "Only one of BGSUB_LEARN or BGSUB_DETECT should be defined."
#endif

int main(int argc, char* argv[])
{
    Torch::CmdLine cmd;

    cmd.addText("\n******************************************************************************");
#ifdef BGSUB_LEARN
    cmd.addText("**   PROGRAM: Learning Background Model from a Set of Images or a Video     **");
#else
    cmd.addText("**   PROGRAM: Foreground Detection from a Set of Images or a Video          **");
#endif
    cmd.addText("**                                                                          **");
    cmd.addText("**             Author:   Dr. Jian Yao                                       **");
    cmd.addText("**                       IDIAP Research Institute                           **");
    cmd.addText("**                       CH-1920, Martigny, Switzerland                     **");
    cmd.addText("**             Emails:   Jian.Yao@idiap.ch                                  **");
    cmd.addText("**                       jianyao75@gmail.com                                **");
    cmd.addText("**             Date:     October 29, 2008                                   **");
    cmd.addText("**                                                                          **");
    cmd.addText("**   For some parameters, the parenthesis denotes the notation used in      **");
    cmd.addText("**   'Multi-Layer Background Subtraction Based on Color and Texture',       **");
    cmd.addText("**   CVPR-VS, Minneapolis, June 2007.                                       **");
    cmd.addText("**                                                                          **");
    cmd.addText("******************************************************************************");

    char msg[2048];

    cmd.addText("\nArguments:");
    char* image_list_or_video_fn;
    char* bg_model_fn;
    cmd.addSCmdArg("images list or video", &image_list_or_video_fn,
                   "the list of images or an video to be used for learning background model");
    cmd.addSCmdArg("background model", &bg_model_fn,
                   "the output file name of learned background model");

    cmd.addText("\nOptions:");

#ifdef BGSUB_LEARN
    char indent_space[] = "                         ";
    cmd.addText("\n *Saved model type:");
    int bg_model_save_type = 2;
    cmd.addICmdOption("-bgmt", &bg_model_save_type, 2,
                      "the saving type of the learned background model");
    sprintf(msg, "%s   0 - background model information (pixel by pixel)", indent_space);
    cmd.addText(msg);
    sprintf(msg, "%s   1 - background model parameters", indent_space);
    cmd.addText(msg);
    sprintf(msg, "%s   2 - both background information (pixel by pixel) and parameters", indent_space);
    cmd.addText(msg);
    char* bg_model_preload;
    cmd.addSCmdOption("-incr", &bg_model_preload, "",
                      "An existing background model to preload (if you need to leanrn on multiple sequences)");
#endif

    cmd.addText("\n *Frames selection:");
    int start_frame_idx = 0;
    int end_frame_idx = 1000000;
    int skipped_frames_num = 0;
    cmd.addICmdOption("-sf", &start_frame_idx, 0,
                      "the start frame index for foreground detection");
    cmd.addICmdOption("-ef", &end_frame_idx, 1000000,
                      "the end frame index for foreground detection");
    cmd.addICmdOption("-sfn", &skipped_frames_num, 0,
                      "the number of skipped frames between two used frames for detection");

    cmd.addText("\n *Pre-processing:");
    real width_scale = 1.0;
    real height_scale = 1.0;
    real noise_sigma = 0;
    cmd.addRCmdOption("-xs", &width_scale, 1.0,
                      "the image scaling factor at the x-direction (width)");
    cmd.addRCmdOption("-ys", &height_scale, 1.0,
                      "the image scaling factor at the y-direction (height)");
    cmd.addRCmdOption("-ngaus", &noise_sigma, 0,
                      "the Gaussian sigma to remove noise on original images");

    cmd.addText("\n *Learning rates:");
#ifdef BGSUB_DETECT
    bool disable_learning = false;
    cmd.addBCmdOption("-nolearn", &disable_learning, disable_learning,
                      "fully disable online learning of background [false]");
#endif
    real mode_learn_rate_per_second = 0.02;
    real weight_learn_rate_per_second = 0.02;
    real init_mode_weight = 0.001;
    real weight_updating_constant = 5.0;
    real frame_duration = 1.0 / 25.0;

#ifdef BGSUB_LEARN
    mode_learn_rate_per_second = 0.5;
    weight_learn_rate_per_second = 0.5;
    init_mode_weight = 0.05;
#endif

    cmd.addRCmdOption("-mlr", &mode_learn_rate_per_second, mode_learn_rate_per_second,
                      "the learning rate per second for the background modes in the detection processing (\\alpha and \\beta)");
    cmd.addRCmdOption("-wlr", &weight_learn_rate_per_second, weight_learn_rate_per_second,
                      "the learning rate per second for the background mode weights in the detection processing (\\alpha_w)");
    cmd.addRCmdOption("-initw", &init_mode_weight, init_mode_weight,
                      "the initial mode weight for a newly added background mode in the detection processing(w_{init})");
    cmd.addRCmdOption("-wuc", &weight_updating_constant, weight_updating_constant,
                      "the mode weight hysteresis (increase/decrease) updating constant (\\tau)");
    cmd.addRCmdOption("-fdur", &frame_duration, frame_duration,
                      "the duration between two used frames");

    cmd.addText("\n *Model parameters:");
    real texture_weight = 0.5;
    real bg_mode_percent = 0.6;
    real bg_prob_updating_threshold = 0.2;
    int max_mode_num = 5;
    real shadow_rate = 0.7;
    real highlight_rate = 1.2;
    int robust_LBP_constant = 15;
    int robust_color_offset = 15;
    real min_noised_angle = 10.0 / 180.0 * PI;
    cmd.addICmdOption("-maxm", &max_mode_num, 5,
                      "maximal mode number for each pixel (K)");
    cmd.addRCmdOption("-tw", &texture_weight, 0.5,
                      "distance weight by using LBP texture feature(\\lambda)");
    cmd.addRCmdOption("-bmp", &bg_mode_percent, 0.6,
                      "reliable background mode percent for foreground detection(T_{bw})");
    cmd.addRCmdOption("-bgut", &bg_prob_updating_threshold, 0.2,
                      "background/foreground distance threshold in the modeling process(T_{bg})");
    cmd.addRCmdOption("-sr", &shadow_rate, 0.7,
                      "robust shadow rate for color illumination changes(\\mu)");
    cmd.addRCmdOption("-hr", &highlight_rate, 1.2,
                      "robust highlight rate for color illumination changes(\\nu)");
    cmd.addICmdOption("-rc", &robust_color_offset, 15,
                      "noise offset for intensity changes(n_{c})");
    cmd.addICmdOption("-rl", &robust_LBP_constant, 15,
                      "robust color offset in each channel for LBP (n)");
    cmd.addRCmdOption("-mna", &min_noised_angle, min_noised_angle,
                      "minimal noised angle between the mode color and the observed color (\\theta_n)");

    cmd.addText("\n *Post-processing:");
    int pattern_neig_half_size = 4;
    real pattern_neig_gaus_sigma = 3.0;
    real bilater_filter_sigma_s = 3.0;
    real bilater_filter_sigma_r = 0.1;
    real bg_prob_threshold = 0.2;
    cmd.addICmdOption("-phs", &pattern_neig_half_size, 5,
                      "neighboring half size of pattern window for gaussian filtering on distance map to remove noise");
    cmd.addRCmdOption("-pgs", &pattern_neig_gaus_sigma, 3.0,
                      "gaussian sigma used to remove noise applied on distance map");
    cmd.addRCmdOption("-bfss", &bilater_filter_sigma_s, 3.0,
                      "spatial sigma for cross bilateral filter used to remove noise on distance map (\\sigma_s)");
    cmd.addRCmdOption("-bfsr", &bilater_filter_sigma_r, 0.1,
                      "normalized color radius sigma for cross bilateral filter used to remove noise on distance map (\\sigma_r)");
    cmd.addRCmdOption("-bgt", &bg_prob_threshold, 0.2,
                      "background/foreground distance threshold in the detection process (T_{bg})");

    cmd.addText("\n Warping with vertical vanishing point:");
    char* cam_paras_fn = NULL;
    char* updated_cam_paras_fn = NULL;
    cmd.addSCmdOption("-cam", &cam_paras_fn, "",
                      "the camera parameters file (3x4 projection matrix + distortion coefficients) (YAML format)");
    cmd.addSCmdOption("-ucam", &updated_cam_paras_fn, "",
                      "the updated camera parameters file (YAML format)");

    cmd.addText("\n Outputs:");

    char* output_warped_mask_fn = NULL;
    cmd.addSCmdOption("-owm", &output_warped_mask_fn, "",
                      "the output warped mask file");
#ifdef BGSUB_DETECT
    char* output_dir = NULL;
    bool export_org_img = false;
    bool export_fg_img = false;
    bool export_fg_mask_img = false;
    bool export_bg_img = false;
    bool export_fg_prob_img = false;
    bool export_merged_img = false;
    bool no_bilateral = false; // REMI: why only for detection? (intended?)
    bool changedetection = false;
    char* mask_img_fn = NULL;
    cmd.addText("");
    cmd.addSCmdOption("-od", &output_dir, "",
                      "the detection output directory");
    cmd.addBCmdOption("-oog", &export_org_img, false,
                      "export the original image [false]");
    cmd.addBCmdOption("-ofi", &export_fg_img, false,
                      "export the foreground image [false]");
    cmd.addBCmdOption("-ofmi", &export_fg_mask_img, false,
                      "export the foreground mask image [false]");
    cmd.addBCmdOption("-ofpi", &export_fg_prob_img, false,
                      "export the foreground probability image [false]");
    cmd.addBCmdOption("-obi", &export_bg_img, false,
                      "export the background image [false]");
    cmd.addBCmdOption("-omerged", &export_merged_img, false,
                      "export the merged image [false]");
    cmd.addSCmdOption("-mask", &mask_img_fn, "",
                      "the mask image file");
    cmd.addBCmdOption("-nobilateral", &no_bilateral, false,
    		          "disable the bilateral filtering on the distance map (speed up the process) [false]");
    /* end: added by cc */
    /*Added option for ChangeDetection output format (http://www.changedetection.net/)
     *
     */
    cmd.addBCmdOption("-changedetection", &changedetection, false,
					"export output frames in the format required by Change Detection [false]");
#endif

    cmd.addText("\n *Others:");
    bool display_results = false;
    cmd.addBCmdOption("-os", &display_results, false,
                      "displaying the learning results");

    cmd.addText("");

    /* Read the command line */
    cmd.read(argc, argv);

    CHECK_STRING(cam_paras_fn);
    CHECK_STRING(updated_cam_paras_fn);
    CHECK_STRING(output_warped_mask_fn);
#ifdef BGSUB_DETECT
    CHECK_STRING(output_dir);
    /* added by CC */
    CHECK_STRING(mask_img_fn);
    /* end: added by cc */
    /*alopez: if change detection is chosen, then automatically set export fg to true
     */
    if (changedetection) export_fg_mask_img=true;

#endif

    const char* disp_win_name = "LT: Original | RT: Background | LB: Distance | RB: Foreground";

    /* check the input data: image list or avi video */
    CvCapture* VIDEO = NULL; /* the opencv video capture handle */
    std::string convertedFile(image_list_or_video_fn);
    if (convertedFile.find(".txt")==std::string::npos)
    {
        VIDEO = cvCaptureFromFile(image_list_or_video_fn);
    }
    CFileList* LIST = NULL;
    if (VIDEO == NULL) {
        LIST = new CFileList(image_list_or_video_fn);
    } else {
        skip_first_frames(VIDEO, start_frame_idx);
    }

    CTimer Timer;

    /* load the images */
    IplImage* org_img = VIDEO ? cvQueryFrame(VIDEO) : cvLoadImage(LIST->GetFileName(0));

    Camera* CAM = NULL;
    CvSize img_size;
    if (cam_paras_fn) {
        CAM = new Camera(org_img);
        CAM->loadCalibration(cam_paras_fn);
        CAM->computeTransforms(width_scale, height_scale);
        img_size = cvGetSize(CAM->undistorted_image);
        if (updated_cam_paras_fn)
            CAM->SaveUpdatedCalibration(updated_cam_paras_fn);

        if (output_warped_mask_fn)
            cvSaveImage(output_warped_mask_fn, CAM->mask);
    } else
        img_size = cvSize(cvCeil((double) org_img->width * width_scale),
                          cvCeil((double) org_img->height * height_scale));

    /* allocate memories */
    IplImage* img = cvCreateImage(img_size, org_img->depth, org_img->nChannels);
    IplImage* fg_img = cvCreateImage(img_size, org_img->depth, org_img->nChannels);
    IplImage* fg_prob_img3 = cvCreateImage(img_size, org_img->depth, org_img->nChannels);
    IplImage* fg_prob_img = cvCreateImage(img_size, org_img->depth, 1);
    IplImage* bg_img = cvCreateImage(img_size, org_img->depth, org_img->nChannels);
    IplImage* fg_mask_img = cvCreateImage(img_size, org_img->depth, 1);

    IplImage* merged_img = NULL;
#ifdef BGSUB_DETECT
    if (display_results || (output_dir && export_merged_img))
#else
        if (display_results)
#endif
            merged_img = cvCreateImage(cvSize(img_size.width * 2, img_size.height * 2), org_img->depth, org_img->nChannels);

    if (!VIDEO)
        cvReleaseImage(&org_img);

    /************************************************************************/
    /* INITALIZATION                                                        */
    /************************************************************************/

    CMultiLayerBGS* BGS = new CMultiLayerBGS();

    /* BGS initialization */
    BGS->Init(img_size.width, img_size.height);

    /* set the foreground mask image pointer */
    BGS->SetForegroundMaskImage(fg_mask_img);

    /* set the foreground probability image pointer */
    BGS->SetForegroundProbImage(fg_prob_img);

#ifdef BGSUB_DETECT
    BGS->m_disableLearning = disable_learning;
    if ( mask_img_fn ) {
            IplImage* mask_img = cvLoadImage(mask_img_fn, 0);
            BGS->SetBkMaskImage(mask_img);
            cvReleaseImage(&mask_img);
    }
    BGS->m_disableBilateral = no_bilateral;
    /* end: added by CC */
    /* Load background model: parameters and modeling information */
    BGS->Load(bg_model_fn);
#else
    if (strlen(bg_model_preload)>0) BGS->Load(bg_model_preload);
#endif

    /* parameter setting */
    BGS->SetMaximalLBPModeNumber(max_mode_num);
    BGS->m_fWeightUpdatingConstant = weight_updating_constant;
    BGS->m_fTextureWeight = texture_weight;
    BGS->m_fBackgroundModelPercent = bg_mode_percent;
    BGS->m_nPatternDistSmoothNeigHalfSize = pattern_neig_half_size;
    BGS->m_fPatternDistConvGaussianSigma = pattern_neig_gaus_sigma;
    BGS->m_fPatternColorDistBgThreshold = bg_prob_threshold;
    BGS->m_fPatternColorDistBgUpdatedThreshold = bg_prob_updating_threshold;
    BGS->m_fRobustColorOffset = ((robust_color_offset > 0) ? robust_color_offset : 0);
    BGS->m_fRobustLBPOffset = ((robust_LBP_constant > 0 )? robust_LBP_constant : 0);
    BGS->m_fMinNoisedAngle = min_noised_angle;
    BGS->m_fRobustShadowRate = shadow_rate;
    BGS->m_fRobustHighlightRate = highlight_rate;
    BGS->m_fSigmaS = bilater_filter_sigma_s;
    BGS->m_fSigmaR = bilater_filter_sigma_r;


    /************************************************************************/
    /*  BACKGROUND DETECTION PROCESSING                                     */
    /************************************************************************/

    /* set frame rate using frame duration, for example, 25 frames per second, the frame duration = 1/25 */
    BGS->SetFrameRate(frame_duration);

    /* set main background learning parameters */
    BGS->SetParameters(max_mode_num, mode_learn_rate_per_second, weight_learn_rate_per_second, init_mode_weight);

    // background learning process
    if (display_results)
        cvNamedWindow(disp_win_name);

    for (int frame_idx = start_frame_idx; frame_idx < end_frame_idx; frame_idx += skipped_frames_num + 1) {
        if (!get_next_image(VIDEO, LIST, CAM, img, noise_sigma, frame_idx, skipped_frames_num))
            break;

        /* record the start time of background subtraction */
        Timer.Start();

        /* set the new image data */
        BGS->SetRGBInputImage(img);

        /* do background subtraction process */
        BGS->Process();

        /* record the stop time of background subtraction */
        Timer.Stop();
        Timer.PrintElapsedTimeMsg(msg);

        /* print the time */
        printf("Learning - Frame : %d\tTime : %s\n", frame_idx, msg);

        if (display_results) {
            /* get background image */
            BGS->GetBackgroundImage(bg_img);

            /* get foreground mask image */
            //BGS->GetForegroundMaskImage(m_pFgMaskImg);

            /* get foreground image */
            BGS->GetForegroundImage(fg_img);

            /* get foreground probability image */
            BGS->GetForegroundProbabilityImage(fg_prob_img3);

            /* merge the above 4 images into one image for display */
            BGS->MergeImages(4, img, bg_img, fg_prob_img3, fg_img, merged_img);

            /* show merged image */
            cvShowImage(disp_win_name, merged_img);

            cvWaitKey(10);
        }

#ifdef BGSUB_DETECT
        if (output_dir) {
            if (export_org_img) {
                char* file_name = get_file_name(output_dir, "org_img", "jpg", frame_idx);
                cvSaveImage(file_name, img);
                delete [] file_name;
            }
            if (export_fg_img) {
                if (!display_results)
                    BGS->GetForegroundImage(fg_img);
                char* file_name = get_file_name(output_dir, "fg_img", "jpg", frame_idx);
                cvSaveImage(file_name, fg_img);
                delete [] file_name;
            }
            if (export_fg_mask_img) {
            	 //alopez: Export bin%06d.png for Change detection
                if (!display_results)
                    BGS->GetForegroundMaskImage(fg_mask_img);
                char* file_name;
                if (changedetection)
                {
                	file_name = get_file_name(output_dir, "bin", "png", frame_idx+1, true);
                }
                else
                	file_name = get_file_name(output_dir, "fg_mask_img", "jpg", frame_idx);
                cvSaveImage(file_name, fg_mask_img);
                delete [] file_name;
            }
            if (export_fg_prob_img) {
                char* file_name = get_file_name(output_dir, "fg_prob_img", "jpg", frame_idx);
                cvSaveImage(file_name, fg_prob_img);
                delete [] file_name;
            }
            if (export_bg_img) {
                if (!display_results)
                    BGS->GetBackgroundImage(bg_img);
                char* file_name = get_file_name(output_dir, "bg_img", "jpg", frame_idx);
                cvSaveImage(file_name, bg_img);
                delete [] file_name;
            }
            if (export_merged_img) {
                if (!display_results) {
                    if (!export_bg_img)
                        BGS->GetBackgroundImage(bg_img);
                    BGS->GetForegroundProbabilityImage(fg_prob_img3);
                    if (!export_fg_img)
                        BGS->GetForegroundImage(fg_img);
                    BGS->MergeImages(4, img, bg_img, fg_prob_img3, fg_img, merged_img);
                }
                char* file_name = get_file_name(output_dir, "merged_img", "jpg", frame_idx);
                cvSaveImage(file_name, merged_img);
                delete [] file_name;
            }

        }
#endif
    }
#ifdef BGSUB_LEARN
    /* save the learned background model */
    printf("\nSaving the background model: %s\n", bg_model_fn);
    BGS->Save(bg_model_fn, bg_model_save_type);
#endif
    /* release memories */
    if (merged_img)
        cvReleaseImage(&merged_img);
    delete BGS;

    cvReleaseImage(&img);
    cvReleaseImage(&fg_img);
    cvReleaseImage(&fg_prob_img);
    cvReleaseImage(&fg_prob_img3);
    cvReleaseImage(&bg_img);

    if (LIST) delete LIST;
    if (CAM) delete CAM;

    return 0;
}
