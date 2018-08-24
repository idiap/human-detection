# Getting started with Idiap's Human Detection code

## Dependencies and compilation

The human detection code requires at least the OpenCV 2.3 library. CMake is also required
to build the human detector.

### Installing dependencies

If you have OpenCV and CMake installed, then everything should be ok. In other cases,
OpenCV and CMake can be installed in different ways depending on your operating system
and you environment. You can refer to the [OpenCV wiki](http://opencv.willowgarage.com/wiki/)
and the [CMake webpage](http://www.cmake.org/) for download and installation instructions.

If you have installed OpenCV manually, please make sure the `OpenCV_DIR` environment
variable points to the opencv installation directory (contains subdirectories bin,
include, ...) so CMake will be able to find it.

If your are using Ubuntu and CMake fails to find OpenCV, then you can install better
OpenCV packages. Please refer to
[the dedicated page on the OpenCV wiki](http://opencv.willowgarage.com/wiki/Ubuntu_Packages)
for details.

### Compiling the background subtraction and the human detector

Once OpenCV is properly installed, everything can be compiled using CMake:

```
mkdir build
cd build
cmake ../src -DCMAKE_BUILD_TYPE=Release
make
cd ..
ls build/bin
```

You should obtain 3 executables: `bgsub_detect`, `bgsub_learn` and `human_detect`.


## Running the programs

### Required video data

To test the program, you need a video. This video must be recorded from a static camera
and should contain persons (if you want to do person detections). For testing purpose,
we provide a video (OneStopNoEnter1cor.mpg) that can also be downloaded from the
[index of the CAVIAR dataset](http://groups.inf.ed.ac.uk/vision/CAVIAR/CAVIARDATA1/).

### Learning background model

The first step is to learn a background model. This model should be learnt using a video
sequence containing as few static objects as possible. You can generate a
`models/bgmodel.yml` output model like this: 

*(you can add the `-os` option to look at the learning process)*

```
build/bin/bgsub_learn data/OneStopNoEnter1cor.mpg models/bgmodel.yml -sfn 4
```

### Doing background subtraction (no human detection)

You will use the generated background model to do background detection. Here, we will work
on the same video sequence. The process will generate foreground mask probability images
in the `results/demo` folder. 

*(you can add the `-os` option to look at the learning process)*

```
mkdir results/demo
build/bin/bgsub_detect data/OneStopNoEnter1cor.mpg models/bgmodel.yml   \
                       -nolearn -od results/demo -ofpi
```

### Doing human detection

Here again, we need a background model and an input video. We also need a configuration
for the detector, we use `models/human.yml`. You can run human detection with interactive
visualization like this:

```
build/bin/human_detect models/human.yml data/OneStopNoEnter1cor.mpg     \
                       -bgm models/bgmodel.yml -ddet -sfn 4
```

## Understanding the parameters

### General parameters

You can get the list of parameters of any executable by just running it with no
parameters or with the `--help` option. Each parameter is accompanied by a description.
Some categories of parameters can be found for most of the executable and are good to
know, these are:

  * **Image sequence**: you can control at which frame to start and end the processing.
  It is also possible to control how many frame should be additionnaly skipped at each
  processing step.
  * **Preprocessing**: you can apply a gaussian filter and/or resize the input images if
  you wish. Note that the gaussian smooting is applied after the possible image resizing.
  * **Display output**: you can enable the display of the results on your screen.
  Depending on the algorithm, different display options are available.

Next sections given more informations about parameters of specific executables. These are
not intended to supersede the executable built-in help but rather it should give another
point of view on these parameters.

### Parameters for background subtraction (learning and detection)

The generic background subtraction algorithm can generate, refine or simply use a
background model. The `bgsub_learn` executable will only learn a model from scratch.
The `bgsub_detect` executable uses a learnt model to segment foreground/background in
images. During detection, the background model is refined by default, you can disable it
(and improve speed) with the `-nolearn` option.

### Parameters for human detection

You can use the `--help` option to get the list of all parameters. For further
explanation, please refer to the corresponding papers.
