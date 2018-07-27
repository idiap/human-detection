project(HUMAN_DETECTOR)
cmake_minimum_required(VERSION 2.6)

# set library and executable locations
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/libs)

## OpenCV library
# defines where to find the built opencv library 
# set(OpenCV_DIR "opencv_build")
# this also sets the include directories
find_package(OpenCV REQUIRED)
message(STATUS "Found OpenCV library version ${OpenCV_VERSION}")

# set include directories 
include_directories(cmdline common blob)
include_directories(cross_bf)
include_directories(bg_sub)
#include_directories(sp_tracker)
include_directories(human_det)

# add each library
add_subdirectory(common)
add_subdirectory(blob)
add_subdirectory(cmdline)
add_subdirectory(cross_bf)
add_subdirectory(bg_sub)
#add_subdirectory(sp_tracker)
add_subdirectory(human_det)

#add_subdirectory(mp_tracker)


# add examples
add_subdirectory(main)
