
add_executable(bgsub_learn bgsub_learn.cpp)
target_link_libraries(bgsub_learn ${OpenCV_LIBS} blob cmdline common background)

add_executable(bgsub_detect bgsub_detect.cpp)
target_link_libraries(bgsub_detect ${OpenCV_LIBS} blob cmdline common background)

add_executable(human_detect human_detect.cpp)
target_link_libraries(human_detect human_detection ${OpenCV_LIBS} blob cmdline common background)
