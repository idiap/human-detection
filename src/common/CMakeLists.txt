include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SRCS "*.cpp")

add_library(common ${SRCS})
target_link_libraries(common cmdline ${OpenCV_LIBS})
