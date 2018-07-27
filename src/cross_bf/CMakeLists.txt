include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SRCS "*.cpp")

add_library(cross_bf ${SRCS})
target_link_libraries(cross_bf ${OpenCV_LIBS})
