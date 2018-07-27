include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SRCS "*.cpp")

add_library(blob ${SRCS})
target_link_libraries(blob ${OpenCV_LIBS})
