include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SRCS "*.cpp")

add_library(human_detection ${SRCS})
target_link_libraries(human_detection ${OpenCV_LIBS})
