include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SRCS "*.cpp")

add_library(background ${SRCS})
target_link_libraries(background ${OpenCV_LIBS} blob)

add_library(bgsub SHARED ${SRCS})
add_dependencies(bgsub ${OpenCV_LIBS} blob)
