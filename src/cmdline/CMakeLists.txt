include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SRCS "*.cc")

add_library(cmdline ${SRCS})
target_link_libraries(cmdline)
