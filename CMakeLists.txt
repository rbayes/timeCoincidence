cmake_minimum_required (VERSION 2.8.11)
set (CMAKE_CXX_STANDARD 11)

project (timeCoincidence)

find_package(ROOT)
include(${ROOT_USE_FILE})

add_library(HistContainer HistContainer.cpp)
target_link_libraries(HistContainer ${ROOT_LIBRARIES})

Include_directories ("${PROJECT_SOURCE_DIR}")
add_executable(timeCoincidence timeCoincidence.cpp)

target_link_libraries(timeCoincidence HistContainer)
target_link_libraries(timeCoincidence ${ROOT_LIRARIES})
