cmake_minimum_required(VERSION 3.13)
project(NOvATB_TOF)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR})

add_executable(NOvATB_TOF main.cpp HitFinderAlg.cpp Run.cpp TreeAccess.cpp TimeMatchingAlg.cpp TimeMatchingAlg.h)
target_link_libraries(NOvATB_TOF ${ROOT_LIBRARIES})