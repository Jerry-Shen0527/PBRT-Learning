# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.5)

# Set the project name
project (RayTracing_in1weekend)

#include head_file
include_directories(include)

FILE(GLOB_RECURSE CURRENT_INCLUDE include/*.h)
FILE(GLOB_RECURSE CURRENT_SOURCES src/*.cpp)
# Add an executable with the above sources
add_executable(RayTracing_in1weekend ${CURRENT_SOURCES} ${CURRENT_INCLUDE})

#OpenMp
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
	message("Found OpenMP")
    target_link_libraries(RayTracing_in1weekend PUBLIC OpenMP::OpenMP_CXX)
endif()

