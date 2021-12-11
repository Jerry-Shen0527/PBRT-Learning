# CMakeList.txt: PBRT 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

project ("PBRT-Learning")

include_directories(include)

# set the directory of executable files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PBRT-Learning_SOURCE_DIR}/bin)

# find required opencv
find_package(OpenCV REQUIRED)
# directory of opencv headers
include_directories(${OpenCV_INCLUDE_DIRS})

FILE(GLOB SC_FILES "*.cpp" "include/*.cpp" "include/*.h")

# name of executable file and path of source file
add_executable(${PROJECT_NAME} ${SC_FILES})

# directory of opencv library
link_directories(${OpenCV_LIBRARY_DIRS})
# opencv libraries

find_package(Assimp REQUIRED)

message("Assimp lib has ${ASSIMP_LIBS}")
include_directories(${ASSIMP_INCLUDE_DIRS})


target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS} assimp::assimp)
