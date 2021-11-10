# CMakeList.txt: PBRT 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

project ("PBRT-Learning")

FILE(GLOB SC_FILES "*.cpp" "*.h")

# 将源代码添加到此项目的可执行文件。
add_executable(${PROJECT_NAME} ${SC_FILES})