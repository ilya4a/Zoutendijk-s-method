cmake_minimum_required(VERSION 3.16)
project(ZoitendijkMethod)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EIGEN_INCLUDE_DIR "/usr/include/eigen3")

find_package(PkgConfig REQUIRED)
pkg_check_modules(CLP REQUIRED clp)

find_package(autodiff REQUIRED)

add_executable(ZoitendijkMethod
        src/main.cpp
        src/AvailDirs.cpp
)

target_include_directories(ZoitendijkMethod PRIVATE
        include
        ${EIGEN_INCLUDE_DIR}
        ${CLP_INCLUDE_DIRS}
)

target_link_libraries(ZoitendijkMethod PRIVATE
        ${CLP_LIBRARIES}
        autodiff::autodiff
)

target_compile_options(ZoitendijkMethod PRIVATE ${CLP_CFLAGS_OTHER})