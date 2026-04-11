cmake_minimum_required(VERSION 3.16)
project(ZoitendijkMethod)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Eigen ---
set(EIGEN_INCLUDE_DIR "/usr/include/eigen3")

# --- CLP ---
find_package(PkgConfig REQUIRED)
pkg_check_modules(CLP REQUIRED clp)

#-----
find_package(autodiff REQUIRED)

# --- Исполняемый файл ---
add_executable(ZoitendijkMethod main.cpp AvailDirs.cpp AvailDirs.h diff.h)

target_include_directories(ZoitendijkMethod PRIVATE
        ${EIGEN_INCLUDE_DIR}
        ${CLP_INCLUDE_DIRS}
)

target_link_libraries(ZoitendijkMethod PRIVATE ${CLP_LIBRARIES} autodiff::autodiff)

target_compile_options(ZoitendijkMethod PRIVATE ${CLP_CFLAGS_OTHER})
