cmake_minimum_required(VERSION 3.16)
project(ZoitendijkMethod
        VERSION 1.0.0
        LANGUAGES CXX
)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(CLP REQUIRED clp)

find_package(autodiff REQUIRED)

add_library(availdirs STATIC
        src/AvailDirs.cpp
)

target_include_directories(availdirs PUBLIC
        include
)

target_link_libraries(availdirs PUBLIC
        Eigen3::Eigen
        ${CLP_LIBRARIES}
        autodiff::autodiff
)

target_include_directories(availdirs SYSTEM PUBLIC
        ${CLP_INCLUDE_DIRS}
)

target_compile_options(availdirs PRIVATE ${CLP_CFLAGS_OTHER})

add_executable(run_examples examples/main.cpp)
target_link_libraries(run_examples PRIVATE availdirs)
