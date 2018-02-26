cmake_minimum_required(VERSION 2.8.3)
project(cartpend_sac)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(main main.cpp)
target_link_libraries(main ${ARMADILLO_LIBRARIES})
