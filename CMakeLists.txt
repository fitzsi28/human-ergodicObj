cmake_minimum_required(VERSION 2.8.3)
project(cartpend_sac)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(sac_cartpend sac_cartpend.cpp)
target_link_libraries(sac_cartpend ${ARMADILLO_LIBRARIES})

#add_executable(mda_cartpend mda_cartpend.cpp)
#target_link_libraries(mda_cartpend ${ARMADILLO_LIBRARIES})

#add_executable(mig_mda_cartpend mig_mda_cartpend.cpp)
#target_link_libraries(mig_mda_cartpend ${ARMADILLO_LIBRARIES})

add_executable(ergsac_cartpend ergsac_cartpend.cpp)
target_link_libraries(ergsac_cartpend ${ARMADILLO_LIBRARIES})

add_executable(sac_doubleint sac_doubleint.cpp)
target_link_libraries(sac_doubleint ${ARMADILLO_LIBRARIES})

add_executable(ergsac_doubleint ergsac_doubleint.cpp)
target_link_libraries(ergsac_doubleint ${ARMADILLO_LIBRARIES})

