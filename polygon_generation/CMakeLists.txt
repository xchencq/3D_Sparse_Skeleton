cmake_minimum_required(VERSION 3.0.2)
project(polygon_generation)

set(CMAKE_VERBOSE_MAKEFILE "false")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -DEIGEN_DONT_PARALLELIZE ${CMAKE_CXX_FLAGS} -g -O3 -Wall")

find_package(catkin REQUIRED COMPONENTS
    roscpp
    std_msgs
    nav_msgs
    sensor_msgs
    quadrotor_msgs
    tf
    kdtree
    plan_env)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)

catkin_package(
  INCLUDE_DIRS include
)

include_directories(
    SYSTEM 
    include
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${PCL_INCLUDE_DIRS})

link_directories(
    ${PCL_LIBRARY_DIRS})

add_executable (skeleton_finder_3D_test
    src/skeleton_finder_3D_test.cpp
    src/skeleton_finder_3D.cpp
    utils/quickhull/QuickHull.cpp
    utils/A_star/A_star.cpp
    )

target_link_libraries (skeleton_finder_3D_test
    ${catkin_LIBRARIES}
    ${PCL_LIBRARIES})

add_executable (map 
    src/map_generator.cpp)

target_link_libraries (map
    ${catkin_LIBRARIES}
    ${PCL_LIBRARIES})

    


