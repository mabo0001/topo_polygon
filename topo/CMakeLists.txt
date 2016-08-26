cmake_minimum_required(VERSION 3.5)
project(topo)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    test_topo/test.c
    test_topo/test_topo.sln
    test_topo/test_topo.vcproj
    list.c
    list.h
    rtree_2d.c
    rtree_2d.h
    topo.c
    topo.h
    topo.lin
    unistd.h)

add_executable(topo ${SOURCE_FILES})