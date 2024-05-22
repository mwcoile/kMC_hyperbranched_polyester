cmake_minimum_required(VERSION 3.6)
set(CMAKE_BUILD_TYPE Debug)
project(hyperbranchedpolyesters)

set (CMAKE_CXX_STANDARD 11)
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.75.0) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    add_executable(hyperbranchedpolyesters homo.cpp molecularWeight.cpp chainUpdate.cpp tests.cpp analysis.cpp inputs.cpp)
    target_link_libraries(hyperbranchedpolyesters ${Boost_LIBRARIES})
endif()

