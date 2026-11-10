cmake_minimum_required(VERSION 3.18)
project(examm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
find_package(OpenMP)

pybind11_add_module(_core cpp/core.cpp)
target_include_directories(_core PRIVATE cpp)
if(OpenMP_CXX_FOUND)
  target_link_libraries(_core PRIVATE OpenMP::OpenMP_CXX)
else()
  message(WARNING "OpenMP not found: examm will run single-threaded")
endif()

install(TARGETS _core LIBRARY DESTINATION examm)
