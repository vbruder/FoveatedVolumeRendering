#####################################
#####  OpenCL Volume Raycaster  #####
##### (c) 2018 Valentin Bruder  #####
#####################################

cmake_minimum_required(VERSION 3.9)

project(VolumeRaycasterCL VERSION 1.0 DESCRIPTION "Volume raycaster based on OpenCL" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")
## TODO: install target
# if (NOT DEFINED VolumeRaycasterCL_INSTALL_BIN_DIR)
#   set(VolumeRaycasterCL_INSTALL_BIN_DIR "bin")
# endif()

### OpenCL
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})
find_package(OpenCL REQUIRED)
include_directories(${OPENCL_INCLUDE_DIR})

find_package(OpenGL REQUIRED)
find_package(OpenMP)

### Qt
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

# Find Qt libraries
find_package(Qt5Core CONFIG REQUIRED)
find_package(Qt5Widgets CONFIG REQUIRED)
find_package(Qt5Gui CONFIG REQUIRED)
find_package(Qt5Concurrent CONFIG REQUIRED)
# Check Qt minor version
if (Qt5Core_FOUND)
    if (Qt5Core_VERSION VERSION_LESS 5.6.0)
        message(FATAL_ERROR "Minimum supported Qt5 version is 5.6.0!")
    endif()
else()
    message(SEND_ERROR "The Qt5 library could not be found!")
endif(Qt5Core_FOUND)

# set headers
set(raycast_headers
  src/io/datrawreader.h
  src/oclutil/openclutilities.h
  src/oclutil/openclglutilities.h
  src/qt/mainwindow.h
  src/qt/transferfunctionwidget.h
  src/qt/volumerenderwidget.h
  src/qt/colorutils.h
  src/qt/colorwheel.h
  src/qt/hoverpoints.h
  src/core/volumerendercl.h
  inc/CL/cl2.hpp
  )

# set sources
set(raycast_sources
  src/io/datrawreader.cpp
  src/oclutil/openclutilities.cpp
  src/oclutil/openclglutilities.cpp
  src/qt/main.cpp
  src/qt/mainwindow.cpp
  src/qt/mainwindow.ui
  src/qt/transferfunctionwidget.cpp
  src/qt/volumerenderwidget.cpp
  src/qt/colorutils.cpp
  src/qt/colorwheel.cpp
  src/qt/hoverpoints.cpp
  src/core/volumerendercl.cpp
  )

set(EXECUTABLE VolumeRaycasterCL)
add_executable(${EXECUTABLE} ${raycast_sources} ${raycast_headers})

# link Qt libraries
target_link_libraries(VolumeRaycasterCL PRIVATE Qt5::Core)
target_link_libraries(VolumeRaycasterCL PRIVATE Qt5::Widgets)
target_link_libraries(VolumeRaycasterCL PRIVATE Qt5::Gui)
target_link_libraries(VolumeRaycasterCL PRIVATE Qt5::Concurrent)

# link OpenCL/OpenGL
target_link_libraries(VolumeRaycasterCL PRIVATE OpenCL::OpenCL)
target_link_libraries(VolumeRaycasterCL PRIVATE OpenGL::GL)
# optional: link OpenMP
if(OPENMP_FOUND)
    target_link_libraries(VolumeRaycasterCL PRIVATE OpenMP::OpenMP_CXX)
endif()

# copy OpenCL kernel source file to build directory (compiled @ runtime)
configure_file(./src/kernel/volumeraycast.cl ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kernels/volumeraycast.cl COPYONLY)
