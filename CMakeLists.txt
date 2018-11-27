#####################################
#####  OpenCL Volume Raycaster  #####
##### (c) 2018 Valentin Bruder  #####
#####################################

cmake_minimum_required(VERSION 3.9)

set(PROJECT "VolumeRaycasterCL")
project(${PROJECT} VERSION 1.0 DESCRIPTION "Volume raycaster based on OpenCL" LANGUAGES CXX)

# require c++14
set(CMAKE_CXX_STANDARD 14)
add_definitions(-DNOMINMAX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")
## TODO: set install target
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
# msvc only; mingw needs different logic
IF(MSVC)
    SET(QT_MISSING True)
    # look for user-registry pointing to qtcreator
    GET_FILENAME_COMPONENT(QT_BIN [HKEY_CURRENT_USER\\Software\\Classes\\Applications\\QtProject.QtCreator.cpp\\shell\\Open\\Command] PATH)

    # get root path so we can search for 5.3, 5.4, 5.5, etc
    STRING(REPLACE "/Tools" ";" QT_BIN "${QT_BIN}")
    LIST(GET QT_BIN 0 QT_BIN)
    #FILE(GLOB QT_VERSIONS "${QT_BIN}/5.*")
    #LIST(SORT QT_VERSIONS)

    # assume the latest version will be last alphabetically
    #LIST(REVERSE QT_VERSIONS)
    #LIST(GET QT_VERSIONS 0 QT_VERSION)

    # fix any double slashes which seem to be common
    #STRING(REPLACE "//" "/"  QT_VERSION "${QT_VERSION}")

	SET(QT_VERSION "C:/Qt/5.11.2")
	
    # infer the folder from Cmake msvc_version string
    # - qt uses (e.g.) "msvc2013"
    # - cmake uses (e.g.) "1800"
    # - see also https://cmake.org/cmake/help/v3.0/variable/MSVC_VERSION.html
    if(${MSVC_VERSION} VERSION_EQUAL "1800")
        SET(QT_MSVC 2013)
    elseif(${MSVC_VERSION} VERSION_EQUAL "1900")
        SET(QT_MSVC 2015)
    elseif(${MSVC_VERSION} VERSION_GREATER "1900")
        SET(QT_MSVC 2017)
    endif()

    # check for 64-bit os
    # may need to be removed for older compilers as it wasn't always offered
    IF(CMAKE_SYSTEM_PROCESSOR MATCHES 64)
        SET(QT_MSVC "${QT_MSVC}_64")
    ENDIF()
    SET(QT_PATH "${QT_VERSION}/msvc${QT_MSVC}")
    SET(QT_MISSING False)

	IF(NOT QT_MISSING)
		MESSAGE("-- Qt found: ${QT_PATH}")
		SET(Qt5_DIR "${QT_PATH}/lib/cmake/Qt5/")
		SET(Qt5Test_DIR "${QT_PATH}/lib/cmake/Qt5Test")
	ENDIF()
ENDIF()

# Find Qt libraries
FIND_PACKAGE(Qt5 COMPONENTS Core Gui Widgets Concurrent REQUIRED)

# Check Qt minor version
if (Qt5Core_FOUND)
    if (Qt5Core_VERSION VERSION_LESS 5.6.0)
        message(FATAL_ERROR "Minimum supported Qt5 version is 5.6.0!")
    endif()
else()
    message(SEND_ERROR "The Qt5 library could not be found!")
endif(Qt5Core_FOUND)

IF(${Qt5Gui_OPENGL_IMPLEMENTATION} MATCHES "GLES")
     MESSAGE(STATUS "Building an OpenGL ES build (${Qt5Gui_OPENGL_IMPLEMENTATION})")
     SET(STEL_GLES_LIBS Qt5::Gui_EGL Qt5::Gui_GLESv2)
ENDIF()

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

add_executable(${PROJECT} ${raycast_sources} ${raycast_headers})

# link Qt libraries
target_link_libraries(${PROJECT} PRIVATE Qt5::Core Qt5::Widgets Qt5::Gui Qt5::Concurrent)
# link OpenCL/OpenGL
target_link_libraries(${PROJECT} PRIVATE OpenCL::OpenCL)
target_link_libraries(${PROJECT} PRIVATE OpenGL::GL)
# optional: link OpenMP
if(OPENMP_FOUND)
    target_link_libraries(${PROJECT} PRIVATE OpenMP::OpenMP_CXX)
endif()

# include tobii research sdk
find_library(TOBII_LIBRARY NAME "tobii_research" PATHS "${CMAKE_CURRENT_SOURCE_DIR}/lib/")
target_link_libraries(${PROJECT} PRIVATE ${TOBII_LIBRARY})
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/lib/tobii_research.dll" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release/tobii_research.dll" COPYONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/lib/tobii_research.dll" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug/tobii_research.dll" COPYONLY)

# copy runtime files
IF(MSVC)
	# copy OpenCL kernel source file to build directory to support start from within VS (compiled @ runtime)
	configure_file(./src/kernel/volumeraycast.cl ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/../kernels/volumeraycast.cl COPYONLY)
	## copy Qt dlls and kernel in runtime directory
	# release
	configure_file("${QT_PATH}/bin/Qt5Core.dll" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release/Qt5Core.dll" COPYONLY)
	configure_file("${QT_PATH}/bin/Qt5Gui.dll" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release/Qt5Gui.dll" COPYONLY)
	configure_file("${QT_PATH}/bin/Qt5Widgets.dll" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release/Qt5Widgets.dll" COPYONLY)
	configure_file("${QT_PATH}/bin/Qt5Concurrent.dll" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release/Qt5Concurrent.dll" COPYONLY)
	configure_file(./src/kernel/volumeraycast.cl ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug/kernels/volumeraycast.cl COPYONLY)
	# debug
	configure_file("${QT_PATH}/bin/Qt5Cored.dll" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug/Qt5Cored.dll" COPYONLY)
	configure_file("${QT_PATH}/bin/Qt5Guid.dll" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug/Qt5Guid.dll" COPYONLY)
	configure_file("${QT_PATH}/bin/Qt5Widgetsd.dll" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug/Qt5Widgetsd.dll" COPYONLY)
	configure_file("${QT_PATH}/bin/Qt5Concurrentd.dll" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug/Qt5Concurrentd.dll" COPYONLY)
	configure_file(./src/kernel/volumeraycast.cl ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release/kernels/volumeraycast.cl COPYONLY)
ELSE()
	# copy OpenCL kernel source file to build directory (compiled @ runtime)
	configure_file(./src/kernel/volumeraycast.cl ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/kernels/volumeraycast.cl COPYONLY)
ENDIF()
