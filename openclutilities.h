#ifndef OPENCLUTILITIES_H
#define OPENCLUTILITIES_H


#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS 0
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#if defined(__APPLE__) || defined(__MACOSX)
    #include "OpenCL/cl2.hpp"
#else
    #include <inc/CL/cl2.hpp>
#endif
#pragma GCC diagnostic pop

#include <string>
#include <iostream>
#include <fstream>
#include <set>

enum cl_vendor
{
    VENDOR_ANY,
    VENDOR_NVIDIA,
    VENDOR_AMD,
    VENDOR_INTEL
};

typedef struct OpenCL
{
    cl::Context context;
    cl::CommandQueue queue;
    cl::Program program;
    cl::Device device;
    cl::Platform platform;
} OpenCL;

cl::Context createCLContextFromArguments(int argc, char ** argv);

cl::Context createCLContext(cl_device_type type = CL_DEVICE_TYPE_ALL, cl_vendor vendor = VENDOR_ANY);

cl::Platform getPlatform(cl_device_type = CL_DEVICE_TYPE_ALL, cl_vendor vendor = VENDOR_ANY);

cl::Program buildProgramFromSource(cl::Context context, std::string filename, std::string buildOptions = "");

std::string getCLErrorString(cl_int err);


#endif // OPENCLUTILITIES_H
