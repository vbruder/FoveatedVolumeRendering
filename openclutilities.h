/**
 * This code is based on Erik Smithad's OpenCLUtilityLibrary
 * (https://github.com/smistad/OpenCLUtilityLibrary) with the following license:
 *
 * Copyright (c) SINTEF, 2014
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of SINTEF.
 *
 */

#pragma once

#define CL_HPP_ENABLE_EXCEPTIONS    // cl2.h
#define __CL_ENABLE_EXCEPTIONS      // cl.h
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

cl::Program buildProgramFromSource(cl::Context context, const std::string &filename, 
                                   const std::string &buildOptions = "");

std::string getCLErrorString(cl_int err);
