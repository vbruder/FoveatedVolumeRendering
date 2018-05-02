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

#include "openclutilities.h"

#ifdef WIN32
#else
#include <sys/stat.h>
#include <time.h>
#endif

cl::Platform getPlatform(cl_device_type type, cl_vendor vendor)
{
    // Get available platforms
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    if(platforms.size() == 0)
        throw cl::Error(1, "No OpenCL platforms were found");

    int platformID = -1;
    if(vendor != VENDOR_ANY)
    {
        std::string find;
        switch(vendor)
        {
            case VENDOR_NVIDIA:
                find = "NVIDIA";
            break;
            case VENDOR_AMD:
                find = "Advanced Micro Devices";
            break;
            case VENDOR_INTEL:
                find = "Intel";
            break;
            default:
                throw cl::Error(1, "Invalid vendor specified");
            break;
        }

        std::cout << platforms.size() << " OpenCL platforms found: " << std::endl;
        for(unsigned int i = 1; i <= platforms.size(); i++)
            std::cout << "    " << i << ") "
                      << platforms[i-1].getInfo<CL_PLATFORM_NAME>() << std::endl;

        for(unsigned int i = 0; i < platforms.size(); i++)
        {
            if(platforms[i].getInfo<CL_PLATFORM_VENDOR>().find(find) != std::string::npos)
            {
                try {
                    std::vector<cl::Device> devices;
                    platforms[i].getDevices(type, &devices);
                    platformID = i;
                    break;
                } catch(cl::Error e) {
                   continue;
                }
            }
        }
    } else
    {
        for(unsigned int i = 0; i < platforms.size(); i++)
        {
            try
            {
                std::vector<cl::Device> devices;
                platforms[i].getDevices(type, &devices);
                platformID = i;
                break;
            } catch(cl::Error e) {
               continue;
            }
        }
    }

    if(platformID == -1)
        throw cl::Error(1, "No compatible OpenCL platform found");

    cl::Platform platform = platforms[platformID];
    std::cout << "Using platform: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
    std::cout << "From vendor: " << platform.getInfo<CL_PLATFORM_VENDOR>() << std::endl;
    return platform;
}


cl::Context createCLContextFromArguments(int argc, char ** argv)
{
    cl_device_type type = CL_DEVICE_TYPE_ALL;
    cl_vendor vendor = VENDOR_ANY;

    for(int i = 0; i < argc; i++) {
        if(strcmp(argv[i], "--device") == 0) {
            if(strcmp(argv[i+1], "cpu") == 0) {
                type = CL_DEVICE_TYPE_CPU;
            } else if(strcmp(argv[i+1], "gpu") == 0) {
                type = CL_DEVICE_TYPE_GPU;
            }
            i++;
        } else if(strcmp(argv[i], "--vendor") == 0) {
            if(strcmp(argv[i+1], "amd") == 0) {
                vendor = VENDOR_AMD;
            } else if(strcmp(argv[i+1], "intel") == 0) {
                vendor = VENDOR_INTEL;
            } else if(strcmp(argv[i+1], "nvidia") == 0) {
                vendor = VENDOR_NVIDIA;
            }
            i++;
        }
    }

    return createCLContext(type, vendor);
}


cl::Context createCLContext(cl_device_type type, cl_vendor vendor)
{
    cl::Platform platform = getPlatform(type, vendor);

    // Use the preferred platform and create a context
    cl_context_properties cps[] = {
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)(platform)(),
        0
    };

    try {
        cl::Context context = cl::Context(type, cps);
        return context;
    } catch(cl::Error error) {
        throw cl::Error(1, "Failed to create an OpenCL context!");
    }
}


cl::Program buildProgramFromSource(cl::Context context, const std::string &filename,
                                   const std::string &buildOptions)
{
        // Read source file
        std::ifstream sourceFile(filename.c_str());
        if(sourceFile.fail())
            throw std::invalid_argument("Failed to open OpenCL kernel file " + filename);
        std::string sourceCode(std::istreambuf_iterator<char>(sourceFile),
                               (std::istreambuf_iterator<char>()));
        //cl::Program::Sources source(1, std::make_pair(sourceCode.c_str(), sourceCode.length()+1));

        // Make program of the source code in the context
        cl::Program program = cl::Program(context, sourceCode);

        std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

        // Build program for these specific devices
        try {
            program.build(devices, buildOptions.c_str());
        } catch(cl::Error error)
        {
            if(error.err() == CL_BUILD_PROGRAM_FAILURE)
            {
                std::cout << "Build log:" << std::endl
                          << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << std::endl;
            }
            throw error;
        }
        return program;
}


std::string getCLErrorString(cl_int err)
{
    switch (err) {
        case CL_SUCCESS:                          return "Success!";
        case CL_DEVICE_NOT_FOUND:                 return "Device not found.";
        case CL_DEVICE_NOT_AVAILABLE:             return "Device not available";
        case CL_COMPILER_NOT_AVAILABLE:           return "Compiler not available";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:    return "Memory object allocation failure";
        case CL_OUT_OF_RESOURCES:                 return "Out of resources";
        case CL_OUT_OF_HOST_MEMORY:               return "Out of host memory";
        case CL_PROFILING_INFO_NOT_AVAILABLE:     return "Profiling information not available";
        case CL_MEM_COPY_OVERLAP:                 return "Memory copy overlap";
        case CL_IMAGE_FORMAT_MISMATCH:            return "Image format mismatch";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:       return "Image format not supported";
        case CL_BUILD_PROGRAM_FAILURE:            return "Program build failure";
        case CL_MAP_FAILURE:                      return "Map failure";
        case CL_INVALID_VALUE:                    return "Invalid value";
        case CL_INVALID_DEVICE_TYPE:              return "Invalid device type";
        case CL_INVALID_PLATFORM:                 return "Invalid platform";
        case CL_INVALID_DEVICE:                   return "Invalid device";
        case CL_INVALID_CONTEXT:                  return "Invalid context";
        case CL_INVALID_QUEUE_PROPERTIES:         return "Invalid queue properties";
        case CL_INVALID_COMMAND_QUEUE:            return "Invalid command queue";
        case CL_INVALID_HOST_PTR:                 return "Invalid host pointer";
        case CL_INVALID_MEM_OBJECT:               return "Invalid memory object";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:  return "Invalid image format descriptor";
        case CL_INVALID_IMAGE_SIZE:               return "Invalid image size";
        case CL_INVALID_SAMPLER:                  return "Invalid sampler";
        case CL_INVALID_BINARY:                   return "Invalid binary";
        case CL_INVALID_BUILD_OPTIONS:            return "Invalid build options";
        case CL_INVALID_PROGRAM:                  return "Invalid program";
        case CL_INVALID_PROGRAM_EXECUTABLE:       return "Invalid program executable";
        case CL_INVALID_KERNEL_NAME:              return "Invalid kernel name";
        case CL_INVALID_KERNEL_DEFINITION:        return "Invalid kernel definition";
        case CL_INVALID_KERNEL:                   return "Invalid kernel";
        case CL_INVALID_ARG_INDEX:                return "Invalid argument index";
        case CL_INVALID_ARG_VALUE:                return "Invalid argument value";
        case CL_INVALID_ARG_SIZE:                 return "Invalid argument size";
        case CL_INVALID_KERNEL_ARGS:              return "Invalid kernel arguments";
        case CL_INVALID_WORK_DIMENSION:           return "Invalid work dimension";
        case CL_INVALID_WORK_GROUP_SIZE:          return "Invalid work group size";
        case CL_INVALID_WORK_ITEM_SIZE:           return "Invalid work item size";
        case CL_INVALID_GLOBAL_OFFSET:            return "Invalid global offset";
        case CL_INVALID_EVENT_WAIT_LIST:          return "Invalid event wait list";
        case CL_INVALID_EVENT:                    return "Invalid event";
        case CL_INVALID_OPERATION:                return "Invalid operation";
        case CL_INVALID_GL_OBJECT:                return "Invalid OpenGL object";
        case CL_INVALID_BUFFER_SIZE:              return "Invalid buffer size";
        case CL_INVALID_MIP_LEVEL:                return "Invalid mip-map level";
        default:                                  return "Unknown";
    }
}

