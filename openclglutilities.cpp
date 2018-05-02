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

#include "openclglutilities.h"

cl::Context createCLGLContext(cl_device_type type, cl_vendor vendor)
{
    cl::Platform platform = getPlatform(type, vendor);

    //Creating the context
#if defined(__APPLE__) || defined(__MACOSX)
    // Apple (untested)
    cl_context_properties cps[] = {
       CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE,
       (cl_context_properties)CGLGetShareGroup(CGLGetCurrentContext()),
       0};


#else
  #ifdef _WIN32
      // Windows
      cl_context_properties cps[] = {
          CL_GL_CONTEXT_KHR,
          (cl_context_properties)wglGetCurrentContext(),
          CL_WGL_HDC_KHR,
          (cl_context_properties)wglGetCurrentDC(),
          CL_CONTEXT_PLATFORM,
          (cl_context_properties)(platform)(),
          0
      };
  #else
      // Linux
    cl_context_properties cps[] = {
        CL_GL_CONTEXT_KHR,
        (cl_context_properties)glXGetCurrentContext(),
        CL_GLX_DISPLAY_KHR,
        (cl_context_properties)glXGetCurrentDisplay(),
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)(platform)(),
        0
    };
#endif
#endif

    try
    {
        // We need to check if there is more than one device first
        std::vector<cl::Device> devices;
        std::vector<cl::Device> singleDevice;
        platform.getDevices(type, &devices);
        cl::Context context;

        // If more than one CL device, find out which one is associated with GL context
        if(devices.size() > 1)
        {
#if !(defined(__APPLE__) || defined(__MACOSX))
            cl::Device interopDevice = getValidGLCLInteropDevice(platform, cps);
            singleDevice.push_back(interopDevice);
            context = cl::Context(singleDevice, cps);
#else
            context = cl::Context(type,cps);
#endif
        }
        else
        {
            context = cl::Context(type, cps);
        }
        return context;
    }
    catch(cl::Error error)
    {
        std::cerr << "Error in " << error.what() << "(" << getCLErrorString(error.err()) << ")"
                  << std::endl;
        throw std::runtime_error("Failed to create an OpenCL context from the OpenGL context. Make \
                                  sure the displaying device is your selected compute device, \
                                  otherwise the context sharing won't work.");
    }
}

#if !(defined(__APPLE__) || defined(__MACOSX))
cl::Device getValidGLCLInteropDevice(cl::Platform platform, cl_context_properties* properties) {
    // Function for finding a valid device for CL-GL context.
    // Thanks to Jim Vaughn for this contribution
//    cl::Device displayDevice;
    cl_device_id interopDeviceId;

    int status;
    size_t deviceSize = 0;

    // TODO select desired platform (not just the first one)
    cl_platform_id platform_id;
    clGetPlatformIDs(1, &platform_id, NULL);
    // Load extension function call
    clGetGLContextInfoKHR_fn glGetGLContextInfo_func =
            (clGetGLContextInfoKHR_fn)clGetExtensionFunctionAddressForPlatform(platform_id,
                                                                               "clGetGLContextInfoKHR");

    // Ask for the CL device associated with the GL context
    status = glGetGLContextInfo_func(properties,
                                     CL_CURRENT_DEVICE_FOR_GL_CONTEXT_KHR,
                                     sizeof(cl_device_id),
                                     &interopDeviceId,
                                     &deviceSize);

    if(deviceSize == 0)
    {
        throw cl::Error(1, "No OpenGL capable devices found for current platform.");
    }
    if(status != CL_SUCCESS)
    {
        throw cl::Error(1, "Could not get CL-GL interop device for the current platform. \
                            Failure occured during call to clGetGLContextInfoKHR.");
    }

    return cl::Device(interopDeviceId);
}
#endif
