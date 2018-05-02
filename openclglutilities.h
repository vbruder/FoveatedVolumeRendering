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

#include "openclutilities.h"

#if defined (__APPLE__) || defined(MACOSX)
    #define GL_SHARING_EXTENSION "cl_APPLE_gl_sharing"
#else
    #define GL_SHARING_EXTENSION "cl_khr_gl_sharing"
#endif

#if defined(__APPLE__) || defined(__MACOSX)
    #include <OpenCL/cl_gl.h>
    #include <OpenGL/OpenGL.h>
#elif _WIN32
    #include <windows.h>
    #include <GL\GL.h>
    #include <GL\GLU.h>
    #include <CL\cl_gl.h>
    
    #pragma comment(lib, "opengl32.lib")
    #pragma comment(lib, "glu32.lib")
    #pragma warning( disable : 4996)
#else   // Linux
    #include <GL/glx.h>
    #include <CL/cl_gl.h>
#endif


cl::Context createCLGLContext(cl_device_type type = CL_DEVICE_TYPE_ALL,
                              cl_vendor vendor = VENDOR_ANY);

#if !(defined(__APPLE__) || defined(__MACOSX))
cl::Device getValidGLCLInteropDevice(cl::Platform platform, cl_context_properties* properties);
#endif

