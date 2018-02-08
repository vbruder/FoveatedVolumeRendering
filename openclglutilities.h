#ifndef OPENCLGLUTILITIES_H
#define OPENCLGLUTILITIES_H

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

#endif // OPENCLGLUTILITIES_H
