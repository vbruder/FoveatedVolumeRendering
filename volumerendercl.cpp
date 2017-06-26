#include "volumerendercl.h"

#include <algorithm>
#include <numeric>
#include "omp.h"

/**
 * @brief VolumeRenderCL::VolumeRenderCL
 */
VolumeRenderCL::VolumeRenderCL()
{
}


/**
 * @brief VolumeRenderCL::~VolumeRenderCL
 */
VolumeRenderCL::~VolumeRenderCL()
{
}


/**
 * @brief VolumeRenderCL::logCLerror
 * @param error
 */
void VolumeRenderCL::logCLerror(cl::Error error)
{
    std::cerr << "Error in " << error.what() << ": "
              << getCLErrorString(error.err()) << std::endl;
    // TODO: logging
    throw std::runtime_error( "ERROR: " + std::string(error.what()) + "("
                              + getCLErrorString(error.err()) + ")");
}


/**
 * @brief VolumeRenderCL::initialize
 */
void VolumeRenderCL::initialize()
{
    try // opencl scope
    {
        // TODO: replace if no NVIDIA GPU
        _contextCL = createCLGLContext(CL_DEVICE_TYPE_GPU, VENDOR_NVIDIA); // VENDOR_ANY
        _queueCL = cl::CommandQueue(_contextCL);
    }
    catch (cl::Error err)
    {
        logCLerror(err);
    }

    initKernel("../RaycastLight/kernels/volumeraycast.cl", "-DIMAGE_SUPPORT=1 -DCL_STD=CL1.2");
}


/**
 * @brief VolumeRenderCL::initKernel
 * @param fileName
 * @param buildFlags
 */
void VolumeRenderCL::initKernel(const std::string fileName, const std::string buildFlags)
{
    try
    {
        cl::Program program = buildProgramFromSource(_contextCL, fileName, buildFlags);
        _raycastKernel = cl::Kernel(program, "volumeRender");
        cl_float16 view = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
        _raycastKernel.setArg(VIEW, view);
        _raycastKernel.setArg(SAMPLING_RATE, 0.5f);     // default step size 0.5*voxel size
        _raycastKernel.setArg(ORTHO, 0);            // perspective cam by default
        _raycastKernel.setArg(ILLUMINATION, 1);     // illumination on by default
        _raycastKernel.setArg(BOX, 1);
        _raycastKernel.setArg(LINEAR, 1);
        cl_float4 bgColor = {{1.f, 1.f, 1.f, 1.f}};
        _raycastKernel.setArg(BACKGROUND, bgColor);
    }
    catch (cl::Error err)
    {
        logCLerror(err);
    }
}


/**
 * @brief VolumeRenderCL::setMemObjects
 */
void VolumeRenderCL::setMemObjects()
{
    _raycastKernel.setArg(VOLUME, _volumeMem);
    _raycastKernel.setArg(OUTPUT, _outputMem);
    _raycastKernel.setArg(TFF, _tffMem);
}


/**
 * @brief VolumeRenderCL::calcScaling
 */
void VolumeRenderCL::calcScaling()
{
    if (!_dr.has_data())
        return;

    _modelScale = std::valarray<double>(_dr.properties().volume_res.size());
    for (size_t i = 0; i < _dr.properties().volume_res.size(); ++i)
    {
        _modelScale[i] = _dr.properties().volume_res.at(i); // TODO
    }
    std::valarray<double> thickness(_dr.properties().slice_thickness.data(),
                                    _dr.properties().slice_thickness.size());
    _modelScale *= thickness*(1.0/thickness[0]);
#undef max  // error here if I don't undef max
    _modelScale = _modelScale.max() / _modelScale;
//    std::cout << "Scaling volume: (" << _scale[0] << ", " << _scale[1] << ", "
//              << _scale[2] << ")" << std::endl;
}


/**
 * @brief VolumeRenderCL::scaleVolume
 * @param scale
 */
void VolumeRenderCL::scaleVolume(std::valarray<double> scale)
{
    calcScaling();
    _modelScale *= scale;
}


/**
 * @brief VolumeRenderCL::updateKernelArgs
 * @param viewMat
 */
void VolumeRenderCL::updateView(const std::array<float, 16> viewMat)
{
    if (!_dr.has_data() || _modelScale.size() < 3)
        return;

    cl_float16 modelViewMat;
    for (size_t i = 0; i < 16; ++i)
    {
        if (i < 4)
            modelViewMat.s[i] = viewMat[i] * _modelScale[0];
        else if (i < 8)
            modelViewMat.s[i] = viewMat[i] * _modelScale[1];
        else if (i < 12)
            modelViewMat.s[i] = viewMat[i] * _modelScale[2];
        else
            modelViewMat.s[i] = viewMat[i];
    }

    try{
        _raycastKernel.setArg(VIEW, modelViewMat);
    } catch (cl::Error err) {
        logCLerror(err);
    }
}


/**
 * @brief VolumeRenderCL::updateStepSize
 * @param stepSize
 */
void VolumeRenderCL::updateSamplingRate(const double samplingRate)
{
    try{
        _raycastKernel.setArg(SAMPLING_RATE, static_cast<cl_float>(samplingRate));
    } catch (cl::Error err) {
        logCLerror(err);
    }
}


/**
 * @brief VolumeRenderCL::updateOutputImg
 * @param width
 * @param height
 */
void VolumeRenderCL::updateOutputImg(const size_t width, const size_t height, cl_GLuint texId)
{
    cl::ImageFormat format;
    format.image_channel_order = CL_RGBA;
    format.image_channel_data_type = CL_FLOAT;
    try
    {
        _outputMem = cl::ImageGL(_contextCL,
                                 CL_MEM_WRITE_ONLY,
                                 GL_TEXTURE_2D,
                                 0,
                                 texId);
        _raycastKernel.setArg(OUTPUT, _outputMem);
        _outputData.resize(width * height * 4, 0);
    }
    catch (cl::Error err)
    {
        logCLerror(err);
    }
}


/**
 * @brief VolumeRenderCL::runRaycast
 * @param imgSize
 */
void VolumeRenderCL::runRaycast(const size_t width, const size_t height)
{
    if (!this->_volLoaded)
        return;
    try // opencl scope
    {
        setMemObjects();
        cl::NDRange globalThreads(width, height);
        cl::Event ndrEvt;

        std::vector<cl::Memory> memObj;
        memObj.push_back(_outputMem);
        _queueCL.enqueueAcquireGLObjects(&memObj);
        _queueCL.enqueueNDRangeKernel(
                    _raycastKernel, cl::NullRange, globalThreads, cl::NullRange, NULL, &ndrEvt);
        _queueCL.enqueueReleaseGLObjects(&memObj);
        _queueCL.finish();    // global sync

        // Profiling: enable before usage
//        cl_ulong start = 0;
//        cl_ulong end = 0;
//        ndrEvt.getProfilingInfo(CL_PROFILING_COMMAND_START, &start);
//        ndrEvt.getProfilingInfo(CL_PROFILING_COMMAND_END, &end);
//        double time = static_cast<double>(end - start)*1e-9;
//        std::cout << "Kernel time: " << time << std::endl << std::endl;

        // Only use this method if not employing CL-GL sharing.
        // Downloading the output buffer from device (GPU) to Host memory.
//        cl::Event readEvt;
//        std::array<size_t, 3> origin = {{0, 0, 0}};
//        std::array<size_t, 3> region = {{width, height, 1}};
//        _queueCL.enqueueReadImage(_outputMem,
//                                  CL_TRUE,
//                                  origin,
//                                  region,
//                                  0,
//                                  0,
//                                  _outputData.data(),
//                                  NULL,
//                                  &readEvt);
//        _queueCL.flush();    // global sync
    }
    catch (cl::Error err)
    {
        logCLerror(err);
    }
}


/**
 * @brief VolumeRenderCL::volDataToCLmem
 * @param volumeData
 */
void VolumeRenderCL::volDataToCLmem(const std::vector<char> &volumeData)
{
    if (!_dr.has_data())
        return;
    try
    {
        cl::ImageFormat format;
        format.image_channel_order = CL_R;

        if (_dr.properties().format == "UCHAR")
            format.image_channel_data_type = CL_UNORM_INT8;
        else if (_dr.properties().format == "USHORT")
            format.image_channel_data_type = CL_UNORM_INT16;
        else if (_dr.properties().format == "FLOAT")
            format.image_channel_data_type = CL_FLOAT;
        else
            throw std::invalid_argument("Unknown or invalid volume data format.");

        _volumeMem = cl::Image3D(_contextCL,
                                 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                 format,
                                 _dr.properties().volume_res[0],
                                 _dr.properties().volume_res[1],
                                 _dr.properties().volume_res[2],
                                 0, 0,
                                 (void *)volumeData.data());
    }
    catch (cl::Error err)
    {
        throw std::runtime_error( "ERROR: " + std::string(err.what()) + "("
                                  + getCLErrorString(err.err()) + ")");
    }
}

/**
 * @brief VolumeRenderCL::loadVolumeData
 * @param fileName
 */
void VolumeRenderCL::loadVolumeData(const std::string fileName)
{
    this->_volLoaded = false;
    std::cout << "Loading volume data defined in " << fileName << std::endl;
    try
    {
        _dr.read_files(fileName);
        std::cout << _dr.data().size() << " bytes have been read." << std::endl;
        std::cout << _dr.properties().to_string() << std::endl;
        volDataToCLmem(_dr.data());
        calcScaling();
    }
    catch (std::runtime_error e)
    {
        std::cerr << e.what() << std::endl;
        throw e;
    }
    // set initally a simple linear transfer function
    std::vector<unsigned char> tff(8*4, 0);
    std::iota(tff.begin() + 2, tff.end(), 0);
    setTransferFunction(tff);

    this->_volLoaded = true;
}


/**
 * @brief VolumeRenderCL::hasData
 * @return
 */
bool VolumeRenderCL::hasData()
{
    return this->_volLoaded;
}


/**
 * @brief VolumeRenderCL::getResolution
 * @return
 */
const std::array<unsigned int, 3> VolumeRenderCL::getResolution() const
{
    if (!_dr.has_data())
        return std::array<unsigned int, 3> {{0, 0, 0}};
    return _dr.properties().volume_res;
}


/**
 * @brief VolumeRenderCL::setTransferFunction
 * @param tff
 */
void VolumeRenderCL::setTransferFunction(std::vector<unsigned char> &tff)
{
    if (!_dr.has_data())
        return;

    try
    {
        cl::ImageFormat format;
        format.image_channel_order = CL_RGBA;
        format.image_channel_data_type = CL_UNORM_INT8;

        cl_mem_flags flags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;
        _tffMem = cl::Image1D(_contextCL, flags, format, tff.size(), tff.data());
    }
    catch (cl::Error err)
    {
        logCLerror(err);
    }
}


/**
 * @brief VolumeRenderCL::setCamOrtho
 * @param setCamOrtho
 */
void VolumeRenderCL::setCamOrtho(bool setCamOrtho)
{
    if (!this->hasData())
        return;

    try {
        _raycastKernel.setArg(ORTHO, (cl_uint)setCamOrtho);
    } catch (cl::Error err) { logCLerror(err); }
}


/**
 * @brief VolumeRenderCL::setIllumination
 * @param illum
 */
void VolumeRenderCL::setIllumination(bool illum)
{
    try {
        _raycastKernel.setArg(ILLUMINATION, (cl_uint)illum);
    } catch (cl::Error err) { logCLerror(err); }
}


/**
 * @brief VolumeRenderCL::setBoundingBox
 * @param boundingBox
 */
void VolumeRenderCL::setBoundingBox(bool boundingBox)
{
    try {
        _raycastKernel.setArg(BOX, (cl_uint)boundingBox);
    } catch (cl::Error err) { logCLerror(err); }
}


/**
 * @brief VolumeRenderCL::setLinearSampling
 * @param linearSampling
 */
void VolumeRenderCL::setLinearInterpolation(bool linearSampling)
{
    try {
        _raycastKernel.setArg(LINEAR, (cl_uint)linearSampling);
    } catch (cl::Error err) { logCLerror(err); }
}


/**
 * @brief VolumeRenderCL::setBackground
 * @param color
 */
void VolumeRenderCL::setBackground(std::array<float, 4> color)
{
    cl_float3 bgColor = {{color[0], color[1], color[2], color[3]}};
    try {
        _raycastKernel.setArg(BACKGROUND, bgColor);
    } catch (cl::Error err) { logCLerror(err); }
}
