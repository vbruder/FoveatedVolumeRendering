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

        std::array<float, 16> zeroMat;
        zeroMat.fill(0);
        _viewMem = cl::Buffer(_contextCL,
                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                              zeroMat.size() * sizeof(cl_float),
                              zeroMat.data());
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
        _raycastKernel.setArg(VIEW, _viewMem);
        _raycastKernel.setArg(STEP_SIZE, 0.5f);     // default step size 0.5*voxel size
        _raycastKernel.setArg(ORTHO, 0);            // perspective cam by default
        _raycastKernel.setArg(ILLUMINATION, 1);     // illumination on by default

        // TESTING
//        std::vector<cl_uint> gCnt(1, 0);
//        _gCountMem = cl::Buffer(_contextCL, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//                                    sizeof(cl_int), (void*)gCnt.data());
//        _raycastKernel.setArg(__TEST_GCOUNT, _gCountMem);

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
    _raycastKernel.setArg(VIEW, _viewMem);

//    std::vector<cl_uint> gCnt(1, 0);
//    _gCountMem = cl::Buffer(_contextCL, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//                                sizeof(cl_int), (void*)gCnt.data());
//    _raycastKernel.setArg(__TEST_GCOUNT, _gCountMem);
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

    std::array<float, 16> modelViewMat;
    for (size_t i = 0; i < modelViewMat.size(); ++i)
    {
        if (i < 4)
            modelViewMat[i] = viewMat[i] * _modelScale[0];
        else if (i < 8)
            modelViewMat[i] = viewMat[i] * _modelScale[1];
        else if (i < 12)
            modelViewMat[i] = viewMat[i] * _modelScale[2];
        else
            modelViewMat[i] = viewMat[i];
    }

    try{
        _queueCL.enqueueWriteBuffer(_viewMem, CL_TRUE, 0, 16*sizeof(cl_float), modelViewMat.data());
    } catch (cl::Error err) {
        logCLerror(err);
    }
}


/**
 * @brief VolumeRenderCL::updateStepSize
 * @param stepSize
 */
void VolumeRenderCL::updateStepSize(const double stepSize)
{
    try{
        _raycastKernel.setArg(STEP_SIZE, static_cast<cl_float>(stepSize));
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

        cl_uint4 resolution = {{_dr.properties().volume_res[0],
                                _dr.properties().volume_res[1],
                                _dr.properties().volume_res[2], 0}};
        _raycastKernel.setArg(RESOLUTION, resolution);
        calcScaling();
    }
    catch (std::runtime_error e)
    {
        std::cerr << e.what() << std::endl;
        throw e;
    }
    // set initally a simple linear transfer function of gray values
    std::vector<unsigned char> tff(4*4, 0);
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
        _raycastKernel.setArg(ORTHO, setCamOrtho);
    } catch (cl::Error err) { logCLerror(err); }
}


/**
 * @brief VolumeRenderCL::setIllumination
 * @param illum
 * @param viewId
 */
void VolumeRenderCL::setIllumination(bool illum)
{
    try {
        _raycastKernel.setArg(ILLUMINATION, illum);
    } catch (cl::Error err) { logCLerror(err); }
}
