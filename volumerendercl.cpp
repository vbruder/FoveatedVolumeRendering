#include "volumerendercl.h"

#include <functional>
#include <algorithm>
#include <numeric>
#include "omp.h"

/**
 * @brief VolumeRenderCL::VolumeRenderCL
 */
VolumeRenderCL::VolumeRenderCL() :
    _volLoaded(false)
  , _lastExecTime(0.0)
  , _modelScale{1.0, 1.0, 1.0}
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
        // NOTE: replace if no NVIDIA GPU
        _contextCL = createCLGLContext(CL_DEVICE_TYPE_GPU, VENDOR_NVIDIA); // VENDOR_ANY
        cl_command_queue_properties cqp = 0;
#ifdef CL_QUEUE_PROFILING_ENABLE
        cqp = CL_QUEUE_PROFILING_ENABLE;
#endif
        _queueCL = cl::CommandQueue(_contextCL, cqp);
    }
    catch (cl::Error err)
    {
        logCLerror(err);
    }

    initKernel("../RaycastLight/kernels/volumeraycast.cl",
               "-DIMAGE_SUPPORT=1 -DCL_STD=CL1.2 -DESS");
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
        _raycastKernel.setArg(SAMPLING_RATE, 1.5f);     // default step size 0.5*voxel size
        _raycastKernel.setArg(ORTHO, 0);                // perspective cam by default
        _raycastKernel.setArg(ILLUMINATION, 1);         // illumination on by default
        _raycastKernel.setArg(BOX, 0);
        _raycastKernel.setArg(LINEAR, 1);
        cl_float4 bgColor = {{1.f, 1.f, 1.f, 1.f}};
        _raycastKernel.setArg(BACKGROUND, bgColor);
        _raycastKernel.setArg(AO, 1);                   // ambient occlusion on by default

        _genBricksKernel = cl::Kernel(program, "generateBricks");
    }
    catch (cl::Error err)
    {
        logCLerror(err);
    }
}


/**
 * @brief VolumeRenderCL::setMemObjects
 */
void VolumeRenderCL::setMemObjectsRaycast(const int t)
{
    _raycastKernel.setArg(VOLUME, _volumesMem.at(t));
    _raycastKernel.setArg(TFF, _tffMem);
    _raycastKernel.setArg(BRICKS, _bricksMem.at(t));
    _raycastKernel.setArg(OUTPUT, _outputMem);
    _raycastKernel.setArg(TFF_PREFIX, _tffPrefixMem);
    cl_float3 modelScale = {_modelScale[0], _modelScale[1], _modelScale[2]};
    _raycastKernel.setArg(MODEL_SCALE, modelScale);
}


/**
 * @brief VolumeRenderCL::setMemObjectsBrickGen
 */
void VolumeRenderCL::setMemObjectsBrickGen(const int t)
{
    _genBricksKernel.setArg(VOLUME, _volumesMem.at(t));
    _genBricksKernel.setArg(TFF, _tffMem);
    _genBricksKernel.setArg(BRICKS, _bricksMem.at(t));
}


/**
 * @brief VolumeRenderCL::calcScaling
 */
void VolumeRenderCL::calcScaling()
{
    if (!_dr.has_data())
        return;

    _modelScale = { static_cast<float>(_dr.properties().volume_res.at(0)),
                    static_cast<float>(_dr.properties().volume_res.at(1)),
                    static_cast<float>(_dr.properties().volume_res.at(2)) };

    std::valarray<float> thickness = { static_cast<float>(_dr.properties().slice_thickness.at(0)),
                                       static_cast<float>(_dr.properties().slice_thickness.at(1)),
                                       static_cast<float>(_dr.properties().slice_thickness.at(2)) };
    _modelScale *= thickness*(1.f/thickness[0]);
#undef max  // error here if we don't undef max
    _modelScale = _modelScale.max() / _modelScale;
}


/**
 * @brief VolumeRenderCL::scaleVolume
 * @param scale
 */
void VolumeRenderCL::scaleVolume(std::valarray<float> scale)
{
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

    cl_float16 view;
    for (size_t i = 0; i < 16; ++i)
        view.s[i] = viewMat[i];
    try{
        _raycastKernel.setArg(VIEW, view);
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
void VolumeRenderCL::updateOutputImg(const size_t width, const size_t height, GLuint texId)
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
#ifdef NO_GL
        _outputData.resize(width * height * 4, 0);
#endif
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
void VolumeRenderCL::runRaycast(const size_t width, const size_t height, const int t)
{
    if (!this->_volLoaded)
        return;
    try // opencl scope
    {
        setMemObjectsRaycast(t);
        cl::NDRange globalThreads(width, height);
        cl::Event ndrEvt;

        std::vector<cl::Memory> memObj;
        memObj.push_back(_outputMem);
        _queueCL.enqueueAcquireGLObjects(&memObj);
        _queueCL.enqueueNDRangeKernel(
                    _raycastKernel, cl::NullRange, globalThreads, cl::NullRange, NULL, &ndrEvt);
        _queueCL.enqueueReleaseGLObjects(&memObj);
        _queueCL.finish();    // global sync

#ifdef NO_GL
        cl::Event readEvt;
        std::array<size_t, 3> origin = {{0, 0, 0}};
        std::array<size_t, 3> region = {{width, height, 1}};
        _queueCL.enqueueReadImage(_outputMem,
                                  CL_TRUE,
                                  origin,
                                  region,
                                  0,
                                  0,
                                  _outputData.data(),
                                  NULL,
                                  &readEvt);
        _queueCL.flush();    // global sync
#endif

#ifdef CL_QUEUE_PROFILING_ENABLE
        cl_ulong start = 0;
        cl_ulong end = 0;
        ndrEvt.getProfilingInfo(CL_PROFILING_COMMAND_START, &start);
        ndrEvt.getProfilingInfo(CL_PROFILING_COMMAND_END, &end);
        _lastExecTime = static_cast<double>(end - start)*1e-9;
//        std::cout << "Kernel time: " << _lastExecTime << std::endl << std::endl;
#endif
    }
    catch (cl::Error err)
    {
        logCLerror(err);
    }
}


/**
 * @brief RoundPow2
 * @param iNumber
 * @return
 */
static uint RoundPow2(uint n)
{
    // next highest power of 2
    // (cf: http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2)
    uint val = n - 1u;
    val |= val >> 1;
    val |= val >> 2;
    val |= val >> 4;
    val |= val >> 8;
    val |= val >> 16;
    val++;
    // previous power of 2
    uint x = val >> 1;
    // round to nearest of the two
    return (val - n) > (n - x) ? x : val;
}

/**
 * @brief VolumeRenderCL::generateBricks
 * @param volumeData
 */
void VolumeRenderCL::generateBricks()
{
    if (!_dr.has_data())
        return;
    try
    {
        // calculate brick size
        const uint numBricks = 64u;
        std::array<uint, 3> brickRes = {1u, 1u, 1u};
        brickRes.at(0) = RoundPow2(_dr.properties().volume_res.at(0)/numBricks);
        brickRes.at(1) = RoundPow2(_dr.properties().volume_res.at(1)/numBricks);
        brickRes.at(2) = RoundPow2(_dr.properties().volume_res.at(2)/numBricks);
        std::array<uint, 3> bricksTexSize = {1u, 1u, 1u};
        bricksTexSize.at(0) = ceil(_dr.properties().volume_res.at(0)/(double)brickRes.at(0));
        bricksTexSize.at(1) = ceil(_dr.properties().volume_res.at(1)/(double)brickRes.at(1));
        bricksTexSize.at(2) = ceil(_dr.properties().volume_res.at(2)/(double)brickRes.at(2));

        // set memory object
        cl::ImageFormat format;
        format.image_channel_order = CL_RG;  // NOTE: CL_RG for min+max

        if (_dr.properties().format == "UCHAR")
            format.image_channel_data_type = CL_UNORM_INT8;
        else if (_dr.properties().format == "USHORT")
            format.image_channel_data_type = CL_UNORM_INT16;
        else if (_dr.properties().format == "FLOAT")
            format.image_channel_data_type = CL_FLOAT;
        else
            throw std::invalid_argument("Unknown or invalid volume data format.");

        _bricksMem.clear();
        for (size_t i = 0; i < _dr.properties().raw_file_names.size(); ++i)
        {
            _bricksMem.push_back(cl::Image3D(_contextCL,
                                             CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS,
                                             format,
                                             bricksTexSize.at(0),
                                             bricksTexSize.at(1),
                                             bricksTexSize.at(2),
                                             0, 0,
                                             NULL));
            // run aggregation kernel
            setMemObjectsBrickGen(i);
            cl::NDRange globalThreads(bricksTexSize.at(0), bricksTexSize.at(1), bricksTexSize.at(2));
            cl::Event ndrEvt;
            _queueCL.enqueueNDRangeKernel(
                        _genBricksKernel, cl::NullRange, globalThreads, cl::NullRange, NULL, &ndrEvt);
            _queueCL.finish();    // global sync
        }
    }
    catch (cl::Error err)
    {
        throw std::runtime_error( "ERROR: " + std::string(err.what()) + "("
                                  + getCLErrorString(err.err()) + ")");
    }
}

/**
 * @brief VolumeRenderCL::volDataToCLmem
 * @param volumeData
 */
void VolumeRenderCL::volDataToCLmem(const std::vector<std::vector<char>> &volumeData)
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

        _volumesMem.clear();
        for (const auto &v : volumeData)
        {
            _volumesMem.push_back(cl::Image3D(_contextCL,
                                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              format,
                                              _dr.properties().volume_res[0],
                                              _dr.properties().volume_res[1],
                                              _dr.properties().volume_res[2],
                                              0, 0,
                                              (void *)v.data()));
        }
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
int VolumeRenderCL::loadVolumeData(const std::string fileName)
{
    this->_volLoaded = false;
    std::cout << "Loading volume data defined in " << fileName << std::endl;
    try
    {
        _dr.read_files(fileName);
        std::cout << _dr.data().front().size()*_dr.data().size() << " bytes have been read from "
                  << _dr.data().size() << " file(s)." << std::endl;
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
    std::vector<unsigned char> tff(256*4, 0);
    std::iota(tff.begin() + 3, tff.end(), 0);
    setTransferFunction(tff);

    std::vector<ushort> prefixSum(256, 0);
#pragma omp for
    for (int i = 0; i < (int)prefixSum.size(); ++i)
        prefixSum.at(i) = i*4;

    std::partial_sum(prefixSum.begin(), prefixSum.end(), prefixSum.begin());
    setTffPrefixSum(prefixSum);

    this->_volLoaded = true;
    return _dr.data().size();
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
        // divide size by 4 because of RGBA
        _tffMem = cl::Image1D(_contextCL, flags, format, tff.size() / 4, tff.data());
        generateBricks();

        std::vector<ushort> prefixSum;
        // copy only alpha values (every fourth element)
        for (int i = 3; i < static_cast<int>(tff.size()); i += 4)
            prefixSum.push_back(static_cast<ushort>(tff.at(i)));
        std::partial_sum(prefixSum.begin(), prefixSum.end(), prefixSum.begin());
        setTffPrefixSum(prefixSum);
    }
    catch (cl::Error err)
    {
        logCLerror(err);
    }
}


void VolumeRenderCL::setTffPrefixSum(std::vector<unsigned short> &tffPrefixSum)
{
    if (!_dr.has_data())
        return;

    try
    {
        cl::ImageFormat format;
        format.image_channel_order = CL_R;
        format.image_channel_data_type = CL_UNSIGNED_INT16;

        cl_mem_flags flags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;
        // divide size by 4 because of RGBA
        _tffPrefixMem = cl::Image1D(_contextCL, flags, format, tffPrefixSum.size(), tffPrefixSum.data());
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
 * @brief VolumeRenderCL::setAmbientOcclusion
 * @param illum
 */
void VolumeRenderCL::setAmbientOcclusion(bool ao)
{
    try {
        _raycastKernel.setArg(AO, (cl_uint)ao);
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


/**
 * @brief VolumeRenderCL::getLastExecTime
 * @return
 */
double VolumeRenderCL::getLastExecTime()
{
    return _lastExecTime;
}
