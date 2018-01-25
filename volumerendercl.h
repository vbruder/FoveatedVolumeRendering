#ifndef VOLUMERENDERCL_H
#define VOLUMERENDERCL_H

#define CL_QUEUE_PROFILING_ENABLE
#define CL_HPP_ENABLE_EXCEPTIONS
#include "openclglutilities.h"

#include "datrawreader.h"

#include <valarray>

class VolumeRenderCL
{
public:
    enum kernel_arg
    {
          VOLUME     = 0 // volume data set                          image3d_t
        , TFF        = 1 // transfer function array                  image1d_t
        , BRICKS         // low resolution brick volume              image3d_t
        , OUTPUT         // output image                             image2d_t
        , SAMPLING_RATE  // step size factor                         cl_float
        , VIEW           // view matrix                              float16
        , ORTHO          // use orthographic camera                  cl_uint (bool)
        , ILLUMINATION   // use illumination (per view)              cl_uint (bool)
        , BOX            // show bounding box aroud volume           cl_uint (bool)
        , LINEAR         // use linear interpolation, not nearest    cl_uint (bool)
        , BACKGROUND     // background color RGBA                    cl_float4
        , TFF_PREFIX     // prefix sum of transfer function          image1d_t
        , AO             // use ambient occlusion                    cl_uint (bool)
        , MODEL_SCALE    // model scaling factor                     cl_float3
    };

    // mipmap down-scaling metric
    enum scaling_metric
    {
        MIN = 0,
        MAX,
        AVG,
        DENSITY,
    };

    VolumeRenderCL();

    ~VolumeRenderCL();

    /**
     * @brief Initialize the volume raycaster, i.e. the OpenCL context, queue and kernel.
     */
    void initialize();

    /**
     * @brief Update the view matrix argument for the raycasting kernel.
     * @param viewMat the 4x4 transposed view matrix.
     */
    void updateView(const std::array<float, 16> viewMat);

    /**
     * @brief Update the integration step size factor kernel argument for the volume raycast.
     * @param samplingRate the sampling rate relative per voxel.
     */
    void updateSamplingRate(const double samplingRate);

    /**
     * @brief Update the output image kernel argument and vector size.
     * @param width The image width in pixels.
     * @param height The image height in pixels.
     */
    void updateOutputImg(const size_t width, const size_t height, cl_GLuint texId);

    /**
     * @brief Run the actual OpenCL volume raycasting kernel.
     * @param width The image width in pixels, used as one dimension of the global thread size.
     * @param height The image height in pixels, used as one dimension of the global thread size.
     * @param t time series id, defaults to 0 if no time series
     */
     void runRaycast(const size_t width, const size_t height, const int t = 0);

    /**
     * @brief Load volume data from a given .dat file name.
     * @param fileName The full path to the volume data file.
     * @return number of loaded volume time steps
     */
    int loadVolumeData(const std::string fileName);

    /**
     * @brief Answers if volume data has been loaded.
     * @return true, if volume data has been loaded, false otherwise.
     */
    bool hasData();

    /**
     * @brief getResolution
     * @return
     */
    const std::array<unsigned int, 3> getResolution() const;

    /**
     * @brief Set the transfer function for the volume raycast as a kernel argument.
     * @param tff a vector of the RGBA transfer function values.
     * @param rangeMin clamp range to minimum
     * @param rangeMax clamp range to maximum
     */
    void setTransferFunction(std::vector<unsigned char> &tff);

    /**
     * @brief Set the prefix sum of the transfer function.
     * @param tffPrefixSum The prefix sum as vector of unsigned ints.
     */
    void setTffPrefixSum(std::vector<unsigned short> &tffPrefixSum);

    /**
     * @brief VolumeRenderCL::scaleVolume
     * @param scale
     */
    void scaleVolume(std::valarray<double> scale);

    /**
     * @brief buildScaledVol
     * @param factor
     * @param metric
     * @param useTff
     * @return
     */
    cl::Image3D buildScaledVol(std::array<unsigned int, 3> factor,
                               scaling_metric metric, bool useTff, cl::Image3D volData);

    // TODO: set background / toggle black/white

    /**
     * @brief setCamOrtho
     * @param setCamOrtho
     */
    void setCamOrtho(bool setCamOrtho);
    /**
     * @brief setIllumination
     * @param illum
     */
    void setIllumination(bool illum);
    /**
     * @brief setBoundingBox
     * @param boundingBox
     */
    void setBoundingBox(bool boundingBox);
    /**
     * @brief setLinearSampling
     * @param linearSampling
     */
    void setLinearInterpolation(bool linearSampling);
    /**
     * @brief setBackground
     * @param color
     */
    void setBackground(std::array<float, 4> color);

    /**
     * @brief Get the execution time of the last kernel run.
     * @return The kernel runtime in seconds.
     */
    double getLastExecTime();

    /**
     * @brief setAmbientOcclusion
     * @param ao
     */
    void setAmbientOcclusion(bool ao);

private:
    /**
     * @brief Generate coarse grained volume bricks that can be used for ESS.
     * @param volumeData
     */
    void generateBricks();

    /**
     * @brief Calculate the scaling vector for the volume data.
     */
    void calcScaling();

    /**
     * @brief volDataToCLmem
     * @param volumeData
     */
    void volDataToCLmem(const std::vector<std::vector<char> > &volumeData);

    /**
     *
     */
    template<class T>
    void volDataToCLmem(const std::vector<char> &volumeData)
    {
        // reinterpret raw data (char) to input format
        auto s = reinterpret_cast<const T *>(volumeData.data());
        auto e = reinterpret_cast<const T *>(volumeData.data() + volumeData.size());
        // convert imput vector to the desired output precision
        std::vector<unsigned char> convertedData(s, e);
        try
        {
            cl::ImageFormat format;
            format.image_channel_order = CL_R;
            format.image_channel_data_type = CL_UNORM_INT8;

            _volumesMem = cl::Image3D(_contextCL,
                                     CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                     format,
                                     _dr.properties().volume_res[0],
                                     _dr.properties().volume_res[1],
                                     _dr.properties().volume_res[2],
                                     0, 0,
                                     (void *)convertedData.data());
        }
        catch (cl::Error err)
        {
            throw std::runtime_error( "ERROR: " + std::string(err.what()) + "("
                                      + getCLErrorString(err.err()) + ")");
        }
    }

    void setMemObjectsRaycast(const int t);
    void setMemObjectsBrickGen(const int t);

    void initKernel(const std::string fileName, const std::string buildFlags = "");

    void logCLerror(cl::Error err);

    cl::Context _contextCL;
    cl::CommandQueue _queueCL;
    cl::Kernel _raycastKernel;
    cl::Kernel _genBricksKernel;

    std::vector<cl::Image3D> _volumesMem;
    std::vector<cl::Image3D> _bricksMem;
    cl::ImageGL _outputMem;
    cl::ImageGL _overlayMem;
    cl::Image1D _tffMem;
    cl::Image1D _tffPrefixMem;
#ifdef NO_GL
    std::vector<float> _outputData;
#endif
    bool _volLoaded;
    double _lastExecTime;
    std::valarray<double> _modelScale;
    DatRawReader _dr;
};

#endif // VOLUMERENDERCL_H
