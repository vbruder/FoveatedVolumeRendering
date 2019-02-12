/**
 * \file
 *
 * \author Valentin Bruder
 *
 * \copyright Copyright (C) 2018 Valentin Bruder
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#pragma once

#define CL_QUEUE_PROFILING_ENABLE
#define CL_HPP_ENABLE_EXCEPTIONS

#include <QImage>
#include <QFile>

#include "src/oclutil/openclglutilities.h"
#include "src/oclutil/openclutilities.h"

#include "src/io/datrawreader.h"

#include <valarray>

/**
 * @brief The volume renderer class based on OpenCL.
 */
class VolumeRenderCL
{
public:
    /**
     * @brief The OpenCL kernel argument enum.
     */
    enum kernel_arg
    {
          VOLUME     = 0 // volume data set                         image3d_t
        , BRICKS     = 1 // low resolution brick volume             image3d_t
        , TFF        = 2 // transfer function array                 image1d_t
        , OUTPUT         // output image                            image2d_t
        , SAMPLING_RATE  // step size factor                        cl_float
        , VIEW           // view matrix                             float16
        , ORTHO          // use orthographic camera                 cl_uint (bool)
        , ILLUMINATION   // use illumination                        cl_uint
        , SHOW_ESS       // visualize empty space skipping          cl_uint (bool)
        , LINEAR         // use linear interpolation, not nearest   cl_uint (bool)
        , BACKGROUND     // background color RGBA                   cl_float4
        , TFF_PREFIX     // prefix sum of transfer function         image1d_t
        , AO             // use ambient occlusion                   cl_uint (bool)
        , MODEL_SCALE    // model scaling factor                    cl_float3
        , CONTOURS       // show contour lines                      cl_uint (bool)
        , AERIAL         // use aerial perspective                  cl_uint (bool)
        , IN_HIT_IMG     // input image for image order ESS         image2d_t (UINT)
        , OUT_HIT_IMG    // output image for image order ESS        image2d_t (UINT)
        , IMG_ESS        // image order empty space skipping        cl_uint (bool)
		, RMODE			 // rendering mode							cl_uint
		, GPOINT		 // gaze point								cl_int2
		, SDSAMPLES		 // amount of samples						cl_uint
        , IMAP			 // Index map extends						cl_uint2    // image2d_t
		, SDATA			 // Sampling map buffer						(buffer)
        , MIP_1
        , MIP_2
        , MIP_3
        , MIP_4
	};

	enum ip_kernel_arg
	{
		  IP_INIMG = 0
		, IP_IMAP
		, IP_OUTIMG
        , IP_LAST_FRAMES
		, IP_GPOINT
        , IP_SDSAMPLES
        , IP_FRAME_ID
		, IP_SDATA
        , IP_ID
        , IP_WEIGHT
        , IP_THIS_FRAME
        , IP_FRAME_CNT
        , IP_VIEW_CHANGED
        , IP_GAZE_CHANGED
	};

    // mipmap down-scaling metric
    enum scaling_metric
    {
        MIN = 0,
        MAX,
        AVG,
        DENSITY,
    };

    /**
     * @brief Ctor
     */
    VolumeRenderCL();
    /**
     * @brief Dtor
     */
    ~VolumeRenderCL();

    /**
     * @brief Initialize the volume raycaster, i.e. the OpenCL context, queue and kernel.
     * @param useGL Use OpenGL context sharing.
     * @param useCPU Use CPU device instead of GPU device.
     * @param vendor OpenCL platform vendor.
     * @param deviceName Name of the OpenCL device to use.
     */
    void initialize(const bool useGL = false, const bool useCPU = false, 
					const cl_vendor vendor = VENDOR_ANY,
                    const std::string deviceName = "", const int platformId = -1);

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
     * @brief updateOutputTex
     * @param texId
     */
    void updateOutputTex(const GLuint texId);

    /**
     * @brief Update the output image kernel argument and vector size.
     * @param width The image width in pixels.
     * @param height The image height in pixels.
     */
    void updateOutputImg(const size_t width, const size_t height, const cl_GLuint texId);

    /**
     * @brief Run the actual OpenCL volume raycasting kernel.
     * @param width The image width in pixels, used as one dimension of the global thread size.
     * @param height The image height in pixels, used as one dimension of the global thread size.
     * @param t time series id, defaults to 0 if no time series
     */
     void runRaycast(const size_t width, const size_t height, const size_t t = 0);

     /**
      * @brief Run the actual OpenCL volume raycasting kernel without OpenGL context shring.
      * @param width The image width in pixels, used as one dimension of the global thread size.
      * @param height The image height in pixels, used as one dimension of the global thread size.
      * @param t time series id, defaults to 0 if no time series
      * @param pixel color output data of the frame
      */
     void runRaycastNoGL(const size_t width, const size_t height, const size_t t,
                         std::vector<float> &output);

	 /**
	 * @brief Runs the Raycast but sets the kernel parameters to perform lbg-sampling.
	 * @param width The image width in pixels, used to determine the section of the lbg sampling texture 
	   and also the size of the output texture.
	   Has to be less or equal the one third of the width of the lbg sampling texture.
	 * @param height The image height in pixels, used to determine the section of the lbg sampling texture 
	   and also the size of the output texture.
	   Has to be less or equal the one third of the height of the lbg sampling texture.
	 * @param t time series id, defaults to 0 if no time series
	 */
     void runRaycastLBG(const size_t t);

	 /**
	  * @brief Run the actual OpenCL volume raycasting kernel without OpenGL context shring.
	  * @param width The image width in pixels, used to determine the section of the lbg sampling texture 
	   and also the size of the output texture.
	   Has to be less or equal the one third of the width of the lbg sampling texture.
	  * @param height The image height in pixels, used to determine the section of the lbg sampling texture 
	   and also the size of the output texture.
	   Has to be less or equal the one third of the height of the lbg sampling texture.
	  * @param t time series id, defaults to 0 if no time series
	  * @param pixel color output data of the frame
	  */
	 void runRaycastLBGNoGL(const size_t width, const size_t height, const size_t t,
		 std::vector<float> &output);

	 /**
	 * @brief Interpolate the image previously computed by runRaycastLBG.
	 */
	 void interpolateLBG(const size_t width, const size_t height, GLuint inTexId, GLuint outTexId);

    /**
     * @brief Load volume data from a given .dat file name.
     * @param fileName The full path to the volume data file.
     * @return number of loaded volume time steps
     */
    size_t loadVolumeData(const std::string fileName);

	/**
	 * @brief Load an index map and a sampling map, both are .png files.
	 * @param fileNameIndexMap The full path to the index map file.
	 * @param fileNameSamplingMap The full path to the sampling map file.
	 * @return void
	 */
    void loadIndexAndSamplingMap(const std::string &fileNameIndexMap,
                                 const std::string &fileNameSamplingMap,
                                 const QString &fileNameNeighborIndex,
                                 const QString &fileNameNeighborWeights);

    /**
     * @brief Answers if volume data has been loaded.
     * @return true, if volume data has been loaded, false otherwise.
     */
	bool hasData() const;

    /**
     * @brief Return spacial and temporal resolution of loaded volume data set.
     * @return 4D array containing the resolution in x,y,z direction and number of time steps.
     */
    const std::array<size_t, 4> getResolution() const;

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
    void setTffPrefixSum(std::vector<unsigned int> &tffPrefixSum);

    /**
     * @brief VolumeRenderCL::scaleVolume
     * @param scale
     */
    void scaleVolume(const std::valarray<float> scale);

    /**
     * @brief buildScaledVol
     * @param factor
     * @param metric
     * @param useTff
     * @return
     */
    cl::Image3D buildScaledVol(std::array<unsigned int, 3> factor,
                               scaling_metric metric, bool useTff, cl::Image3D volData);

    /**
     * @brief Set ortographic camera kernel parameter.
     * @param setCamOrtho
     */
    void setCamOrtho(const bool setCamOrtho);
    /**
     * @brief Set show illumination kernel paramter.
     * @param illum
     */
    void setIllumination(const unsigned int illum);
    /**
     * @brief Set show empty space skipping kernel parameter.
     * @param boundingBox
     */
    void setShowESS(const bool showESS);
    /**
     * @brief Set linear sampling kernel parameter.
     * @param linearSampling
     */
    void setLinearInterpolation(const bool linearSampling);
    /**
     * @brief Set show contours kernel paramter.
     * @param contours
     */
    void setContours(const bool contours);
    /**
     * @brief Set aerial perspective kernel parameter.
     * @param aerial
     */
    void setAerial(const bool aerial);
    /**
     * @brief Set image order empty space skipping kernel parameter.
     * @param useEss
     */
    void setImgEss(const bool useEss);
    /**
     * @brief Set object order empty space skipping kernel parameter.
     * @param useEss
     */
    void setObjEss(const bool useEss);
    /**
     * @brief Set background color kernel parameter.
     * @param color
     */
    void setBackground(const std::array<float, 4> color);
	/**
	 * @brief Sets the gaze point to be used by lbg sampling.
	 * @param gaze_point
	 */
	void setGazePoint(QPoint gaze_point);

	/**
	 * @brief Sets the gaze point to be used by lbg sampling.
	 * @param gaze_point
	 */
	void setGazePoint(cl_float2 gaze_point);

    /**
     * @brief Get the execution time of the last kernel run.
     * @return The kernel runtime in seconds.
     */
    double getLastExecTime() const;

	/*
	Updates the parameters of the raycast kernel according to the rendering method.
	*/
	void updateRenderingParameters(unsigned int renderingMethod);

    /**
     * @brief getPlatformNames
     * @return platform names
     */
    const std::vector<std::string> getPlatformNames() const;

    /**
     * @brief getDeviceNames
     * @param platformId OpenCL platform id
     * @param type string describing the OpenCL device type ('GPU' or 'CPU')
     * @return device names
     */
    const std::vector<std::string> getDeviceNames(size_t platformId, const std::string &type);

    /**
     * @brief Get the OpenCL device that is currently in use.
     * @return The name of the current OpenCL device in use.
     *         Returns an empty string if no device is currently used.
     */
    const std::string & getCurrentDeviceName() const;

    /**
     * @brief setAmbientOcclusion
     * @param ao
     */
    void setAmbientOcclusion(const bool ao);

    /**
     * @brief Generate a downsampling of the currently loaded volume file.
     * @param t Timestep to be downsampled
     * @param factor downsampling factor, uniform for all 3 dimensions
     */
    const std::string generateVolumeDownsampling(const size_t t, const int factor);

	// Returns the extends of the index map. {0,0} if it hasn't been loaded yet.
	QPoint getIndexMapExtends();

    /**
     * @brief Generate a full mipmap stack of the currently loaded volume texture.
     * @param levelCnt Number of mipmap levels.
     */
    void generateMipmaps(size_t levelCnt);

private:
    /**
     * @brief Generate coarse grained volume bricks that can be used for ESS.
     * @param volumeData
     */
    void generateBricks();

    /**
     * @brief Downsample a volume data set.
     */
    cl::Image3D downsampleVolume(const cl::ImageFormat &format,
                                       const std::array<size_t, 3> &newSize,
                                       const cl::Image3D &volumeMem);

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
     * @brief Convert volume data to UCHAR format and generate OpenCL image textue memory object.
     * @param volumeData The raw volume data.
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

            _volumesMem.push_back(cl::Image3D(_contextCL,
                                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                              format,
                                              _dr.properties().volume_res[0],
                                              _dr.properties().volume_res[1],
                                              _dr.properties().volume_res[2],
                                              0, 0,
                                              const_cast<unsigned char*>(convertedData.data())));
        }
        catch (cl::Error err)
        {
            throw std::runtime_error( "ERROR: " + std::string(err.what()) + "("
                                      + getCLErrorString(err.err()) + ")");
        }
    }

    /**
     * @brief Set OpenCL memory objects for volume raycast kernel.
     * @param t number of volume timesteps.
     */
    void setMemObjectsRaycast(const size_t t);

	/**
	 * @brief Set OpenCL memory objects for interpolateLBG Kernel.
	 */
	void setMemObjectsInterpolationLBG(GLuint inTexId, GLuint outTexId);

    /**
     * @brief Set OpenCL memory objects for brick generation kernel.
     * @param t number of volume timesteps.
     * @throws runtime_error when the number of volume timesteps
     *         is not equal to number of brick timesteps
     */
    void setMemObjectsBrickGen(const size_t t);

    /**
     * @brief Initialize OpenCL kernel with default paramters.
     * @param fileName File name of the kernel source file.
     * @param buildFlags Comiler flags for building the kernel.
     */
    void initKernel(const std::string fileName, const std::string buildFlags = "");

    /**
     * @brief Log a OpenCL error message.
     * @param err The OpenCL error object to be logged.
     * @throws Runtime error
     **/
    [[ noreturn ]] void logCLerror(const cl::Error err) const;

    // *** member variables ***
    cl::Context _contextCL;
    cl::CommandQueue _queueCL;
    cl::Kernel _raycastKernel;
    cl::Kernel _genBricksKernel;
    cl::Kernel _downsamplingKernel;
	cl::Kernel _interpolateLBGKernel;

    std::vector<cl::Image3D> _volumesMem;
    std::vector<cl::Image3D> _bricksMem;
    cl::ImageGL _outputMem;
	cl::ImageGL _inputMem;	// Image Object to hold temporary Images like pre interpolation for LBG
    cl::Image1D _tffMem;
    cl::Image1D _tffPrefixMem;
    cl::Image2D _outputMemNoGL;
    cl::Image2D _outputHitMem;
    cl::Image2D _inputHitMem;
	cl::Buffer _place_holder_smd;
	cl::Image2D _place_holder_imap;
	cl::Buffer _samplingMapData; // The width of the sample map defines the total number of work items to be started for lbg sampling raycast
	cl::Image2D _indexMap;
    cl::Image2DArray _lastFramesMem;
    cl::Image2D _thisFrameMem;
    cl::Buffer _neighborIdMap;
    cl::Buffer _neighborWeightMap;
    std::vector<cl::Image3D> _volMipmapsMem;

    QPoint _indexMapExtends; // The width and height of the Index Map
	size_t _amountOfSamples;	// The indexes of the samplingMap (scanlLine width).

    bool _imsmLoaded = false;	// is set to true iff _samplingMapData and _indexMap have been loaded successfully.
    bool _volLoaded = false;
    double _lastExecTime = 1.0;
    std::valarray<float> _modelScale;
    bool _useGL = true;
    bool _useImgESS = false;
    std::string _currentDevice;
    cl_uint _frameId = 0;
    cl_uint _frameIpCnt = 5;    // max 8
    cl_uint _viewChanged = false;
    cl_uint _gazeChanged = false;
    cl_float2 _gazePoint = {{0,0}};

    std::vector<float> _output;

    DatRawReader _dr;
};
