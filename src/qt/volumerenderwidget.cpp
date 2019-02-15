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

#include "src/qt/volumerenderwidget.h"

#include <QPainter>
#include <QGradient>
#include <QImage>
#include <QCoreApplication>
#include <QScreen>
#include <QInputDialog>
#include <QJsonObject>
#include <QErrorMessage>
#include <QLoggingCategory>
#include <QMessageBox>
#include <QFileDialog>

#include <thread>
#include <chrono>

const static double Z_NEAR = 1.0;
const static double Z_FAR = 500.0;

static const char *pVsScreenQuadSource =
    "#version 330\n"
    "layout(location = 0) in highp vec3 vertex;\n"
    "out highp vec2 texCoord;\n"
    "uniform mat4 projMatrix;\n"
    "uniform mat4 mvMatrix;\n"
    "void main() {\n"
    "   texCoord = vec2(0.5f) + 0.5f * vertex.xy;\n"
    "   gl_Position = projMatrix * mvMatrix * vec4(vertex.xy, 1.0f, 1.0f);\n"
    "}\n";

static const char *pFsScreenQuadSource =
    "#version 330\n"
    "in highp vec2 texCoord;\n"
    "out highp vec4 fragColor;\n"
    "uniform highp sampler2D outTex;\n"
    "void main() {\n"
    "   fragColor = texture(outTex, texCoord);\n"
    "   fragColor.a = 1.0f;\n"
    "}\n";


/**
 * @brief VolumeRenderWidget::VolumeRenderWidget
 * @param parent
 */
VolumeRenderWidget::VolumeRenderWidget(QWidget *parent)
    : QOpenGLWidget(parent)
    , _tffRange(QPoint(0, 255))
    , _timestep(0)
    , _eyetracker(nullptr)
    , _monitor_offset(QPoint(0,0))
    , _curr_monitor_width(0)
    , _curr_monitor_height(0)
    , _lastLocalCursorPos(QPoint(0,0))
    , _rotQuat(QQuaternion(1, 0, 0, 0))
    , _translation(QVector3D(0, 0, 2.0))
    , _useEyetracking(false)
    , _noUpdate(true)
    , _loadingFinished(false)
    , _writeImage(false)
    , _recordVideo(false)
    , _imgCount(0)
    , _imgSamplingRate(1)
    , _useGL(true)
    , _showOverlay(true)
    , _logView(false)
    , _logInteraction(false)
    , _contRendering(false)
    , _renderingMethod(Standard)
{
    this->setMouseTracking(true);
}


/**
 * @brief VolumeRenderWidget::~VolumeRenderWidget
 */
VolumeRenderWidget::~VolumeRenderWidget()
{
}


/**
 * @brief VolumeRenderWidget::paintOrientationAxis
 */
void VolumeRenderWidget::paintOrientationAxis(QPainter &p)
{
    QMatrix4x4 proj;
    proj.perspective(53.14f, 1.0f, 0.1f, 1.0);
    QMatrix4x4 viewProj(proj * _coordViewMX);
    QVector4D x =         viewProj * QVector4D( 20,  0,  0, 0);
    QVector4D xArrLeft =  viewProj * QVector4D( 16, -2,  0, 0);
    QVector4D xArrRight = viewProj * QVector4D( 16, +2,  0, 0);
    QVector4D y =         viewProj * QVector4D(  0, 20,  0, 0);
    QVector4D yArrLeft =  viewProj * QVector4D( -2, 16,  0, 0);
    QVector4D yArrRight = viewProj * QVector4D( +2, 16,  0, 0);
    QVector4D z =         viewProj * QVector4D(  0,  0, 20, 0);
    QVector4D zArrLeft =  viewProj * QVector4D( -2,  0, 16, 0);
    QVector4D zArrRight = viewProj * QVector4D( +2,  0, 16, 0);

    p.resetTransform();
    p.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    p.translate(66, height() - 66);
    int textOffset = 5;
    // x axis
    p.setPen(Qt::red);
    p.drawLine(0, 0, x.x(), x.y());
    p.drawLine(xArrLeft.x(), xArrLeft.y(), x.x(), x.y());
    p.drawLine(xArrRight.x(), xArrRight.y(), x.x(), x.y());
    p.drawText(x.x() + textOffset, x.y() + textOffset, "x");
    // y axis
    p.setPen(Qt::green);
    p.drawLine(0, 0, y.x(), y.y());
    p.drawLine(yArrLeft.x(), yArrLeft.y(), y.x(), y.y());
    p.drawLine(yArrRight.x(), yArrRight.y(), y.x(), y.y());
    p.drawText(y.x() + textOffset, y.y() + textOffset, "y");
    // z axis
    p.setPen(Qt::blue);
    p.drawLine(0, 0, z.x(), z.y());
    p.drawLine(zArrLeft.x(), zArrLeft.y(), z.x(), z.y());
    p.drawLine(zArrRight.x(), zArrRight.y(), z.x(), z.y());
    p.drawText(z.x() + textOffset, z.y() + textOffset, "z");
}


/**
 * @brief paintFPS
 */
void VolumeRenderWidget::paintFps(QPainter &p, const double fps, const double lastTime)
{
    p.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    p.setPen(Qt::darkGreen);
    p.setFont(QFont("Helvetica", 11));
    QString s = "FPS: " + QString::number(fps);
    p.drawText(10, 20, s);
    s = "Last: " + QString::number(lastTime);
    p.drawText(10, 36, s);
    s = QString(_volumerender.getCurrentDeviceName().c_str());
    p.drawText(10, 52, s);
}


/**
 * @brief VolumeRenderWidget::initializeGL
 */
void VolumeRenderWidget::initializeGL()
{
    connect( context(), &QOpenGLContext::aboutToBeDestroyed, this,
             &VolumeRenderWidget::cleanup );

    initializeOpenGLFunctions();
    makeCurrent();

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    _spScreenQuad.addShaderFromSourceCode(QOpenGLShader::Vertex, pVsScreenQuadSource);
    _spScreenQuad.addShaderFromSourceCode(QOpenGLShader::Fragment, pFsScreenQuadSource);
    _spScreenQuad.bindAttributeLocation("vertex", 0);
    _spScreenQuad.link();

    _spScreenQuad.bind();
    _screenQuadVao.create();
    _screenQuadVao.bind();

    const int numQuadVertices = 8;
    GLfloat quadVertices[numQuadVertices] =
    {
        -1.0f, -1.0f,
        +1.0f, -1.0f,
        -1.0f, +1.0f,
        +1.0f, +1.0f
    };

    // Setup vertex buffer object.
    _quadVbo.create();
    _quadVbo.bind();
    _quadVbo.allocate( quadVertices, numQuadVertices * sizeof( GLfloat ) );
    _quadVbo.release();
    // Store the vertex attribute bindings for the program.
    setupVertexAttribs();
    _viewMX.setToIdentity();
    _viewMX.translate( 0.0f, 0.0f, -1.0f );
    // set quad model matrix
    _modelMX.setToIdentity();
    _modelMX.rotate( 180.0f, 1, 0, 0 );
    _spScreenQuad.release();
    _screenQuadVao.release();

    initVolumeRenderer();
	// FIXME: dual gpu setup
	//initVolumeRenderer(false, false); 
}

/**
 * @brief VolumeRenderWidget::initVolumeRenderer
 */
void VolumeRenderWidget::initVolumeRenderer(bool useGL, const bool useCPU)
{
    try
    {
        // try to use GPU with GL context sharing
        _volumerender.initialize(useGL, useCPU);
    }
    catch (std::invalid_argument e)
    {
        qCritical() << e.what();
    }
    catch (std::runtime_error e)
    {
        _useGL = false;
		try
		{
			qCritical() << e.what() << "\nSDisabling OpenGL context sharing.";
			_volumerender.initialize(_useGL, false);
		}
		catch (std::runtime_error e)
		{
			qCritical() << e.what() << "\nSwitching to CPU rendering mode.";
			_volumerender.initialize(_useGL, true);
		}
    }
    catch (...)
    {
        qCritical() << "An unknown error occured initializing OpenCL/OpenGL.";
    }
}

/**
 * @brief VolumeRenderWidget::saveFrame
 */
void VolumeRenderWidget::saveFrame()
{
    _writeImage = true;
    update();
}

/**
 * @brief VolumeRenderWidget::toggleVideoRecording
 */
void VolumeRenderWidget::toggleVideoRecording()
{
    qInfo() << (_recordVideo ? "Stopped recording." : "Started recording.");

    _recordVideo = !_recordVideo;
    _writeImage = true;
    update();
}

/**
 * @brief VolumeRenderWidget::toggleViewRecording
 * Toggle writing camera configuration (rotation as quaternion and translation as a vecor)
 * into two files every time camera parameters change.
 */
void VolumeRenderWidget::toggleViewRecording()
{
    qInfo() << (_logView ? "Stopped view config recording." : "Started view config recording.");

    _logView = !_logView;
    if (_logView)
    {
        QFileDialog dialog;
        _viewLogFile= dialog.getSaveFileName(this, tr("Save camera path"),
                                                QDir::currentPath(), tr("All files"));
        if (_viewLogFile.isEmpty())
            return;
    }
    updateView();
}

/**
 * @brief VolumeRenderWidget::toggleInteractionLogging
 */
void VolumeRenderWidget::toggleInteractionLogging()
{
	qInfo() << (_logInteraction ? "Stopped view config recording." : "Started view config recording.");

	_logInteraction = !_logInteraction;
	if (_logInteraction)
	{
		QFileDialog dialog;
		_interactionLogFile = dialog.getSaveFileName(this, tr("Save camera path"),
			QDir::currentPath(), tr("All files"));
		if (_interactionLogFile.isEmpty())
			return;
		_timer.restart();

		// log initial configuration
		QString s;
		s += QString::number(_timer.elapsed());
		s += "; tffInterpolation; ";
		s += _tffInterpol == QEasingCurve::InOutQuad ? "quad" : "linear";
		s += "\n";
		// tff
		s += QString::number(_timer.elapsed());
		s += "; transferFunction; ";
		const std::vector<unsigned char> tff = getRawTransferFunction(_tffStops);
		foreach(unsigned char c, tff)
		{
			s += QString::number(static_cast<int>(c)) + " ";
		}
		s += "\n";
		// camera
		s += QString::number(_timer.elapsed());
		s += "; camera; ";
		s += QString::number(_rotQuat.toVector4D().w()) + " ";
		s += QString::number(_rotQuat.x()) + " " + QString::number(_rotQuat.y()) + " ";
		s += QString::number(_rotQuat.z()) + ", ";
		s += QString::number(_translation.x()) + " " + QString::number(_translation.y()) + " ";
		s += QString::number(_translation.z()) + "\n";
		// timestep
		s += QString::number(_timer.elapsed());
		s += "; timestep; ";
		s += QString::number(_timestep) + "\n";
		logInteraction(s);
	}
}


/**
 * @brief VolumeRenderWidget::setTimeStep
 * @param timestep
 */
void VolumeRenderWidget::setTimeStep(int timestep)
{
    _timestep = timestep;
    update();

	if (_logInteraction)
	{
		QString s;
		s += QString::number(_timer.elapsed());
		s += "; timestep; ";
		s += QString::number(_timestep) + "\n";

		logInteraction(s);
	}
}

/**
 * @brief VolumeRenderWidget::setImageSamplingRate
 * @param samplingRate
 */
void VolumeRenderWidget::setImageSamplingRate(const double samplingRate)
{
    _imgSamplingRate = samplingRate;
    this->resizeGL(this->width(), this->height());
}

void VolumeRenderWidget::setRenderingMethod(int rm)
{
	this->_renderingMethod = static_cast<VolumeRenderWidget::RenderingMethod>(rm);
	_volumerender.updateRenderingParameters(_renderingMethod);
	this->resizeGL(this->size().width(), this->size().height()); // calling resize to create new textures
}

#ifdef _WIN32
void VolumeRenderWidget::setEyetracking(bool eyetracking)
{
	if (check_eyetracker_availability()) {
		if (eyetracking) {
			// start using eyetracking: subscribe to data
			TobiiResearchStatus status = tobii_research_subscribe_to_gaze_data(_eyetracker, &VolumeRenderWidget::gaze_data_callback, &_gaze_data);
			if (status != TOBII_RESEARCH_STATUS_OK) {
				qCritical() << "Something went wrong when trying to subscribe to eyetracker data.\n";
			}
		}
		else {
			// stop using eyetracking: unsubscribe to data
			TobiiResearchStatus status = tobii_research_unsubscribe_from_gaze_data(_eyetracker, &VolumeRenderWidget::gaze_data_callback);
			if (status != TOBII_RESEARCH_STATUS_OK) {
				qCritical() << "Something went wrong when trying to unsubscribe to eyetracker data.\n";
			}

		}
		_useEyetracking = eyetracking;
	}
	else {
		_useEyetracking = false;
		QWidget* ppW = parentWidget()->parentWidget();
		QCheckBox* eyeCB = ppW->findChild<QCheckBox*>("chbEyetracking");
		eyeCB->setChecked(false);
	}
	update();
	// std::cout << "use eyetracking: " << _useEyetracking << std::endl;
}

/*
Shows the available eyetracking devices in a drop down menu and lets the user select one of them.
*/
void VolumeRenderWidget::showSelectEyetrackingDevice()
{
	TobiiResearchEyeTrackers* eyetrackers = NULL;
	std::vector<std::string> eyetracker_device_names;

	TobiiResearchStatus result;
	size_t i = 0;
	result = tobii_research_find_all_eyetrackers(&eyetrackers);
	if (result != TOBII_RESEARCH_STATUS_OK) {
		qCritical() << "Finding trackers failed. Error: " << result << "\n";
		return;
	}
	else {
		if (eyetrackers->count == 0) {
			qDebug() << "No Eyetrackers found!\n";
		}
	}

	for (i = 0; i < eyetrackers->count; i++) {
		TobiiResearchEyeTracker* eyetracker = eyetrackers->eyetrackers[i];
		char* device_name;
		char* serial_number;
		tobii_research_get_device_name(eyetracker, &device_name);
		tobii_research_get_serial_number(eyetracker, &serial_number);
		eyetracker_device_names.push_back(std::string(device_name).append(", serial number: ").append(std::string(serial_number)));
		tobii_research_free_string(device_name);
		tobii_research_free_string(serial_number);
	}

	QStringList platforms;
	bool ok = false;
	bool only_one = false;
	QString platform;
	int eyetracker_index = 0;

	for (std::string &s : eyetracker_device_names)
		platforms.append(QString::fromStdString(s));

	if (eyetracker_device_names.size() > 1) {
		platform = QInputDialog::getItem(this, tr("Select Eyetracker"),
			tr("Select Eyetracking Device:"),
			platforms, 0, false, &ok);
		eyetracker_index = platforms.indexOf(platform);
	}
	else {
		only_one = true;
		eyetracker_index = 0;
	}
	if (ok && !platform.isEmpty() || only_one && eyetrackers->count > 0)
	{
		_eyetracker = eyetrackers->eyetrackers[eyetracker_index];
	}

	tobii_research_free_eyetrackers(eyetrackers);
}

/*
Shows the available monitors in a drop down menu and lets the user select one of them.
The monitor to be selected should be the monitor on which the current eyetracking device is calibrated to.
*/
void VolumeRenderWidget::actionSelectMonitor()
{
	DISPLAY_DEVICE dd;
	dd.cb = sizeof(dd);

	std::vector<std::string> device_names;
	std::vector<std::string> device_strings;

	for (int deviceIndex = 0; EnumDisplayDevices(0, deviceIndex, &dd, 0); deviceIndex++) {
		std::string dn = dd.DeviceName;
		for (int monitorIndex = 0; EnumDisplayDevices(dn.c_str(), monitorIndex, &dd, 0); monitorIndex++) {
			device_names.push_back(dd.DeviceName);
			device_strings.push_back(dd.DeviceString);
		}
	}

	std::string result;

	for (int i = 0; i < device_names.size(); i++) {
		result.append("Device ").append(std::to_string(i)).append(": ").append(device_names[i]).append(", ").append(device_strings[i]).append("\n");

		// only need the substring of th device name to be able to compare it with the name returned by GetMonitorInfo
		device_names[i] = device_names[i].substr(0, 12);
	}

	if (device_names.size() != device_strings.size()) {
		qCritical() << "Something went wrong! device_names and device_strings have different sizes!\n";
		return;
	}

	QStringList platforms;
	for (int i = 0; i < device_names.size(); i++)
		platforms.append(QString::fromStdString(device_names[i]).append(QString::fromStdString(", ")).append(QString::fromStdString(device_strings[i])));

	bool ok = false;
	bool only_one = false;
	QString platform;
	int monitor_index;

	if (device_names.size() > 1) {
		platform = QInputDialog::getItem(this, tr("Select Monitor"),
			tr("Select the Monitor the Eyetracker is calibrated to:"),
			platforms, 0, false, &ok, Qt::MSWindowsFixedSizeDialogHint);
		monitor_index = platforms.indexOf(platform);
	}
	else {
		only_one = true;
		monitor_index = 0;
	}

	if (ok && !platform.isEmpty() || only_one)
	{
		struct return_data_struct {
			bool success;
			long left;
			long top;
			long right;
			long bottom;
			std::string device_name;
		} return_data;

		return_data.success = false;
		return_data.device_name = device_names[monitor_index];

		// enumerate display monitors to retrive handles and information
//		EnumDisplayMonitors(NULL, NULL, reinterpret_cast<MONITORENUMPROC>(&VolumeRenderWidget::MonitorEnumProc), reinterpret_cast<LPARAM>(&return_data));

		if (return_data.success) {
			_monitor_offset = QPoint(return_data.left, return_data.top);
			_curr_monitor_width = return_data.right - return_data.left;
			_curr_monitor_height = return_data.bottom - return_data.top;
		}
		else {
			qCritical() << "Could not retrive data for selected Monitor!\n";
		}
	}
}

/*bool VolumeRenderWidget::MonitorEnumProc(HMONITOR monitor, HDC hdcMnitor, LPRECT rect, LPARAM param)
{
	struct return_data_struct {
		bool success;
		long left;
		long top;
		long right;
		long bottom;
		std::string device_name;
	};

	return_data_struct* return_data_ptr = reinterpret_cast<return_data_struct*>(param);

	MONITORINFOEX qmi;
	qmi.cbSize = sizeof(MONITORINFOEX);

	bool success = GetMonitorInfo(monitor, &qmi);

	if (success) {
		if (std::string(qmi.szDevice).compare(return_data_ptr->device_name) == 0) {
			return_data_ptr->success = true;
			return_data_ptr->bottom = qmi.rcMonitor.bottom;
			return_data_ptr->left = qmi.rcMonitor.left;
			return_data_ptr->right = qmi.rcMonitor.right;
			return_data_ptr->top = qmi.rcMonitor.top;
		}
	}

	return success;
}*/

bool VolumeRenderWidget::check_eyetracker_availability()
{
	if (_eyetracker == nullptr) {
		return false;
	}
	else {
		TobiiResearchEyeTrackers* eyetrackers = NULL;

		TobiiResearchStatus result;
		result = tobii_research_find_all_eyetrackers(&eyetrackers);
		if (result != TOBII_RESEARCH_STATUS_OK) {
			qCritical() << "Finding trackers to check status failed. Error: " << result << "\n";
			return false;
		}

		bool eyetracker_exists = false;
		for (int i = 0; i < eyetrackers->count; i++) {
			if (_eyetracker == eyetrackers->eyetrackers[i]) {
				eyetracker_exists = true;
			}
		}
		tobii_research_free_eyetrackers(eyetrackers);
		return eyetracker_exists;
	}

}
#endif

void VolumeRenderWidget::gaze_data_callback(TobiiResearchGazeData * gaze_data, void * user_data)
{
	memcpy(user_data, gaze_data, sizeof(*gaze_data));
}

/**
 * @brief VolumeRenderWidget::paintGL
 */
void VolumeRenderWidget::paintGL()
{
	switch (_renderingMethod) {
	case LBG_Sampling:
		cl_float2 lcpf;
        if (_bench.active)
        {
            if (_bench.active && _bench.isCameraIteration())
                updateView();
            else
            {
                lcpf.x = static_cast<float>(std::generate_canonical<float, std::numeric_limits<float>::digits>(_prng));
                lcpf.y = static_cast<float>(std::generate_canonical<float, std::numeric_limits<float>::digits>(_prng));
            }
            _bench.iteration++;
        }
        else
        {
            lcpf.x = static_cast<cl_float>(_lastLocalCursorPos.x() / static_cast<cl_float>(this->size().width()));
            lcpf.y = static_cast<cl_float>(_lastLocalCursorPos.y() / static_cast<cl_float>(this->size().height()));
        }
		_volumerender.setGazePoint(lcpf);
		paintGL_LBG_sampling();

        if (_bench.active)
        {
            QString outString;
            QTextStream out(&outString);
            out << _bench.iteration << "; ";
            out << _rotQuat.toVector4D().w() << " " << _rotQuat.x() << " " << _rotQuat.y() << " "
                       << _rotQuat.z() << "; ";
            out << _translation.x() << " " << _translation.y() << " " << _translation.z() << "; ";
            out << lcpf.x << " " << lcpf.y << "; ";
            out << _volumerender.getLastExecTime() << "\n";
            _bench.writeState(out.readAll());
        }
		break;
	default:
		paintGL_standard();
		break;
	}
    if (_contRendering)
        update();
}


void VolumeRenderWidget::paintGL_standard() {
	double fps = 0.0;
	if (this->_loadingFinished && _volumerender.hasData() && !_noUpdate)
	{
		// OpenCL raycast
		try
		{
			if (_useGL)
				_volumerender.runRaycast(floor(this->size().width() * _imgSamplingRate),
					floor(this->size().height()* _imgSamplingRate), _timestep);
			else
			{
				std::vector<float> d;
				_volumerender.runRaycastNoGL(floor(this->size().width() * _imgSamplingRate),
					floor(this->size().height()* _imgSamplingRate),
					_timestep, d);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F,
					floor(this->size().width() * _imgSamplingRate),
					floor(this->size().height()* _imgSamplingRate),
					0, GL_RGBA, GL_FLOAT,
					d.data());
				glGenerateMipmap(GL_TEXTURE_2D);
				_volumerender.updateOutputImg(static_cast<size_t>(width()),
					static_cast<size_t>(height()), _outTexId);
			}
		}
		catch (std::runtime_error e)
		{
			qCritical() << e.what();
		}
		fps = getFps();
	}

	QPainter p(this);
	p.beginNativePainting();
	{
		// render the ray casting output
		// clear to white to avoid getting colored borders outside the quad
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// draw screen quad
		//
		_screenQuadVao.bind();
		_quadVbo.bind();
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, nullptr);

		// render screen quad
		//
		_spScreenQuad.bind();
		_spScreenQuad.setUniformValue(_spScreenQuad.uniformLocation("projMatrix"),
			_screenQuadProjMX);
		_spScreenQuad.setUniformValue(_spScreenQuad.uniformLocation("mvMatrix"),
			_viewMX * _modelMX);

		_spScreenQuad.setUniformValue(_spScreenQuad.uniformLocation("outTex"), GL_TEXTURE0);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
		_screenQuadVao.release();
		_quadVbo.release();
		_spScreenQuad.release();

		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		if (_volumerender.hasData() && _writeImage)
		{
			QImage img = this->grabFramebuffer();
			QString number = QString("%1").arg(_imgCount++, 6, 10, QChar('0'));
			if (!_recordVideo)
			{
				QLoggingCategory category("screenshot");
				qCInfo(category, "Writing current frame img/frame_%s.png",
					number.toStdString().c_str());
				_writeImage = false;
			}
			if (!QDir("img").exists())
				QDir().mkdir("img");
			img.save("img/frame_" + number + "_"
				+ QString::number(_volumerender.getLastExecTime()) + ".png");
		}
	}
	p.endNativePainting();

	// render overlays
	if (_showOverlay)
	{
        paintFps(p, fps, _volumerender.getLastExecTime());
		paintOrientationAxis(p);
	}

	// recover opengl texture
	p.beginNativePainting();
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, _outTexId);
	}
	p.endNativePainting();
	p.end();
}

void VolumeRenderWidget::paintGL_LBG_sampling() {
	double fps = 0.0;
	if (this->_loadingFinished && _volumerender.hasData() && !_noUpdate)
	{
		// OpenCL raycast
		try
		{
            if (_useGL)
            {
				// set first Texture to extends of index map
                _volumerender.updateOutputTex(_tmpTexId);

                _volumerender.runRaycastLBG(floor(this->size().width() * _imgSamplingRate),
                                            floor(this->size().height()* _imgSamplingRate),
                                            _timestep);

				fps = _volumerender.getLastExecTime();

				// second texture needs to have one third in each dimension of the index map
				//_volumerender.updateOutputImg(floor(_volumerender.getIndexMapExtends().x() / 2.0), floor(_volumerender.getIndexMapExtends().y() / 2.0),
				//	_outTexId);

                _volumerender.interpolateLBG(floor(_volumerender.getIndexMapExtends().x() / 2.0),
                                             floor(_volumerender.getIndexMapExtends().y() / 2.0),
                                             _tmpTexId, _outTexId);
				

				/*generateOutputTextures(floor(_volumerender.getIndexMapExtends().x()), floor(_volumerender.getIndexMapExtends().y()),
					&_tmpTexId, GL_TEXTURE1);

				_volumerender.runRaycastLBG(_timestep);

				fps = _volumerender.getLastExecTime();

				// second texture needs to have one third in each dimension of the index map
				//_volumerender.updateOutputImg(floor(_volumerender.getIndexMapExtends().x() / 2.0), floor(_volumerender.getIndexMapExtends().y() / 2.0),
				//	_outTexId);

				_volumerender.interpolateLBG(floor(_volumerender.getIndexMapExtends().x() / 2.0), floor(_volumerender.getIndexMapExtends().y() / 2.0), _tmpTexId, _outTexId);
				*/
			}
			else
			{
				std::vector<float> d;
				_volumerender.runRaycastLBGNoGL(floor(this->size().width() * _imgSamplingRate),
					floor(this->size().height()* _imgSamplingRate),
					_timestep, d);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F,
					floor(this->size().width() * _imgSamplingRate),
					floor(this->size().height()* _imgSamplingRate),
					0, GL_RGBA, GL_FLOAT,
					d.data());
                //glGenerateMipmap(GL_TEXTURE_2D);
				_volumerender.updateOutputImg(static_cast<size_t>(width()),
					static_cast<size_t>(height()), _outTexId);
			}
		}
		catch (std::runtime_error e)
		{
			qCritical() << e.what();
		}
        // offset for render pass + interpolation
        fps = getFps(fps);
	}

	QPainter p(this);
	p.beginNativePainting();
	{
		// render the ray casting output
		// clear to white to avoid getting colored borders outside the quad
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// draw screen quad
		//
		_screenQuadVao.bind();
		_quadVbo.bind();
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, nullptr);

		// render screen quad
		//
		_spScreenQuad.bind();
		_spScreenQuad.setUniformValue(_spScreenQuad.uniformLocation("projMatrix"),
			_screenQuadProjMX);
		_spScreenQuad.setUniformValue(_spScreenQuad.uniformLocation("mvMatrix"),
			_viewMX * _modelMX);

		_spScreenQuad.setUniformValue(_spScreenQuad.uniformLocation("outTex"), GL_TEXTURE0);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
		_screenQuadVao.release();
		_quadVbo.release();
		_spScreenQuad.release();

		glDisable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		if (_volumerender.hasData() && _writeImage)
		{
			QImage img = this->grabFramebuffer();
			QString number = QString("%1").arg(_imgCount++, 6, 10, QChar('0'));
			if (!_recordVideo)
			{
				QLoggingCategory category("screenshot");
				qCInfo(category, "Writing current frame img/frame_%s.png",
					number.toStdString().c_str());
				_writeImage = false;
			}
			if (!QDir("img").exists())
				QDir().mkdir("img");
			img.save("img/frame_" + number + "_"
				+ QString::number(_volumerender.getLastExecTime()) + ".png");
		}
	}
	p.endNativePainting();

	// render overlays
	if (_showOverlay)
	{
        paintFps(p, fps, _volumerender.getLastExecTime());
		paintOrientationAxis(p);
	}

	// recover opengl texture
	p.beginNativePainting();
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, _outTexId);
	}
	p.endNativePainting();
	p.end();
}

/**
 * @brief VolumeRenderWidget::resizeGL
 * @param w
 * @param h
 */
void VolumeRenderWidget::resizeGL(const int w, const int h)
{
    _screenQuadProjMX.setToIdentity();
    _screenQuadProjMX.perspective(53.14f, 1.0f, Z_NEAR, Z_FAR);

    _overlayProjMX.setToIdentity();
    _overlayProjMX.perspective(53.14f, qreal(w)/qreal(h ? h : 1), Z_NEAR, Z_FAR);

    try
    {
		//generateOutputTextures(floor(_volumerender.getIndexMapExtends().x()), floor(_volumerender.getIndexMapExtends().y()),
		//	&_outTexId, GL_TEXTURE0);
        generateOutputTextures(floor(_volumerender.getIndexMapExtends().x()),
                               floor(_volumerender.getIndexMapExtends().y()),
                               &_tmpTexId, GL_TEXTURE1);
//        if(_renderingMethod == LBG_Sampling)
//            generateOutputTextures(floor(_volumerender.getIndexMapExtends().x() / 2.0),
//                                   floor(_volumerender.getIndexMapExtends().y() / 2.0),
//                                   &_outTexId, GL_TEXTURE0);
//		else
            generateOutputTextures(floor(w*_imgSamplingRate), floor(h*_imgSamplingRate),
                                   &_outTexId, GL_TEXTURE0);

		/*
		// Debug
		std::cout << "Resize _tmpImgId: (" << floor(_volumerender.getIndexMapExtends().x()) << ", " << floor(_volumerender.getIndexMapExtends().y()) << ")" << std::endl;
		std::cout << "size / 2.0: (" << floor(_volumerender.getIndexMapExtends().x()) / 2.0 << ", " << floor(_volumerender.getIndexMapExtends().y()) / 2.0 << ")\n" << std::endl;
		*/
	}
    catch (std::runtime_error e)
    {
        qCritical() << "An error occured while generating output texture." << e.what();
    }

    emit frameSizeChanged(this->size());
}


/**
 * @brief VolumeRenderWidget::generateOutputTextures
 */
void VolumeRenderWidget::generateOutputTextures(const int width, const int height,
                                                GLuint* texture, GLuint tex_unit)
{
	glDeleteTextures(1, texture);
	glGenTextures(1, texture);

    glActiveTexture(tex_unit);
    glBindTexture(GL_TEXTURE_2D, *texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    if (!_volumerender.hasData())
    {
        QImage img(width, height, QImage::Format_RGBA8888);
        img.fill(Qt::white);
        QPainter p(&img);
        p.setFont(QFont("Helvetia", 12));
        p.drawText(width/2 - 110, height/2, "Drop your volume data file here.");
        p.end();
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,
                     width, height, 0,
                     GL_RGBA, GL_UNSIGNED_BYTE,
                     img.bits());
    }
    else
    {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,
                     width, height, 0,
                     GL_RGBA, GL_UNSIGNED_BYTE,
                     nullptr);
    }
    glGenerateMipmap(GL_TEXTURE_2D);

    _volumerender.updateOutputImg(static_cast<size_t>(width), static_cast<size_t>(height),
		*texture);

    updateView(0, 0);
}

void VolumeRenderWidget::setShowOverlay(bool showOverlay)
{
    _showOverlay = showOverlay;
    updateView();
}

QQuaternion VolumeRenderWidget::getCamRotation() const
{
    return _rotQuat;
}

void VolumeRenderWidget::setCamRotation(const QQuaternion &rotQuat)
{
    _rotQuat = rotQuat;
}

QVector3D VolumeRenderWidget::getCamTranslation() const
{
    return _translation;
}

void VolumeRenderWidget::setCamTranslation(const QVector3D &translation)
{
    _translation = translation;
}

/**
 * @brief VolumeRenderWidget::showSelectOpenCL
 */
void VolumeRenderWidget::showSelectOpenCL()
{
    std::vector<std::string> names;
    try
    {
        names = _volumerender.getPlatformNames();
    }
    catch (std::runtime_error e)
    {
        qCritical() << "An error occured while trying to retrieve OpenCL platform information."
                    << e.what();
        return;
    }

    QStringList platforms;
    for (std::string &s : names)
        platforms.append(QString::fromStdString(s));

    bool ok = false;
    QString platform = QInputDialog::getItem(this, tr("Select platform"),
                                             tr("Select OpenCL platform:"),
                                             platforms, 0, false, &ok);
    if (ok && !platform.isEmpty())
    {
        cl_vendor vendor = VENDOR_ANY;
        QString type = "GPU";
        _useGL = false;

        // filter GPU only platforms
        if (platform.contains("NVIDIA"))
            vendor = VENDOR_NVIDIA;
        else
        {
            if (!platform.contains("Graphics"))
            {
                QStringList types;
                types << "GPU" << "CPU";
                type = QInputDialog::getItem(this, tr("Select type"),
                                             tr("Select device type:"), types, 0, false, &ok);
            }
            if (platform.contains("Advanced Micro Devices", Qt::CaseInsensitive))
                vendor = VENDOR_AMD;
            else if (platform.contains("Intel", Qt::CaseInsensitive))
                vendor = VENDOR_INTEL;

        }
        if (!type.isEmpty())
        {
            if (type == "GPU")
            {
                QMessageBox msgBox;
                msgBox.setText("Do you wish to try OpenGL context sharing using this platform?");
                msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
                msgBox.setDefaultButton(QMessageBox::Yes);
                int ret = msgBox.exec();
                _useGL = ret == QMessageBox::Yes;
            }

            size_t platformId = static_cast<size_t>(platforms.indexOf(platform));
            try
            {
                names = _volumerender.getDeviceNames(platformId, type.toStdString());
            }
            catch (std::runtime_error e)
            {
                qCritical() << "No capable device found on using the selected platform and type."
                            << e.what();
                return;
            }
            QStringList devices;
            for (std::string &s : names)
                devices.append(QString::fromStdString(s));

            QString device;
            if (devices.empty())
                throw std::runtime_error("No capable device found on the selected platform.");
            else if (devices.size() == 1)
                device = devices.front();
            else
                device = QInputDialog::getItem(this, tr("Select device"),
                                                   tr("Select OpenCL device:"),
                                                   devices, 0, false, &ok);
            if (!device.isEmpty())
            {
                try {
                    _volumerender.initialize(_useGL, type == "CPU", vendor, device.toStdString(),
                                             static_cast<int>(platformId));
                } catch (std::runtime_error e) {
                    qCritical() << e.what() << "\nSwitching to CPU fallback mode.";
                    _useGL = false;
                    _volumerender.initialize(false, true);
                }
                catch (...) {
                    qCritical() << "An unknown error occured initializing OpenCL/OpenGL.";
                }
                updateTransferFunction(_tffStops);
                this->resizeGL(this->width(), this->height());
            }
        }
    }
}


/**
 * @brief VolumeRenderWidget::setupVertexAttribs
 */
void VolumeRenderWidget::setupVertexAttribs()
{
    // screen quad
    _quadVbo.bind();
    glEnableVertexAttribArray( 0 );
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, nullptr);
    _quadVbo.release();
}


/**
 * @brief VolumeRenderWidget::setVolumeData
 * @param fileName
 */
void VolumeRenderWidget::setVolumeData(const QString &fileName)
{
    this->_noUpdate = true;
    size_t timesteps = 0;
    try
    {
        timesteps = _volumerender.loadVolumeData(fileName.toStdString());
    }
    catch (std::invalid_argument e)
    {
        qCritical() << e.what();
    }
    catch (std::runtime_error e)
    {
        qCritical() << e.what();
    }
    if (timesteps > 1)
        emit timeSeriesLoaded(static_cast<int>(timesteps - 1));

    _overlayModelMX.setToIdentity();
    QVector3D res = getVolumeResolution().toVector3D();
    _overlayModelMX.scale(res / qMax(res.x(), qMax(res.y(), res.z())));
    this->_noUpdate = false;
    update();
}

void VolumeRenderWidget::setIndexandSamplingMap(const QString &fileNameIndexMap,
                                                const QString &fileNameSamplingMap,
                                                const QString &fileNameNeighborIndex,
                                                const QString &fileNameNeighborWeights)
{
	this->_noUpdate = true;
	// qDebug() << fileNameIndexMap << " " << fileNameSamplingMap;

	try
	{
        _volumerender.loadIndexAndSamplingMap(fileNameIndexMap.toStdString(),
                                              fileNameSamplingMap.toStdString(),
                                              fileNameNeighborIndex,
                                              fileNameNeighborWeights);
	}
	catch (std::invalid_argument e)
	{
		qCritical() << e.what();
	}
	catch (std::runtime_error e)
	{
		qCritical() << e.what();
	}

	this->_noUpdate = false;
	update();
}


/**
 * @brief VolumeRenderWidget::hasData
 * @return
 */
bool VolumeRenderWidget::hasData() const 
{
    return _volumerender.hasData();
}


/**
 * @brief VolumeRenderWidget::getVolumeResolution
 * @return
 */
const QVector4D VolumeRenderWidget::getVolumeResolution() const
{
    if (_volumerender.hasData() == false)
        return QVector4D();

    return QVector4D(_volumerender.getResolution().at(0),
                     _volumerender.getResolution().at(1),
                     _volumerender.getResolution().at(2),
                     _volumerender.getResolution().at(3));
}


/**
 * @brief VolumeRenderWidget::updateStepSize
 * @param stepSize
 */
void VolumeRenderWidget::updateSamplingRate(double samplingRate)
{
    _volumerender.updateSamplingRate(samplingRate);
    update();
}


/**
 * @brief VolumeRenderWidget::setInterpolation
 * @param method
 */
void VolumeRenderWidget::setTffInterpolation(const QString method)
{
    if (method.contains("Quad"))
        _tffInterpol = QEasingCurve::InOutQuad;
    else if (method.contains("Linear"))
        _tffInterpol = QEasingCurve::Linear;

	if (_logInteraction)
	{
		QString s;
		s += QString::number(_timer.elapsed());
		s += "; tffInterpolation; ";
		s += _tffInterpol == QEasingCurve::InOutQuad ? "quad" : "linear";
		s += "\n";

		logInteraction(s);
	}
}

/**
 * @brief VolumeRenderWidget::setRawTransferFunction
 * @param tff
 */
void VolumeRenderWidget::setRawTransferFunction(std::vector<unsigned char> tff)
{
    try
    {
        _volumerender.setTransferFunction(tff);
        update();
    }
    catch (std::runtime_error e)
    {
        qCritical() << e.what();
    }
}

/**
 * @brief VolumeRenderWidget::updateTransferFunction
 * @param stops
 */
void VolumeRenderWidget::updateTransferFunction(QGradientStops stops)
{
    const int tffSize = 256;
    const qreal granularity = 4096.0;
    std::vector<uchar> tff(tffSize*4);
    std::vector<unsigned int> prefixSum(tffSize);

    QPropertyAnimation interpolator;
    interpolator.setEasingCurve(_tffInterpol);
    interpolator.setDuration(static_cast<int>(granularity));
    foreach (QGradientStop stop, stops)
    {
        interpolator.setKeyValueAt(stop.first, stop.second);
    }
//    tff.at(0) = (uchar)0;
//    tff.at(1) = (uchar)0;
//    tff.at(2) = (uchar)0;
//    tff.at(3) = (uchar)0;
#pragma omp for
    for (int i = 0; i < tffSize; ++i)
    {
        interpolator.setCurrentTime((i/static_cast<double>(tffSize)) * granularity);
        tff.at(i*4 + 0) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().red()   - 3);
        tff.at(i*4 + 1) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().green() - 3);
        tff.at(i*4 + 2) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().blue()  - 3);
        tff.at(i*4 + 3) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().alpha() - 3);
        prefixSum.at(i) = tff.at(i*4 + 3);
    }
    try
    {
        _volumerender.setTransferFunction(tff);
        // TODO: replace with std::exclusicve_scan(std::par, ... )
        std::partial_sum(prefixSum.begin(), prefixSum.end(), prefixSum.begin());
        _volumerender.setTffPrefixSum(prefixSum);
    }
    catch (std::runtime_error e)
    {
        qCritical() << e.what();
    }
    update();
    _tffStops = stops;

	if (_logInteraction)
	{
		QString s;
		s += QString::number(_timer.elapsed());
		s += "; transferFunction; ";
		const std::vector<unsigned char> tff = getRawTransferFunction(stops);
		foreach (unsigned char c, tff)
		{
			s += QString::number(static_cast<int>(c)) + " ";
		}
		s += "\n";
	}
}

std::vector<unsigned char> VolumeRenderWidget::getRawTransferFunction(QGradientStops stops) const
{
    const size_t tffSize = 256;
    const qreal granularity = 4096.0;
    std::vector<uchar> tff(tffSize*4);

    QPropertyAnimation interpolator;
    interpolator.setEasingCurve(_tffInterpol);
    interpolator.setDuration(static_cast<int>(granularity));
    foreach (QGradientStop stop, stops)
    {
        interpolator.setKeyValueAt(stop.first, stop.second);
    }

    for (size_t i = 0; i < tffSize; ++i)
    {
        interpolator.setCurrentTime((i/static_cast<double>(tffSize)) * granularity);
        tff.at(i*4 + 0) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().red()   - 3);
        tff.at(i*4 + 1) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().green() - 3);
        tff.at(i*4 + 2) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().blue()  - 3);
        tff.at(i*4 + 3) = (uchar)qMax(0, interpolator.currentValue().value<QColor>().alpha() - 3);
    }
    return tff;
}


/**
 * @brief VolumeRenderWidget::cleanup
 */
void VolumeRenderWidget::cleanup()
{
//    makeCurrent();
//    if (_quadVbo.isCreated())
//        _quadVbo.destroy();
}



/**
 * @brief VolumeRenderWidget::mousePressEvent
 * @param event
 */
void VolumeRenderWidget::mousePressEvent(QMouseEvent *event)
{
    // nothing yet
    _lastLocalCursorPos = event->pos();
}


/**
 * @brief VolumeRenderWidget::mouseReleaseEvent
 * @param event
 */
void VolumeRenderWidget::mouseReleaseEvent(QMouseEvent *event)
{
    event->accept();
}


/**
 * @brief VolumeRenderWidget::recordView
 */
void VolumeRenderWidget::recordViewConfig() const
{
    QFile saveQuat(_viewLogFile + "_quat.txt");
    QFile saveTrans(_viewLogFile + "_trans.txt");
    if (!saveQuat.open(QFile::WriteOnly | QFile::Append)
            || !saveTrans.open(QFile::WriteOnly | QFile::Append))
    {
        qWarning() << "Couldn't open file for saving the camera configurations: "
                   << _viewLogFile;
        return;
    }

    QTextStream quatStream(&saveQuat);
    quatStream << _rotQuat.toVector4D().w() << " " << _rotQuat.x() << " " << _rotQuat.y() << " "
               << _rotQuat.z() << "; ";
    QTextStream transStream(&saveTrans);
    transStream << _translation.x() << " " << _translation.y() << " " << _translation.z() << "; ";
}


/**
 *
 */
void VolumeRenderWidget::logInteraction(const QString &str) const
{
	QFile f(_interactionLogFile);
	if (!f.open(QFile::WriteOnly | QFile::Append))
	{
		qWarning() << "Couldn't open file for saving the camera configurations: "
				   << _interactionLogFile;
		return;
	}
	QTextStream fileStream(&f);
	fileStream << str;
}

/**
 * @brief VolumeRenderWidget::resetCam
 */
void VolumeRenderWidget::resetCam()
{
    _rotQuat = QQuaternion(1, 0, 0, 0);
    _translation = QVector3D(0, 0, 2.0);
    updateView();
}

/**
 * @brief update camera view
 * @param dx
 * @param dy
 */
void VolumeRenderWidget::updateView(float dx, float dy)
{
    if (_bench.active && _bench.isCameraIteration())
    {
        _bench.iteration++;
        dx = static_cast<float>(std::generate_canonical<float, std::numeric_limits<float>::digits>(_prng));
        dy = static_cast<float>(std::generate_canonical<float, std::numeric_limits<float>::digits>(_prng));
        _translation.setZ(static_cast<float>(std::generate_canonical<float, std::numeric_limits<float>::digits>(_prng)) * 4);
    }

    QVector3D rotAxis = QVector3D(dy, dx, 0.0f).normalized();
    float angle = QVector2D(dx, dy).length()*360.f;
    _rotQuat = _rotQuat * QQuaternion::fromAxisAndAngle(rotAxis, -angle);

    QMatrix4x4 viewMat;
//    viewMat.scale(QVector3D(1,1,1./_zScale));
    viewMat.rotate(_rotQuat);

    _coordViewMX.setToIdentity();
    _coordViewMX.scale(1, -1, 1);
    _coordViewMX.translate(_translation * -1.0);
    _coordViewMX *= QMatrix4x4(_rotQuat.toRotationMatrix().transposed());

    viewMat.translate(_translation);
    viewMat.scale(_translation.z());

    std::array<float, 16> viewArray;
    for (size_t i = 0; i < viewArray.size(); ++i)
        viewArray.at(i) = viewMat.transposed().constData()[i];

    try
    {
        _volumerender.updateView(viewArray);
    }
    catch (std::runtime_error e)
    {
        qCritical() << e.what();
    }
    update();

    if (_logView)
        recordViewConfig();
	if (_logInteraction)
	{
		QString s;
		s += QString::number(_timer.elapsed());
		s += "; camera; ";
		s += QString::number(_rotQuat.toVector4D().w()) + " ";
		s += QString::number(_rotQuat.x()) + " " + QString::number(_rotQuat.y()) + " ";
		s += QString::number(_rotQuat.z()) + ", ";

		s += QString::number(_translation.x()) + " " + QString::number(_translation.y()) + " ";
		s += QString::number(_translation.z()) + "\n";

		logInteraction(s);
	}
}


/**
 * @brief VolumeRenderWidget::mouseMoveEvent
 * @param event
 */
void VolumeRenderWidget::mouseMoveEvent(QMouseEvent *event)
{
    float dx = (float)(event->pos().x() - _lastLocalCursorPos.x()) / width();
    float dy = (float)(event->pos().y() - _lastLocalCursorPos.y()) / height();

    // rotate object
    if (event->buttons() & Qt::LeftButton)
    {
        if (event->modifiers() & Qt::ShiftModifier)
        {
            dx *= 0.1f;
            dy *= 0.1f;
        }
        updateView(dx, dy);
    }
    // translate object
    if (event->buttons() & Qt::MiddleButton)
    {
        int sensitivity = 6;
        if (event->modifiers() & Qt::ShiftModifier)
            sensitivity = 1;

        _translation.setX(_translation.x() - dx*sensitivity);
        _translation.setY(_translation.y() + dy*sensitivity);
        updateView();
    }

    _lastLocalCursorPos = event->pos();
    event->accept();
}


/**
 * @brief VolumeRenderWidget::wheelEvent
 * @param event
 */
void VolumeRenderWidget::wheelEvent(QWheelEvent *event)
{
    double t = 1600.0;
    if (event->modifiers() & Qt::ShiftModifier)
        t *= 6.0;

    // limit translation to origin, otherwise camera setup breaks (flips)
    _translation.setZ(qMax(0.01, _translation.z() - event->angleDelta().y() / t));
    updateView();
    event->accept();
}


/**
 * @brief VolumeRenderWidget::mouseDoubleClickEvent
 * @param event
 */
void VolumeRenderWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
    // nothing yet
    event->accept();
}



/**
 * @brief VolumeRenderWidget::keyReleaseEvent
 * @param event
 */
void VolumeRenderWidget::keyReleaseEvent(QKeyEvent *event)
{
    // nothing yet
    event->accept();
}


/**
 * @brief VolumeRenderWidget::getLoadingFinished
 * @return
 */
bool VolumeRenderWidget::getLoadingFinished() const
{
    return _loadingFinished;
}


/**
 * @brief VolumeRenderWidget::setLoadingFinished
 * @param loadingFinished
 */
void VolumeRenderWidget::setLoadingFinished(const bool loadingFinished)
{
    this->setTimeStep(0);
    _loadingFinished = loadingFinished;
}


/**
 * @brief VolumeRenderWidget::setCamOrtho
 * @param camOrtho
 */
void VolumeRenderWidget::setCamOrtho(const bool camOrtho)
{
    _volumerender.setCamOrtho(camOrtho);
    _overlayProjMX.setToIdentity();
    if (camOrtho)
        _overlayProjMX.ortho(QRect(0, 0, width(), height()));
    else
        _overlayProjMX.perspective(53.14f, qreal(width())/qreal(height() ? height() : 1), Z_NEAR, Z_FAR);
    this->updateView();
}

void VolumeRenderWidget::setContRendering(const bool contRendering)
{
    _contRendering = contRendering;
    this->updateView();
}

/**
 * @brief VolumeRenderWidget::setIllumination
 * @param illum
 */
void VolumeRenderWidget::setIllumination(const int illum)
{
    _volumerender.setIllumination(static_cast<unsigned int>(illum));
    this->updateView();
}


/**
 * @brief VolumeRenderWidget::setAmbientOcclusion
 * @param ao
 */
void VolumeRenderWidget::setAmbientOcclusion(const bool ao)
{
    _volumerender.setAmbientOcclusion(ao);
    this->updateView();
}

/**
 * @brief VolumeRenderWidget::setLinearInterpolation
 * @param linear
 */
void VolumeRenderWidget::setLinearInterpolation(const bool linear)
{
    _volumerender.setLinearInterpolation(linear);
    this->updateView();
}

/**
 * @brief VolumeRenderWidget::setContours
 * @param contours
 */
void VolumeRenderWidget::setContours(const bool contours)
{
    _volumerender.setContours(contours);
    this->updateView();
}

/**
 * @brief VolumeRenderWidget::setAerial
 * @param aerial
 */
void VolumeRenderWidget::setAerial(const bool aerial)
{
    _volumerender.setAerial(aerial);
    this->updateView();
}

/**
 * @brief VolumeRenderWidget::setImgEss
 * @param useEss
 */
void VolumeRenderWidget::setImgEss(const bool useEss)
{
    if (useEss)
        _volumerender.updateOutputImg(static_cast<size_t>(width() * _imgSamplingRate),
                                      static_cast<size_t>(height() * _imgSamplingRate), _outTexId);
    _volumerender.setImgEss(useEss);
    this->updateView();
}

/**
 * @brief VolumeRenderWidget::setImgEss
 * @param useEss
 */
void VolumeRenderWidget::setObjEss(const bool useEss)
{
    _volumerender.setObjEss(useEss);
    this->updateView();
}

/**
 * @brief VolumeRenderWidget::setDrawBox
 * @param box
 */
void VolumeRenderWidget::setDrawBox(const bool box)
{
    _volumerender.setShowESS(box);
    this->updateView();
}


/**
 * @brief VolumeRenderWidget::setBackgroundColor
 * @param col
 */
void VolumeRenderWidget::setBackgroundColor(const QColor col)
{
    std::array<float, 4> color = {{ static_cast<float>(col.redF()),
                                    static_cast<float>(col.greenF()),
                                    static_cast<float>(col.blueF()),
                                    static_cast<float>(col.alphaF()) }};
    _volumerender.setBackground(color);
    this->updateView();
}


/**
 * @brief VolumeRenderWidget::calcFPS
 * @return
 */
double VolumeRenderWidget::getFps(double offset)
{
    _times.push_back(_volumerender.getLastExecTime() + offset);
    if (_times.length() > 64)
        _times.pop_front();

    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < _times.size(); ++i)
        sum += _times.at(i);

    double fps = 1.0/(sum/_times.length());

//    qDebug() << fps;
    return fps;
}


/**
 * @brief VolumeRenderWidget::generateLowResVolume
 * @param factor
 */
void VolumeRenderWidget::generateLowResVolume()
{
    bool ok = false;
    int factor = QInputDialog::getInt(this, tr("Factor"),
                                         tr("Select downsampling factor:"), 2, 2, 64, 1, &ok);
    if (ok)
    {
        try
        {
            std::string name = _volumerender.generateVolumeDownsampling(
                        static_cast<size_t>(_timestep),factor);
            QLoggingCategory category("volumeDownSampling");
            qCInfo(category, "Successfully created down-sampled volume data set: '%s.raw'",
                   name.c_str());
        }
        catch (std::invalid_argument e)
        {
            qCritical() << e.what();
        }
        catch (std::runtime_error e)
        {
            qCritical() << e.what();
        }
    }
}

/**
 * @brief VolumeRenderWidget::read
 * @param json
 */
void VolumeRenderWidget::read(const QJsonObject &json)
{
    if (json.contains("camRotation"))
    {
        QStringList sl = json["camRotation"].toVariant().toString().split(' ');
        if (sl.length() >= 4)
        {
            _rotQuat.setScalar(sl.at(0).toFloat());
            _rotQuat.setX(sl.at(1).toFloat());
            _rotQuat.setY(sl.at(2).toFloat());
            _rotQuat.setZ(sl.at(3).toFloat());
        }
    }
    if (json.contains("camTranslation"))
    {
        QStringList sl = json["camTranslation"].toVariant().toString().split(' ');
        if (sl.length() >= 3)
        {
            _translation.setX(sl.at(0).toFloat());
            _translation.setY(sl.at(1).toFloat());
            _translation.setZ(sl.at(2).toFloat());
        }
    }
    updateView();
}

/**
 * @brief VolumeRenderWidget::write
 * @param json
 */
void VolumeRenderWidget::write(QJsonObject &json) const
{
    QString sTmp = QString::number(_rotQuat.scalar()) + " " + QString::number(_rotQuat.x())
                   + " " + QString::number(_rotQuat.y()) + " " + QString::number(_rotQuat.z());
    json["camRotation"] = sTmp;
    sTmp = QString::number(_translation.x()) + " " + QString::number(_translation.y())
            + " " + QString::number(_translation.z());
    json["camTranslation"] = sTmp;
}

/**
 * @brief VolumeRenderWidget::reloadKernels
 */
void VolumeRenderWidget::reloadKernels()
{
    // NOTE: this reload resets all previously defined rendering settings to default values
    initVolumeRenderer();
    resizeGL(width(), height());
}

/**
 * @brief VolumeRenderWidget::toggleBenchmark
 */
void VolumeRenderWidget::toggleBenchmark()
{
    qInfo() << (_bench.active ? "Stopped benchmark run." : "Started benchmark run.");

    _bench.active = !_bench.active;
    if (_bench.active)
    {
        QFileDialog dialog;
        _bench.logFileName = dialog.getSaveFileName(this, tr("Save benchmark results"),
                                                   QDir::currentPath(), tr("All files"));
        if (_bench.logFileName.isEmpty())
        {
            _bench.active = false;
            return;
        }
        bool ok = false;
        _bench.gaze_iterations = QInputDialog::getInt(this, tr("Gaze iterations"),
                                             tr("Select gaze iterations:"), 100, 1, 10000, 1, &ok);
        _prng = QRandomGenerator64(42);
        _bench.iteration = 0;
    }
    updateView();
}
