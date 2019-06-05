# Description #

An Voronoi-based foveated volume renderer. Check out the paper for more details: https://vbruder.github.io/publication/bruder-sbfwe-19/bruder-sbfwe-19.pdf

The Weighted-Linde-Buzo-Gray algorithm can be used to generate the sampling mask in a preprocessing step.

Besides foveated rendering, the raycaster features early ray termination, object and image order empty space skipping, local illumination, and various gradient based shading techniques.
It can display volume (timeseries) data sets and uses the [Qt](https://www.qt.io) framework for the GUI. 

# Setup and build #

To compile the code you need:

* An OpenCL 1.2 capable device and drivers/libraries with image support. It is recommended to update your GPU driver before building/running.
* Qt version 5.10 or higher.
* A compiler capable of c++14.
* CMake version 3.9 or higher.

Use CMake to build the volume rasycaster:
```
git clone https://theVall@bitbucket.org/theVall/basicvolumeraycaster.git
cd basicvolumeraycaster
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=/path/to/Qt/install
make -j `nproc`
```
Make sure to replace the CMAKE_PREFIX_PATH with the path to your Qt install directory, e.g. ```/home/username/Qt/5.11.2/gcc_64/```

For Tobii eye tracking support, you also need to link the respective library.

# Confirmed to build/run on the following configurations #

* NVIDIA Maxwell & Pascal, AMD Fiji & Vega, Intel Gen9 GPU & Skylake CPU
* GCC 5.3.1 & 7.3.0, Visual Studio 2015 (v140), Clang 6.0
* Qt 5.11.2
* CMake 3.10.2 & 3.12.2
* Tobii Pro Spectrum eye tracker

# Screenshots #

![eye-tracker](hhttps://github.com/vbruder/FoveatedVolumeRendering/blob/master/resources/eyeTracker.jpg)

![richtmyer](https://github.com/vbruder/FoveatedVolumeRendering/blob/master/resources/richtmyer.png)

![vortex](https://github.com/vbruder/FoveatedVolumeRendering/blob/master/resources/vortex.png)

# Planned changes/extensions #

* Support for head mounted devices with integrated eye tracking.

# License #

Copyright (C) 2017-2019 Valentin Bruder vbruder@gmail.com

This software is licensed under [LGPLv3+](https://www.gnu.org/licenses/lgpl-3.0.en.html).

# Credits #
	
  * Color wheel from Mattia Basaglia's Qt-Color-Widgets: https://github.com/mbasaglia/Qt-Color-Widgets
  * OpenCL utils based on Erik Smistad's OpenCLUtilityLibrary: https://github.com/smistad/OpenCLUtilityLibrary
  * Transfer function editor based on Qt sample code.
  * Stippling is based on O. Deussen, M. Spicker, Q. Zeng's Weighted Linde-Buzo-Gray Algorithm: http://graphics.uni-konstanz.de/publikationen/Deussen2017LindeBuzoGray/index.html
