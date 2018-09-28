# Description #

An interactive, cross platform volume raycaster based on the OpenCL compute API.
It features early ray termination, object and image order empty space skipping, local illumination, and various gradient based shading techniques.
It can display volume (timeseries) data sets and uses the [Qt](https://www.qt.io) framework for the GUI. 

# Setup #

To compile/run the code you need:

* An OpenCL 1.2 capable device and drivers/libraries with image support. It is recommended to update your GPU driver before building/running.
* Qt, at leat version 5.6 is required.
* A compiler capable of c++14
* CMake version 3.9 or higher 

# Confirmed to build/run on the following configuration #

* NVIDIA Kepler & Pascal, AMD Fiji & Vega, Intel Gen9 GPU & Skylake CPU
* GCC 5.3.1 & 7.3.0, Visual Studio 2015 (v140), Clang 6.0
* Qt 5.11.2
* CMake 3.10.2 & 3.12.2

# Screenshots #

![2017-12-19-nova](https://bytebucket.org/theVall/basicvolumeraycaster/raw/6d3ef5483cd67d8a6416620887a19d36ca6e4d67/screenshots/2018-04-23-nova.png)

![2017-12-19-cham](https://bytebucket.org/theVall/basicvolumeraycaster/raw/6d3ef5483cd67d8a6416620887a19d36ca6e4d67/screenshots/2018-04-23-cham.png)

# Planned changes/extensions #

*  Out of core rendering for timeseries data on dGPUs.

# License #

Copyright (C) 2017-2018 Valentin Bruder vbruder@gmail.com

This software is licensed under [LGPLv3+](https://www.gnu.org/licenses/lgpl-3.0.en.html).

# Credits #
	
  * Color wheel from Mattia Basaglia's Qt-Color-Widgets: https://github.com/mbasaglia/Qt-Color-Widgets
  * OpenCL utils based on Erik Smistad's OpenCLUtilityLibrary: https://github.com/smistad/OpenCLUtilityLibrary
  * Transfer function editor based on Qt sample code.