# Description #

This is an interactive, cross platform volume raycaster based on the OpenCL compute API.
It features early ray termination, image order empty space skipping, local illumination and ambient occlusion based on central differences.
It can display volume (timeseries) data sets in .dat/.raw format and uses the [Qt](https://www.qt.io) framework for the GUI. 

# Setup #

To compile/run the code you need:

*  An OpenCL 1.2 capable device and drivers with image support. If you are compiling under Windows, make sure to link against the OpenCL/OpenGL libraries provided by your vendor. See the description at [StreamHPC](https://streamhpc.com/blog/2015-03-16/how-to-install-opencl-on-windows/) for details on how to set up Visual Studio. 
*  A reasonably new Qt version is required, tested with Qt 5.9.x and Qt 5.10. Under Windows, make sure to copy the required Qt .dlls in the directory next to your executable.  
*  A compiler capable of C++14, tested with gcc 7.1.0 and Visual Studio 2017 (MSVC v141).

# Screenshots #

![2017-12-19-nova](https://bytebucket.org/theVall/basicvolumeraycaster/raw/b29bb112fdde3784923e22f35ef56d7d9408b6f6/screenshots/2017-12-19-nova.png)

![2017-12-19-cham](https://bytebucket.org/theVall/basicvolumeraycaster/raw/b29bb112fdde3784923e22f35ef56d7d9408b6f6/screenshots/2017-12-19-cham.png)

# Planned changes/extensions #

*  Switch to CMake as build system

# License #

Copyright (C) 2017-2018 Valentin Bruder vbruder@gmail.com

This software is licensed under [LGPLv3+](https://www.gnu.org/licenses/lgpl-3.0.en.html).

# Credits #
	
  * Color wheel from Mattia Basaglia's Qt-Color-Widgets: https://github.com/mbasaglia/Qt-Color-Widgets
  * OpenCL utils based on Erik Smistad's OpenCLUtilityLibrary: https://github.com/smistad/OpenCLUtilityLibrary
  * Transfer function editor based on Qt sample code.