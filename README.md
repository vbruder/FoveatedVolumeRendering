# Description #

An Voronoi-based foveated volume renderer. Check out our EuroVis 2019 short paper for more details: https://vbruder.github.io/publication/bruder-sbfwe-19/bruder-sbfwe-19.pdf

The volume renderer is based on my OpenCL volume renderer **VolumeRendererCL**.
Check out the [repo](https://github.com/vbruder/VolumeRendererCL) for more details on the code, supported data sets and how to build. 

The Weighted-Linde-Buzo-Gray algorithm can be used to generate the sampling mask in a preprocessing step.
Use the `lbg-stippling` to generate a custom sampling mask or use the pre-created one in `resources/mc1024_92702.zip` (fixed to an output size of 1024x1024px).

For Tobii eye tracking support, you also need to link the respective library.

## Screenshots ##

![eye-tracker](https://github.com/vbruder/FoveatedVolumeRendering/blob/master/resources/eyeTracker.jpg)
![richtmyer](https://github.com/vbruder/FoveatedVolumeRendering/blob/master/resources/richtmyer.png)
![vortex](https://github.com/vbruder/FoveatedVolumeRendering/blob/master/resources/vortex.png)

*Copyright belongs to the Eurographics Association.*

## License ##

Copyright (C) 2017-2019 Valentin Bruder vbruder@gmail.com

This software is licensed under [LGPLv3+](https://www.gnu.org/licenses/lgpl-3.0.en.html).

## Credits ##
	
  * Color wheel from Mattia Basaglia's Qt-Color-Widgets: https://github.com/mbasaglia/Qt-Color-Widgets
  * OpenCL utils based on Erik Smistad's OpenCLUtilityLibrary: https://github.com/smistad/OpenCLUtilityLibrary
  * Transfer function editor based on Qt sample code.
  * Stippling is based on O. Deussen, M. Spicker, Q. Zeng's [Weighted Linde-Buzo-Gray Algorithm](http://graphics.uni-konstanz.de/publikationen/Deussen2017LindeBuzoGray/index.html)
