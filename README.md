# README #

This is a volume raycaster based on OpenCL featuring early ray termination, image order empty space skipping and local illumination based on central differences.
It can display (timeseries) volumes in dat-raw format and uses the Qt framework for the GUI. 

To compile/run the code you need an OpenCL 1.2 capable device and drivers.
A reasonably new Qt version is required (tested with Qt 5.9.0) as well as a compiler capable of C++14 (tested with gcc 7.1.0).

![2017-12-19-nova](https://bytebucket.org/theVall/basicvolumeraycaster/raw/b29bb112fdde3784923e22f35ef56d7d9408b6f6/screenshots/2017-12-19-nova.png)

![2017-12-19-cham](https://bytebucket.org/theVall/basicvolumeraycaster/raw/b29bb112fdde3784923e22f35ef56d7d9408b6f6/screenshots/2017-12-19-cham.png)
