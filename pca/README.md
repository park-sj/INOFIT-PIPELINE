# Data-Driven Nonrigid ICP
This is a C++ implementation of mesh registration based on data-driven nonrigid ICP scheme.
Given face meshes with the same topology, called template meshes, per-region deformation prior is computed by PCA.
The deformation prior is then used to register the template meshes onto a target mesh.
#### Reference paper
Schneider, David C., and Peter Eisert. "Fast nonrigid mesh registration with a data-driven deformation prior." 2009 IEEE 12th International Conference on Computer Vision Workshops, ICCV Workshops. IEEE, 2009.
## Dependencies
1. CGAL 5.0.2
2. Eigen 3.3.7
3. Boost 1.71.0
4. jsoncpp 1.7.4
5. CMake 3.16.3
#### Optional
1. CUDA 10.2
2. OpenCL 1.2, and ViennaCL 1.7.1
3. CGAL 5.1, OpenGR, and libpointmatcher 1.1.0
4. OpenMesh 8.1
## Compliation
The CMakeLists.txt file is still under development.
#### For Ubuntu and Mac
```
cd build
cmake ..
make
cd ..
```
You can also specify CMAKE_BUILD_TYPE.
```
cmake .. -DCMAKE_BUILD_TYPE=Debug
```
```
cmake .. -DCMAKE_BUILD_TYPE=Release
```
Default build type is RelWithDebInfo.   
   
If CMake cannot find jsoncpp or ViennaCL automatically, jsoncpp_INCLUDE_DIR or ViennaCL_INCLUDE_DIR can be set explicitly.
```
cmake .. -Djsoncpp_INCLUDE_DIR=/usr/local/include -DViennaCL_INCLUDE_DIR=/usr/include
```
If CMake cannot find other packages or libraries automatically, please see CMakeLists.txt to modify certain lines.
#### For Windows
```
cd build
cmake ..
msbuild Data_Driven_Nonrigid_ICP.sln
cd ..
```
You can also specify CMAKE_BUILD_TYPE.
```
cmake .. -DCMAKE_BUILD_TYPE=Debug
```
```
cmake .. -DCMAKE_BUILD_TYPE=Release
```
Default build type is RelWithDebInfo.

The command "msbuild" might need additional parameters as follows.
```
msbuild Data_Driven_Nonrigid_ICP.sln /property:Configuration=Release /property:Platform=x64
```
If CMake cannot find other packages or libraries automatically, please see CMakeLists.txt to modify certain lines.
Also note that none of additional dependencies were tested on Windows environment.
#### Tested environment
- macOS 10.14.5, packages installed via homebrew, with OpenGR and libpointmatcher installed manually
- Ubuntu 20.04, packages installed via apt, with CUDA installed manually
- Windows 10, packages installed via vcpkg
## Usage
- preprocess
```
build/bin/app/preprocess params/preprocess/A.json
```
This program computes per-region PCA from template meshes.   
You can adjust the configuration for preprocess by modifying params/preprocess/~~.json file.
- registration
```
build/bin/app/registration params/registration/A.json
```
This program registers the template meshes onto a target mesh and generates result mesh file.   
You can adjust the configuration for registration by modifying params/registration/~~.json file.
