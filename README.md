This repo contains some code that use the vtk-m to implement the uncertainty algorithm.

### Build

install example on ubuntu machine:

```
clone the repo
$ cd install_scripts
$ /bin/bash ubuntu_cpu.sh 
...
try to add library path by executing:
export LD_LIBRARY_PATH=/home/zw/cworkspace/UCV/install_scripts/ubuntu_cpu/install/vtk-m/lib:${LD_LIBRARY_PATH}]
$ export LD_LIBRARY_PATH=/home/zw/cworkspace/UCV/install_scripts/ubuntu_cpu/install/vtk-m/lib:${LD_LIBRARY_PATH}]
```

### Example

extracting the ensemble dataset:

```
$ ./ubuntu_cpu/install/UCV/ucv_extract ../dataset/raw_data_123_208_208.vtk instance
```

computing the probability marching cube:

```
$ ./ubuntu_cpu/install/UCV/ucv_umc ../dataset/raw_data_123_208_208_Derived.vtk instance
```