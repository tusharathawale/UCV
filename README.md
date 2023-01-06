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
$ ./ubuntu_cpu/install/UCV/ucv_umc ../dataset/raw_data_123_208_208_Derived.vtk instance 900
```

or combining two steps together (size of hixel block is 4 and isovalue is 900)

```
$ ./ubuntu_cpu/install/UCV/ucv_reduce_umc ../dataset/raw_data_123_208_208.vtk instance 4 900
```

The generated dataset raw_data_123_208_208_Prob.vtk contains the field such as cross_prob and entropy.