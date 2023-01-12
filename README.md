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


using the uniform distribution
```
$ ./ucv_reduce_umc ../../../../dataset/raw_data_128_208_208.vtk instance uni 4 900
```

using the indepednet gaussian distribution

```
$ ./ucv_reduce_umc ../../../../dataset/raw_data_128_208_208.vtk instance ig 4 900
```

using the multivariant gaussian distribution (only works for hixlesize=4 currently)

```
$ ./ucv_reduce_umc ../../../../dataset/raw_data_128_208_208.vtk instance mg 4 900
```

