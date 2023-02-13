This repo contains some code that use the vtk-m to implement the uncertainty algorithm.

The original name is the UCV, the updated name might be UVM (uncertianty visualizatio based on vkt-m), which might be more make sense.

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
$ ./ucv_reduce_umc ../../../../dataset/beetle_496_832_832.vtk ground_truth uni 4 900
```

using the indepednet gaussian distribution

```
$ ./ucv_reduce_umc ../../../../dataset/beetle_496_832_832.vtk ground_truth ig 4 900
```

using the multivariant gaussian distribution (only works for #vertexies=4 currently)

```
$ ./ucv_reduce_umc ../../../../dataset/raw_data_128_208_208.vtk instance mg 4 900
```

### TODO

support multi blocks version