This branch contains the code that uses the vtk-m to implement the uncertainty algorithm to find the critical points.

The details can be found in the following paper:

Uncertainty Visualization of Critical Points of 2D Scalar Fields for Parametric and Nonparametric Probabilistic Models

### Build testing files

The installing scripts are listed under the `UCV/install_scripts`.

For example, on the mac platform, we need to execute 

`/bin/bash mac.sh`

under the install_scripts folder. This script will create a folder called mac. Under this folder, the dependency, namely the vtk-m, will be downloaded under the `src` folder and installed into a specific folder. The executable file will be installed at `UCV/install_scripts/mac/install/UCV`

The `exp_scripts/frontier_run_weather.sh` is the script to run experiments on the Frontier supercomputer.

### Build ParaView plugin

We first need to build the ParaView on a specific platform. The user can refer to [this document](https://gitlab.kitware.com/paraview/paraview/blob/master/Documentation/dev/build.md) to build the ParaView manually.

When installing the uncertainty critical point filter through ParaView, we need to set

`DParaView_DIR` and `DVTKm_DIR` as the proper path. In this case, we use the vtk-m embedded as the third-party dependency of the ParaView.

The current `ubuntu_cpu.sh` provides an example of how to set the ParaView path and install the filter through the ParaView plugin.

After installing the ParaView plugin, the user could load the associated plugin through the ParaView filter and run the loaded filter.

The [demo video](https://drive.google.com/file/d/1GS0OJW_HQWHP5HyS8xV0cxbDHKK_sRgR/view?usp=sharing) shows how to run the filter through the ParaView plugin.
