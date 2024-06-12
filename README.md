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

Please using [this branch](https://github.com/wangzhezhe/UCV/tree/exp_critical_point) to build the plugin.

We first need to build the ParaView on a specific platform. The user can refer to [this document](https://gitlab.kitware.com/paraview/paraview/blob/master/Documentation/dev/build.md) to build the ParaView manually.

When installing the uncertainty critical point filter through ParaView, we need to set

`DParaView_DIR` and `DVTKm_DIR` as the proper path. In this case, we use the vtk-m embedded as the third-party dependency of the ParaView.

The current `ubuntu_cpu.sh` provides an example of how to set the ParaView path and install the filter through the ParaView plugin.

After installing the ParaView plugin, the user could load the associated plugin through the ParaView filter and run the loaded filter.

The [demo video](https://drive.google.com/file/d/1GS0OJW_HQWHP5HyS8xV0cxbDHKK_sRgR/view?usp=sharing) shows how to run the filter through the ParaView plugin.

### Using different approaches to compute critical points

In the paper, we describes multiple appraoches to compute the critical point, including monte carlo sampling, closed form appraoch based on indepedent uniform distribution, multi-variant gaussian distribution, histogram and kernel density estimation based on Epanechnikov. 
 - The `LoadEnsAndProcessByUniform` can load the ensemble data and uses the uniform kernel to compute the critical point (minimal value).
 - The `LoadEnsAndProcessByMC` can load the ensmeble data and uses Monte Carlo Sampling approach to compute the critical point (minimal value).
 - The `LoadEnsAndProcessByMVG` can load the ensmeble data and uses multi-variant gauusian approach to compute the critical point (minimal value).
 - The `LoadWeatherProcessByEpanech` can load the ensmeble data and uses Epanechnikov compute the critical point (minimal value).

These testing files can be used for both redsea data and weather data described in paper, in particular, [this script](https://github.com/wangzhezhe/UCV/blob/exp_critical_point_noplugin/exp_scripts/frontier_run_redsea.sh) and [this script](https://github.com/wangzhezhe/UCV/blob/exp_critical_point_noplugin/exp_scripts/frontier_run_weather.sh) shows details and associated parameters of how to run testing fiels on Frontier platform for redsea data and weather data, respectively.

