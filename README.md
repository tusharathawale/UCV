The original code comes from 
https://gitlab.kitware.com/nrushad2001/vtk-m/-/blob/FiberUncertaintyVis/vtkm/filter/uncertainty/Fiber.h?ref_type=heads

We compile the code separately to make it work based on the latest vtkm

Run example:

Extract ens members from raw data
./mac/install/UCV/uncertainty/testing/DecomposeSuperNovaToEns /Users/zw1/Documents/Uncertainty/src/dataset/supernova_visit_100_100_100.vtk 4 supernova_25_data


Run fiber exp
./mac/install/UCV/uncertainty/testing/TestSuperNova ./supernova_25_data 64

TODO:
Double check why the ouput is zero
Use the supernova data for testing (Nickle and Iron)
