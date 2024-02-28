The original code comes from 
https://gitlab.kitware.com/nrushad2001/vtk-m/-/blob/FiberUncertaintyVis/vtkm/filter/uncertainty/Fiber.h?ref_type=heads

We compile the code separately to make it work based on the latest vtkm

Run example:

Extract ens members from raw data
./mac/install/UCV/uncertainty/testing/DecomposeSuperNovaToEns /Users/zw1/Documents/Uncertainty/src/dataset/supernova_visit_100_100_100.vtk 4 supernova_25_data


Run fiber exp
./mac/install/UCV/uncertainty/testing/TestSuperNova ./supernova_25_data 64

TODO:

Divide the Fiber into two filter, one is the closedForm another is MonteCarlo

Run on readsea data, using mean+-stdev to get min and max

Associated paper is this one
https://diglib.eg.org/bitstream/handle/10.2312/evs20211053/043-047.pdf?sequence=1&isAllowed=y

Get the performance number for GPU on andes machine maybe, show that it can run on different backend (andes has a low gpu version)

Fix the issues, when ens_one min max equals to ens_two min max, the output between closed form and monte carlo are different (one is zero (closed form) and another is 1 (montecarlo))