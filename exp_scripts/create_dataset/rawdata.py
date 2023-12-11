import numpy as np
import vtk
import os
from vtk.util import numpy_support
import matplotlib.pyplot as plt
from vtkmodules.vtkCommonDataModel import vtkStructuredPoints


def writeStructuredDs(fname, ds):
    writer = vtk.vtkStructuredPointsWriter()
    writer.SetFileName(fname)
    writer.SetFileVersion(42)
    writer.SetInputData(ds)
    writer.Update()
    writer.Write() 


xdim=800
ydim=800
zdim=1

x,y = np.meshgrid(np.linspace(-1.0,1.0,xdim), np.linspace(-1.0,1.0,ydim))
d=np.sqrt(x*x+y*y)
sigma, mu = 1.0, 0.0

g=np.exp(-((d-mu)**2/(2.0*sigma**2)))
#print("2D gaussian like array")
#print(g)
#plt.matshow(g)
#plt.savefig("raw_data.png")

structured_dataset = vtkStructuredPoints()
structured_dataset.SetDimensions(xdim, ydim, zdim)
structured_dataset.SetOrigin(0, 0, 0)
#print(np.array(g).shape)
vtkArray = numpy_support.numpy_to_vtk(np.array(g).flatten())
vtkArray.SetNumberOfComponents(1)
vtkArray.SetName("TestField")

structured_dataset.GetPointData().AddArray(vtkArray)
structured_dataset.GetPointData().SetActiveScalars("TestField")

writeStructuredDs("RawdataPointScalar.vtk",structured_dataset)

#print(structured_dataset)
