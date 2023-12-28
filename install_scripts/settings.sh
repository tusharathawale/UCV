#!/bin/bash

VTKM_REPO=https://gitlab.kitware.com/vtk/vtk-m.git
#this one is the updated in deviceadaptor of kokkos
#which remove these private command
#VTKM_REPO=https://gitlab.kitware.com/zhe.wang/vtk-m.git
#VTKM_VERSION=vtkm_ucv
VTKM_VERSION=master


#using binary file
KOKKOS_REPO=https://github.com/kokkos/kokkos.git
KOKKOS_VERSION=3.7.01

#EasyLinalg repo
EASY_LINALG_REPO=https://github.com/wangzhezhe/EasyLinalg.git
EASY_LINALG_VERSION=main