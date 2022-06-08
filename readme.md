LSM3D

The repository contains code for an implementation of the particle level set method. The code is developed as part of a Master's thesis by Sommersel (2022). The main computations are programmed in C++, while the post-processing is programmed in Python.

main.cpp
This file contains the main structure of the particle level set method and reinitialization. The computation of
the different measures of error are also included.

initialization
These files contain code to create the initial signed distance field with reference from a sphere surface and a
function to save the signed distance field to a .txt-file.

particleLSM
These files contain code to set up the particles in the particle level set method, use these particles to correct the interface, and a function to save the particle coordinates to a .txt-file. particleLSM.h is given first and contains the header file, while particleLSM.cpp contains the whole implementation.

schemes
These files contain all the numerical schemes used in the implemented particle level set method.

testCases
These files contain code to generate the velocity fields and the computation of the error measures.

vectorUtilities
These files contain code used for various vector operations. vectorUtilities.h is given first and contains the header file, while vectorUtilities.cpp contains the whole implementation.

plotter.py
This file contains the Python code to find the interface with the marching cubes algorithm using the scikit-image Python library, and also code to plot the interface and the particles used in the particle level set method.
